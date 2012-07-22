module topog

  use util, only: streq, s2s
  use ncutils, only: nc_open, nc_get_dim, nc_create, nc_def_dim, handle_err, nc_get_text
  use ncutils, only: nc_close, nc_enddef, nc_def_double, nc_get_double, nc_put_double
  use intsubs, only: xintp

  use box, only: box_type
  use grid, only: grid_type

  implicit none

  private

  type topog_type

     double precision, allocatable :: dtop(:,:)
     double precision:: davg
     double precision, allocatable :: ddyn(:,:)
     integer :: k_topo

     ! dtop is the (nxp,nyp) array of topography at the
     ! bottom of the domain, tabulated at p points (m)
     ! davg is its average value (m)
     ! ddyn is the dynamic (rescaled) topography (s^-1)

  end type topog_type

  public topout_nc
  public topog_type
  public load_topog
  public check_atmos_topog_over_ocean
  public print_topog_summary

contains

  subroutine load_topog(b, foo, topo, indir)
    type(box_type), intent(in) :: b
    character, allocatable, intent(in) :: foo(:)
    type(topog_type), intent(inout) :: topo
    character (len=*), intent(in) :: indir

    integer :: tempid, idim, jdim
    double precision :: davg
    integer :: j

    character subnam*(*)
    parameter ( subnam = 'load_topog' )

    allocate(topo%ddyn(b%nxp,b%nyp))
    allocate(topo%dtop(b%nxp,b%nyp))
    if (b%dz_sign == 1) then       
       topo%k_topo = 1
    else
       topo%k_topo = b%nl
    endif

    if (streq(foo,'flat')) then
       topo%dtop(:,:) = 0.0d0
    else if (streq(foo, 'mid_atlantic')) then
       call set_mid_atlantic_ridge(b, topo)
    else if (streq(foo, 'rocky')) then
       call rocky_mountains(topo, b)
    else
       ! Case of topography to be read from a file
       ! -----------------------------------------
       ! Read ocean topography from an existing netCDF file
       tempid = nc_open(trim(indir)//"/"//s2s(foo), subnam)
       ! First check that the dimensions are correct
       idim = nc_get_dim(tempid, 'xp', subnam)
       jdim = nc_get_dim(tempid, 'yp', subnam)
       if (idim /= b%nxp) then
          print *,' topography netCDF file error'
          print *,' xp dimension not equal to nxp:'
          print *,' dimensions are: ',idim,b%nxp
          print *,' Program terminates in topset'
          stop
       endif
       if (jdim /= b%nyp) then
          print *,' topography netCDF file error'
          print *,' yp dimension not equal to nyp:'
          print *,' dimensions are: ',jdim,b%nyp
          print *,' Program terminates in topset'
          stop
       endif
       ! Now read ocean topography data
       topo%dtop(:,:) = nc_get_double(tempid, 'dtop', b%nxp, b%nyp, subnam)
       call nc_close(tempid, subnam)
    endif

    ! If cyclic and not flat-bottom, check cyclicity of topography
    ! ------------------------------------------------------------
    if (b%cyclic) then
       do j=1,b%nyp
          if (topo%dtop(1,j) /= topo%dtop(b%nxp,j)) then
             print *,' *** WARNING *** problem with specified topography'
             print *,' Topography not exactly cyclic for j = ',j
             print *,' dtop values are: ',topo%dtop(1,j),topo%dtop(b%nxp,j)
             stop
          endif
       enddo
    endif

    davg = xintp(topo%dtop, b%nxp, b%nyp)
    topo%davg = davg*b%norm
    topo%ddyn(:,:) = (b%fnot/b%h(topo%k_topo))*topo%dtop(:,:)

  end subroutine load_topog

  pure subroutine rocky_mountains(topat, ga)
    type(topog_type), intent(inout) :: topat
    type(box_type), intent(in) :: ga

    double precision :: dxlo, dxhi, dcent, dwidth
    integer :: itlo, ithi, i, j

    dxlo = 1800.0d3
    dxhi = 5400.0d3
    itlo = 1 + nint( dxlo/ga%dx )
    ithi = 1 + nint( dxhi/ga%dx )
    dcent = 0.5d0*( dxlo + dxhi )
    dwidth = 0.5d0*( dxhi - dxlo )
    do j=1,ga%nyp
       do i=1,ga%nxp
          topat%dtop(i,j) = 1000.0d0*( 1.0d0 - abs(ga%xp(i)-dcent)/dwidth )
          topat%dtop(i,j) = max( 0.0d0, topat%dtop(i,j) )
       enddo
    enddo

  end subroutine rocky_mountains


  pure subroutine set_mid_atlantic_ridge(go, topoc)

    type(box_type), intent(in) :: go
    type(topog_type), intent(inout) :: topoc

    double precision :: dxlo, dxhi, dcent, dwidth, xrel
    integer :: i, j, ithi, itlo

    dxlo = 2000.0d3
    dxhi = 2600.0d3
    itlo = 1 + nint( dxlo/go%dx )
    ithi = 1 + nint( dxhi/go%dx )
    dcent = 0.5d0*( dxlo + dxhi )
    dwidth = 0.5d0*( dxhi - dxlo )
    topoc%dtop(:,:) = 0.0d0
    do j=1,go%nyp
       do i=1,go%nxp
          xrel = go%xp(i) - go%xp(1)
          topoc%dtop(i,j) = 1000.0d0*( 1.0d0 - abs(xrel-dcent)/dwidth )
       enddo
    enddo
    topoc%dtop(:,:) = max( 0.0d0, topoc%dtop(:,:) )

  end subroutine set_mid_atlantic_ridge

  subroutine topout_nc (b, topo, outdir, filename)

    type(box_type), intent(in) :: b
    type(topog_type), intent(in) :: topo
    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename

    character :: subnam*(*)
    parameter ( subnam = 'topout_nc' )

    integer :: xpdim,ypdim
    integer :: xp_id,yp_id,dtop_id
    integer :: topog_id
    
    topog_id = nc_create(outdir, filename, subnam)

    !! Initialise netcdf file
    xpdim = nc_def_dim(topog_id, 'xp', b%nxp, subnam)
    ypdim = nc_def_dim(topog_id, 'yp', b%nyp, subnam)

    xp_id = nc_def_double(topog_id, 'xp', xpdim, 'km', subnam)
    yp_id = nc_def_double(topog_id, 'yp', ypdim, 'km', subnam)
    dtop_id = nc_def_double(topog_id, 'dtop', (/ xpdim, ypdim /), 'm', subnam)

    !! Leave definition mode: entering data mode.
    call nc_enddef(topog_id, subnam)

    !! Convert from m -> km
    call nc_put_double(topog_id, xp_id, 1.0d-3*( b%xp(:) - b%x0 ), subnam)
    call nc_put_double(topog_id, yp_id, 1.0d-3*( b%yp(:) - b%y0 ), subnam)

    !! Write out the topography
    call nc_put_double(topog_id, dtop_id, topo%dtop, subnam)

    call nc_close(topog_id, subnam)

  end subroutine topout_nc

  subroutine check_atmos_topog_over_ocean(g, ga, topat)

    type(grid_type), intent(in) :: g
    type(box_type), intent(in) :: ga
    type(topog_type), intent(in) :: topat

    integer :: i, j

    ! Check there is no atmosphere topography over
    ! the ocean (including ocean boundary points)
    ! --------------------------------------------
    ! This check may need to be switched off during
    ! the phase of producing suitable topography,
    ! but should be turned on during normal running

    do j=g%ny1,g%ny1+g%nyaooc
       do i=g%nx1,g%nx1+g%nxaooc
          if (topat%dtop(i,j) /= 0.0d0) then
             print *,' Nonzero atmosphere topography over ocean'
             print *,' Problem occurs at i, j, xpa(km), ypa(km) = ', &
                  i,j,1.0d-3*ga%xp(i),1.0d-3*ga%yp(j)
             print '(a,2f11.2)', '  Ocean limits in x (km) are: ', &
                  1.0d-3*ga%xp(g%nx1),1.0d-3*ga%xp(g%nx1+g%nxaooc)
             print '(a,2f11.2)', '  Ocean limits in y (km) are: ', &
                  1.0d-3*ga%yp(g%ny1),1.0d-3*ga%yp(g%ny1+g%nyaooc)
             print *,' Program terminates in topset'
             stop
          endif
       enddo
    enddo

  end subroutine check_atmos_topog_over_ocean

  subroutine print_topog_summary(b, topog, title)
    type(box_type), intent(in) :: b
    type(topog_type), intent(in) :: topog
    character (len=*), intent(in) :: title

    print *, "Topog Summary: ", title
    print *, "Min/Max/Avg (m): ", minval(topog%dtop), maxval(topog%dtop), topog%davg
    print *, "Min/max frac. topog: ", minval(topog%dtop)/b%h(topog%k_topo), maxval(topog%dtop)/b%h(topog%k_topo)

  end subroutine print_topog_summary

end module topog
