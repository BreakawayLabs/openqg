module tavsubs

  use box, only: box_type
  use mixed, only: mixed_type
  use qg, only: qg_type
  use ekman, only: ekman_type
  use ncutils, only: nc_enddef, nc_def_dim, nc_def_double, nc_put_double
  use ncutils, only: nc_create, nc_close
  use clock, only: clock_type, days_to_steps

  implicit none

  type basin_averages_type

     logical :: active = .false.

     double precision, allocatable :: txav(:,:), tyav(:,:)
     double precision, allocatable :: wtav(:,:)
     double precision, allocatable :: fmav(:,:)
     double precision, allocatable :: stav(:,:)
     double precision, allocatable :: pav(:,:,:),qav(:,:,:)
     ! txatav, tyatav are dynamic stress components          (m^2 s^-2)
     ! wtatav is the Ekman velocity at atmospheric T points  (m s^-1)
     ! fmatav is the net atmos. mixed layer forcing at T pts (W m^-2)
     ! astav  is the atmos. mix. layer temperature (anomaly) (K)
     ! patav  is the dynamic pressure in each layer          (m^2 s^-1)
     ! qatav  is the vorticity in each layer                 (s^-1)

     double precision, allocatable :: uuf(:,:),tuf(:,:),utuf(:,:)
     double precision, allocatable :: vvf(:,:),tvf(:,:),vtvf(:,:)

     integer :: ntav, nmid

     ! to the running totals for the atmosphere
     integer :: nsum

     character (len=64) :: filename

  end type basin_averages_type

  public basin_averages_type
  public init_basin_average
  public update_avg
  public avg_step
  public tavout

contains

  type(basin_averages_type) function init_basin_average(b, dtav, clk, filename)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: dtav
    type(clock_type), intent(in) :: clk
    character (len=*), intent(in) :: filename

    init_basin_average%active = (dtav > 0.0d0)
    call finish_loading_avg(b, dtav, clk, filename, init_basin_average)
    if ( init_basin_average%active ) then
       write(*,204) '  Averaging int. (day) = ',dtav
       write(*,201) '  No. of steps in avg. inter. = ',init_basin_average%ntav
       if ( mod(clk%ntsrun,init_basin_average%ntav).ne.0 .or. mod(init_basin_average%ntav,2).ne.0 ) then
          print *,' Unsuitable choice of dtav; program stops'
          stop 1
       endif
    endif

201 format(a,9i13)
204 format(a,9f13.4)

  end function init_basin_average

  subroutine finish_loading_avg(b, dt, clk, filename, avg)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: dt
    type(clock_type), intent(in) :: clk
    character (len=*), intent(in) :: filename
    type(basin_averages_type), intent(inout) :: avg

    if (avg%active) then
       avg%ntav = days_to_steps(dt, clk)
       avg%nmid = avg%ntav/2

       allocate(avg%txav(b%nxp,b%nyp))
       allocate(avg%tyav(b%nxp,b%nyp))
       allocate(avg%wtav(b%nxt,b%nyt))
       allocate(avg%fmav(b%nxt,b%nyt))
       allocate(avg%stav(b%nxt,b%nyt))
       allocate(avg%pav(b%nxp,b%nyp,b%nl))
       allocate(avg%qav(b%nxp,b%nyp,b%nl))
       allocate(avg%uuf(b%nxp,b%nyt))
       allocate(avg%tuf(b%nxp,b%nyt))
       allocate(avg%utuf(b%nxp,b%nyt))
       allocate(avg%vvf(b%nxt,b%nyp))
       allocate(avg%tvf(b%nxt,b%nyp))
       allocate(avg%vtvf(b%nxt,b%nyp))

       avg%filename = filename

       ! Windstress fields (p-grid fields)
       avg%txav(:,:) = 0.0d0
       avg%tyav(:,:) = 0.0d0

       ! Ekman velocity, forcing fields & mixed
       ! layer temperatures (T-grid fields)
       avg%wtav(:,:) = 0.0d0
       avg%fmav(:,:) = 0.0d0
       avg%stav(:,:) = 0.0d0

       ! Dynamic pressure and vorticity in QG layers
       avg%pav(:,:,:) = 0.0d0
       avg%qav(:,:,:) = 0.0d0

       ! Fields for mixed layer temperature advection calculations
       avg%uuf (:,:) = 0.0d0
       avg%tuf (:,:) = 0.0d0
       avg%utuf(:,:) = 0.0d0
     
       avg%vvf (:,:) = 0.0d0
       avg%tvf (:,:) = 0.0d0
       avg%vtvf(:,:) = 0.0d0

       ! Counter
       avg%nsum = 0
    endif

  end subroutine finish_loading_avg

  logical pure function avg_step(nt, avg)

    integer, intent(in) :: nt
    type(basin_averages_type), intent(in) :: avg

    if (avg%active) then 
       avg_step = mod(nt, avg%ntav) == avg%nmid
    else 
       avg_step = .false.
    endif

  end function avg_step

  subroutine update_avg(qg, b, ek, st, fnet, avg)
    
    type(qg_type), intent(in) :: qg
    type(box_type), intent(in) :: b
    type(ekman_type), intent(in) :: ek
    type(mixed_type), intent(in) :: st
    double precision, intent(in) :: fnet(qg%b%nxt,qg%b%nyt)
    type(basin_averages_type), intent(inout) :: avg

    integer :: i,j
    double precision :: uvgfac,rhf0hm,uuf(b%nxp),tuf(b%nxp)
    double precision :: utuf(b%nxp), vvf,tvf,vtvf

    uvgfac = st%ycexp*b%rdxf0
    rhf0hm = 0.5d0/(b%fnot*b%hm)

    ! Windstress fields
    avg%txav(:,:) = avg%txav(:,:) + ek%taux(:,:)
    avg%tyav(:,:) = avg%tyav(:,:) + ek%tauy(:,:)

    ! Ekman velocity
    avg%wtav(:,:) = avg%wtav(:,:) + ek%wekt(:,:)

    if (st%active) then
       avg%fmav(:,:) = avg%fmav(:,:) + fnet(:,:)
       avg%stav(:,:) = avg%stav(:,:) + st%data(:,:)

    ! Fields for mixed layer temperature advection calculations
    ! C-grid advection, second order accurate,
    ! borrowed from subroutine omladv in omlsubs.f

    ! Zonal advection
    ! In finite box case, no normal flux for temperature, so deem u to
    ! be zero. Deem boundary temperature to be that of internal point.
       do j=1,b%nyt
          ! Western boundary
          if (b%cyclic) then
             ! Zonally cyclic ocean
             uuf (1) = -uvgfac*( qg%p(1,j+1,1)-qg%p(1,j,1) ) &
                  +rhf0hm*(ek%tauy(1,j+1)+ek%tauy(1,j))
             tuf (1) = 0.5d0*( st%data(1,j) + st%data(b%nxt,j) )
             utuf(1) = uuf(1)*tuf(1)
          else
             ! Finite box ocean (no normal heat flux)
             uuf ( 1 ) = 0.0d0
             tuf ( 1 ) = st%data(1,j)
             utuf( 1 ) = 0.0d0
          endif
          ! Inner points.
          do i=2,b%nxp-1
             uuf (i) = -uvgfac*( qg%p(i,j+1,1)-qg%p(i,j,1) ) &
                  +rhf0hm*(ek%tauy(i,j+1)+ek%tauy(i,j))
             tuf (i) = 0.5d0*( st%data(i,j) + st%data(i-1,j) )
             utuf(i) = uuf(i)*tuf(i)
          enddo
          ! Eastern boundary
          if (b%cyclic) then
             ! Zonally cyclic ocean
             uuf (b%nxp) = uuf (1)
             tuf (b%nxp) = tuf (1)
             utuf(b%nxp) = utuf(1)
          else
             ! Finite box ocean (no normal heat flux)
             uuf (b%nxp) = 0.0d0
             tuf (b%nxp) = st%data(b%nxt,j)
             utuf(b%nxp) = 0.0d0
          endif
          do i=1,b%nxp
             avg%uuf (i,j) = avg%uuf (i,j) + uuf (i)
             avg%tuf (i,j) = avg%tuf (i,j) + tuf (i)
             avg%utuf(i,j) = avg%utuf(i,j) + utuf(i)
          enddo
       enddo

       ! Meridional advection
       ! Points excluding zonal boundaries
       do j=2,b%nyt
          do i=1,b%nxt
             vvf  =  uvgfac*( qg%p(i+1,j,1)-qg%p(i,j,1) ) &
                  -rhf0hm*(ek%taux(i+1,j)+ek%taux(i,j))
             tvf  = 0.5d0*( st%data(i,j) + st%data(i,j-1) )
             vtvf = vvf*tvf
             avg%vvf (i,j) = avg%vvf (i,j) + vvf
             avg%tvf (i,j) = avg%tvf (i,j) + tvf
             avg%vtvf(i,j) = avg%vtvf(i,j) + vtvf
          enddo
       enddo

       ! Zonal boundaries including corners; choice of BCs
       do i=1,b%nxt
          ! Southern boundary
          if (st%fluxs) then
             ! Advection consistent with an outflow across the
             ! southern boundary equal to the Ekman transport,
             ! carrying fluid of a specified temperature tsbdy.
             ! p contribution to vm vanishes because p is uniform along bdy.
             vvf  = -rhf0hm*( ek%taux(i+1,1) + ek%taux(i,1) )
             tvf  = 0.5d0*( st%data(i,1) + st%sval )
             vtvf = vvf*tvf
          else
             ! Advection (no normal mass flux => vm = 0)
             vvf  = 0.0d0
             tvf  = st%data(i,1)
             vtvf = 0.0d0
          endif
          avg%vvf (i, 1 ) = avg%vvf (i, 1 ) + vvf
          avg%tvf (i, 1 ) = avg%tvf (i, 1 ) + tvf
          avg%vtvf(i, 1 ) = avg%vtvf(i, 1 ) + vtvf

          ! Northern boundary
          if (st%fluxn) then
             ! Advection consistent with an outflow across the
             ! northern boundary equal to the Ekman transport,
             ! carrying fluid of a specified temperature tnbdy.
             ! p contribution to vp vanishes because p is uniform along bdy.
             vvf  = -rhf0hm*( ek%taux(i+1,b%nyt+1) + ek%taux(i,b%nyt+1) )
             tvf  = 0.5d0*( st%data(i,b%nyt) + st%nval )
             vtvf = vvf*tvf
          else
             ! Advection (no normal mass flux => vp = 0)
             vvf  = 0.0d0
             tvf  = st%data(i,b%nyt)
             vtvf = 0.0d0
          endif
          avg%vvf (i,b%nyp) = avg%vvf (i,b%nyp) + vvf
          avg%tvf (i,b%nyp) = avg%tvf (i,b%nyp) + tvf
          avg%vtvf(i,b%nyp) = avg%vtvf(i,b%nyp) + vtvf
       enddo
    endif

    ! Dynamic pressure and vorticity in QG layers
    avg%pav(:,:,:) = avg%pav(:,:,:) + qg%p(:,:,:)
    avg%qav(:,:,:) = avg%qav(:,:,:) + qg%q(:,:,:)

    ! Update counter
    avg%nsum = avg%nsum + 1

  end subroutine update_avg

  subroutine tavout (outdir, b, avg)

    ! Computes time averages of various quantities (mainly forcing
    ! fields) at the end of a run, and outputs a netCDF dump of them.

    character (len=*), intent(in) :: outdir
    type(box_type), intent(in) :: b
    type(basin_averages_type), intent(inout) :: avg

    character :: subnam*(*)
    parameter ( subnam = 'tavout' )

    integer :: i,j,k
    double precision :: rns,uptp(b%nxp,b%nyt),vptp(b%nxt,b%nyp)

    integer :: scadim,dims(2),pdims(3)
    integer :: xpdim,ypdim,ldim,xtdim,ytdim,xp_id,xt_id
    integer :: yp_id,yt_id,l_id,st_id,wekt_id,fm_id,tx_id
    integer :: ty_id,p_id,q_id,ut_id,vt_id
    double precision :: xx(b%nxp),yy(b%nyp),tmp(b%nl)
    integer :: avncid

    if (avg%active) then
       if ( avg%nsum.eq.0 ) then
          rns = 0.0d0
       else
          rns = 1.0d0/dble(avg%nsum)
       endif

       ! Windstress fields & wekpo (p-grid fields)
       avg%txav(:,:) = rns*avg%txav(:,:)
       avg%tyav(:,:) = rns*avg%tyav(:,:)

       ! Ekman velocity, forcing fields & mixed
       ! layer temperatures (T-grid fields)
       avg%wtav(:,:) = rns*avg%wtav(:,:)
       avg%fmav(:,:) = rns*avg%fmav(:,:)
       avg%stav(:,:) = rns*avg%stav(:,:)

       ! Dynamic pressure and vorticity in QG layers
       avg%pav(:,:,:) = rns*avg%pav(:,:,:)
       avg%qav(:,:,:) = rns*avg%qav(:,:,:)

       ! Fields for mixed layer temperature flux calculations
       ! Compute means, and derive eddy fluxes
       avg%uuf (:,:) = rns*avg%uuf (:,:)
       avg%tuf (:,:) = rns*avg%tuf (:,:)
       avg%utuf(:,:) = rns*avg%utuf(:,:)
       uptp(:,:) = avg%utuf(:,:) - avg%uuf(:,:)*avg%tuf(:,:)

       avg%vvf (:,:) = rns*avg%vvf (:,:)
       avg%tvf (:,:) = rns*avg%tvf (:,:)
       avg%vtvf(:,:) = rns*avg%vtvf(:,:)
       vptp(:,:) = avg%vtvf(:,:) - avg%vvf(:,:)*avg%tvf(:,:)

       avncid = nc_create(outdir, avg%filename, subnam)

       !! now initialise netCDF file
       scadim = nc_def_dim(avncid, 'scalar', 1, subnam)

       xpdim = nc_def_dim(avncid, 'xp', b%nxp, subnam)
       ypdim = nc_def_dim(avncid, 'yp', b%nyp, subnam)
       ldim = nc_def_dim(avncid, 'z', b%nl, subnam)

       xtdim = nc_def_dim(avncid, 'xt', b%nxt, subnam)
       ytdim = nc_def_dim(avncid, 'yt', b%nyt, subnam)

       xp_id = nc_def_double(avncid, 'xp', xpdim, 'km', subnam, 'X axis (p-grid)')
       xt_id = nc_def_double(avncid, 'xt', xtdim, 'km', subnam, 'X axis (T-grid)')
       yp_id = nc_def_double(avncid, 'yp', ypdim, 'km', subnam, 'Y axis (p-grid)')
       yt_id = nc_def_double(avncid, 'yt', ytdim, 'km', subnam, 'Y axis (T-grid)')
       l_id = nc_def_double(avncid, 'z', ldim, 'km', subnam, 'depth axis')

       dims = (/xtdim, ytdim/)
       st_id = nc_def_double(avncid, 'st', dims, 'K', subnam, 'surface temperature')
       wekt_id = nc_def_double(avncid, 'wekt', dims, 'm/s', subnam, 'Ekman velocity (T-grid)')
       fm_id = nc_def_double(avncid, 'fnet', dims, 'W/m^2', subnam, 'surface forcing (T-grid)')

       dims = (/xpdim, ypdim/)
       tx_id = nc_def_double(avncid, 'taux', dims, 'm^2/s^2', subnam, 'surface x stress (p-grid)')
       ty_id = nc_def_double(avncid, 'tauy', dims, 'm^2/s^2', subnam, 'surface y stress (p-grid)')
       
       pdims = (/xpdim, ypdim, ldim/)
       p_id = nc_def_double(avncid, 'p', dims, 'm^2/s^2', subnam, 'dynamic pressure')
       q_id = nc_def_double(avncid, 'q', dims, 's^-1', subnam, 'vorticity')

       dims = (/xpdim, ytdim/)
       ut_id = nc_def_double(avncid, 'uptp', dims, 'K.m/s', subnam, '<upTp>')

       dims = (/xtdim, ypdim/)
       vt_id = nc_def_double(avncid, 'vptp', dims, 'K.m/s', subnam, '<vpTp>')

       !! Leave definition mode: entering data mode.
       call nc_enddef(avncid, subnam)
    
       !!  Convert all x- and y-axes into km
       !!  Calculate x gridpoints and store in 'x' arrays
       !!  p-grid points
       do i=1,b%nxp
          xx(i) = 1.0d-3*( b%xp(i) - b%xp(1) )
       enddo
       call nc_put_double(avncid, xp_id, xx, subnam)
       !!  T-grid points
       do i=1,b%nxt
          xx(i) = 1.0d-3*( b%xt(i) - b%xp(1) )
       enddo
       call nc_put_double(avncid, xt_id, xx(:b%nxt), subnam)
       !!  Calculate y gridpoints and store in 'y' arrays
       !!  p-grid points
       do j=1,b%nyp
          yy(j) = 1.0d-3*( b%yp(j) - b%yp(1) )
       enddo
       call nc_put_double(avncid, yp_id, yy, subnam)
       !!  T-grid points
       do j=1,b%nyt
          yy(j) = 1.0d-3*( b%yt(j) - b%yp(1) )
       enddo
       call nc_put_double(avncid, yt_id, yy(:b%nyt), subnam)

       !!  Convert mid-layer depth into km and store in 'z'
       tmp(1) = 0.5d-3*b%h(1)
       do k=2,b%nl
          tmp(k) = tmp(k-1) + 0.5d-3*( b%h(k-1) + b%h(k) )
       enddo
       call nc_put_double(avncid, l_id, tmp, subnam)

       call nc_put_double(avncid, st_id, avg%stav, subnam)
       call nc_put_double(avncid, wekt_id, avg%wtav, subnam)
       call nc_put_double(avncid, fm_id, avg%fmav, subnam)
       call nc_put_double(avncid, tx_id, avg%txav, subnam)
       call nc_put_double(avncid, ty_id, avg%tyav, subnam)
       call nc_put_double(avncid, p_id, avg%pav, subnam)
       call nc_put_double(avncid, q_id, avg%qav, subnam)
       call nc_put_double(avncid, ut_id, uptp, subnam)
       call nc_put_double(avncid, vt_id, vptp, subnam)

       call nc_close(avncid, subnam)

    endif

  end subroutine tavout

end module tavsubs
