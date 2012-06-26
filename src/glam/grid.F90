module grid

  use ncutils, only: nc_open, nc_close, nc_get_double, nc_get_dim, nc_get_int
  use box, only: box_type

  implicit none

  private

  type grid_type

     ! User settable grid parameters:
     integer :: nxaooc, nyaooc
     integer :: ndxr

     ! Derived grid parameters:
     integer :: nxtaor,nytaor,nxpaor,nypaor, &
          nx1,ny1
     integer :: nxtoar,nytoar,nxpoar,nypoar

     !logical :: cyclic

     !     nxtaor, nytaor are the numbers of atmos. T points at ocean resolution
     !     nxpaor, nypaor are the numbers of atmos. p points at ocean resolution
     !     nxtoar, nytoar are the numbers of ocean T points at atmos. resolution
     !     nxpoar, nypoar are the numbers of ocean p points at atmos. resolution
     !     (required for the new windstress formulation introduced at v1.4.0)
     !     nx1, ny1 are starting indices for the ocean in the atmos. grid.
     !     We choose them to centre the ocean in the atmospheric domain (if possible).
     !     atnorm, ocnorm are normalisation factors for computing mean
     !     values, and are both equal to 1/(number of T gridcells).

     !double precision :: dyo ! dxo, dyo are the oceanic gridlengths (m)
      ! dxa, dya are the atmospheric gridlengths (m)

     !     hto = total ocean depth = hoc(1) + ... + hoc(nlo) (m)
     !     hta = total atmos. depth = hat(1) + ... + hat(nla) (m)

     !     xpo = vector of p-gridpoint positions in ocean
     !           (including offset w.r.t. atmosphere) (m)
     !     ypo = vector of p-gridpoint positions in ocean
     !           (including offset w.r.t. atmosphere) (m)
     !     xto = vector of T-gridpoint positions in ocean
     !           (including offset w.r.t. atmosphere) (m)
     !     yto = vector of T-gridpoint positions in ocean
     !           (including offset w.r.t. atmosphere) (m)
     !     xpa = vector of p-gridpoint positions in atmosphere (m)
     !     ypa = vector of p-gridpoint positions in atmosphere (m)
     !     xta = vector of T-gridpoint positions in atmosphere (m)
     !     yta = vector of T-gridpoint positions in atmosphere (m)


     !     hoc = unperturbed oceanic layer thicknesses (m)
     !     hat = unperturbed atmospheric layer thicknesses (m)

  end type grid_type

  public grid_type
  public load_grid
  public print_grid
  public new_load_grid

contains

  type(grid_type) function new_load_grid(go, ga)
    type(box_type), intent(in) :: go
    type(box_type), intent(in) :: ga

    new_load_grid%ndxr = nint(ga%dx/go%dx)

    new_load_grid%nxtaor = ga%nxt*new_load_grid%ndxr
    new_load_grid%nytaor = ga%nyt*new_load_grid%ndxr
    new_load_grid%nxpaor = new_load_grid%nxtaor + 1
    new_load_grid%nypaor = new_load_grid%nytaor + 1

    new_load_grid%nxtoar = go%nxt/new_load_grid%ndxr
    new_load_grid%nytoar = go%nyt/new_load_grid%ndxr
    new_load_grid%nxpoar = new_load_grid%nxtoar + 1
    new_load_grid%nypoar = new_load_grid%nytoar + 1

    new_load_grid%nxaooc = go%nxt/new_load_grid%ndxr
    new_load_grid%nyaooc = go%nyt/new_load_grid%ndxr

    new_load_grid%nx1 = 1
    new_load_grid%ny1 = 1 + (ga%nyt - new_load_grid%nyaooc)/2

  end function new_load_grid

  subroutine new_print_grid(g)
    type(grid_type), intent(in) :: g

    print *, "Coupling Grid Summary"
    print *, "Atmos/Ocean ratio (ndxr): ", g%ndxr
    print *, "# Atmos T. cells at ocean res (nxtaor,nytaor): ", g%nxtaor, g%nytaor
    print *, "# Atmos cells over ocean (nxaooc, nyaooc): ", g%nxaooc, g%nyaooc
    print *, "# Index on atmos T-grid of first ocean cell (nx1, ny1): ", g%nx1, g%ny1

  end subroutine new_print_grid

  type(grid_type) function load_grid(filename, go, ga)

    character (len=*), intent(in) :: filename
    type(box_type), intent(in) :: go
    type(box_type), intent(in) :: ga

    integer :: grid_id

    ! Namelist specified values
    grid_id = nc_open(filename, 'Load grid')

    load_grid%nxaooc = nc_get_int(grid_id, 'nxaooc', 'load_grid')
    load_grid%nyaooc = nc_get_int(grid_id, 'nyaooc', 'load_grid')
    load_grid%ndxr = nc_get_int(grid_id, 'ndxr', 'load_grid')


    ! Derived quantities
    load_grid%nxtaor = ga%nxt*load_grid%ndxr
    load_grid%nytaor = ga%nyt*load_grid%ndxr
    load_grid%nxpaor = load_grid%nxtaor+1
    load_grid%nypaor = load_grid%nytaor+1
    load_grid%nx1 = 1
    load_grid%ny1 = 1 + (ga%nyt-load_grid%nyaooc)/2

    ! fnot is the Coriolis parameter at the
    ! central latitude of the domain (rad s^-1)
    ! beta is the y-derivative of the Coriolis parameter (rad s^-1 m^-1)


    call nc_close(grid_id, 'load_grid')

    if (go%cyclic) then
       ! Check cyclic ocean case is properly set up
       if ( ga%nxt.ne.load_grid%nxaooc ) then
          print *,' '
          print *,' nxta, nxaooc = ',ga%nxt,load_grid%nxaooc
          print *,' For cyclic ocean nxta should equal nxaooc'
          print *,' Program terminates'
          stop 1
       endif
    else
       ! Check finite box ocean case is properly set up
       if ( ga%nxt.lt.load_grid%nxaooc ) then
          print *,' '
          print *,' nxta, nxaooc = ',ga%nxt,load_grid%nxaooc
          print *,' Inadequate nxta needs to be at least nxaooc'
          print *,' Program terminates'
          stop 1
       endif
    endif
    
    if ( ga%nyt.lt.load_grid%nyaooc ) then
       print *,' '
       print *,' nyta, nyaooc = ',ga%nyt,load_grid%nyaooc
       print *,' Inadequate nyta needs to be at least nyaooc'
       print *,' Program terminates'
       stop 1
    endif

    if ( ga%nl.lt.2 .or. go%nl.lt.2 ) then
       print *,' '
       print *,' nla, nlo = ',ga%nl,go%nl
       print *,' Inadequate nla or nlo, needs to be at least 2'
       print *,' Program terminates'
       stop 1
    endif

  end function load_grid

  subroutine print_grid(g, go, ga)

    type(grid_type), intent(in) :: g
    type(box_type), intent(in) :: go
    type(box_type), intent(in) :: ga

    integer :: k

    print *,' '
    print *,' Grid parameters:'
    print *,' ----------------'
    write(*,201) '  Atmos/ocean grid ratio ndxr = ',g%ndxr
    write(*,201) '  Atmos. gridcells over ocean = ',g%nxaooc,g%nyaooc
    write(*,201) '  Ocn start indices  nx1, ny1 = ',g%nx1,g%ny1
    write(*,214) '  Coriolis par. f0 (rad s^-1) = ',go%fnot
    write(*,214) '  Beta =df/dy (rad s^-1 m^-1) = ',go%beta

    print *,' '
    print *,' Oceanic grid:'
    print *,' -------------'
    write(*,201) '  No. of ocean QG layers  nlo = ',go%nl
    write(*,201) '  No. of gridcells nxto, nyto = ',go%nxt,go%nyt
    write(*,204) '  Gridlength dxo         (km) = ',1.0d-3*go%dx
    write(*,203) '  Domain sizes xlo, ylo  (km) = ', &
         1.0d-3*go%xl,1.0d-3*go%yl
    write(*,205) '  Rossby number   Beta*ylo/f0 = ',go%beta*go%yl/abs(go%fnot)
    write(*,214) '  f range S -> N   (rad s^-1) = ', &
         go%fnot+go%beta*go%yprel(1),go%fnot+go%beta*go%yprel(go%nyp)
    write(*,214) '  Midlatitude Coriolis param  = ', &
         go%fnot+go%beta*0.5d0*( go%yprel(1) + go%yprel(go%nyp) )
    write(*,203) '  Layer thicknesses hoc   (m) = ',(go%h(k),k=1,go%nl)
    write(*,203) '  Total thickness   hto   (m) = ',sum(go%h)

    print *,' '
    print *,' Atmospheric grid:'
    print *,' -----------------'
    write(*,201) '  No. of atmos. QG layers nla = ',ga%nl
    write(*,201) '  No. of gridcells nxta, nyta = ',ga%nxt,ga%nyt
    write(*,204) '  Gridlength dxa         (km) = ',1.0d-3*ga%dx
    write(*,203) '  Domain sizes xla, yla  (km) = ', &
         1.0d-3*ga%xl,1.0d-3*ga%yl
    write(*,205) '  Rossby number   Beta*yla/f0 = ',go%beta*ga%yl/abs(ga%fnot)
    write(*,214) '  f range S -> N   (rad s^-1) = ', &
         ga%fnot+ga%beta*ga%yprel(1),ga%fnot+go%beta*ga%yprel(ga%nyp)
    write(*,214) '  Midlatitude Coriolis param  = ', &
         ga%fnot+ga%beta*0.5d0*( ga%yprel(1) + ga%yprel(ga%nyp) )
    write(*,203) '  Layer thicknesses hat   (m) = ',(ga%h(k),k=1,ga%nl)
    write(*,203) '  Total thickness   hta   (m) = ',sum(ga%h)

  201 format(a,9i13)
  203 format(a,9f13.3)
  204 format(a,9f13.4)
  205 format(a,9f13.5)
  214 format(a,1p,9d13.4)

  end subroutine print_grid

end module grid
