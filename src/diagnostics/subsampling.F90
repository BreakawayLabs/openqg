module subsampling

  use units, only: m_to_km
  use box, only: box_type
  use clock, only: clock_type, days_to_steps
  use mixed, only: mixed_type
  use qg, only: qg_type
  use ekman, only: ekman_type
  use ncutils, only: nc_def_float, nc_def_dim, nc_create, nc_open_w, nc_close, nc_enddef, nc_put_double

  implicit none

  private

  type subsampling_type

     character (len=64):: filename_p
     character (len=64):: filename_t

     integer :: nsk
     logical :: outfl(7)
     integer :: nout
     
     logical :: active = .false.

  end type subsampling_type

  public init_ocn_subsamp
  public init_atm_subsamp
  public ocnc_out
  public atnc_out
  public subsampling_type
  public subsamp_step

contains

  type(subsampling_type) function init_ocn_subsamp(filename_p, filename_t, nsk, outfl, clk, &
       subsamp_days, outdir, b)

    character (len=*), intent(in) :: filename_p
    character (len=*), intent(in) :: filename_t
    integer, intent(in) :: nsk
    logical, intent(in) :: outfl(7)
    type(clock_type), intent(in) :: clk
    double precision, intent(in) :: subsamp_days
    character (len=*), intent(in) :: outdir
    type(box_type), intent(in) :: b

    integer :: nout, noutstep

    nout = days_to_steps(subsamp_days, clk)
    noutstep = clk%ntsrun/nout + 1

    init_ocn_subsamp%filename_p = filename_p
    init_ocn_subsamp%filename_t = filename_t
    init_ocn_subsamp%nsk = nsk
    init_ocn_subsamp%outfl = outfl
    init_ocn_subsamp%nout = nout
    call ocnc_init(outdir, noutstep, b, init_ocn_subsamp)
    
    init_ocn_subsamp%active = .true.

  end function init_ocn_subsamp

  type(subsampling_type) function init_atm_subsamp(filename_p, filename_t, nsk, outfl, clk, &
       subsamp_days, outdir, b)

    character (len=*), intent(in) :: filename_p
    character (len=*), intent(in) :: filename_t
    integer, intent(in) :: nsk
    logical, intent(in) :: outfl(7)
    type(clock_type), intent(in) :: clk
    double precision, intent(in) :: subsamp_days
    character (len=*), intent(in) :: outdir
    type(box_type), intent(in) :: b

    integer :: nout, noutstep

    nout = days_to_steps(subsamp_days, clk)
    noutstep = clk%ntsrun/nout + 1

    init_atm_subsamp%filename_p = filename_p
    init_atm_subsamp%filename_t = filename_t
    init_atm_subsamp%nsk = nsk
    init_atm_subsamp%outfl = outfl
    init_atm_subsamp%nout = nout
    call atnc_init(outdir, noutstep, b, init_atm_subsamp)
    
    init_atm_subsamp%active = .true.

  end function init_atm_subsamp

  logical function subsamp_step(ntdone, subsamp, force)
    
    integer, intent(in) :: ntdone
    type(subsampling_type), intent(in) :: subsamp
    logical, intent(in) :: force

    subsamp_step = subsamp%active .and. (force .or. mod(ntdone,subsamp%nout) == 0)

  end function subsamp_step

  subroutine ocnc_init (outdir, numoutsteps, go, subsamp)

    character (len=*), intent(in) :: outdir
    integer, intent(in) :: numoutsteps
    type(box_type), intent(in) :: go
    type(subsampling_type), intent(in) :: subsamp

    integer :: ocpid, octid

    character :: subnam*(*)
    parameter ( subnam = 'ocnc_init' )

    integer :: timopdim, xopdim, yopdim, lodim, lomdim
    integer :: timotdim, xotdim, yotdim

    double precision :: tmp(go%nl)
    integer :: i

    ! Definition section: define dimensions and variables
    ! Define four dimensions: x, y, z, time
    ocpid = nc_create(outdir, subsamp%filename_p, subnam)
    print *,' ocpo.nc file created'
    octid = nc_create(outdir, subsamp%filename_t, subnam)
    print *,' ocsst.nc file created'

    ! Dimension definitions for p-grid output file
    timopdim = nc_def_dim(ocpid, 'time', numoutsteps, subnam)

    ! x and y dimensions for the (subsampled) p-grid
    xopdim = nc_def_dim(ocpid, 'xp', size(go%xp(::subsamp%nsk)), subnam)
    yopdim = nc_def_dim(ocpid, 'yp', size(go%yp(::subsamp%nsk)), subnam)

    ! And here the z dimension
    lodim = nc_def_dim(ocpid, 'z', go%nl, subnam)
    ! Also define a zi dimension for the interfaces
    lomdim = nc_def_dim(ocpid, 'zi', go%nl-1, subnam)

    ! Dimension definitions for T-grid output file
    timotdim = nc_def_dim(octid, 'time', numoutsteps, subnam)

    ! x and y dimensions for the (subsampled) T-grid
    xotdim = nc_def_dim(octid, 'xt', size(go%xt(::subsamp%nsk)), subnam)
    yotdim = nc_def_dim(octid, 'yt', size(go%yt(::subsamp%nsk)), subnam)

    ! Grid variable definitions for p-grid files
    call nc_def_float(ocpid, 'xp', xopdim, 'km', subnam, 'Ocean X axis (p-grid)')
    call nc_def_float(ocpid, 'yp', yopdim, 'km', subnam, 'Ocean Y axis (p-grid)')
    call nc_def_float(ocpid, 'z', lodim, 'km', subnam, 'Ocean mid-layer depth axis')
    call nc_def_float(ocpid, 'zi', lomdim, 'km', subnam, 'Ocean interface depth axis')
      
    ! Grid variable definitions for T-grid files
    call nc_def_float(octid, 'xt', xotdim, 'km', subnam, 'Ocean X axis (T-grid)')
    call nc_def_float(octid, 'yt', yotdim, 'km', subnam, 'Ocean Y axis (T-grid)')

    call nc_def_float(ocpid, 'time', timopdim, 'years', subnam, 'Time axis')
    call nc_def_float(octid, 'time', timotdim, 'years', subnam,'Time axis')

    if ( subsamp%outfl(1) ) then
       call nc_def_float(octid, 'sst', (/xotdim,yotdim,timotdim/),'K', subnam, 'Ocean surface temperature')
    endif
    if ( subsamp%outfl(2) ) then
       call nc_def_float(ocpid, 'p', (/xopdim,yopdim,lodim,timopdim/), 'm^2/s^2', subnam, 'Ocean dynamic pressure')
    endif
    if ( subsamp%outfl(3) ) then
       call nc_def_float(ocpid, 'q', (/xopdim,yopdim,lodim,timopdim/), 's^-1', subnam, 'Ocean potential vorticity')
    endif
    if ( subsamp%outfl(4) ) then
       call nc_def_float(octid, 'wekt', (/xotdim,yotdim,timotdim/), 'm/s', subnam, 'Ocean Ekman velocity (T-grid)')
    endif
    if ( subsamp%outfl(5) ) then
       call nc_def_float(ocpid, 'h', (/xopdim,yopdim,lomdim,timopdim/), 'm', subnam, 'Ocean interface displacement')
    endif
    if ( subsamp%outfl(6) ) then
       call nc_def_float(ocpid, 'taux', (/xopdim,yopdim,timopdim/), 'm^2/s^2', subnam,'Zonal dynamic wind stress')
       call nc_def_float(ocpid, 'tauy', (/xopdim,yopdim,timopdim/), 'm^2/s^2', subnam,'Meridional dynamic wind stress')
    endif

    ! The oceanic mixed layer thickness is fixed in this
    ! version of the model, so outfloc(7) is irrelevant

    ! Leave definition mode: entering data mode
    call nc_enddef(ocpid, subnam)
    call nc_enddef(octid, subnam)

    ! Calculate x gridpoints and store in 'x' arrays
    call nc_put_double(ocpid, 'xp', m_to_km(go%xp(::subsamp%nsk) - go%xp(1)), subnam)
    call nc_put_double(octid, 'xt', m_to_km(go%xt(::subsamp%nsk) - go%xp(1)), subnam)

    ! Calculate y gridpoints and store in 'y' arrays
    call nc_put_double(ocpid, 'yp', m_to_km(go%yp(::subsamp%nsk) - go%yp(1)), subnam)
    call nc_put_double(octid, 'yt', m_to_km(go%yt(::subsamp%nsk) - go%yp(1)), subnam)

    ! Convert mid-layer depths into km and store in 'z'
    tmp(1) = 0.5d0*go%h(1)
    do i=2,go%nl
       tmp(i) = tmp(i-1) + 0.5d0*(go%h(i-1) + go%h(i))
    enddo
    call nc_put_double(ocpid, 'z', m_to_km(tmp), subnam)

    ! Convert interface depths into km and store in 'zi'
    tmp(1) = go%h(1)
    do i=2,go%nl-1
       tmp(i) = tmp(i-1) + go%h(i)
    enddo
    call nc_put_double(ocpid, 'zi', m_to_km(tmp(:go%nl-1)), subnam)

    print *,' Ocean netCDF files initialised'
    call nc_close(ocpid)
    call nc_close(octid)

  end subroutine ocnc_init

  subroutine atnc_init (outdir, numoutsteps, ga, subsamp)

    character (len=*), intent(in) :: outdir
    integer, intent(in) :: numoutsteps
    type(box_type), intent(in) :: ga
    type(subsampling_type), intent(in) :: subsamp

    integer :: atpid, attid

    character :: subnam*(*)
    parameter ( subnam = 'atnc_init' )

    ! netCDF variables used locally
    integer :: timapdim, xapdim, yapdim, ladim, lamdim
    integer :: timatdim, xatdim, yatdim

    ! Other variables used locally
    double precision tmp(ga%nl)
    integer i

    atpid = nc_create(outdir, subsamp%filename_p, subnam)
    print *,' atpa.nc file created'
    attid = nc_create(outdir, subsamp%filename_t, subnam)
    print *,' atast.nc file created'

    ! Define four dimensions: x, y, z, time
    ! Do pressure fields first
    timapdim = nc_def_dim(atpid, 'time', numoutsteps, subnam)

    ! Note that x and y are designed for subsampling
    ! Here are the x and y dimensions for the p-grid
    xapdim = nc_def_dim(atpid, 'xp', size(ga%xp(::subsamp%nsk)), subnam)
    yapdim = nc_def_dim(atpid, 'yp', size(ga%yp(::subsamp%nsk)), subnam)

    ! And here the z dimension
    ladim = nc_def_dim(atpid, 'z', ga%nl, subnam)
    ! Also define a zi dimension for the interfaces
    lamdim = nc_def_dim(atpid, 'zi', ga%nl-1, subnam)

    ! Now do it for the T-grid
    timatdim = nc_def_dim(attid, 'time', numoutsteps, subnam)

    ! Here are the x and y dimensions for the T-grid
    xatdim = nc_def_dim(attid, 'xt', size(ga%xt(::subsamp%nsk)), subnam)
    yatdim = nc_def_dim(attid, 'yt', size(ga%yt(::subsamp%nsk)), subnam)

    call nc_def_float(atpid, 'xp', xapdim, 'km', subnam, 'Atmosphere X axis (p-grid)')
    call nc_def_float(attid, 'xt', xatdim, 'km', subnam, 'Atmosphere X axis (T-grid)')
    call nc_def_float(atpid, 'yp', yapdim, 'km', subnam, 'Atmosphere Y axis (p-grid)')
    call nc_def_float(attid, 'yt', yatdim, 'km', subnam, 'Atmosphere Y axis (T-grid)')

    call nc_def_float(atpid, 'z', ladim, 'km', subnam, 'Atmosphere mid-layer height axis')
    call nc_def_float(atpid, 'zi', lamdim, 'km', subnam, 'Atmosphere interface height axis')

    call nc_def_float(atpid, 'time', timapdim, 'years', subnam, 'Time axis')
    call nc_def_float(attid, 'time', timatdim, 'years', subnam, 'Time axis')

    if ( subsamp%outfl(1) ) then
       call nc_def_float(attid, 'ast', (/xatdim,yatdim, timatdim/), 'K', subnam, 'Atmosphere surface temperature')
    endif
    if ( subsamp%outfl(2) ) then
       call nc_def_float(atpid, 'p', (/xapdim,yapdim,ladim,timapdim/), 'm^2/s^2', subnam, 'Atmosphere dynamic pressure')
    endif
    if ( subsamp%outfl(3) ) then
       call nc_def_float(atpid, 'q', (/xapdim,yapdim,ladim,timapdim/), 's^-1', subnam, 'Atmosphere potential vorticity')
    endif
    if ( subsamp%outfl(4) ) then
       call nc_def_float(attid, 'wekt', (/xatdim,yatdim,timatdim/), 'm/s', subnam, 'Atmosphere Ekman velocity (T-grid)')
    endif
    if ( subsamp%outfl(5) ) then
       call nc_def_float(atpid, 'h', (/xapdim,yapdim,lamdim,timapdim/), 'm', subnam, 'Atmosphere interface displacement')
    endif
    if ( subsamp%outfl(6) ) then
       call nc_def_float(atpid, 'taux', (/xapdim,yapdim,timapdim/), 'm^2/s^2', subnam, 'Zonal wind stress')
       call nc_def_float(atpid, 'tauy', (/xapdim,yapdim,timapdim/), 'm^2/s^2', subnam, 'Meridional wind stress')
    endif
    if ( subsamp%outfl(7) ) then
       call nc_def_float(attid, 'hmixa', (/xatdim,yatdim,timatdim/), 'm', subnam, 'Mixed layer interface height')
    endif

    ! Leave definition mode: entering data mode
    call nc_enddef(atpid, subnam)
    call nc_enddef(attid, subnam)

    ! Calculate x gridpoints and store in 'x' arrays
    call nc_put_double(atpid, 'xp', m_to_km(ga%xp(::subsamp%nsk)), subnam)
    call nc_put_double(attid, 'xt', m_to_km(ga%xt(::subsamp%nsk)), subnam)

    ! Calculate y gridpoints and store in 'y' arrays
    call nc_put_double(atpid, 'yp', m_to_km(ga%yp(::subsamp%nsk)), subnam)
    call nc_put_double(attid, 'yt', m_to_km(ga%yt(::subsamp%nsk)), subnam)

    ! Convert mid-layer heights into km and store in 'z'
    tmp(1) = 0.5d0*ga%h(1)
    do i=2,ga%nl
       tmp(i) = tmp(i-1) + 0.5d0*(ga%h(i-1) + ga%h(i))
    enddo
    call nc_put_double(atpid, 'z', m_to_km(tmp), subnam)

    ! Convert interface heights into km and store in 'zi'
    tmp(1) = ga%h(1)
    do i=2,ga%nl-1
       tmp(i) = tmp(i-1) + ga%h(i)
    enddo
    call nc_put_double(atpid, 'zi', m_to_km(tmp(:ga%nl-1)), subnam)

    print *,' Atmos. netCDF files initialised'
    call nc_close(atpid)
    call nc_close(attid)

  end subroutine atnc_init

  subroutine ocnc_out (outdir, nt, tyrs, sst, qgo, eko, go, subsamp)

    character (len=*), intent(in) :: outdir
    integer, intent(in) :: nt
    double precision, intent(in) :: tyrs
    type(mixed_type), intent(in) :: sst
    type(qg_type), intent(in) :: qgo
    type(ekman_type), intent(in) :: eko
    type(box_type), intent(in) :: go
    type(subsampling_type), intent(in) :: subsamp

    character :: subnam*(*)
    parameter ( subnam = 'ocnc_out' )

    integer :: start, k, ocpid, octid
    double precision :: eta(go%nxp,go%nyp,go%nl-1)

    ocpid = nc_open_w(outdir, subsamp%filename_p, subnam)
    octid = nc_open_w(outdir, subsamp%filename_t, subnam)

    ! Store current time as part of 'time' vector
    start = nt/subsamp%nout+1
    ! Pressure file:
    call nc_put_double(ocpid, 'time', start, tyrs, subnam)
    ! Temperature file:
    call nc_put_double(octid, 'time', start, tyrs, subnam)

    if (subsamp%outfl(1) .and. sst%active) then
       call nc_put_double(octid, 'sst', start, sst%data(::subsamp%nsk,::subsamp%nsk), subnam)
    endif
    if (subsamp%outfl(2)) then
       call nc_put_double(ocpid, 'p', start, qgo%p(::subsamp%nsk,::subsamp%nsk,:), subnam)
    endif
    if (subsamp%outfl(3)) then
       call nc_put_double(ocpid, 'q', start, qgo%q(::subsamp%nsk,::subsamp%nsk,:), subnam)
    endif
    if (subsamp%outfl(4)) then
       call nc_put_double(octid, 'wekt', start, eko%wekt(::subsamp%nsk,::subsamp%nsk), subnam)
    endif
    if (subsamp%outfl(5)) then
       do k=1,go%nl-1
          eta(:,:,k) = qgo%b%dz_sign*(qgo%p(:,:,k) - qgo%p(:,:,k+1))/qgo%gp(k)
       enddo
       call nc_put_double(ocpid, 'h', start, eta(::subsamp%nsk,::subsamp%nsk,:), subnam)
    endif
    if (subsamp%outfl(6)) then
       call nc_put_double(ocpid, 'taux', start, eko%taux(::subsamp%nsk,::subsamp%nsk), subnam)
       call nc_put_double(ocpid, 'tauy', start, eko%tauy(::subsamp%nsk,::subsamp%nsk), subnam)
    endif

    ! The oceanic mixed layer thickness is fixed in this
    ! version of the model, so outfloc(7) is irrelevant

    call nc_close(ocpid)
    call nc_close(octid)

  end subroutine ocnc_out

  subroutine atnc_out (outdir, nt, tyrs, ast, hmixa, qga, eka, ga, subsamp)

    character (len=*), intent(in) :: outdir
    integer, intent(in) :: nt
    double precision, intent(in) :: tyrs
    type(mixed_type), intent(in) :: ast
    type(mixed_type), intent(in) :: hmixa
    type(qg_type), intent(in) :: qga
    type(ekman_type), intent(in) :: eka
    type(box_type), intent(in) :: ga
    type(subsampling_type), intent(in) :: subsamp

    character :: subnam*(*)
    parameter ( subnam = 'atnc_out' )

    integer :: start, k, atpid, attid
    double precision :: eta(ga%nxp,ga%nyp,ga%nl-1)

    atpid = nc_open_w(outdir, subsamp%filename_p, subnam)
    attid = nc_open_w(outdir, subsamp%filename_t, subnam)

    ! Store current time as part of 'time' vector
    start = nt/subsamp%nout+1
    ! Pressure file:
    call nc_put_double(atpid, 'time', start, tyrs, subnam)
    ! Temperature file:
    call nc_put_double(attid, 'time', start, tyrs, subnam)

    if (subsamp%outfl(1) .and. ast%active) then
       call nc_put_double(attid, 'ast', start, ast%data(::subsamp%nsk,::subsamp%nsk), subnam)
    endif
    if (subsamp%outfl(2)) then
       call nc_put_double(atpid, 'p', start, qga%p(::subsamp%nsk,::subsamp%nsk,:), subnam)
    endif
    if (subsamp%outfl(3)) then
       call nc_put_double(atpid, 'q', start, qga%q(::subsamp%nsk,::subsamp%nsk,:), subnam)
    endif
    if (subsamp%outfl(4)) then
       call nc_put_double(attid, 'wekt', start, eka%wekt(::subsamp%nsk,::subsamp%nsk), subnam)
    endif
    if (subsamp%outfl(5)) then
       do k=1,ga%nl-1
          eta(:,:,k) = qga%b%dz_sign*(qga%p(:,:,k) - qga%p(:,:,k+1))/qga%gp(k)
       enddo
       call nc_put_double(atpid, 'h', start, eta(::subsamp%nsk,::subsamp%nsk,:), subnam)
    endif
    if (subsamp%outfl(6)) then
       call nc_put_double(atpid, 'taux', start, eka%taux(::subsamp%nsk,::subsamp%nsk), subnam)
       call nc_put_double(atpid, 'tauy', start, eka%tauy(::subsamp%nsk,::subsamp%nsk), subnam)
    endif
    if (subsamp%outfl(7) .and. hmixa%active) then
       call nc_put_double(attid, 'hmixa', start, hmixa%data(::subsamp%nsk,::subsamp%nsk), subnam)
    endif

    call nc_close(atpid)
    call nc_close(attid)

  end subroutine atnc_out

end module subsampling
