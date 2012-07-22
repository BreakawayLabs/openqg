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
    integer :: timotdim, xotdim, yotdim, tmp_id
    integer :: xop_id, yop_id, xot_id, yot_id, lo_id, lom_id

    double precision :: xo(go%nxp),yo(go%nyp),tmp(go%nl)
    integer :: i, iwk, mwk

    ! Definition section: define dimensions and variables
    ! Define four dimensions: x, y, z, time
    ocpid = nc_create(outdir, subsamp%filename_p, subnam)
    print *,' ocpo.nc file created'
    octid = nc_create(outdir, subsamp%filename_t, subnam)
    print *,' ocsst.nc file created'

    ! Dimension definitions for p-grid output file
    timopdim = nc_def_dim(ocpid, 'time', numoutsteps, subnam)

    ! x and y dimensions for the (subsampled) p-grid
    mwk = mod(go%nxp,subsamp%nsk)
    iwk = min(mwk,1) + (go%nxp-mwk)/subsamp%nsk
    xopdim = nc_def_dim(ocpid, 'xp', iwk, subnam)
    mwk = mod(go%nyp,subsamp%nsk)
    iwk = min(mwk,1) + (go%nyp-mwk)/subsamp%nsk
    yopdim = nc_def_dim(ocpid, 'yp', iwk, subnam)

    ! And here the z dimension
    lodim = nc_def_dim(ocpid, 'z', go%nl, subnam)
    ! Also define a zi dimension for the interfaces
    lomdim = nc_def_dim(ocpid, 'zi', go%nl-1, subnam)

    ! Dimension definitions for T-grid output file
    timotdim = nc_def_dim(octid, 'time', numoutsteps, subnam)

    ! x and y dimensions for the (subsampled) T-grid
    mwk = mod(go%nxt,subsamp%nsk)
    iwk = min(mwk,1) + (go%nxt-mwk)/subsamp%nsk
    xotdim = nc_def_dim(octid, 'xt', iwk, subnam)
    mwk = mod(go%nyt,subsamp%nsk)
    iwk = min(mwk,1) + (go%nyt-mwk)/subsamp%nsk
    yotdim = nc_def_dim(octid, 'yt', iwk, subnam)

    ! Grid variable definitions for p-grid files
    xop_id = nc_def_float(ocpid, 'xp', xopdim, 'km', subnam, 'Ocean X axis (p-grid)')
    yop_id = nc_def_float(ocpid, 'yp', yopdim, 'km', subnam, 'Ocean Y axis (p-grid)')
    lo_id = nc_def_float(ocpid, 'z', lodim, 'km', subnam, 'Ocean mid-layer depth axis')
    lom_id = nc_def_float(ocpid, 'zi', lomdim, 'km', subnam, 'Ocean interface depth axis')
      
    ! Grid variable definitions for T-grid files
    xot_id = nc_def_float(octid, 'xt', xotdim, 'km', subnam, 'Ocean X axis (T-grid)')
    yot_id = nc_def_float(octid, 'yt', yotdim, 'km', subnam, 'Ocean Y axis (T-grid)')

    tmp_id = nc_def_float(ocpid, 'time', timopdim, 'years', subnam, 'Time axis')
    tmp_id = nc_def_float(octid, 'time', timotdim, 'years', subnam,'Time axis')

    if ( subsamp%outfl(1) ) then
       tmp_id = nc_def_float(octid, 'sst', (/xotdim,yotdim,timotdim/),'K', subnam, 'Ocean surface temperature')
    endif
    if ( subsamp%outfl(2) ) then
       tmp_id = nc_def_float(ocpid, 'p', (/xopdim,yopdim,lodim,timopdim/), 'm^2/s^2', subnam, 'Ocean dynamic pressure')
    endif
    if ( subsamp%outfl(3) ) then
       tmp_id = nc_def_float(ocpid, 'q', (/xopdim,yopdim,lodim,timopdim/), 's^-1', subnam, 'Ocean potential vorticity')
    endif
    if ( subsamp%outfl(4) ) then
       tmp_id = nc_def_float(octid, 'wekt', (/xotdim,yotdim,timotdim/), 'm/s', subnam, 'Ocean Ekman velocity (T-grid)')
    endif
    if ( subsamp%outfl(5) ) then
       tmp_id = nc_def_float(ocpid, 'h', (/xopdim,yopdim,lomdim,timopdim/), 'm', subnam, 'Ocean interface displacement')
    endif
    if ( subsamp%outfl(6) ) then
       tmp_id = nc_def_float(ocpid, 'taux', (/xopdim,yopdim,timopdim/), 'm^2/s^2', subnam,'Zonal dynamic wind stress')
       tmp_id = nc_def_float(ocpid, 'tauy', (/xopdim,yopdim,timopdim/), 'm^2/s^2', subnam,'Meridional dynamic wind stress')
    endif

    ! The oceanic mixed layer thickness is fixed in this
    ! version of the model, so outfloc(7) is irrelevant

    ! Leave definition mode: entering data mode
    call nc_enddef(ocpid, subnam)
    call nc_enddef(octid, subnam)

    ! Calculate x gridpoints and store in 'x' arrays
    ! p-grid points
    mwk = mod(go%nxp,subsamp%nsk)
    iwk = min(mwk,1) + (go%nxp-mwk)/subsamp%nsk
    do i=1,iwk
       xo(i) = m_to_km(go%xp(1+(i-1)*subsamp%nsk) - go%xp(1))
    enddo
    call nc_put_double(ocpid, xop_id, xo(:iwk), subnam)
    ! T-grid points
    mwk = mod(go%nxt,subsamp%nsk)
    iwk = min(mwk,1) + (go%nxt-mwk)/subsamp%nsk
    do i=1,iwk
       xo(i) = m_to_km(go%xt(1+(i-1)*subsamp%nsk) - go%xp(1))
    enddo
    call nc_put_double(octid, xot_id, xo(:iwk), subnam)

    ! Calculate y gridpoints and store in 'y' arrays
    ! p-grid points
    mwk = mod(go%nyp,subsamp%nsk)
    iwk = min(mwk,1) + (go%nyp-mwk)/subsamp%nsk
    do i=1,iwk
       yo(i) = m_to_km(go%yp(1+(i-1)*subsamp%nsk) - go%yp(1))
    enddo
    call nc_put_double(ocpid, yop_id, yo(:iwk), subnam)
    ! T-grid points
    mwk = mod(go%nyt,subsamp%nsk)
    iwk = min(mwk,1) + (go%nyt-mwk)/subsamp%nsk
    do i=1,iwk
       yo(i) = m_to_km(go%yt(1+(i-1)*subsamp%nsk) - go%yp(1))
    enddo
    call nc_put_double(octid, yot_id, yo(:iwk), subnam)

    ! Convert mid-layer depths into km and store in 'z'
    tmp(1) = 0.5d0*m_to_km(go%h(1))
    do i=2,go%nl
       tmp(i) = tmp(i-1) + 0.5d0*m_to_km(go%h(i-1) + go%h(i))
    enddo
    call nc_put_double(ocpid, lo_id, tmp, subnam)

    ! Convert interface depths into km and store in 'zi'
    tmp(1) = m_to_km(go%h(1))
    do i=2,go%nl-1
       tmp(i) = tmp(i-1) + m_to_km(go%h(i))
    enddo
    call nc_put_double(ocpid, lom_id, tmp(:go%nl-1), subnam)

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
    integer :: xap_id, yap_id, xat_id, yat_id, la_id, lam_id

    ! Other variables used locally
    double precision xa(ga%nxp),ya(ga%nyp),tmp(ga%nl)
    integer i, iwk, mwk
    integer :: tmp_id

    atpid = nc_create(outdir, subsamp%filename_p, subnam)
    print *,' atpa.nc file created'
    attid = nc_create(outdir, subsamp%filename_t, subnam)
    print *,' atast.nc file created'

    ! Define four dimensions: x, y, z, time
    ! Do pressure fields first
    timapdim = nc_def_dim(atpid, 'time', numoutsteps, subnam)

    ! Note that x and y are designed for subsampling
    ! Here are the x and y dimensions for the p-grid
    mwk = mod(ga%nxp,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nxp-mwk)/subsamp%nsk
    xapdim = nc_def_dim(atpid, 'xp', iwk, subnam)
    mwk = mod(ga%nyp,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nyp-mwk)/subsamp%nsk
    yapdim = nc_def_dim(atpid, 'yp', iwk, subnam)

    ! And here the z dimension
    ladim = nc_def_dim(atpid, 'z', ga%nl, subnam)
    ! Also define a zi dimension for the interfaces
    lamdim = nc_def_dim(atpid, 'zi', ga%nl-1, subnam)

    ! Now do it for the T-grid
    timatdim = nc_def_dim(attid, 'time', numoutsteps, subnam)

    ! Here are the x and y dimensions for the T-grid
    mwk = mod(ga%nxt,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nxt-mwk)/subsamp%nsk
    xatdim = nc_def_dim(attid, 'xt', iwk, subnam)
    mwk = mod(ga%nyt,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nyt-mwk)/subsamp%nsk
    yatdim = nc_def_dim(attid, 'yt', iwk, subnam)

    xap_id = nc_def_float(atpid, 'xp', xapdim, 'km', subnam, 'Atmosphere X axis (p-grid)')
    xat_id = nc_def_float(attid, 'xt', xatdim, 'km', subnam, 'Atmosphere X axis (T-grid)')
    yap_id = nc_def_float(atpid, 'yp', yapdim, 'km', subnam, 'Atmosphere Y axis (p-grid)')
    yat_id = nc_def_float(attid, 'yt', yatdim, 'km', subnam, 'Atmosphere Y axis (T-grid)')

    la_id = nc_def_float(atpid, 'z', ladim, 'km', subnam, 'Atmosphere mid-layer height axis')
    lam_id = nc_def_float(atpid, 'zi', lamdim, 'km', subnam, 'Atmosphere interface height axis')

    tmp_id = nc_def_float(atpid, 'time', timapdim, 'years', subnam, 'Time axis')
    tmp_id = nc_def_float(attid, 'time', timatdim, 'years', subnam, 'Time axis')

    if ( subsamp%outfl(1) ) then
       tmp_id = nc_def_float(attid, 'ast', (/xatdim,yatdim, timatdim/), 'K', subnam, 'Atmosphere surface temperature')
    endif
    if ( subsamp%outfl(2) ) then
       tmp_id = nc_def_float(atpid, 'p', (/xapdim,yapdim,ladim,timapdim/), 'm^2/s^2', subnam, 'Atmosphere dynamic pressure')
    endif
    if ( subsamp%outfl(3) ) then
       tmp_id = nc_def_float(atpid, 'q', (/xapdim,yapdim,ladim,timapdim/), 's^-1', subnam, 'Atmosphere potential vorticity')
    endif
    if ( subsamp%outfl(4) ) then
       tmp_id = nc_def_float(attid, 'wekt', (/xatdim,yatdim,timatdim/), 'm/s', subnam, 'Atmosphere Ekman velocity (T-grid)')
    endif
    if ( subsamp%outfl(5) ) then
       tmp_id = nc_def_float(atpid, 'h', (/xapdim,yapdim,lamdim,timapdim/), 'm', subnam, 'Atmosphere interface displacement')
    endif
    if ( subsamp%outfl(6) ) then
       tmp_id = nc_def_float(atpid, 'taux', (/xapdim,yapdim,timapdim/), 'm^2/s^2', subnam, 'Zonal wind stress')
       tmp_id = nc_def_float(atpid, 'tauy', (/xapdim,yapdim,timapdim/), 'm^2/s^2', subnam, 'Meridional wind stress')
    endif
    if ( subsamp%outfl(7) ) then
       tmp_id = nc_def_float(attid, 'hmixa', (/xatdim,yatdim,timatdim/), 'm', subnam, 'Mixed layer interface height')
    endif

    ! Leave definition mode: entering data mode
    call nc_enddef(atpid, subnam)
    call nc_enddef(attid, subnam)

    ! Calculate x gridpoints and store in 'x' arrays
    ! p-grid points
    mwk = mod(ga%nxp,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nxp-mwk)/subsamp%nsk
    do i=1,iwk
       xa(i) = m_to_km(ga%xp(1+(i-1)*subsamp%nsk))
    enddo
    call nc_put_double(atpid, xap_id, xa(:iwk), subnam)
    ! T-grid points
    mwk = mod(ga%nxt,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nxt-mwk)/subsamp%nsk
    do i=1,iwk
       xa(i) = m_to_km(ga%xt(1+(i-1)*subsamp%nsk))
    enddo
    call nc_put_double(attid, xat_id, xa(:iwk), subnam)

    ! Calculate y gridpoints and store in 'y' arrays
    ! p-grid points
    mwk = mod(ga%nyp,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nyp-mwk)/subsamp%nsk
    do i=1,iwk
       ya(i) = m_to_km(ga%yp(1+(i-1)*subsamp%nsk))
    enddo
    call nc_put_double(atpid, yap_id, ya(:iwk), subnam)
    ! T-grid points
    mwk = mod(ga%nyt,subsamp%nsk)
    iwk = min(mwk,1) + (ga%nyt-mwk)/subsamp%nsk
    do i=1,iwk
       ya(i) = m_to_km(ga%yt(1+(i-1)*subsamp%nsk))
    enddo
    call nc_put_double(attid, yat_id, ya(:iwk), subnam)

    ! Convert mid-layer heights into km and store in 'z'
    tmp(1) = 0.5d0*m_to_km(ga%h(1))
    do i=2,ga%nl
       tmp(i) = tmp(i-1) + 0.5d0*m_to_km(ga%h(i-1) + ga%h(i))
    enddo
    call nc_put_double(atpid, la_id, tmp, subnam)

    ! Convert interface heights into km and store in 'zi'
    tmp(1) = m_to_km(ga%h(1))
    do i=2,ga%nl-1
       tmp(i) = tmp(i-1) + m_to_km(ga%h(i))
    enddo
    call nc_put_double(atpid, lam_id, tmp(:ga%nl-1), subnam)

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

    integer :: start
    double precision :: eta(go%nxp,go%nyp,go%nl-1)
    integer :: k,ipwk,jpwk,itwk,jtwk,mwk
    integer :: ocpid, octid

    ocpid = nc_open_w(outdir, subsamp%filename_p, subnam)
    octid = nc_open_w(outdir, subsamp%filename_t, subnam)

    ! Store current time as part of 'time' vector
    start = nt/subsamp%nout+1
    ! Pressure file:
    call nc_put_double(ocpid, 'time', start, tyrs, subnam)
    ! Temperature file:
    call nc_put_double(octid, 'time', start, tyrs, subnam)

    ! Compute subsampling array indices
    mwk = mod(go%nxp,subsamp%nsk)
    ipwk = min(mwk,1) + (go%nxp-mwk)/subsamp%nsk
    mwk = mod(go%nxt,subsamp%nsk)
    itwk = min(mwk,1) + (go%nxt-mwk)/subsamp%nsk
    mwk = mod(go%nyp,subsamp%nsk)
    jpwk = min(mwk,1) + (go%nyp-mwk)/subsamp%nsk
    mwk = mod(go%nyt,subsamp%nsk)
    jtwk = min(mwk,1) + (go%nyt-mwk)/subsamp%nsk

    if (subsamp%outfl(1) .and. sst%active) then
       call write_subsample(octid, 'sst', start, sst%data, &
            itwk, jtwk, subsamp%nsk, subnam)
    endif
    if (subsamp%outfl(2)) then
       call write_subsample_3d(ocpid, 'p', start, qgo%p, &
            ipwk, jpwk, subsamp%nsk, go%nl, subnam)
    endif
    if (subsamp%outfl(3)) then
       call write_subsample_3d(ocpid, 'q', start, qgo%q, &
            ipwk, jpwk, subsamp%nsk, go%nl, subnam)
    endif
    if (subsamp%outfl(4)) then
       call write_subsample(octid, 'wekt', start, eko%wekt, &
            itwk, jtwk, subsamp%nsk, subnam)
    endif
    if (subsamp%outfl(5)) then
       do k=1,go%nl-1
          eta(:,:,k) = qgo%b%dz_sign*(qgo%p(:,:,k) - qgo%p(:,:,k+1))/qgo%gp(k)
       enddo
       call write_subsample_3d(ocpid, 'h', start, eta, &
            ipwk, jpwk, subsamp%nsk, go%nl-1, subnam)
    endif
    if (subsamp%outfl(6)) then
       call write_subsample(ocpid, 'taux', start, eko%taux, &
            ipwk, jpwk, subsamp%nsk, subnam)
       call write_subsample(ocpid, 'tauy', start, eko%tauy, &
            ipwk, jpwk, subsamp%nsk, subnam)
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

    integer :: start
    double precision :: eta(ga%nxp,ga%nyp,ga%nl-1)
    integer :: k,ipwk,jpwk,itwk,jtwk,mwk
    integer :: atpid, attid

    atpid = nc_open_w(outdir, subsamp%filename_p, subnam)
    attid = nc_open_w(outdir, subsamp%filename_t, subnam)

    ! Store current time as part of 'time' vector
    start = nt/subsamp%nout+1
    ! Pressure file:
    call nc_put_double(atpid, 'time', start, tyrs, subnam)
    ! Temperature file:
    call nc_put_double(attid, 'time', start, tyrs, subnam)

    ! Compute subsampling array indices
    mwk = mod(ga%nxp,subsamp%nsk)
    ipwk = min(mwk,1) + (ga%nxp-mwk)/subsamp%nsk
    mwk = mod(ga%nxt,subsamp%nsk)
    itwk = min(mwk,1) + (ga%nxt-mwk)/subsamp%nsk
    mwk = mod(ga%nyp,subsamp%nsk)
    jpwk = min(mwk,1) + (ga%nyp-mwk)/subsamp%nsk
    mwk = mod(ga%nyt,subsamp%nsk)
    jtwk = min(mwk,1) + (ga%nyt-mwk)/subsamp%nsk

    if (subsamp%outfl(1) .and. ast%active) then
       call write_subsample(attid, 'ast', start, ast%data, &
            itwk, jtwk, subsamp%nsk, subnam)
    endif
    if (subsamp%outfl(2)) then
       call write_subsample_3d(atpid, 'p', start, qga%p, &
            ipwk, jpwk, subsamp%nsk, ga%nl, subnam)
    endif
    if (subsamp%outfl(3)) then
       call write_subsample_3d(atpid, 'q', start, qga%q, &
            ipwk, jpwk, subsamp%nsk, ga%nl, subnam)
    endif
    if (subsamp%outfl(4)) then
       call write_subsample(attid, 'wekt', start, eka%wekt, &
            itwk, jtwk, subsamp%nsk, subnam)
    endif
    if (subsamp%outfl(5)) then
       do k=1,ga%nl-1
          eta(:,:,k) = qga%b%dz_sign*(qga%p(:,:,k) - qga%p(:,:,k+1))/qga%gp(k)
       enddo
       call write_subsample_3d(atpid, 'h', start, eta, &
            ipwk, jpwk, subsamp%nsk, ga%nl-1, subnam)
    endif
    if (subsamp%outfl(6)) then
       call write_subsample(atpid, 'taux', start, eka%taux, &
            ipwk, jpwk, subsamp%nsk, subnam)
       call write_subsample(atpid, 'tauy', start, eka%tauy, &
            ipwk, jpwk, subsamp%nsk, subnam)
    endif
    if (subsamp%outfl(7) .and. hmixa%active) then
       call write_subsample(attid, 'hmixa', start, hmixa%data, &
            itwk, jtwk, subsamp%nsk, subnam)
    endif

    call nc_close(atpid)
    call nc_close(attid)

  end subroutine atnc_out

  subroutine write_subsample(ncid, varname, start, data, &
       itwk, jtwk, nsk, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: start
    double precision, intent(in) :: data(:,:)
    integer, intent(in) :: itwk, jtwk, nsk
    character (len=*), intent(in) :: subnam

    integer :: i, j
    double precision :: wrk(itwk,jtwk)

    do j=1,jtwk
       do i=1,itwk
          wrk(i,j) = data(1+(i-1)*nsk,1+(j-1)*nsk)
       enddo
    enddo
    call nc_put_double(ncid, varname, start, wrk, subnam)
    
  end subroutine write_subsample

  subroutine write_subsample_3d(ncid, varname, start, data, &
       itwk, jtwk, nsk, nl, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: start
    double precision, intent(in) :: data(:,:,:)
    integer, intent(in) :: itwk, jtwk, nsk, nl
    character (len=*), intent(in) :: subnam

    integer :: i, j, k
    double precision :: wrk(itwk,jtwk,nl)

    do k=1,nl
       do j=1,jtwk
          do i=1,itwk
             wrk(i,j,k) = data(1+(i-1)*nsk,1+(j-1)*nsk,k)
          enddo
       enddo
    enddo
    call nc_put_double(ncid, varname, start, wrk, subnam)
    
  end subroutine write_subsample_3d

end module subsampling
