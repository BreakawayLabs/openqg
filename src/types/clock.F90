module clock

  use ncutils, only: nc_open, nc_close, nc_get_double, nc_get_int

  implicit none

  private

  type clock_type

     ! user defined values
     double precision :: dta ! Atmos. timestep in seconds
     integer :: nstr ! Timestep ratio dto/dta
     double precision :: tini
     double precision :: trun

     ! Derived values
     double precision :: dto, tdto, tdta

     !! dto = ocean timestep (s)
     !! tdto = 2.0d0*dto (s)

     integer :: nsteps0, nsteps, ntsrun
     integer :: noutre

     ! 0.013698630D0   !! trun    = Length of run in years   (5 days)
     ! 0.027397260D0   !! trun    = Length of run in years  (10 days)
     ! 0.041095890D0   !! trun    = Length of run in years  (15 days)
     ! 0.054794521D0   !! trun    = Length of run in years  (20 days)
     ! 0.082191781D0   !! trun    = Length of run in years  (30 days)
     ! 0.109589041D0   !! trun    = Length of run in years  (40 days)
     ! 0.164383562D0   !! trun    = Length of run in years  (60 days)
     ! 0.20D0          !! trun    = Length of run in years  (73 days)
     ! 0.246575342D0   !! trun    = Length of run in years  (90 days)
     ! 0.328767123D0   !! trun    = Length of run in years (120 days)
     ! 0.40D0          !! trun    = Length of run in years (146 days)
     ! 0.60D0          !! trun    = Length of run in years (219 days)
     ! 5.00D0         !! trun    = Length of run in years

  end type clock_type

  public clock_type
  public load_clock
  public days_to_steps

contains

  type(clock_type) function load_clock(filename)

    character (len=*), intent(in) :: filename
    double precision :: resday

    integer :: clock_id 
    double precision :: tend

    double precision :: secday,daysyr,secsyr
    parameter ( secday=86400.0d0, daysyr=365.0d0, &
         secsyr=secday*daysyr )

    ! Namelist specified values
    clock_id = nc_open(filename, 'load_clock')
    load_clock%dta = nc_get_double(clock_id, 'dta', 'load_clock')
    load_clock%tini = nc_get_double(clock_id, 'tini', 'load_clock')
    load_clock%trun = nc_get_double(clock_id, 'trun', 'load_clock')
    load_clock%nstr = nc_get_int(clock_id, 'nstr', 'load_clock')
    resday = nc_get_double(clock_id, 'resday', 'load_clock')
    call nc_close(clock_id, 'load_clock')

    ! Derived values
    load_clock%dto = load_clock%nstr*load_clock%dta
    load_clock%tdto = 2.0d0*load_clock%dto
    load_clock%tdta = 2.0d0*load_clock%dta
    tend = load_clock%tini + load_clock%trun
    load_clock%nsteps0 = nint(load_clock%tini*secsyr/load_clock%dta)
    load_clock%nsteps = nint(tend*secsyr/load_clock%dta)
    load_clock%ntsrun = load_clock%nsteps - load_clock%nsteps0

    write(*,201) '  Oc/atm. timestep ratio nstr = ',load_clock%nstr
    write(*,205) '  Timestep dto      (minutes) = ',load_clock%dto/60.0d0
    write(*,205) '  Timestep dta      (minutes) = ',load_clock%dta/60.0d0
    write(*,204) '  Start time tini     (years) = ',load_clock%tini
    write(*,204) '  Run length trun     (years) = ',load_clock%trun
    write(*,204) '  Final time tend     (years) = ',tend
    write(*,201) '  Start no. of (atmos)  steps = ',load_clock%nsteps0
    write(*,201) '  Final no. of (atmos)  steps = ',load_clock%nsteps
    write(*,201) '  Total no. of (atmos)  steps = ',load_clock%ntsrun
    if ( resday.gt.0.0d0 ) then
       write(*,204) '  Restart dump interval (day) = ',resday
       load_clock%noutre = days_to_steps(resday, load_clock)
    else
       load_clock%noutre = load_clock%ntsrun + 1
    endif

201 format(a,9i13)
205 format(a,9f13.5)
204 format(a,9f13.4)

  end function load_clock

  integer pure function days_to_steps(days, clk)

    double precision, intent(in) :: days
    type(clock_type), intent(in) :: clk

    double precision :: secday
    parameter ( secday=86400.0d0 )

    days_to_steps = nint( days*secday/clk%dta )

  end function days_to_steps

end module clock
