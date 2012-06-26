module print

  use qg, only: qg_type
  use omlsubs, only: ocean_mixed_type
  use amlsubs, only: atmos_mixed_type

  implicit none

  type print_type
     
     integer :: nocmon
     logical :: active = .false.

  end type print_type

  private

  public print_type

  public init_print
  public print_step
  public print_ocn
  public print_atm

contains

  type(print_type) function init_print(nocmon)

    integer, intent(in) :: nocmon

    init_print%nocmon = nocmon
    init_print%active = .true.

  end function init_print

  logical function print_step(print, ntdone, force)

    type(print_type), intent(in) :: print
    integer, intent(in) :: ntdone
    logical, intent(in) :: force

    print_step = print%active .and. (force .or. mod(ntdone,print%nocmon).eq.0 )

  end function print_step

  subroutine print_ocn(qgo, oml)

    type(qg_type), intent(in) :: qgo
    type(ocean_mixed_type), intent(in) :: oml

    integer :: nxco,nyco,k
    double precision :: sstmin,sstmax

    ! Print some ocean spot values and extrema
    if (qgo%active) then
       nxco = (qgo%b%nxp+1)/2
       nyco = (qgo%b%nyp+1)/2

       write(*,217) '  po(k) at centre = ',(qgo%p(nxco,nyco,k),k=1,qgo%b%nl)
       write(*,217) '  qo(k) at centre = ',(qgo%q(nxco,nyco,k),k=1,qgo%b%nl)
    endif
    if (oml%active) then
       sstmin = minval(oml%sst%data)
       sstmax = maxval(oml%sst%data)
       write(*,217) '  s.s.t: min, max = ',sstmin,sstmax
    endif

217 format(a,1p,9d15.7)

  end subroutine print_ocn

  subroutine print_atm(qga, aml)

    type(qg_type), intent(in) :: qga
    type(atmos_mixed_type), intent(in):: aml

    integer :: nxca,nyca,k
    double precision :: astmin,astmax,hmxmin,hmxmax

    ! Print some atmos. spot values and extrema
    if (qga%active) then
       nxca = (qga%b%nxp+1)/2
       nyca = (qga%b%nyp+1)/2

       write(*,217) '  pa(k) at centre = ',(qga%p(nxca,nyca,k),k=1,qga%b%nl)
       write(*,217) '  qa(k) at centre = ',(qga%q(nxca,nyca,k),k=1,qga%b%nl)
    endif
    if (aml%active) then
       astmin = minval(aml%ast%data)
       astmax = maxval(aml%ast%data)
       hmxmin = minval(aml%hmixa%data)
       hmxmax = maxval(aml%hmixa%data)
       write(*,217) '  a.s.t: min, max = ',astmin,astmax
       write(*,217) '  hmixa: min, max = ',hmxmin,hmxmax
    endif

217 format(a,1p,9d15.7)

  end subroutine print_atm

end module print
