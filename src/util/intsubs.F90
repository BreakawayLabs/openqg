module intsubs

  implicit none

  private

  public trapin

contains

  pure double precision function trapin(fofz, nz, delz)

    ! Computes the extended trapezoidal rule approximation to
    ! the integral of a tabulated function. fofz is a vector
    ! of length nz containing the function tabulation (including
    ! at the end points), delz is the tabulation interval,
    ! and vinteg is the returned value of the integral.

    integer, intent(in) :: nz
    double precision, intent(in) :: fofz(nz),delz

    trapin = delz*(0.5d0*fofz(1) + sum(fofz(2:nz-1)) + 0.5d0*fofz(nz))

  end function trapin
  
end module intsubs
