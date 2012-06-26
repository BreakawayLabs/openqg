module intsubs

  implicit none

  private

  public trapin
  public xintp
  public genint

contains

  double precision function xintp(valp, nxp, nyp)

    ! Computes area integral of a field valp tabulated at p points
    ! Includes weighted contribution from boundary points
    ! Modified version with reduced accumulator error
    
    double precision, intent(in) :: valp(nxp,nyp)
    integer, intent(in) :: nxp, nyp

    integer i,j
    double precision sump,xxs,xxn

    sump = 0.0d0
    xxs = 0.5d0*valp(1, 1 )
    xxn = 0.5d0*valp(1,nyp)

    ! Inner points + 0.5d0*( W & E boundaries)
    do j=2,nyp-1
       sump = sump + 0.5d0*valp(1,j) + sum(valp(2:nxp-1,j)) + 0.5d0*valp(nxp,j)
    enddo

    ! N & S boundary contributions from inner
    ! points (do corners separately because
    ! they need an additional factor of 0.5)
    do i=2,nxp-1
       xxs = xxs + valp(i, 1 )
       xxn = xxn + valp(i,nyp)
    enddo

    xxs = xxs + 0.5d0*valp(nxp, 1 )
    xxn = xxn + 0.5d0*valp(nxp,nyp)

    xintp = sump + 0.5d0*( xxs + xxn )

  end function xintp


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


  double precision function genint (val, nx, ny, facwe, facsn)

    ! Computes area integral of internal values of field
    ! val(nx,ny). facwe controls contribution from western
    ! & eastern boundaries; facsn controls contribution from
    ! southern & northern boundaries. Integral returned as answer
    ! Modified version with reduced accumulator error

    ! N.B. Parallelised using "orphaned directives";
    ! only intended for use in the parallel regions
    ! established earlier in calling subroutine diagno

    integer, intent(in) :: nx,ny
    double precision, intent(in) :: val(nx,ny),facwe,facsn

    integer :: i,j
    double precision :: sumi,xxs,xxn

    genint = 0.0d0
    ! Inner points + facwe*(W & E boundaries)
    do j=2,ny-1
       sumi = facwe*val(1,j)
       do i=2,nx-1
          sumi = sumi + val(i,j)
       enddo
       sumi = sumi + facwe*val(nx,j)
       genint = genint + sumi
    enddo

    ! N & S boundary contributions
    xxs = facwe*val(1, 1)
    xxn = facwe*val(1,ny)
    ! Inner points
    do i=2,nx-1
       xxs = xxs + val(i, 1)
       xxn = xxn + val(i,ny)
    enddo
    xxs = xxs + facwe*val(nx, 1)
    xxn = xxn + facwe*val(nx,ny)

    genint = genint + facsn*( xxs + xxn )

  end function genint
  
end module intsubs
