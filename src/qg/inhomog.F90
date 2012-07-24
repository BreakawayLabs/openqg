module inhomog

  use box, only: box_type
  use constants, only: PI, TWOPI

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  private

  type inhomog_type

     double precision :: a
     double precision, allocatable :: bd2(:)

     type(C_PTR) :: ocn_plan_f
     type(C_PTR) :: ocn_plan_b
     real(C_DOUBLE), pointer :: fftw_ocn_in(:,:)
     real(C_DOUBLE), pointer :: fftw_ocn_out_r(:,:)
     complex(C_DOUBLE), pointer :: fftw_ocn_out_c(:,:)
     type(C_PTR) :: in_p
     type(C_PTR) :: out_p

  end type inhomog_type

  public inhomog_type
  public init_inhomog

  public hsbx
  public hscy

contains

  type(inhomog_type) function init_inhomog(b)

    type(box_type), intent(in) :: b

    integer :: i

    init_inhomog%a = 1.0d0/( b%dy*b%dy )

    if (b%cyclic) then
       allocate(init_inhomog%bd2(b%nxt/2 + 1))
       do i=1,b%nxt/2 + 1
          init_inhomog%bd2(i) = -2.0d0*init_inhomog%a + 2.0d0*b%dxm2*(cos((i-1)*TWOPI/b%nxt)  - 1.0d0)
       enddo
    else
       allocate(init_inhomog%bd2(b%nxt-1))
       do i=1,b%nxt-1
          init_inhomog%bd2(i) = -2.0d0*init_inhomog%a + 2.0d0*b%dxm2*( cos( i*PI/b%nxt ) - 1.0d0 ) ! FIXME: nxt-1
       enddo
    endif

    if (b%cyclic) then
       init_inhomog%in_p = fftw_alloc_real(int(b%nxt*(b%nyt-1), C_SIZE_T))
       call c_f_pointer(init_inhomog%in_p, init_inhomog%fftw_ocn_in, [b%nxt, b%nyt-1])

       init_inhomog%out_p = fftw_alloc_complex(int((b%nxt/2 + 1)*(b%nyt-1), C_SIZE_T))
       call c_f_pointer(init_inhomog%out_p, init_inhomog%fftw_ocn_out_c, [b%nxt/2 + 1, b%nyt-1])

       init_inhomog%ocn_plan_f = fftw_plan_many_dft_r2c(1, (/b%nxt/), b%nyt-1, &
            init_inhomog%fftw_ocn_in, (/b%nxt/), 1, b%nxt, &
            init_inhomog%fftw_ocn_out_c, (/b%nxt/2 + 1/), 1, b%nxt/2 + 1, &
            ior(FFTW_PATIENT,FFTW_DESTROY_INPUT))
       
       init_inhomog%ocn_plan_b = fftw_plan_many_dft_c2r(1, (/b%nxt/), b%nyt-1, &
            init_inhomog%fftw_ocn_out_c, (/b%nxt/2 + 1/), 1, b%nxt/2 + 1, &
            init_inhomog%fftw_ocn_in, (/b%nxt/), 1, b%nxt, &
            ior(FFTW_PATIENT,FFTW_DESTROY_INPUT))
    else
       init_inhomog%in_p = fftw_alloc_real(int((b%nxt-1)*(b%nyt-1), C_SIZE_T))
       call c_f_pointer(init_inhomog%in_p, init_inhomog%fftw_ocn_in, [b%nxt-1, b%nyt-1])

       init_inhomog%out_p = fftw_alloc_real(int((b%nxt-1)*(b%nyt-1), C_SIZE_T))
       call c_f_pointer(init_inhomog%out_p, init_inhomog%fftw_ocn_out_r, [b%nxt-1, b%nyt-1])
       print *, "PLAN"
       init_inhomog%ocn_plan_f = fftw_plan_many_r2r(1, (/b%nxt-1/), b%nyt-1, &
            init_inhomog%fftw_ocn_in, (/b%nxt-1/), 1, b%nxt-1, &
            init_inhomog%fftw_ocn_out_r, (/b%nxt-1/), 1, b%nxt-1, &
            (/FFTW_RODFT00/), ior(FFTW_EXHAUSTIVE,FFTW_DESTROY_INPUT))
       print *, "DONE"
    endif

  end function init_inhomog

  subroutine hsbx (inhom, b, rhs, bb, inhomog)

    ! Solves the inhomogeneous Helmholtz equation for
    ! given rhs in a domain with meridional boundaries
    ! On entry wrk contains the rhs for a given mode.
    ! On exit wrk contains the modal pressure
    ! solution, including all boundary values.
    ! aa and bb are the coefficients of the
    ! sine transformed version of the equation.
    ! Only needed for box ocean case, hence the #ifndef

    ! This version uses the FFTPACK routine DSINT

    type(inhomog_type), intent(inout) :: inhom
    type(box_type), intent(in) :: b
    double precision, intent(inout) :: rhs(b%nxp,b%nyp)
    double precision, intent(in) :: bb(b%nxt-1)
    double precision, intent(inout) :: inhomog(b%nxp,b%nyp)

    double precision ftnorm

    integer j
   
    double precision gam(b%nxp-2,b%nyp-2),betinv(b%nxp-2),uvec(b%nxp-2,b%nyp-2)

    ftnorm = 0.5d0/b%nxt! FIXME: nxt -1
 
    ! Create a local copy of the fft coeffts + workspace

    ! Compute sine transform of rhs along latitude lines
    ! --------------------------------------------------
    ! N.B. uses extra element wrk(nxp,j) as workspace.
    ! Value in this location does not affect inverse.
    inhom%fftw_ocn_in(:,:) = rhs(2:b%nxp-1,2:b%nyp-1)
    call fftw_execute_r2r(inhom%ocn_plan_f, inhom%fftw_ocn_in, inhom%fftw_ocn_out_r)

    ! For each wavenumber i, solve for sine transform of p
    ! ----------------------------------------------------
    ! Compute solution in vector uvec
    ! Tridiagonal solver routine now inlined.
    betinv(:) = 1.0d0/bb(:)
    uvec(:,1) = inhom%fftw_ocn_out_r(:,1)*betinv(:)
    ! Decomposition and forward substitution.
    do j=2,b%nyp-2
       gam(:,j) = inhom%a*betinv(:)
       betinv(:) = 1.0d0/( bb(:) - inhom%a*gam(:,j) )
       uvec(:,j) = ( inhom%fftw_ocn_out_r(:,j) - inhom%a*uvec(:,j-1) )*betinv(:)
    enddo
    ! Backsubstitution.
    do j=b%nyp-3,1,-1
       uvec(:,j) = uvec(:,j) - gam(:,j+1)*uvec(:,j+1)
    enddo
    ! Copy back solution and rescale
    inhom%fftw_ocn_in(:,:) = ftnorm*uvec(:,:)

    ! Inverse sine transform along latitude lines
    ! -------------------------------------------
    call fftw_execute_r2r(inhom%ocn_plan_f, inhom%fftw_ocn_in, inhom%fftw_ocn_out_r)
    inhomog(2:b%nxp-1,2:b%nyp-1) = inhom%fftw_ocn_out_r(:,:)
    ! Zero pressure on Western & Eastern boundaries
    inhomog(1,:) = 0.0d0
    inhomog(b%nxp,:) = 0.0d0

    ! Impose N & S boundary values, which are
    ! implicit in the tridiagonal formulation.
    inhomog(:,1) = 0.0d0
    inhomog(:,b%nyp) = 0.0d0

  end subroutine hsbx

  subroutine hscy (inhom, b, rhs, bb, inhomog)

    ! Solves the inhomogeneous Helmholtz equation
    ! for given rhs in a zonally periodic domain.
    ! On entry wrk contains the rhs for a given mode.
    ! On exit wrk contains the modal pressure
    ! solution, including zonal boundary values.
    ! aa and bb are the coefficients of the
    ! Fourier transformed form of the equation.
    ! Version specifically optimised for ocean arrays

    type(inhomog_type), intent(inout) :: inhom
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rhs(b%nxp,b%nyp)
    double precision, intent(in) :: bb(b%nxt/2 + 1)
    double precision, intent(out) :: inhomog(b%nxp,b%nyp)

    double precision :: ftnorm

    integer :: j
    double precision :: gam(b%nxt/2+1,b%nyp-2),betinv(b%nxt/2+1)
    double complex :: uvec(b%nxt/2+1,b%nyp-2)

    ftnorm = 1.0d0/b%nxt

    ! Compute FFT of rhs along latitude lines
    ! ---------------------------------------
    inhom%fftw_ocn_in(:,:) = rhs(:b%nxp-1,2:b%nyp-1)
    call fftw_execute_dft_r2c(inhom%ocn_plan_f, inhom%fftw_ocn_in, inhom%fftw_ocn_out_c)

    ! For each wavenumber i, solve for FFT of p
    ! -----------------------------------------
    ! Compute solution in vector uvec
    ! Tridiagonal solver routine now inlined.
    betinv(:) = 1.0d0/bb(:)
    uvec(:,1) = inhom%fftw_ocn_out_c(:,1)*betinv(:)
    ! Decomposition and forward substitution.
    do j=2,b%nyp-2
       gam(:,j) = inhom%a*betinv(:)
       betinv(:) = 1.0d0/( bb(:) - inhom%a*gam(:,j) )
       uvec(:,j) = ( inhom%fftw_ocn_out_c(:,j) - inhom%a*uvec(:,j-1) )*betinv(:)
    enddo
    ! Backsubstitution.
    do j=b%nyp-3,1,-1
       uvec(:,j) = uvec(:,j) - gam(:,j+1)*uvec(:,j+1)
    enddo
    ! Copy back solution and rescale
    inhom%fftw_ocn_out_c(:,:) = ftnorm*uvec(:,:)

    call fftw_execute_dft_c2r(inhom%ocn_plan_b, inhom%fftw_ocn_out_c, inhom%fftw_ocn_in)
    inhomog(:b%nxp-1,2:b%nyp-1) = inhom%fftw_ocn_in(:,:)
    ! Impose N & S boundary values, which are
    ! implicit in the tridiagonal formulation.
    inhomog(:,1) = 0.0d0
    inhomog(:,b%nyp) = 0.0d0
    ! Cyclic condition
    inhomog(b%nxp,:) = inhomog(1,:)

  end subroutine hscy

end module inhomog
