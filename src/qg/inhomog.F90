module inhomog
  ! This module provides a solver for the inhomogeneous Helmholtz equation.
  !
  ! The inhomogenenous Helmholtz equation is given by
  !
  ! (\Delta^2_H - 1/r_m^2)p_m = Q_m.
  !
  ! A solver of type(inhomog_type) can be created for 
  ! a given vector 1/r_m^2 over a domain of type(box_type).
  !
  ! This solver can then be used to compute p_m given Q_m.
  ! Furthermore, given a function L with the property
  ! \Delta^2_H L = 0, a solution to the homogeneous equation
  !
  ! (\Delta^2_H - 1/r_m^2)p_m = 0
  ! 
  ! can be generated.
  !
  ! The general technique for solving the inhomogeneous equation
  ! for a given Q_m is
  ! 1. Fourier transform Q_m
  ! 2. Solve a tridiagonal matrix equation
  ! 3. Inverse fourier transform the result to get the solution p_m.
  !
  ! The library FFTW (fftw.org) is used for FFT computations. The 
  ! tridagonal matrix algorithm is taken from Numerical Recipes (Press et. al. 1992).
  ! 
  ! Full documentation can be found in doc/openqg-guide.

  use box, only: box_type
  use constants, only: TWOPI

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  private

  type inhomog_type

     ! Solver input configuration 
     double precision, allocatable :: rdm2(:)
     type(box_type) :: b

     ! FFTW parameters
     type(C_PTR) :: ocn_plan_f
     type(C_PTR) :: ocn_plan_b
     real(C_DOUBLE), pointer :: fftw_ocn_in(:,:)
     real(C_DOUBLE), pointer :: fftw_ocn_out_r(:,:)
     complex(C_DOUBLE), pointer :: fftw_ocn_out_c(:,:)
     type(C_PTR) :: in_p
     type(C_PTR) :: out_p
     double precision :: ftnorm

     ! Matrix parameters
     integer :: nx
     double precision :: a
     double precision, allocatable :: bb(:,:)

  end type inhomog_type

  public inhomog_type
  public init_inhomog
  public solve_inhomog_eqn
  public generate_homog_soln

contains

  type(inhomog_type) function init_inhomog(b, rdm2)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)

    integer :: i, m, N, n_in, n_out
    double precision, allocatable :: bd2(:)

    ! Input parameters
    init_inhomog%b = b
    allocate(init_inhomog%rdm2(b%nl))
    init_inhomog%rdm2(:) = rdm2(:)

    ! Compute the logical and physical size of the transform arrays
    if (b%cyclic) then
       n_in = b%nxp - 1   ! Physical number of transform input points
       n_out = n_in/2 + 1 ! Physical number of transform output points
       N = n_in           ! The logical size of the transform
    else
       n_in = b%nxp - 2
       n_out = n_in
       N = 2*(n_in + 1)
    endif
    init_inhomog%nx = n_out

    ! Setup FFTW inputs/output/plans
    if (b%cyclic) then
       init_inhomog%in_p = fftw_alloc_real(int(n_in*(b%nyp-2), C_SIZE_T)) ! Real input array
       call c_f_pointer(init_inhomog%in_p, init_inhomog%fftw_ocn_in, [n_in, b%nyp-2])

       init_inhomog%out_p = fftw_alloc_complex(int(n_out*(b%nyp-2), C_SIZE_T)) ! Complex output array
       call c_f_pointer(init_inhomog%out_p, init_inhomog%fftw_ocn_out_c, [n_out, b%nyp-2])

       init_inhomog%ocn_plan_f = fftw_plan_many_dft_r2c(1, (/N/), b%nyp-2, & ! Forward transform
            init_inhomog%fftw_ocn_in,    (/n_in/),  1, n_in,  &  ! Takes n_in real elements as input
            init_inhomog%fftw_ocn_out_c, (/n_out/), 1, n_out, &  ! Makes n_out complex elements as output
            ior(FFTW_PATIENT,FFTW_DESTROY_INPUT))
       
       init_inhomog%ocn_plan_b = fftw_plan_many_dft_c2r(1, (/N/), b%nyp-2, & ! Backwards transform
            init_inhomog%fftw_ocn_out_c, (/n_out/), 1, n_out, & ! Takes n_out complex elements as input
            init_inhomog%fftw_ocn_in,    (/n_in/),  1, n_in,  & ! Makes n_in real elements as input
            ior(FFTW_PATIENT,FFTW_DESTROY_INPUT))
    else
       init_inhomog%in_p = fftw_alloc_real(int(n_in*(b%nyp-2), C_SIZE_T)) ! real input array
       call c_f_pointer(init_inhomog%in_p, init_inhomog%fftw_ocn_in, [n_in, b%nyp-2])

       init_inhomog%out_p = fftw_alloc_real(int(n_out*(b%nyp-2), C_SIZE_T)) ! real output array
       call c_f_pointer(init_inhomog%out_p, init_inhomog%fftw_ocn_out_r, [n_out, b%nyp-2])

       init_inhomog%ocn_plan_f = fftw_plan_many_r2r(1, (/n_in/), b%nyp-2, &
            init_inhomog%fftw_ocn_in,    (/n_in/),  1, n_in,  & ! Takes n_in real elements as input
            init_inhomog%fftw_ocn_out_r, (/n_out/), 1, n_out, & ! Makes n_out real elements as output
            (/FFTW_RODFT00/), ior(FFTW_EXHAUSTIVE,FFTW_DESTROY_INPUT))
    endif
    init_inhomog%ftnorm = 1.0d0/N ! Normalisation factor

    ! Tridiagonal matrix parameters
    init_inhomog%a = 1.0d0/( b%dy*b%dy )
    allocate(bd2(init_inhomog%nx))
    if (b%cyclic) then
       do i=1,init_inhomog%nx
          bd2(i) = -2.0d0*init_inhomog%a + 2.0d0*b%dxm2*(cos((i-1)*TWOPI/N)  - 1.0d0)
       enddo
    else
       do i=1,init_inhomog%nx
          bd2(i) = -2.0d0*init_inhomog%a + 2.0d0*b%dxm2*( cos( i*TWOPI/N ) - 1.0d0 )
       enddo
    endif
    allocate(init_inhomog%bb(init_inhomog%nx,b%nl))
    do m=1,b%nl
       init_inhomog%bb(:,m) = bd2(:) - rdm2(m)
    enddo
    deallocate(bd2)

  end function init_inhomog

  function generate_homog_soln(inhom, m, L)
    type(inhomog_type), intent(inout) :: inhom
    integer, intent(in) :: m
    double precision, intent(in) :: L(inhom%b%nxp,inhom%b%nyp)
    double precision :: generate_homog_soln(inhom%b%nxp,inhom%b%nyp)

    generate_homog_soln(:,:) = L(:,:) + inhom%rdm2(m)*solve_inhomog_eqn(inhom, m, L)

  end function generate_homog_soln

  function solve_inhomog_eqn(inhom, m, rhs)
    type(inhomog_type), intent(inout) :: inhom
    integer, intent(in) :: m
    double precision, intent(in) :: rhs(inhom%b%nxp,inhom%b%nyp)
    double precision :: solve_inhomog_eqn(inhom%b%nxp,inhom%b%nyp)

    if (inhom%b%cyclic) then
       call solve_inhomog_cyclic(inhom, inhom%b, rhs, m, solve_inhomog_eqn)
    else
       call solve_inhomog_box(inhom, inhom%b, rhs, m, solve_inhomog_eqn)
    endif

  end function solve_inhomog_eqn

  subroutine solve_inhomog_box(inhom, b, rhs, m, inhomog)

    ! Solves the inhomogeneous Helmholtz equation for given rhs in a domain with meridional boundaries
    type(inhomog_type), intent(inout) :: inhom
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rhs(b%nxp,b%nyp)

    integer, intent(in) :: m
    double precision, intent(inout) :: inhomog(b%nxp,b%nyp)

    integer :: j   
    double precision :: gam(inhom%nx,b%nyp-2),betinv(inhom%nx),uvec(inhom%nx,b%nyp-2)
 
    ! Compute sine transform of rhs along latitude lines
    inhom%fftw_ocn_in(:,:) = rhs(2:b%nxp-1,2:b%nyp-1)
    call fftw_execute_r2r(inhom%ocn_plan_f, inhom%fftw_ocn_in, inhom%fftw_ocn_out_r)

    ! For each wavenumber i, solve for sine transform of p
    betinv(:) = 1.0d0/inhom%bb(:,m)
    uvec(:,1) = inhom%fftw_ocn_out_r(:,1)*betinv(:)
    ! Decomposition and forward substitution.
    do j=2,b%nyp-2
       gam(:,j) = inhom%a*betinv(:)
       betinv(:) = 1.0d0/( inhom%bb(:,m) - inhom%a*gam(:,j) )
       uvec(:,j) = ( inhom%fftw_ocn_out_r(:,j) - inhom%a*uvec(:,j-1) )*betinv(:)
    enddo
    ! Back substitution.
    do j=b%nyp-3,1,-1
       uvec(:,j) = uvec(:,j) - gam(:,j+1)*uvec(:,j+1)
    enddo

    ! Rescale and perform inverse transform along latitude lines
    inhom%fftw_ocn_in(:,:) = inhom%ftnorm*uvec(:,:)
    call fftw_execute_r2r(inhom%ocn_plan_f, inhom%fftw_ocn_in, inhom%fftw_ocn_out_r)
    inhomog(2:b%nxp-1,2:b%nyp-1) = inhom%fftw_ocn_out_r(:,:)

    ! Zero pressure on all boundaries
    inhomog(1,:) = 0.0d0
    inhomog(b%nxp,:) = 0.0d0
    inhomog(:,1) = 0.0d0
    inhomog(:,b%nyp) = 0.0d0

  end subroutine solve_inhomog_box

  subroutine solve_inhomog_cyclic(inhom, b, rhs, m, inhomog)

    ! Solves the inhomogeneous Helmholtz equation for given rhs in a zonally periodic domain.
    type(inhomog_type), intent(inout) :: inhom
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rhs(b%nxp,b%nyp)
    integer, intent(in) :: m
    double precision, intent(out) :: inhomog(b%nxp,b%nyp)

    integer :: j
    double precision :: gam(inhom%nx,b%nyp-2),betinv(inhom%nx)
    double complex :: uvec(inhom%nx,b%nyp-2)

    ! Compute FFT of rhs along latitude lines
    inhom%fftw_ocn_in(:,:) = rhs(:b%nxp-1,2:b%nyp-1)
    call fftw_execute_dft_r2c(inhom%ocn_plan_f, inhom%fftw_ocn_in, inhom%fftw_ocn_out_c)

    ! For each wavenumber i, solve for FFT of p
    betinv(:) = 1.0d0/inhom%bb(:,m)
    uvec(:,1) = inhom%fftw_ocn_out_c(:,1)*betinv(:)
    ! Decomposition and forward substitution.
    do j=2,b%nyp-2
       gam(:,j) = inhom%a*betinv(:)
       betinv(:) = 1.0d0/( inhom%bb(:,m) - inhom%a*gam(:,j) )
       uvec(:,j) = ( inhom%fftw_ocn_out_c(:,j) - inhom%a*uvec(:,j-1) )*betinv(:)
    enddo
    ! Back substitution.
    do j=b%nyp-3,1,-1
       uvec(:,j) = uvec(:,j) - gam(:,j+1)*uvec(:,j+1)
    enddo

    ! Rescale and perform inverse transform along latitude lines
    inhom%fftw_ocn_out_c(:,:) = inhom%ftnorm*uvec(:,:)
    call fftw_execute_dft_c2r(inhom%ocn_plan_b, inhom%fftw_ocn_out_c, inhom%fftw_ocn_in)
    inhomog(:b%nxp-1,2:b%nyp-1) = inhom%fftw_ocn_in(:,:)

    ! Zero pressure on N/S boundaries. Cyclic E/W
    inhomog(:,1) = 0.0d0
    inhomog(:,b%nyp) = 0.0d0
    inhomog(b%nxp,:) = inhomog(1,:)

  end subroutine solve_inhomog_cyclic

end module inhomog
