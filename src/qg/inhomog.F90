module inhomog

  use box, only: box_type
  use constants, only: PI, TWOPI

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  private

  type inhomog_type

     double precision, allocatable :: rdm2(:)

     double precision :: a
     type(C_PTR) :: ocn_plan_f
     type(C_PTR) :: ocn_plan_b
     real(C_DOUBLE), pointer :: fftw_ocn_in(:,:)
     real(C_DOUBLE), pointer :: fftw_ocn_out_r(:,:)
     complex(C_DOUBLE), pointer :: fftw_ocn_out_c(:,:)
     type(C_PTR) :: in_p
     type(C_PTR) :: out_p

     integer :: nx
     double precision :: ftnorm
     double precision, allocatable :: bb(:,:)

     type(box_type) :: b

  end type inhomog_type

  public inhomog_type
  public init_inhomog
  public solve_inhomog_eqn
  public generate_homog_soln

contains

  type(inhomog_type) function init_inhomog(b, rdm2)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)

    integer :: i, m, nx
    double precision, allocatable :: bd2(:)

    init_inhomog%b = b

    init_inhomog%a = 1.0d0/( b%dy*b%dy )

    allocate(init_inhomog%rdm2(b%nl))
    init_inhomog%rdm2(:) = rdm2(:)

    if (b%cyclic) then
       nx = b%nxt/2 + 1
       allocate(bd2(nx))
       do i=1,nx
          bd2(i) = -2.0d0*init_inhomog%a + 2.0d0*b%dxm2*(cos((i-1)*TWOPI/b%nxt)  - 1.0d0)
       enddo
       init_inhomog%ftnorm = 1.0d0/b%nxt
    else
       nx = b%nxt - 1
       allocate(bd2(nx))
       do i=1,nx
          bd2(i) = -2.0d0*init_inhomog%a + 2.0d0*b%dxm2*( cos( i*PI/b%nxt ) - 1.0d0 ) ! FIXME: nxt-1
       enddo
       init_inhomog%ftnorm = 0.5d0/b%nxt! FIXME: nxt -1
    endif
    init_inhomog%nx = nx

    allocate(init_inhomog%bb(nx,b%nl))
    do m=1,b%nl
       init_inhomog%bb(:,m) = bd2(:) - rdm2(m)
    enddo
    deallocate(bd2)

    if (b%cyclic) then
       init_inhomog%in_p = fftw_alloc_real(int(b%nxt*(b%nyp-2), C_SIZE_T))
       call c_f_pointer(init_inhomog%in_p, init_inhomog%fftw_ocn_in, [b%nxt, b%nyp-2])

       init_inhomog%out_p = fftw_alloc_complex(int(nx*(b%nyp-2), C_SIZE_T))
       call c_f_pointer(init_inhomog%out_p, init_inhomog%fftw_ocn_out_c, [nx, b%nyp-2])

       init_inhomog%ocn_plan_f = fftw_plan_many_dft_r2c(1, (/b%nxt/), b%nyp-2, &
            init_inhomog%fftw_ocn_in, (/b%nxt/), 1, b%nxt, &
            init_inhomog%fftw_ocn_out_c, (/nx/), 1, nx, &
            ior(FFTW_PATIENT,FFTW_DESTROY_INPUT))
       
       init_inhomog%ocn_plan_b = fftw_plan_many_dft_c2r(1, (/b%nxt/), b%nyp-2, &
            init_inhomog%fftw_ocn_out_c, (/nx/), 1, nx, &
            init_inhomog%fftw_ocn_in, (/b%nxt/), 1, b%nxt, &
            ior(FFTW_PATIENT,FFTW_DESTROY_INPUT))
    else
       init_inhomog%in_p = fftw_alloc_real(int(nx*(b%nyp-2), C_SIZE_T))
       call c_f_pointer(init_inhomog%in_p, init_inhomog%fftw_ocn_in, [nx, b%nyp-2])

       init_inhomog%out_p = fftw_alloc_real(int(nx*(b%nyp-2), C_SIZE_T))
       call c_f_pointer(init_inhomog%out_p, init_inhomog%fftw_ocn_out_r, [nx, b%nyp-2])
       print *, "PLAN"
       init_inhomog%ocn_plan_f = fftw_plan_many_r2r(1, (/nx/), b%nyp-2, &
            init_inhomog%fftw_ocn_in, (/nx/), 1, nx, &
            init_inhomog%fftw_ocn_out_r, (/nx/), 1, nx, &
            (/FFTW_RODFT00/), ior(FFTW_EXHAUSTIVE,FFTW_DESTROY_INPUT))
       print *, "DONE"
    endif

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
       call hscy(inhom, inhom%b, rhs, m, solve_inhomog_eqn)
    else
       call hsbx(inhom, inhom%b, rhs, m, solve_inhomog_eqn)
    endif

  end function solve_inhomog_eqn

  subroutine hsbx (inhom, b, rhs, m, inhomog)

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

  end subroutine hsbx

  subroutine hscy (inhom, b, rhs, m, inhomog)

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

  end subroutine hscy

end module inhomog
