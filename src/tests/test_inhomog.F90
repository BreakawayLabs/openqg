module test_inhomog_mod

  use mesh, only: mesh_type, init_mesh
  use box, only: box_type, init_box_from_mesh
  use inhomog, only: inhomog_type, init_inhomog, solve_inhomog_eqn, generate_homog_soln
  use numerics, only: dP2dx2_bc, dP2dy2_bc

  use testlib, only: start_suite, end_suite, start_test, end_test, check_threshold

  implicit none

  private

  public test_inhomog

contains

  subroutine test_inhomog()
    type(mesh_type) :: mesh
    type(box_type) :: b
    type(inhomog_type) :: inhom
    double precision :: rdm2(3) ! 1/r_m^2 for 3 modes

    call start_suite('inhomog')

    ! 100x100 10km cyclic mesh at 55 degrees north
    mesh = init_mesh(100, 100, 10000.0d0, 10000.0d0, 0.0d0, 0.0d0, .true., 0.0d0, 55.0d0)
    
    ! Create a 3 layer box over the entire domain
    b = init_box_from_mesh(mesh, 100, 100, 1, 1, 1, 1, 3, (/1000.0d0, 2000.0d0, 3000.0d0/), 100.0d0, 1)

    ! Inverse squared deformation radius of modes.
    ! Use 0km (barotropic), 50km, 150km.
    rdm2(:) = (/0.0d0, 1.0d0/(50000.0d0**2), 1.0d0/(150000.0d0**2)/)

    ! Create a solver for this system
    inhom = init_inhomog(b, rdm2)

    ! Test solving the inhomogeneous equation
    call test_constant(b, rdm2, inhom, 'cyclic')
    call test_sine(b, rdm2, inhom, 'cyclic')
    call test_random(b, rdm2, inhom, 'cyclic')

    ! Test generating homogeneous solutions
    call test_constant_homog(b, rdm2, inhom, 'cyclic')
    call test_sine_homog(b, rdm2, inhom, 'cyclic')
    call test_random_homog(b, rdm2, inhom, 'cyclic')

    ! 100x100 10km non-cyclic mesh at 55 degrees north
    mesh = init_mesh(100, 100, 10000.0d0, 10000.0d0, 0.0d0, 0.0d0, .false., 0.0d0, 55.0d0)

    ! Create a 3 layer box over the a subset of the domain (80x100)
    b = init_box_from_mesh(mesh, 80, 100, 1, 1, 1, 1, 3, (/1000.0d0, 2000.0d0, 3000.0d0/), 100.0d0, 1)

    ! Create a solver for this system
    inhom = init_inhomog(b, rdm2)

    ! Test solving the inhomogeneous equation
    call test_constant(b, rdm2, inhom, 'box')
    call test_sine(b, rdm2, inhom, 'box')
    call test_random(b, rdm2, inhom, 'box')

    ! Test generating homogeneous solutions
    call test_constant_homog(b, rdm2, inhom, 'box')
    call test_sine_homog(b, rdm2, inhom, 'box')
    call test_random_homog(b, rdm2, inhom, 'box')

    call end_suite('inhomog')

  end subroutine test_inhomog

  subroutine test_constant(b, rdm2, inhom, test_name)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes
    character (len=*), intent(in) :: test_name

    double precision :: rhs_in(b%nxp,b%nyp)
    
    rhs_in(:,:) = 1.0d0
    
    call test_rhs(b, rdm2, inhom, rhs_in, 'constant_'//test_name)
    
  end subroutine test_constant

  subroutine test_sine(b, rdm2, inhom, test_name)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes
    character (len=*), intent(in) :: test_name

    double precision :: rhs_in(b%nxp,b%nyp)    
    integer :: i, j

    do i=1,b%nxp
       do j=1,b%nyp
          rhs_in(:,:) =  2.0d0*sin(5.0d0*i/b%nxp) + 3.0*sin(2.0d0*j/b%nyp)
       enddo
    enddo
    
    call test_rhs(b, rdm2, inhom, rhs_in, 'sine_'//test_name)
    
  end subroutine test_sine

  subroutine test_random(b, rdm2, inhom, test_name)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes
    character (len=*), intent(in) :: test_name

    double precision :: rhs_in(b%nxp,b%nyp)    

    call random_seed()
    call random_number(rhs_in)
    if (b%cyclic) then
       rhs_in(b%nxp,:) = rhs_in(1,:) ! ensure cyclic
    endif
    call test_rhs(b, rdm2, inhom, rhs_in, 'random_'//test_name)
    
  end subroutine test_random


  subroutine test_rhs(b, rdm2, inhom, rhs_in, test_name)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes
    double precision, intent(in) :: rhs_in(b%nxp,b%nyp)
    character (len=*), intent(in) :: test_name

    double precision :: soln(b%nxp,b%nyp)
    double precision :: rhs_out(b%nxp,b%nyp)
    double precision :: alpha, result
    integer :: m

    call start_test('test_solve_inhomog_eqn.'//test_name)

    ! Solve the system for each mode
    do m=1,3
       soln(:,:) = solve_inhomog_eqn(inhom, m, rhs_in)
       alpha = 0.0d0
       rhs_out(:,:) = dP2dx2_bc(soln(:,:), b, alpha) + dP2dy2_bc(soln(:,:), b, alpha) - rdm2(m)*soln(:,:)
       ! Check that (\Delta^2 - 1/r_m^2)soln = rhs_in
       if (b%cyclic) then
          result = maxval(abs(rhs_out(:,2:b%nyp-1) - rhs_in(:,2:b%nyp-1))) / maxval(abs(rhs_in(:,2:b%nyp-1)))
       else
          result = maxval(abs(rhs_out(2:b%nxp-1,2:b%nyp-1) - rhs_in(2:b%nxp-1,2:b%nyp-1))) / &
               maxval(abs(rhs_in(2:b%nxp-1,2:b%nyp-1)))
       end if
       call check_threshold(result, 1.0d-11, test_name)
    enddo

    call end_test('test_solve_inhomog_eqn.'//test_name)

  end subroutine test_rhs

  subroutine test_constant_homog(b, rdm2, inhom, test_name)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)
    type(inhomog_type), intent(inout) :: inhom
    character (len=*), intent(in) :: test_name

    double precision :: L(b%nxp,b%nyp)

    L(:,:) = 1.0d0

    call test_L(b, rdm2, inhom, L, 'constant_'//test_name)

  end subroutine test_constant_homog

  subroutine test_sine_homog(b, rdm2, inhom, test_name)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)
    type(inhomog_type), intent(inout) :: inhom
    character (len=*), intent(in) :: test_name

    double precision :: L(b%nxp,b%nyp)
    integer :: i, j

    do i=1,b%nxp
       do j=1,b%nyp
          L(i,j) = 2.0d0*sin(5.0d0*i/b%nxp) + 3.0*sin(2.0d0*j/b%nyp)
       enddo
    enddo
    if (b%cyclic) then
       L(b%nxp,:) = L(1,:) ! ensure cyclic
    endif
    call test_L(b, rdm2, inhom, L, 'sine_'//test_name)

  end subroutine test_sine_homog

  subroutine test_random_homog(b, rdm2, inhom, test_name)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)
    type(inhomog_type), intent(inout) :: inhom
    character (len=*), intent(in) :: test_name

    double precision :: L(b%nxp,b%nyp)

    call random_seed()
    call random_number(L)
    if (b%cyclic) then
       L(b%nxp,:) = L(1,:) ! ensure cyclic
    endif
    call test_L(b, rdm2, inhom, L, 'random_'//test_name)

  end subroutine test_random_homog

  subroutine test_L(b, rdm2, inhom, L, test_name)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: L(b%nxp,b%nyp)
    character (len=*), intent(in) :: test_name

    double precision :: soln(b%nxp,b%nyp)
    double precision :: rhs_out1(b%nxp,b%nyp), rhs_out2(b%nxp,b%nyp)
    integer :: m
    double precision :: alpha, result

    call start_test('test_generate_homog_soln.'//test_name)

    do m=1,3
       soln(:,:) = generate_homog_soln(inhom, m, L)
       alpha = 0.0d0
       ! (\Delta^2 - 1/r_m^2)soln
       rhs_out1(:,:) = dP2dx2_bc(soln(:,:), b, alpha) + dP2dy2_bc(soln(:,:), b, alpha) - rdm2(m)*soln(:,:)       
       ! \Delta^2 L
       rhs_out2(:,:) = dP2dx2_bc(L(:,:), b, alpha) + dP2dy2_bc(L(:,:), b, alpha)
       ! Check that (\Delta^2 - 1/r_m^2)soln = \Delta^2 L
       if (b%cyclic) then          
          result = maxval(abs(rhs_out1(:,2:b%nyp-1) - rhs_out2(:,2:b%nyp-1)))/maxval(abs(L(:,2:b%nyp-1)))
       else
          result = maxval(abs(rhs_out1(2:b%nxp-1,2:b%nyp-1) - rhs_out2(2:b%nxp-1,2:b%nyp-1)))/maxval(abs(L(2:b%nxp-1,2:b%nyp-1)))
       end if
       call check_threshold(result, 1.0d-22, test_name)
    enddo

    call end_test('test_generate_homog_soln.'//test_name)

  end subroutine test_L
 
end module test_inhomog_mod
