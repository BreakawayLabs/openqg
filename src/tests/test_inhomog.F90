program test_inhomog

  use mesh, only: mesh_type, init_mesh
  use box, only: box_type, init_box_from_mesh
  use inhomog, only: inhomog_type, init_inhomog, solve_inhomog_eqn, generate_homog_soln
  use numerics, only: dP2dx2_bc, dP2dy2_bc
  use linalg, only: LU_factor, solve_Ax_b

  implicit none

  call main()

contains

  subroutine main()

    type(mesh_type) :: mesh
    type(box_type) :: b
    type(inhomog_type) :: inhom
    double precision :: rdm2(3) ! 1/r_m^2 for 3 modes

    ! 100x100 10km cyclic mesh at 55 degrees north
    mesh = init_mesh(100, 100, 10000.0d0, 10000.0d0, 0.0d0, 0.0d0, .true., 0.0d0, 55.0d0)
    
    ! Create a 3 layer box over the entire domain
    b = init_box_from_mesh(mesh, 100, 100, 1, 1, 1, 1, 3, (/1000.0d0, 2000.0d0, 3000.0d0/), 100.0d0, 1)

    ! Inverse squared deformation radius of modes.
    ! Use 0km (barotropic), 50km, 150km.
    rdm2(:) = (/0.0d0, 1.0d0/(50000.0d0**2), 1.0d0/(150000.0d0**2)/)

    ! Create a solver for this system
    inhom = init_inhomog(b, rdm2)

    call start_suite('inhomog.cyclic')

    ! Test solving the inhomogeneous equation
    call test_constant(b, rdm2, inhom)
    call test_sine(b, rdm2, inhom)
    call test_random(b, rdm2, inhom)

    ! Test generating homogeneous solutions
    call test_constant_homog(b, rdm2, inhom)
    call test_sine_homog(b, rdm2, inhom)
    call test_random_homog(b, rdm2, inhom)

    call end_suite('inhomog.cyclic')
    
    ! 100x100 10km non-cyclic mesh at 55 degrees north
    mesh = init_mesh(100, 100, 10000.0d0, 10000.0d0, 0.0d0, 0.0d0, .false., 0.0d0, 55.0d0)

    ! Create a 3 layer box over the a subset of the domain (80x100)
    b = init_box_from_mesh(mesh, 80, 100, 1, 1, 1, 1, 3, (/1000.0d0, 2000.0d0, 3000.0d0/), 100.0d0, 1)

    ! Create a solver for this system
    inhom = init_inhomog(b, rdm2)

    call start_suite('inhomog.box')    

    ! Test solving the inhomogeneous equation
    call test_constant(b, rdm2, inhom)
    call test_sine(b, rdm2, inhom)
    call test_random(b, rdm2, inhom)

    ! Test generating homogeneous solutions
    call test_constant_homog(b, rdm2, inhom)
    call test_sine_homog(b, rdm2, inhom)
    call test_random_homog(b, rdm2, inhom)

    call end_suite('inhomog.box')

    call test_linalg()

  end subroutine main

  subroutine test_constant(b, rdm2, inhom)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes

    double precision :: rhs_in(b%nxp,b%nyp)
    
    rhs_in(:,:) = 1.0d0
    
    call test_rhs(b, rdm2, inhom, rhs_in, 'solve_homog_constant')
    
  end subroutine test_constant

  subroutine test_sine(b, rdm2, inhom)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes

    double precision :: rhs_in(b%nxp,b%nyp)    
    integer :: i, j

    do i=1,b%nxp
       do j=1,b%nyp
          rhs_in(:,:) =  2.0d0*sin(5.0d0*i/b%nxp) + 3.0*sin(2.0d0*j/b%nyp)
       enddo
    enddo
    
    call test_rhs(b, rdm2, inhom, rhs_in, 'solve_homog_sine')
    
  end subroutine test_sine

  subroutine test_random(b, rdm2, inhom)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes

    double precision :: rhs_in(b%nxp,b%nyp)    

    call random_seed()
    call random_number(rhs_in)
    if (b%cyclic) then
       rhs_in(b%nxp,:) = rhs_in(1,:) ! ensure cyclic
    endif
    call test_rhs(b, rdm2, inhom, rhs_in, 'solve_homog_random')
    
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

    call start_test(test_name)

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

    call end_test(test_name)

  end subroutine test_rhs

  subroutine test_constant_homog(b, rdm2, inhom)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)
    type(inhomog_type), intent(inout) :: inhom

    double precision :: L(b%nxp,b%nyp)

    L(:,:) = 1.0d0

    call test_L(b, rdm2, inhom, L, 'generate_homog_constant')

  end subroutine test_constant_homog

  subroutine test_sine_homog(b, rdm2, inhom)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)
    type(inhomog_type), intent(inout) :: inhom

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
    call test_L(b, rdm2, inhom, L, 'generate_homog_sine')

  end subroutine test_sine_homog

  subroutine test_random_homog(b, rdm2, inhom)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: rdm2(b%nl)
    type(inhomog_type), intent(inout) :: inhom

    double precision :: L(b%nxp,b%nyp)

    call random_seed()
    call random_number(L)
    if (b%cyclic) then
       L(b%nxp,:) = L(1,:) ! ensure cyclic
    endif
    call test_L(b, rdm2, inhom, L, 'generate_homog_random')

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

    call start_test(test_name)

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

    call end_test(test_name)

  end subroutine test_L

  subroutine test_linalg()

    double precision :: A3(3,3)
    double precision :: L3(3,3), U3(3,3)
    double precision :: A10(10,10)
    double precision :: L10(10,10), U10(10,10)
    integer :: d, i, j

    call start_suite('linalg')

    ! Create a matrix with a non-singular U to ensure successful factorisation
    L3 = 0.0d0
    do d=1,3
       L3(d,d) = 1.0d0
    enddo
    do i=2,3
       do j=1,i-1
          L3(i,j) = 1.0d0*(i + j)
       enddo
    enddo
    U3 = 0.0d0
    do i=1,3
       do j=i,3
          U3(i,j) = 1.0d0*(i + j)
       enddo
    enddo
    A3 = matmul(L3, U3)
    call test_LU_factor(A3, 3, 'test_LU_factor.3x3')

    ! Create a matrix with a non-singular U to ensure successful factorisation
    L10 = 0.0d0
    do d=1,10
       L10(d,d) = 1.0d0
    enddo
    do i=2,10
       do j=1,i-1
          L10(i,j) = 1.0d0*(i + j)
       enddo
    enddo
    U10 = 0.0d0
    do i=1,10
       do j=i,10
          U10(i,j) = 1.0d0*(i + j)
       enddo
    enddo
    A10 = matmul(L10, U10)
    call test_LU_factor(A10, 10, 'test_LU_factor.10x10')

    call test_solve_Ax_b(A3, (/1.0d0, 2.0d0, 3.0d0/), 3, 'test_solve_AX_b.3x3')

    call end_suite('linalg')

  end subroutine test_linalg

  subroutine test_solve_Ax_b(A, b, n, test_name)
    double precision, intent(in) :: A(n,n)
    double precision, intent(in) :: b(n)
    integer, intent(in) :: n
    character (len=*), intent(in) :: test_name

    double precision :: LU(n,n), x(n), b_out(n), result
    integer :: ipiv(n)

    call start_test(test_name)
    
    call LU_factor(A, LU, ipiv)
    call solve_Ax_b(A, LU, ipiv, b, x)

    b_out = matmul(A, x)
    result = sum(abs(b - b_out))
    call check_threshold(result, 1.0d-14, test_name)

    call end_test(test_name)    

  end subroutine test_solve_Ax_b

  subroutine test_LU_factor(A, n, test_name)

    double precision, intent(in) :: A(n,n)
    integer, intent(in) :: n
    character (len=*), intent(in) :: test_name

    double precision :: LU(n,n), L(n,n), U(n,n), A_out(n,n), tmp_row(n)
    integer :: ipiv(n)
    integer :: d, i, j
    double precision :: result

    call start_test(test_name)

    call LU_factor(A, LU, ipiv)

    L = 0.0d0
    do d=1,n
       L(d,d) = 1.0d0
    enddo
    do i=2,n
       do j=1,i-1
          L(i,j) = LU(i,j)
       enddo
    enddo
    U = 0.0d0
    do i=1,n
       do j=i,n
          U(i,j) = LU(i,j)
       enddo
    enddo
    A_out = matmul(L, U)
    do i=1,n
       tmp_row = A_out(i,:)
       A_out(i,:) = A_out(ipiv(i),:)
       A_out(ipiv(i),:) = tmp_row
    enddo

    result = abs(sum(A_out - A))
    call check_threshold(result, 1.0d-12, test_name)

    call end_test(test_name)

  end subroutine test_LU_factor

  subroutine start_suite(suite_name)
    character (len=*), intent(in) :: suite_name
    print *, "##teamcity[testSuiteStarted name='", suite_name, "']"
  end subroutine start_suite

  subroutine end_suite(suite_name)
    character (len=*), intent(in) :: suite_name
    print *, "##teamcity[testSuiteFinished name='", suite_name, "']"
  end subroutine end_suite

  subroutine start_test(test_name)
    character (len=*), intent(in) :: test_name
    print *, "##teamcity[testStarted name='", test_name, "' captureStandardOutput='true']"
  end subroutine start_test

  subroutine end_test(test_name)
    character (len=*), intent(in) :: test_name
    print *, "##teamcity[testFinished name='", test_name, "']"
  end subroutine end_test

  subroutine check_threshold(result, threshold, test_name)
    double precision, intent(in) :: result
    double precision, intent(in) :: threshold
    character (len=*), intent(in) :: test_name

    if (result > threshold) then
       print *, "##teamcity[testFailed type='comparisonFailure' name='", test_name, &
            "' message='result > threshold' expected='", threshold , &
            "' actual='", result , "']]"
    end if
  end subroutine check_threshold

end program test_inhomog
