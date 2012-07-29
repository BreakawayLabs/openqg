program openqg

  use mesh, only: mesh_type, init_mesh
  use box, only: box_type, init_box_from_mesh
  use inhomog, only: inhomog_type, init_inhomog, solve_inhomog_eqn, generate_homog_soln
  use numerics, only: dP2dx2_bc, dP2dy2_bc

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

    print *, "##teamcity[testSuiteStarted name='inhomog']"

    call test_constant(b, rdm2, inhom)
    call test_sine(b, rdm2, inhom)

    print *, "##teamcity[testSuiteFinished name='inhomog']"

  end subroutine main

  subroutine test_constant(b, rdm2, inhom)
    type(box_type), intent(in) :: b
    type(inhomog_type), intent(inout) :: inhom
    double precision, intent(in) :: rdm2(b%nl) ! 1/r_m^2 for 3 modes

    double precision :: rhs_in(b%nxp,b%nyp)
    
    rhs_in(:,:) = 1.0d0
    
    call test_rhs(b, rdm2, inhom, rhs_in, 'constant')
    
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
    
    call test_rhs(b, rdm2, inhom, rhs_in, 'sine')
    
  end subroutine test_sine


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


    print *, "##teamcity[testStarted name='", test_name, "' captureStandardOutput='true']"

    ! Solve the system for each mode
    do m=1,3
       soln(:,:) = solve_inhomog_eqn(inhom, m, rhs_in)
       alpha = 0.0d0
       rhs_out(:,:) = dP2dx2_bc(soln(:,:), b, alpha) + dP2dy2_bc(soln(:,:), b, alpha) - rdm2(m)*soln(:,:)
       result = maxval(abs(rhs_out(:,2:b%nyp-1) - rhs_in(:,2:b%nyp-1))) / maxval(abs(rhs_in(:,2:b%nyp-1)))
    enddo
    if (result > 1.0d-12) then
       print *, "##teamcity[testFailed type='comparisonFailure' name='", test_name, &
            "' message='solution error > 1.0d-12' expected='", 1.0d-12 , &
            "' actual='", result , "']]"
    end if

    print *, "##teamcity[testFinished name='", test_name, "']"

  end subroutine test_rhs

end program openqg
