module test_linalg_mod

  use linalg, only: LU_factor, solve_Ax_b, solve_eigenproblem

  use testlib, only: start_suite, end_suite, start_test, end_test, check_threshold

  implicit none

  private

  public test_linalg

contains

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

    call test_solve_eigenproblem(A3, 3, 'test_solve_eigenproblem.3x3')

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

  subroutine test_solve_eigenproblem(A, n, test_name)
    double precision, intent(in) :: A(n,n)
    integer, intent(in):: n
    character (len=*), intent(in) :: test_name

    double precision :: eigval_real(n), eigvec_left(n,n), eigvec_right(n,n)
    double precision :: AR(n), eR(n), LA(n), eL(n)
    integer :: i

    call start_test(test_name)
    
    call solve_eigenproblem(A, n, eigval_real, eigvec_left, eigvec_right)

    do i=1,n
       AR = matmul(A, eigvec_right(:,i))
       eR = eigvec_right(:,i)*eigval_real(i)

       LA = matmul(eigvec_left(:,i), A)
       eL = eigvec_left(:,i)*eigval_real(i)

       call check_threshold(sum(AR - eR), 1.0d-13, test_name)
       call check_threshold(sum(LA - eL), 1.0d-13, test_name)
    enddo

    call end_test(test_name)

  end subroutine test_solve_eigenproblem

end module test_linalg_mod
