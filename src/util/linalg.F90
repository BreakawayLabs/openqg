module linalg
  ! 
  ! This module contains basic linear algebra routines.
  !
  ! These subroutines provide wrappers around LAPACK subroutines
  ! to simplify the interface and handle errors.
  ! 
  ! Subroutines:
  !
  ! LU_factor: factorise square matrix A into P * L * U
  ! solve_Ax_b: solve Ax = b given (square matrix A, L, U, P) and b
  ! solve_eigenproblem: Solve A*R = eR, L*A = eL given symmetric matrix A

  implicit none

  private

  public solve_AX_b
  public LU_factor
  public solve_eigenproblem

contains

  subroutine LU_factor(A, LU, ipiv)
    ! Factorise a square matrix A into
    !
    ! A = P * L * U
    !
    ! where L is lower triangular with unit diagonal elements,
    ! U is upper triangular, and P is a permutation matrix.
    ! 
    ! L and U are returned as the matrix LU, with the unit diagonal
    ! of L discarded.
    !
    ! P is returned as a vector of pivot indices, ipiv, where row i of the
    ! matrix was interchanged with row ipiv(i).
    !
    ! LU and ipiv are suitable inputs for solve_Ax_b.
    !
    ! If the factorisation fails then the program terminates with an error
    ! code of 1.
    double precision, intent(in) :: A(:,:)
    double precision, intent(out) :: LU(:,:)
    integer, intent(out) :: ipiv(:)
    
    integer :: n(2), info

    LU(:,:) = A(:,:)
    n = shape(A)

    ! DGETRF = NAG routine F07ADF
    call DGETRF(n(1), n(2), LU, n(1), ipiv, info)
    call check_result(info, 'DGETRF', 'LU_factor')

  end subroutine LU_factor

  subroutine solve_Ax_b(A, LU, ipiv, b, x)
    ! Solve Ax = b for a square matrix A and its factorisation P*L*U
    ! 
    ! LU and ipiv represent the factorisation of A and can be obtained
    ! from LU_factor
    ! 
    ! b is the (input) RHS vector and x is the (output) solution vector.
    !
    ! If the solution finding fails the program terminates with an error
    ! code of 1
    double precision, intent(in) :: A(:,:)
    double precision, intent(in) :: LU(:,:)
    integer, intent(in) :: ipiv(:)
    double precision, intent(in) :: b(:)
    double precision, intent(out) :: x(:)

    integer :: n(2), info
    double precision :: ferr, berr
    double precision :: work(3*size(x))
    integer :: iwork(size(x))

    n = shape(A)

    ! Matrix equation is A*x = P*LU*x = b
    ! Solve the linear system using the LU factorised matrix LU
    ! DGETRS = NAG routine F07AEF
    x(:) = b(:)
    call DGETRS('Norm', n(1), 1, LU, n(1), ipiv, x, n(1), info)
    call check_result(info, 'DGETRS', 'solve_Ax_b')

    ! Improve the solution by iterative refinement
    ! DGERFS = NAG routine F07AHF
    call DGERFS('Norm', n(1), 1, A, n(1), LU, n(1), ipiv, b, n(1), x, n(1), &
         ferr, berr, work, iwork, info)
    call check_result(info, 'DGERFS', 'solve_Ax_b')
    
  end subroutine solve_Ax_b

  subroutine solve_eigenproblem(A_in, n, eigval_real, eigvec_left, eigvec_right)
    ! For a symmetric matrix A, solve the eigenproblem
    !
    ! A*R = eR
    ! L*A = eL.
    !
    ! Returns the real eigenvalues 'e' as eigval_real and the left and right 
    ! eigenvectors L, and R as eigvec_left and eigvec_right.
    !
    ! If the solver fails the program terminates with an error code of 1
    double precision, intent(in) :: A_in(n,n)
    integer, intent(in) :: n
    double precision, intent(out) :: eigval_real(n)
    double precision, intent(out) :: eigvec_left(n,n), eigvec_right(n,n)

    double precision :: A(n,n)
    double precision :: eigval_imag(n)
    double precision :: work(20*n), scale(n), abnrm, rconde(n), rcondv(n)
    integer :: info, iwork(2*n - 2)
    integer :: ihi, ilo

    A(:,:) = A_in(:,:)
    call DGEEVX('B', 'V', 'V', 'B', n, A, n, eigval_real, eigval_imag, &
         eigvec_left, n, eigvec_right, n, ilo, ihi, scale, abnrm, &
         rconde, rcondv, work, size(work), iwork, info )
    call check_result(info, 'DGEEVX', 'solve_eigenproblem')

  end subroutine solve_eigenproblem

  subroutine check_result(info, lapack_name, sub_name)
    ! Check the return value of a call to an LAPACK subroutine
    ! Terminates the program with an error code of 1 if info is non-zero.
    integer, intent(in) :: info
    character (len=*), intent(in) :: lapack_name
    character (len=*), intent(in) :: sub_name
    if (info /= 0) then
       print *,' Problem in ', lapack_name, ', INFO = ', info
       print *,' program terminates in ', sub_name
       stop 1
    endif

  end subroutine check_result

end module linalg
