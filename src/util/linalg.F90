module linalg

  implicit none

  private

  public solve_AX_b
  public LU_factor

contains

  subroutine LU_factor(A, LU, ipiv)
    double precision, intent(in) :: A(:,:)
    double precision, intent(out) :: LU(:,:)
    integer, intent(out) :: ipiv(:)
    
    integer :: n(2), info

    LU(:,:) = A(:,:)
    n = shape(A)

    ! DGETRF = NAG routine F07ADF
    call DGETRF(n(1), n(2), LU, n(1), ipiv, info)
    if (info /= 0) then
       print *,'  DGETRF for returns info = ', info
       print *,'  program terminates in LU_factor'
       stop
    endif

  end subroutine LU_factor

  subroutine solve_Ax_b(A, LU, ipiv, b, x)
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
    call DGETRS ('Norm', n(1), 1, LU, n(1), ipiv, x, n(1), info)
    if (info /= 0) then
       print *,'  DGETRS returns info = ',info
       print *,'  program terminates.'
       stop
    endif
    ! Improve the solution by iterative refinement
    ! DGERFS = NAG routine F07AHF
    call DGERFS ('Norm', n(1), 1, A, n(1), LU, n(1), ipiv, b, n(1), x, n(1), &
         ferr, berr, work, iwork, info)
    
  end subroutine solve_Ax_b

end module linalg
