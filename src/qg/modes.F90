module modes

  use box, only: box_type
  use topog, only: topog_type

  implicit none

  private

  type modes_type

     double precision, allocatable :: amat(:,:) ! A matrix linking pressures and eta
     double precision, allocatable :: c_phase(:)! phase speeds of modes (m s^-1)
     double precision, allocatable :: rdef(:)   ! deformation radii of modes (m)
     double precision, allocatable :: rdm2(:)   ! 1.0/(deformation radii of modes)**2 (m^-2)
     double precision, allocatable :: ctl2m(:,:)! coefficients for conversion from layers to modes
     double precision, allocatable :: ctm2l(:,:)! transpose of matrix for conversion from layers to modes
     ! transposes are held for reasons of computational efficiency
     ! All the above are got by solving an eigenvalue/vector
     ! equation using subroutine eigmod, called from the main program

  end type modes_type

  public modes_type
  public init_modes
  public eigmod
  public p_modes_rhs
  public modes_to_layers

contains

  type(modes_type) pure function init_modes(nl)
  
    integer, intent(in) :: nl

    allocate(init_modes%amat(nl,nl))
    allocate(init_modes%c_phase(nl))
    allocate(init_modes%rdef(nl))
    allocate(init_modes%rdm2(nl))
    allocate(init_modes%ctl2m(nl,nl))
    allocate(init_modes%ctm2l(nl,nl))

  end function init_modes

  subroutine eigmod (b, gpr, mod)
! *
! *     Derive eigenmodes: eigenvalues and left & right
! *     eigenvectors, which determine the coefficients
! *     for transforming between layer and modal pressures
! *     Works for both oceanic and atmospheric cases

! *     Input arguments:
! *     gpr  : vector of reduced gravities across QG layer interfaces (m s^-2)
! *     h    : vector of unperturbed thicknesses of QG layers (m)
! *     (all the above are unchanged on exit)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: gpr(b%nl-1)
    type(modes_type), intent(inout) :: mod

    logical :: prtval
    parameter ( prtval = .true. )

    double precision :: wre(b%nl),elder(b%nl,b%nl)
    double precision :: evecl(b%nl,b%nl),evecr(b%nl,b%nl)

    call compute_A(b, gpr, mod)

    call compute_eigs(mod%amat, b%nl, wre, evecl, evecr, elder)

    call derive_qg_modes(b, wre, evecl, evecr, elder, mod)

  end subroutine eigmod

  subroutine compute_eigs(amat, nl, eigval_real, eigvec_left, eigvec_right, elder)
    double precision, intent(in) :: amat(nl,nl)
    integer, intent(in) :: nl
    double precision, intent(out) :: eigval_real(nl)
    double precision, intent(out) :: eigvec_left(nl,nl), eigvec_right(nl,nl)
    double precision, intent(out) :: elder(nl,nl)

    double precision :: A(nl,nl)
    double precision :: eigval_imag(nl)
    double precision :: work(20*nl), scale(nl), abnrm, rconde(nl), rcondv(nl)
    integer :: info, iwork(20*nl - 2)
    integer :: n, m, ihi, ilo

    A(:,:) = amat(:,:)
    call DGEEVX('B', 'V', 'V', 'B', nl, A, nl, eigval_real, eigval_imag, &
         eigvec_left, nl, eigvec_right, nl, ilo, ihi, scale, abnrm, &
         rconde, rcondv, work, size(work), iwork, info )
    if (info /= 0) then
       print *,' Problem in DGEEVX, INFO = ', info
       print *,' program terminates in compute_eigs'
       stop 1
    endif

    do n=1,nl
       do m=1,nl
          elder(m,n) = sum(eigvec_left(:,m)*eigvec_right(:,n))
       enddo
    enddo

  end subroutine compute_eigs

  subroutine derive_qg_modes(b, wre, evecl, evecr, elder, mod)
    
    type(box_type), intent(in) :: b
    double precision, intent(in) :: wre(:)
    double precision, intent(in) :: evecl(:,:)
    double precision, intent(in) :: evecr(:,:)
    double precision, intent(in) :: elder(:,:)
    type(modes_type), intent(inout) :: mod

    double precision :: c2rabs(b%nl), eig_val(b%nl), c2temp, cl2m(b%nl,b%nl), cm2l(b%nl,b%nl),&
         ccprod(b%nl,b%nl)
    double precision :: aevec(b%nl),eevec(b%nl)
    integer :: index(b%nl)

    integer :: i, j, k, l, m, inm, indtmp

    logical :: eigchk
    parameter ( eigchk=.true. )

    ! Derive quantities required for OpenQG
    ! ------------------------------------
    ! Eigenvalues are 1/c^2; we want these in increasing order
    ! Barotropic mode is then first in list
    c2rabs(:) = abs( wre(:) )
    do m=1,b%nl
       index(m) = m
    enddo
    ! Index the vector c2rabs, i.e. output the vector index(1:nl) such
    ! that cr2abs( index(m) ) is in ascending order for m = 1,2,...,nl
    do m=2,b%nl
       indtmp = index(m)
       c2temp = c2rabs(indtmp)
       ! Order m-th entry w.r.t. previous (already sorted) ones
       do i=m-1,1,-1
          if (c2rabs(index(i)) <=  c2temp) goto 100
          index(i+1) = index(i)
       enddo
       i = 0
100    index(i+1) = indtmp
    enddo

    ! Modal eigenvalues (barotropic is ~0; replace with exact 0)
    ! and wavespeeds (barotropic is infinite; replace with 0)
    eig_val(1) = 0.0d0
    mod%c_phase(1) = 0.0d0
    do m=2,b%nl
       eig_val(m) = c2rabs(index(m))
       mod%c_phase(m) = 1.0d0/sqrt( c2rabs(index(m)) )
    enddo
    ! Deformation radii
    mod%rdef(1) = 0.0d0
    mod%rdm2(1) = 0.0d0
    do m=2,b%nl
       mod%rdef(m) = 1.0d0/sqrt( c2rabs(index(m)) )/abs(b%fnot)
       mod%rdm2(m) = b%fnot*b%fnot*c2rabs(index(m))
    enddo
    ! Mode/layer conversion coefficients
    ! m = mode number; k = layer number
    do m=1,b%nl
       inm = index(m)
       cl2m(m,:) = evecl(:,inm)/elder(inm,inm)
       mod%ctl2m(:,m) = cl2m(m,:)

       cm2l(:,m) = evecr(:,inm)
       mod%ctm2l(m,:) = cm2l(:,m)
    enddo
    ! Check product of conversion matrices
    do l=1,b%nl
       do k=1,b%nl
          ccprod(k,l) = sum(cm2l(k,:)*cl2m(:,l))
          if (l == k .and. abs(ccprod(k,l) - 1.0d0) > 5.0d-15) then
             print *, "ccprod - Not identity", k, l, ccprod(k,l)
             stop 1
          endif
          if (l /= k .and. abs(ccprod(k,l) - 0.0d0) > 5.0d-15) then
             print *, "ccprod - Not identity", k, l, ccprod(k,l)
             stop 1
          endif
       enddo
    enddo

    ! Write formatted results to standard output
    ! ------------------------------------------
    print *
    print '(2x,a)', 'Eigenmode solver;'
    print '(2x,a)', 'A matrix:'
    print '(2x,a,i9,8i16)', ' j:  ',(j,j=1,b%nl)
    do i=1,b%nl
       print '(2x,i3,1x,1p,9d16.7)', i,(mod%amat(i,j),j=1,b%nl)
    enddo
    print '(2x,a)', 'Modes:'
    print '(2x,a,i9,8i16)', ' m:  ',(m,m=1,b%nl)
    print '(2x,a,1p,9d16.7)', 'eig_val',(eig_val(m),m=1,b%nl)
    print '(2x,a,1p,9d16.7)', 'C_phase',(mod%c_phase(m),m=1,b%nl)
    print '(2x,a,1p,9d16.7)', 'Rdef',(mod%rdef(m),m=1,b%nl)
    print '(2x,a,1p,9d16.7)', 'Rdm2',(mod%rdm2(m),m=1,b%nl)
    print '(2x,a)', 'right eigenvectors Rm(k):'
    do k=1,b%nl
       print '(2x,i3,1x,1p,9d16.7)', k,(evecr(k,index(m)),m=1,b%nl)
    enddo
    print '(2x,a)', 'left  eigenvectors Lm(k):'
    do k=1,b%nl
       print '(2x,i3,1x,1p,9d16.7)', k,(evecl(k,index(m)),m=1,b%nl)
    enddo
    print '(2x,a)', 'layer to mode coefficients cl2m(k,m):'
    do k=1,b%nl
       print '(2x,i3,1x,1p,9d16.7)', k,(cl2m(k,m),m=1,b%nl)
    enddo
    print '(2x,a)', 'mode to layer coefficients cm2l(m,k):'
    do k=1,b%nl
       print '(2x,i3,1x,1p,9d16.7)', k,(cm2l(m,k),m=1,b%nl)
    enddo
    print '(2x,a)', 'matrix product cm2l(m,k)*cl2m(l,m):'
    print '(2x,a,i9,8i16)', ' l:  ',(l,l=1,b%nl)
    do k=1,b%nl
       print '(2x,i3,1x,1p,9d16.7)', k,(ccprod(k,l),l=1,b%nl)
    enddo
    print '(2x,a)', '(should be the identity matrix)'

    ! Optionally verify that the (ordered) eigenvectors are indeed correct
    ! --------------------------------------------------------------------
    if ( eigchk ) then
       print *
       print '(2x,a)', 'Verify ordered eigenmodes:'
       do m=1,b%nl
          inm = index(m)
          do i=1,b%nl
             aevec(i) = sum(mod%amat(i,:)*evecr(:,inm))
             eevec(i) = eig_val(m)*evecr(i,inm)
          enddo
          print '(2x,a,i2)', 'For mode m = ',m
          print '(2x,a,1p,9d17.9)', '  Rm   = ',(evecr(i,inm),i=1,b%nl)
          print '(2x,a,1p,9d17.9)', ' A *Rm = ',(aevec(i),i=1,b%nl)
          print '(2x,a,1p,9d17.9)', 'lam*Rm = ',(eevec(i),i=1,b%nl)
       enddo
    endif

  end subroutine derive_qg_modes

  pure subroutine compute_A(b, gpr, mod)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: gpr(b%nl-1)
    type(modes_type), intent(inout) :: mod

    integer :: k

    ! Setup A matrix
    mod%amat(:,:) = 0.0d0
    ! k = 1 layer (top(ocean) or bottom(atmos.))
    k = 1
    mod%amat(k,k+1) = -1.0d0/( gpr( k )*b%h(k) )
    mod%amat(k, k ) = - mod%amat(k,k+1)
    ! Intermediate layers
    do k=2,b%nl-1
       mod%amat(k,k-1) = -1.0d0/( gpr(k-1)*b%h(k) )
       mod%amat(k,k+1) = -1.0d0/( gpr( k )*b%h(k) )
       mod%amat(k, k ) = - mod%amat(k,k-1) - mod%amat(k,k+1)
    enddo
    ! k = nl layer (bottom(ocean) or top(atmos.))
    k = b%nl
    mod%amat(k,k-1) = -1.0d0/( gpr(k-1)*b%h(k) )
    mod%amat(k, k ) = - mod%amat(k,k-1)

  end subroutine compute_A

  function p_modes_rhs(q, topo, mod, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: q(b%nxp,b%nyp,b%nl)
    type(modes_type), intent(in) :: mod
    type(topog_type), intent(in) :: topo

    double precision :: p_modes_rhs(b%nxp,b%nyp,b%nl)

    double precision :: ql(b%nxp,b%nyp,b%nl)
    integer :: j
    ! Compute vorticity RHS for each mode - equation (8.13)
    do j=1,b%nyp
       ql(:,j,:) = q(:,j,:) - b%beta*b%yprel(j)
       ql(:,j,topo%k_topo) = ql(:,j,topo%k_topo) - topo%ddyn(:,j)
       p_modes_rhs(:,j,:) = b%fnot*matmul(ql(:,j,:), mod%ctl2m(:,:))
    enddo

  end function p_modes_rhs

  pure function modes_to_layers(pm, mod, b)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: pm(b%nxp,b%nyp,b%nl)
    type(modes_type), intent(in) :: mod

    double precision :: modes_to_layers(b%nxp,b%nyp,b%nl)
    integer :: i,j,k

    do j=1,b%nyp
       do i=1,b%nxp
          do k=1,b%nl
             modes_to_layers(i,j,k) = sum(pm(i,j,:)*mod%ctm2l(:,k))
          enddo
       enddo
    enddo

  end function modes_to_layers

end module modes
