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

    call compute_eigs(mod, prtval, b%nl, wre, evecl, evecr, elder)

    call derive_qg_modes(b, wre, evecl, evecr, elder, mod)

  end subroutine eigmod

  subroutine compute_eigs(mod, prtval, nl, wre, evecl, evecr, elder)

    type(modes_type), intent(in) :: mod
    integer, intent(in) :: nl
    logical, intent(in) :: prtval
    double precision, intent(out) :: wre(:)
    double precision, intent(out) :: evecl(:,:)
    double precision, intent(out) :: evecr(:,:)
    double precision, intent(out) :: elder(:,:)

    integer :: nblk,lwork
    parameter ( nblk=64 )

    double precision :: scale(nl), work(nl*nblk), wim(nl), tauvec(nl)
    double precision :: aba(nl,nl),hup(nl,nl),qqq(nl,nl),zzz(nl,nl),uqt(nl,nl)
    logical :: select(nl)

    integer :: i, j, k, n, m, mwant, mgot, info, ilo, ihi

    character :: baljob*(*)
    parameter ( baljob = 'Both' )

    lwork = nl*nblk

    ! Solve for the eigenvalues and L & R eigenvectors of A
    ! -----------------------------------------------------
    ! (real nonsymmetric eigenproblem). Use LAPACK routines
    ! DGEBAL = F08NHF
    ! DGEHRD = F08NEF
    ! DORGHR = F08NFF
    ! DHSEQR = F08PEF
    ! DTREVC = F08QKF
    ! DGEBAK = F08NJF
    ! Set default values of ILO, IHI.
    ilo = 1
    ihi = nl
    ! ------------------------------------------------------
    ! Balance a real general matrix to improve accuracy
    ! Use a copy of A in case the original is required
    ! BALJOB passed to DGEBAK to transform eigenvectors suitably
    aba(:,:) = mod%amat(:,:)
    call DGEBAL (baljob, nl, aba, nl, ilo, ihi, scale, info)
    call check_result(info, 'DGEBAL')

    ! ABA is now overwritten with the balanced matrix
    ! ------------------------------------------------------
    ! Reduce a real general matrix to Hessenberg form
    ! N.B. overwrites matrix with its upper Hessenberg form H
    ! and details of the orthogonal transformation matrix Q
    hup(:,:) = aba(:,:)
    call DGEHRD (nl, ilo, ihi, hup, nl, tauvec, work, lwork, info)
    call check_result(info, 'DHEHRD')

    !------------------------------------------------------
    ! Generate the real orthogonal matrix Q
    ! that reduces A to upper Hessenberg form
    ! Provide as input a copy of the upper Hessenberg matrix
    ! previously generated, which also contains details of
    ! the vectors which define the elementary reflectors
    qqq(:,:) = hup(:,:)
    call DORGHR (nl, ilo, ihi, qqq, nl, tauvec, work, lwork, info)
    call check_result(info, 'DORGHR')

    ! ------------------------------------------------------
    ! Compute all the eigenvalues, and optionally the Schur
    ! factorization, of a real upper Hessenberg matrix
    ! We require the Schur form in order to get the eigenvectors
    ! This means that UQT is overwritten with the upper
    ! quasi-triangular matrix T from the Schur decomposition
    ! Set ZZZ to contain the matrix Q on entry.
    uqt(:,:) = hup(:,:)
    zzz(:,:) = qqq(:,:)
    call DHSEQR ('Schur', 'V', nl, ilo, ihi, uqt, nl, wre, wim, zzz, nl, work, lwork, info)
    call check_result(info, 'DHSEQR')
    ! Check for reality of eigenvalues.
    if (any(wim(:) /= 0.0d0)) then
       print *,' DHSEQR returns complex eigenvalues'
       print *,' program terminates in eigmod'
       stop 1
    endif

    ! ------------------------------------------------------
    ! Routine computes selected left and/or right eigen-
    ! vectors of a real upper quasi-triangular matrix
    ! We require ALL left AND right eigenvectors
    ! EVECL and EVECR need to be set on entry to the
    ! matrix of Schur vectors returned by DHSEQR,
    ! which is not necessarily the same as Q
    evecl(:,:) = zzz(:,:)
    evecr(:,:) = zzz(:,:)
    mwant = nl
    call DTREVC ('Both', 'Back', select, nl, uqt, nl, evecl, nl, evecr, nl, mwant, mgot, work, info)
    call check_result(info, 'DTREVC')
    ! Left and right eigenvectors are stored in the columns of
    ! EVECL, EVECR, i.e. 1st subscript = eigenvector component,
    ! 2nd subscript = eigenmode number

    ! ------------------------------------------------------
    ! Transform the eigenvectors of a balanced matrix to
    ! those of the original real nonsymmetric matrix A
    call DGEBAK (baljob, 'Left', nl, ilo, ihi, scale, mgot, evecl, nl, info)
    call DGEBAK (baljob, 'Righ', nl, ilo, ihi, scale, mgot, evecr, nl, info)
    ! ------------------------------------------------------

    ! At this point we should have a full set of eigen-
    ! values and eigenvectors of the original matrix A

    ! Check orthogonality of eigenvectors
    ! -----------------------------------
    do n=1,nl
       do m=1,nl
          elder(m,n) = sum(evecl(:,m)*evecr(:,n))
       enddo
    enddo

    ! Optionally print eigenvalues and vectors as a check
    ! ---------------------------------------------------
    if ( prtval ) then
       print *
       print '(2x,a)', 'A matrix:'
       do i=1,nl
          print '(2x,i3,1x,1p,9d16.7)', i,(mod%amat(i,j),j=1,nl)
       enddo
       print '(2x,a,i9,8i16)', ' m:  ',(m,m=1,nl)
       print '(2x,a,1p,9d16.7)', ' Wre',(wre(m),m=1,nl)
       print '(2x,a,1p,9d16.7)', ' Wim',(wim(m),m=1,nl)
       print '(2x,a)', 'Right eigenvectors:'
       do k=1,nl
          print '(2x,i3,1x,1p,9d16.7)', k,(evecr(k,m),m=1,nl)
       enddo
       print '(2x,a)', 'Left  eigenvectors:'
       do k=1,nl
          print '(2x,i3,1x,1p,9d16.7)', k,(evecl(k,m),m=1,nl)
       enddo
       print '(2x,a)', 'Dot products evecl(m).evecr(n):'
       print '(2x,a,i9,8i16)', ' n:  ',(n,n=1,nl)
       do m=1,nl
          print '(2x,i3,1x,1p,9d16.7)', m,(elder(m,n),n=1,nl)
       enddo
    endif

    end subroutine

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

    ! *     Check matrix triple product occuring in cyclic constraint equations
    ! *     -------------------------------------------------------------------
    ! *     Matrix is Cl2m*A*Cm2l
    ! **    do j=1,nl
    ! *       Work out A*Cm2l for current j
    ! **      do k=1,nl
    ! **        ac(k) = 0.0d0
    ! **        do m=1,nl
    ! **          ac(k) = ac(k) + mod%amat(k,m)*cm2l(m,j)
    ! **        enddo
    ! **      enddo
    ! *       Work out Cl2m*(A*Cm2l) for current i
    ! **      do i=1,nl
    ! **        ca = 0.0d0
    ! **        do k=1,nl
    ! **          ca = ca + cl2m(i,k)*ac(k)
    ! **        enddo
    ! **        ccprod(i,j) = ca
    ! **      enddo
    ! **    enddo
    ! **    print *
    ! **    print '(2x,a,i9,8i16)', ' j:  ',(j,j=1,nl)
    ! **    print '(2x,a,1p,9d16.7)', 'eig_val',(eig_val(j),j=1,nl)
    ! **    print '(2x,a,a)', 'matrix triple product ',
    ! **   &                     'cl2m(i,k)*A(k,m)*cm2l(m,j):'
    ! **    do i=1,nl
    ! **      print '(2x,i3,1x,1p,9d16.7)', i,(ccprod(i,j),j=1,nl)
    ! **    enddo
    ! *     Confirms that the matrix triple product is a diagonal
    ! *     matrix whose entries are the eigenvalues Lambda

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

  subroutine check_result(info, routine)

    integer, intent(in) :: info
    character, intent(in) :: routine*(*)

    if (info /= 0) then
       print *,' Problem in ', routine, ', INFO = ', info
       print *,' program terminates in eigmod'
       stop 1
    endif

  end subroutine check_result

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

    !do j=1,b%nyp
    !   modes_to_layers(:,j,:) = matmul(pm(:,j,:), mod%ctm2l(:,:))
    !enddo
    do j=1,b%nyp
       do i=1,b%nxp
          do k=1,b%nl
             modes_to_layers(i,j,k) = sum(pm(i,j,:)*mod%ctm2l(:,k))
          enddo
       enddo
    enddo


  end function modes_to_layers

end module modes
