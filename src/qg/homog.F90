module homog

  use box, only: box_type
  use modes, only: modes_type
  use intsubs, only: xintp, trapin
  use numerics, only: int_P_dA
  use constraint, only: core_constr_type, constraint_type
  use constraint, only: step_constraints
  use linalg, only: LU_factor, solve_Ax_b

  use inhomog, only: inhomog_type, solve_inhomog_eqn

  implicit none

  private

  type homog_box_type

     double precision, allocatable :: hom_sol_bc(:,:,:)
     !*     hom_sol_bc contains the homogeneous baroclinic
     !*     solutions, computed in the main program

     double precision, allocatable :: cdhoc(:,:),LU(:,:)
     integer, allocatable :: ipiv(:)
     double precision, allocatable :: cdiffo(:,:)
     !*     cdiffo and cdhoc contain (constant) matrices used in the
     !*     mass constraint equation applied in subroutine ocinvq
     !*     cdiffo(m,k) contains the coefficient of the mode m
     !*     contribution to the interface k constraint equation
     !*
     !*     LU contains the LU factorisation of cdhoc
     !*     ipiv contains the pivot indices for LU

  end type homog_box_type

  type homog_cyclic_type

     double precision, allocatable :: hom_sol_bc1(:,:), hom_sol_bc2(:,:)
     double precision, allocatable :: hom_sol_bt(:)

     double precision, allocatable :: hc1s(:), hc2s(:), hc1n(:), hc2n(:)

     double precision :: hbsi, int_sol_bt
     double precision, allocatable :: int_sol_bc(:)

     ! *     hom_sol_bc1, hom_sol_bc2 are the two homogeneous solutions for each baroclinic
     ! *     mode, computed in subroutine homsol. They are functions of y only
     ! *     (tabulated at p points). int_sol_bc is the area integral of these
     ! *     solutions for each mode (both modes have the same area integral)
     ! *
     ! *     hom_sol_bt is the homogeneous barotropic mode solution, computed in
     ! *     subroutine homsol. There is only one barotropic homogeneous
     ! *     mode, which is a function of y only. int_sol_bt is its area integral
     ! *
     ! *     hc1sat, hc2sat, hc1nat, hc2nat are boundary integrals of the
     ! *     homogeneous baroclinic solutions above, which arise in the momentum
     ! *     constraint equations in a channel model
     ! *
     ! *     hbsiat is the inverse of a simpler equivalent of the above
     ! *     quantities which occurs for the barotropic mode

  end type homog_cyclic_type

  type homog_type

     type(homog_cyclic_type) :: cyc
     type(homog_box_type) :: box
     
  end type homog_type

  
  public homog_type
  public init_homog
  public homsol

  public cyclic_homog
  public box_homog

contains

  type(homog_type) function init_homog(b)

    implicit none
    type(box_type), intent(in) :: b

    if (b%cyclic) then
       allocate(init_homog%cyc%hom_sol_bc1(b%nyp,2:b%nl))
       allocate(init_homog%cyc%hom_sol_bc2(b%nyp,2:b%nl))

       allocate(init_homog%cyc%hom_sol_bt(b%nyp))

       allocate(init_homog%cyc%hc1s(2:b%nl))
       allocate(init_homog%cyc%hc2s(2:b%nl))
       allocate(init_homog%cyc%hc1n(2:b%nl))
       allocate(init_homog%cyc%hc2n(2:b%nl))

       allocate(init_homog%cyc%int_sol_bc(2:b%nl))
    else
       allocate(init_homog%box%hom_sol_bc(b%nxp, b%nyp, 2:b%nl))

       allocate(init_homog%box%cdhoc(b%nl-1,b%nl-1))
       allocate(init_homog%box%LU(b%nl-1,b%nl-1))
       allocate(init_homog%box%ipiv(b%nl-1))
       allocate(init_homog%box%cdiffo(b%nl,b%nl-1))
    endif

    ! Compute tridiagonal coefficients for Helmholtz equations
    ! ========================================================
    ! Compute Del-sqd part once and for all; add def. rad. later.
    ! In all cases the tridiagonal system is diagonally dominant and
    ! thus numerically well-behaved. The Del-sqd operator to be
    ! inverted is the finite-difference form in both x- and y-directions.
    ! Ordering of wavenumber components in bd2oc,
    ! bd2at is appropriate to FFTPACK.

  end function init_homog

  subroutine homsol (b, mod, hom, inhom)

    ! Computes homogeneous baroclinic modal solutions for both ocean
    ! and atmosphere. These are stored in common blocks ochomog.cmn
    ! and athomog.cmn respectively, and used in subroutines ocinvq
    ! and atinvq to satisfy mass and momentum constraints.

    type(box_type), intent(in) :: b
    type(modes_type), intent(in) :: mod
    type(homog_type), intent(inout) :: hom
    type(inhomog_type), intent(inout) :: inhom

    if (b%cyclic) then
       call homsol_cyclic(b, hom%cyc, inhom)
    else
       call homsol_box(b, mod, hom%box, inhom)
    endif

  end subroutine homsol

  subroutine homsol_cyclic(b, hom_cyc, inhom)

    type(box_type), intent(in) :: b
    type(homog_cyclic_type), intent(inout) :: hom_cyc
    type(inhomog_type), intent(inout) :: inhom

    double precision :: wk1(b%nxp,b%nyp),wk2(b%nxp,b%nyp)
    double precision :: rhs1(b%nxp,b%nyp),rhs2(b%nxp,b%nyp)
    double precision pch1yn,pch2yn,pch1ys,pch2ys, &
         pchdet,aipch1,aipch2

    integer :: i, j, m
    ! Compute homogeneous channel solutions
    ! =================================================
    ! Barotropic mode: homogeneous solutions would need
    ! to satisfy d2p/dy2 = 0, but the coefficient of
    ! integral of p vanishes, so only the homogeneous
    ! component with a derivative can be inferred.
    ! W.L.O.G. we choose this to be 1 on the Southern
    ! boundary, and to vanish on the Northern boundary.
    ! We then apply the constraint at the Southern boundary
    do j=1,b%nyp
       hom_cyc%hom_sol_bt(j) = dble(b%nyp-j)/dble(b%nyp-1)
    enddo
    ! hbs = integral of y-derivative along Southern bdy
    hom_cyc%hbsi = b%yl/b%xl
    hom_cyc%int_sol_bt = 0.5d0*b%xl*b%yl
    print *,' '
    print *, ' Homogeneous barotropic solution:'
    print 240, '  int_sol_bt    = ',hom_cyc%int_sol_bt
    print 240, '  hbsi          = ',hom_cyc%hbsi
    
    ! Now compute baroclinic solutions for use in atinvq
    ! Solutions are functions of y only, but use the usual
    ! 2-D Helmholtz solver for simplicity and consistency

    print *
    print *, ' Homogeneous (baroclinic) solutions:'
    do m=2,b%nl
       print *
       print '(a,i2)', '  Mode: ',m

       ! pch = L(y) + rdm2*sol0, where L(y) is linear in y,
       ! and sol0 satisfies Del-sqd(sol0) - rdm2*sol0 = L(y)
       ! with the usual solid boundary condition p = 0.
       ! For hom_sol_bc1, L(y) = 1.0 on S bdy; = 0.0 on N bdy
       ! For hom_sol_bc2, L(y) = 0.0 on S bdy; = 1.0 on N bdy
       ! Specify RHSs
       do j=1,b%nyp
          hom_cyc%hom_sol_bc1(j,m) = ( b%yp(b%nyp)-b%yp(j) )/b%yl
          hom_cyc%hom_sol_bc2(j,m) = ( b%yp(  j )-b%yp(1) )/b%yl
          rhs1(:,j) = hom_cyc%hom_sol_bc1(j,m)
          rhs2(:,j) = hom_cyc%hom_sol_bc2(j,m)
       enddo

       ! Invert these RHSs for baroclinic homog. solutions (sol0 above)
       wk1(:,:) = solve_inhomog_eqn(inhom, b, m, rhs1)
       wk2(:,:) = solve_inhomog_eqn(inhom, b, m, rhs2)
       ! Add Helmholtz solution to L(y) to get full solutions
       ! Solutions in wk1, wk2 are functions of y only, i.e.
       ! independent of i, so just save solution for one i value
       do j=1,b%nyp
          do i=1,b%nxp
             wk1(i,j) = hom_cyc%hom_sol_bc1(j,m) + inhom%rdm2(m)*wk1(i,j)
             wk2(i,j) = hom_cyc%hom_sol_bc2(j,m) + inhom%rdm2(m)*wk2(i,j)
          enddo
          hom_cyc%hom_sol_bc1(j,m) = wk1(1,j)
          hom_cyc%hom_sol_bc2(j,m) = wk2(1,j)
       enddo
       ! Compute area integrals of hom_sol_bc1 and hom_sol_bc2
       aipch1 = xintp(wk1, b%nxp, b%nyp)
       aipch2 = xintp(wk2, b%nxp, b%nyp)
       ! Both solutions should have the same area integral
       hom_cyc%int_sol_bc(m) = 0.5d0*(aipch1+aipch2)*b%dx*b%dy
       
       ! Compute dp/dy half a gridpoint in from the north
       ! and south boundaries, and "integrate" in x
       ! Since these solutions are independent of x,
       ! x integration means just multiply by xla
       pch1ys = ( hom_cyc%hom_sol_bc1(  2 ,m) - hom_cyc%hom_sol_bc1(   1  ,m) )/b%dy
       pch2ys = ( hom_cyc%hom_sol_bc2(  2 ,m) - hom_cyc%hom_sol_bc2(   1  ,m) )/b%dy
       pch1yn = ( hom_cyc%hom_sol_bc1(b%nyp,m) - hom_cyc%hom_sol_bc1(b%nyp-1,m) )/b%dy
       pch2yn = ( hom_cyc%hom_sol_bc2(b%nyp,m) - hom_cyc%hom_sol_bc2(b%nyp-1,m) )/b%dy
       ! Correction for baroclinic modes
       pch1ys = -pch1ys + 0.5d0*b%dy*inhom%rdm2(m)*hom_cyc%hom_sol_bc1(  1 ,m)
       pch2ys = -pch2ys + 0.5d0*b%dy*inhom%rdm2(m)*hom_cyc%hom_sol_bc2(  1 ,m)
       pch1yn =  pch1yn + 0.5d0*b%dy*inhom%rdm2(m)*hom_cyc%hom_sol_bc1(b%nyp,m)
       pch2yn =  pch2yn + 0.5d0*b%dy*inhom%rdm2(m)*hom_cyc%hom_sol_bc2(b%nyp,m)
       ! Convert to line integrals
       pch1ys = b%xl*pch1ys
       pch2ys = b%xl*pch2ys
       pch1yn = b%xl*pch1yn
       pch2yn = b%xl*pch2yn
       ! The above are (for each mode m) the quantities in square
       ! brackets on the RHS of (B.14) and (B.15)
       pchdet = pch1ys*pch2yn - pch2ys*pch1yn
       hom_cyc%hc1s(m) = pch1ys/pchdet
       hom_cyc%hc2s(m) = pch2ys/pchdet
       hom_cyc%hc1n(m) = pch1yn/pchdet
       hom_cyc%hc2n(m) = pch2yn/pchdet
       print *
       print 240, '  aipch1, aipch2 = ',aipch1,aipch2
       print 240, '  int_sol_bc     = ',hom_cyc%int_sol_bc(m)
       print 240, '  pch1ys, pch1yn = ',pch1ys,pch1yn
       print 240, '  pch2ys, pch2yn = ',pch2ys,pch2yn
       print 240, '  pchdet         = ',pchdet
       print 240, '  hc1s, hc2s = ',hom_cyc%hc1s(m),hom_cyc%hc2s(m)
       print 240, '  hc1n, hc2n = ',hom_cyc%hc1n(m),hom_cyc%hc2n(m)
    enddo

240 format(a,1p,9d21.13)

  end subroutine homsol_cyclic

  subroutine homsol_box(b, mod, hom_box, inhom)

    type(box_type), intent(in) :: b
    type(modes_type), intent(in) :: mod
    type(homog_box_type), intent(inout) :: hom_box
    type(inhomog_type), intent(inout) :: inhom

    double precision :: rhs(b%nxp,b%nyp)
    integer :: k, m

    ! the area integrals of the homogeneous baroclinic solutions
    double precision :: aipohs(2:b%nl)

    ! Finite box ocean
    ! ----------------
    print *
    print *, ' Homogeneous (baroclinic) solutions:'
    do m=2,b%nl
       print *
       print '(a,i2)', '  Mode: ',m

       ! Compute new homogeneous solution = (1 + rdm2*sol0)
       ! sol0 satisfies Del-sqd(sol0) - rdm2*sol0 = 1
       ! with the usual boundary condition p = 0
       ! These are baroclinic solutions for use in ocinvq
       ! Setup RHS.
       rhs(:,:) = 1.0d0

       ! Solve for sol0 in ochom.
       hom_box%hom_sol_bc(:,:,m) = solve_inhomog_eqn(inhom, b, m, rhs(:,:))

       ! Add constant offset
       hom_box%hom_sol_bc(:,:,m) = 1.0d0 + inhom%rdm2(m)*hom_box%hom_sol_bc(:,:,m)
       ! Area integral of full homogeneous solution
       aipohs(m) = xintp(hom_box%hom_sol_bc(:,:,m), b%nxp, b%nyp)
       aipohs(m) = aipohs(m)*b%dx*b%dy
       print 240, '  aipohs = ',aipohs(m)
    enddo
    ! Compute the matrices used in the mass constraint equation
    ! dpioc(k) = Area integral of pressure diff ( po(k+1) - po(k) )
    ! Choose sign of dpioc so that +ve dpioc(k) -> +ve eta(k)
    ! cdiffo(m,k) is  coefficient which multiplies the mode
    ! m amplitude to give its contribution to dpioc(k)
    ! cdhoc(k,m) is  coefficient which multiplies a homogeneous
    ! baroclinic mode coefficient to give its contribution to dpioc(k)
    do k=1,b%nl-1
       do m=1,b%nl
          hom_box%cdiffo(m,k) = mod%ctm2l(m,k+1) - mod%ctm2l(m,k)
       enddo
       do m=1,b%nl-1
          hom_box%cdhoc(k,m) = ( mod%ctm2l(m+1,k+1) - mod%ctm2l(m+1,k) )*aipohs(m+1)
       enddo
    enddo
    ! Compute the LU factorization of cdhoc
    call LU_factor(hom_box%cdhoc, hom_box%LU, hom_box%ipiv)

    print *
    print *, ' Mass constraint matrices:'
    print *
    print *, ' cdiffo:'
    do k=1,b%nl-1
       print '(2x,i2,1x,1p,9d17.9)', k,(hom_box%cdiffo(m,k),m=1,b%nl)
    enddo
    print *
    print *, ' cdhoc:'
    do k=1,b%nl-1
       print '(2x,i2,1x,1p,9d17.9)', k,(hom_box%cdhoc(k,m),m=1,b%nl-1)
    enddo
    print *
    print *, ' LU:'
    do k=1,b%nl-1
       print '(2x,i2,1x,1p,9d17.9)', k,(hom_box%LU(k,m),m=1,b%nl-1)
    enddo
    print *
    print *, ' ipiv:'
    print '(4x,9i4)', (hom_box%ipiv(k),k=1,b%nl-1)

240 format(a,1p,9d21.13)

  end subroutine homsol_box

  subroutine cyclic_homog(b, tau_sign, tdt, constr, inhomog, &
       hom_cyc, mod, con, gp, &
       ent_xn, homcor)
    ! Solve constraint equations and add homogeneous solutions

    type(box_type), intent(in) :: b
    integer, intent(in) :: tau_sign
    double precision, intent(in) :: tdt
    type(constraint_type), intent(inout) :: constr
    double precision, intent(in) :: inhomog(b%nxp,b%nyp,b%nl)
    type(homog_cyclic_type), intent(in) :: hom_cyc
    type(modes_type), intent(in) :: mod
    type(core_constr_type), intent(inout) :: con
    double precision, intent(in) :: gp(b%nl-1)
    double precision, intent(in) :: ent_xn(b%nl-1)
    double precision, intent(out) :: homcor(b%nxp,b%nyp,b%nl)

    double precision :: aiplay(b%nl)
    double precision :: c_bc1(2:b%nl), c_bc2(2:b%nl), c_bt
    double precision :: homcor_2d(b%nyp,b%nl)

    ! Compute homogeneous solution coefficients
    ! Accumulate RHSs for the constraint equations for each layer
    ! Compute boundary integrals of entrainment
    call step_constraints(b, tau_sign, tdt, constr)

    ! Compute LHSs for the c_bc1, c_bc2 equations
    call compute_correction_coeffs(c_bc1, c_bc2, c_bt, inhomog, hom_cyc, &
         mod, constr%cs, constr%cn, b)

    call compute_layer_pressure(aiplay, b, inhomog, c_bc1, c_bc2, c_bt, hom_cyc, mod)

    call check_continuity(con, gp, b, tdt, aiplay, ent_xn)

    call homogeneous_correction(homcor_2d(:,:), hom_cyc, c_bc1, c_bc2, c_bt, b%nyp, b%nl)
    homcor(:,:,:) = spread(homcor_2d, 1, b%nxp)

  end subroutine cyclic_homog

  subroutine compute_correction_coeffs(c_bc1, c_bc2, c_bt, inhomog, hom_cyc, mod, cs, cn, b)
    type(box_type), intent(in) :: b
    double precision, intent(out) :: c_bc1(2:b%nl)
    double precision, intent(out) :: c_bc2(2:b%nl)
    double precision, intent(out) :: c_bt
    double precision, intent(in) :: inhomog(b%nxp,b%nyp,b%nl)
    type(homog_cyclic_type), intent(in) :: hom_cyc
    type(modes_type), intent(in) :: mod
    double precision, intent(in) :: cs(b%nl), cn(b%nl)

    integer :: m
    double precision :: clhss(b%nl)
    double precision :: clhsn(b%nl)
    double precision :: ayis, ayin

    do m=1,b%nl
       ! Compute line integrals of p_y for new modal solutions
       ! Integrate along south & north boundaries for all modes
       ! -point formulation, but values on bdy are exactly zero
       ayis = trapin(inhomog(:,2,m), b%nxp, b%dx)
       ayin = trapin(inhomog(:,b%nyp-1,m), b%nxp, b%dx)
       ! Compute LHSs for the c_bc1, c_bc2 equations
       clhss(m) = sum(mod%ctl2m(:,m)*cs(:)) + ayis/b%dy
       clhsn(m) = sum(mod%ctl2m(:,m)*cn(:)) + ayin/b%dy
    enddo

    ! Get coefft for barotropic mode
    c_bt = clhss(1)*hom_cyc%hbsi
    ! Derive c_bc1, c_bc2 for baroclinic modes
    do m=2,b%nl
       c_bc1(m) = hom_cyc%hc2n(m)*clhss(m) - hom_cyc%hc2s(m)*clhsn(m)
       c_bc2(m) = hom_cyc%hc1s(m)*clhsn(m) - hom_cyc%hc1n(m)*clhss(m)
    enddo
  end subroutine compute_correction_coeffs

  subroutine compute_layer_pressure(aiplay, b, inhomog, c_bc1, c_bc2, c_bt, hom_cyc, mod)

    type(box_type), intent(in) :: b
    double precision, intent(out) :: aiplay(b%nl)
    double precision, intent(in) ::  inhomog(b%nxp,b%nyp,b%nl)
    double precision, intent(in) :: c_bc1(2:b%nl)
    double precision, intent(in) :: c_bc2(2:b%nl)
    double precision, intent(in) :: c_bt

    type(homog_cyclic_type), intent(in) :: hom_cyc
    type(modes_type), intent(in) :: mod

    integer :: m, k
    double precision :: xinhom(b%nl), aipmod(b%nl)
    ! Compute area integrals of pressures
    ! -----------------------------------
    ! Integrals of modal pressures
    do m=1,b%nl
       ! Compute area integral of new inhomogeneous solution
       xinhom(m) = int_P_dA(inhomog(:,:,m), b)
    enddo

    aipmod(1) = xinhom(1) + c_bt*hom_cyc%int_sol_bt
    do m=2,b%nl
       aipmod(m) = xinhom(m) + ( c_bc1(m) + c_bc2(m) )*hom_cyc%int_sol_bc(m)
    enddo
    ! Integrals of layer pressures
    do k=1,b%nl
       aiplay(k) = sum(mod%ctm2l(:,k)*aipmod(:))
    enddo

  end subroutine compute_layer_pressure

  pure subroutine check_continuity(con, gp, b, tdt, aiplay, ent_xn)
    ! Check continuity is satisfied at each interface
    ! Update continuity measures at each interface
    type(core_constr_type), intent(inout) :: con
    type(box_type), intent(in) :: b
    double precision, intent(in) :: gp(b%nl-1)
    double precision, intent(in) :: tdt
    double precision, intent(in) :: aiplay(b%nl)
    double precision, intent(in) :: ent_xn(b%nl-1)

    integer :: k
    double precision :: est1, est2, edif, esum

    double precision ecrit
    parameter ( ecrit=1.0d-13 )

    do k=1,b%nl-1
       ! Choose sign of dpioc so that +ve dpioc -> +ve eta
       ! Check continuity is satisfied at each interface
       ! MONITORING - extra section for ermaso, emfroc
       ! Compute alternative estimates of new dpioc
       est1 = b%dz_sign*(aiplay(k) - aiplay(k+1))
       est2 = con%dpip(k) - tdt*gp(k)*ent_xn(k)
       edif = est1 - est2
       esum = abs(est1) + abs(est2)
       con%ermas(k) = edif
       ! Compute fractional error if entrainment is significant;
       ! fraction is meaningless if est1, est2 just noisy zeros
       if (esum > (ecrit*b%xl*b%yl*tdt*gp(k))) then
          con%emfr(k) = 2.0d0*edif/esum
       else
          con%emfr(k) = 0.0d0
       endif
       ! Update continuity constants
       con%dpip(k) = con%dpi(k)
       con%dpi(k) = b%dz_sign*(aiplay(k) - aiplay(k+1))
    enddo

  end subroutine check_continuity

  pure subroutine homogeneous_correction(homcor, hom_cyc, c_bc1, c_bc2, c_bt, nyp, nl)

    double precision, intent(out) :: homcor(nyp,nl)
    type(homog_cyclic_type), intent(in) :: hom_cyc
    double precision, intent(in) :: c_bc1(2:nl)
    double precision, intent(in) :: c_bc2(2:nl)
    double precision, intent(in) :: c_bt
    integer, intent(in) :: nyp, nl

    integer :: m

    ! Compute homogeneous corrections (indep. of i)
    ! Barotropic mode
    homcor(:,1) = c_bt*hom_cyc%hom_sol_bt(:)
    ! Baroclinic modes
    do m=2,nl
       homcor(:,m) = c_bc1(m)*hom_cyc%hom_sol_bc1(:,m) + c_bc2(m)*hom_cyc%hom_sol_bc2(:,m)
    enddo

  end subroutine homogeneous_correction

  subroutine box_homog(b, tdt, inhomog, hom_box, con, gp, ent_xn, homcor)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: tdt
    double precision, intent(in) :: inhomog(b%nxp,b%nyp,b%nl)
    type(homog_box_type), intent(in) :: hom_box
    type(core_constr_type), intent(inout) :: con
    double precision, intent(in) :: gp(b%nl-1)
    double precision, intent(in) :: ent_xn(b%nl-1)
    double precision, intent(out) :: homcor(b%nxp,b%nyp,b%nl)

    double precision :: aient(b%nl-1), xinhom(b%nl), rhs(b%nl-1)
    double precision :: hclco(b%nl-1), aitmp

    integer :: k, m

    ! Get multiples of homogeneous solutions to conserve thickness
    ! aient(k) = Area integral of oceanic entrainment e(k)
    aient(1) = ent_xn(1)
    ! N.B. xon(1) is now zero by construction of entoc in oml
    ! All other entrainments assumed exactly zero in the ocean
    do k=2,b%nl-1
       aient(k) = 0.0d0
    enddo
    ! Area integral of d(eta(k))/dt = - Area integral of entrainment e(k)
    do m=1,b%nl
       ! Compute area integral of new inhomogeneous solution
       xinhom(m) = int_P_dA(inhomog(:,:,m), b)
    enddo

    do k=1,b%nl-1
       aitmp = con%dpi(k)
       con%dpi(k) = con%dpip(k) - tdt*gp(k)*aient(k)
       con%dpip(k) = aitmp
       rhs(k) = con%dpi(k) - sum(hom_box%cdiffo(:,k)*xinhom(:))
    enddo

    ! Matrix equation is cdhoc*hclco = rhs
    ! Solve equation for homogeneous solution coeffts using LAPACK
    ! Solve the linear system using the LU factorised matrix
    call solve_Ax_b(hom_box%cdhoc, hom_box%LU, hom_box%ipiv, rhs, hclco)
         
    ! Barotropic mode               
    homcor(:,:,1) = 0.0d0
    ! Baroclinic modes including homogeneous contribution
    do m=2,b%nl
       homcor(:,:,m) = hclco(m-1)*hom_box%hom_sol_bc(:,:,m)
    enddo

  end subroutine box_homog

end module homog
