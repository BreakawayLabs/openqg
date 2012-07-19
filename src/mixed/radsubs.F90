module radsubs

  use ncutils, only: nc_open, nc_get_double
  use grid, only: grid_type
  use box, only: box_type
  use intsubs, only: trapin
  use mixed, only: temp_type, init_temp_from_rad

  implicit none

  private

  type radiate_type
     ! User defined values
     double precision :: fsbar, fspamp
     double precision, allocatable :: zopt(:)
     double precision :: gamma, xlamda

     double precision, allocatable :: tabsat(:)

     ! derived values
     double precision, allocatable :: Aup(:,:)
     double precision, allocatable :: Adown(:,:)
     double precision, allocatable :: Dup(:)
     double precision :: D0up
     double precision :: Dmdown

     double precision :: fspco
     double precision, allocatable :: rbetat(:)

     !      fsbar = background forcing amplitude (<0)       (W m^-2)
     !      fspamp = perturb.  forcing amplitude (>0)       (W m^-2)
     !      fspco = signed perturbation coefft (sign as f0) (W m^-2)
     !      zopt(0) = optical depth in atmos. mixed layer        (m)
     !      gamma = adiabatic lapse rate                    (K m^-1)
     !      zopt(k) = optical depth in atmos. layer k       (m)
     !      xlamda = (sensible + latent) transfer coefft.   (W m^-2 K^-1)

     !      Coefficients for internal QG interface etas:
     !      Aup(k,l) = layer k upward radiation coefft.     (W m^-3)
     !      (multiplying eta of l-th interface)
     !      Adown(k,l) = layer k downward radiation coefft. (W m^-3)
     !      (multiplying eta of l-th interface)
     !      Coefficients for m.l. temperature perturbation:
     !      D0up = o.m.l. upward radiation coefft.          (W m^-2 K^-1)
     !      Dup(k) = layer k upward radiation coefft.       (W m^-2 K^-1)
     !     Dmdown = a.m.l. downward radiation coefft.      (W m^-2 K^-1)
     !     rbetat(k) = multiplier of fs' for radiative balance
     !     initialisation of interface displacement eta(k)
     !

  end type radiate_type

  public radiate_type

  public radiat
  public load_rad
  public fsprim

contains

  type(radiate_type) function load_rad(ga, fsbar, fspamp, zopt, gamma, xlamda, T_abs)

    type(box_type), intent(in) :: ga

    double precision, intent(in) :: T_abs(ga%nl)
    double precision, intent(in) :: fsbar, fspamp
    double precision, intent(in) :: zopt(0:ga%nl)
    double precision, intent(in) :: gamma, xlamda

    allocate(load_rad%Aup(0:ga%nl,0:ga%nl-1))
    allocate(load_rad%Adown(ga%nl,0:ga%nl-1))
    allocate(load_rad%Dup(0:ga%nl))
    allocate(load_rad%rbetat(ga%nl-1))
    allocate(load_rad%zopt(0:ga%nl))
    allocate(load_rad%tabsat(ga%nl))

    load_rad%fsbar = fsbar
    load_rad%fspamp = fspamp
    load_rad%zopt = zopt
    load_rad%gamma = gamma
    load_rad%xlamda = xlamda 
    load_rad%tabsat = T_abs

    ! Derive signed perturbation coefficient
    if ( load_rad%fspamp.lt.0.0d0 ) then
       print *,' Sign of fspamp appears incorrect'
       print *,' Program terminates in radiat'
       stop 1
    endif
    load_rad%fspco = sign( load_rad%fspamp, ga%fnot)

  end function load_rad


  double precision elemental function fsprim (ga, fspco, yrel)
    
    ! Computes perturbative radiation forcing of system.
    ! fspco is peak-to-trough amplitude of perturbative forcing.
    ! fspco is > 0 for Northern hemisphere, < 0 for Southern.
    ! yrel is y relative to central latitude
    ! Function should be chosen to have zero integral
    ! over y = [0, yla], i.e. over yrel = [-yla/2, yla/2].
    ! Fs is +ve upwards, so function should be -ve near
    ! the equator, and become positive poleward.
        
    type(box_type), intent(in) :: ga
    double precision, intent(in) :: fspco
    double precision, intent(in) ::  yrel
    
    double precision :: PI
    parameter ( PI=3.14159265358979324D0 )

    fsprim = fspco*0.5d0*sin( PI*yrel/ga%yl )

  end function fsprim

  subroutine radiat(ga, rad, rbtmat, rbtmoc, ocn_tmbar, atm_tmbar)

    ! Computes the mean state radiative balance, including the
    ! atmosphere and ocean mixed layer temperatures that ensure
    ! equilibrium. Also computes the perturbation (linearised)
    ! radiation coefficients A, B, C and D, the entrainment
    ! coefficients derived from them, and the radiative balance
    ! initialisation coefficients.

    type(box_type), intent(in) :: ga
    type(radiate_type), intent(inout) :: rad
    double precision, intent(out) :: rbtmat, rbtmoc
    double precision, intent(out) :: ocn_tmbar
    double precision, intent(out) :: atm_tmbar

      
    double precision :: tauk(0:ga%nl), uprad(ga%nl),dnrad(ga%nl)
    double precision :: F0upbar, Fupbar(0:ga%nl),Fmdnbar,Fdnbar(ga%nl)
    double precision :: uprad_m

    call compute_layer_trans(ga, rad%zopt, tauk)

    call compute_mean_radiation(ga, rad%tabsat, rad%gamma, rad%zopt, uprad, dnrad)

    call compute_ml_mean_temp(ga, rad, uprad, tauk, uprad_m, atm_tmbar)

    ocn_tmbar = compute_oml_mean_temp(rad, atm_tmbar)

    call compute_mean_flux(ga, ocn_tmbar, atm_tmbar, uprad_m, uprad, &
         dnrad, tauk, F0upbar, Fupbar, Fmdnbar, Fdnbar)

    call write_mean_rad_params(ga, rad, atm_tmbar, &
         F0upbar, Fmdnbar, tauk, uprad, dnrad, Fupbar, Fdnbar)

    call compute_perturb_rad_params(ga, Fupbar, Fdnbar, uprad, &
         dnrad, tauk, ocn_tmbar, atm_tmbar, rad%tabsat, rad)

    call write_rad_coeffs(ga, rad)

    call compute_rad_bal_coeffs(ga, rbtmoc, rbtmat, rad)

  end subroutine radiat

  subroutine compute_layer_trans(ga, zopt, tauk)
    ! Compute transmissivities through each layer.
    ! The 1st layer is taken only down to the mixed layer, which
    ! has its transmissivity calculated separately
    
    type(box_type), intent(in) :: ga
    double precision, intent(in) :: zopt(0:ga%nl)
    double precision, intent(out) :: tauk(0:ga%nl)
          
    integer :: k

    tauk(0) = trans(ga%hm, zopt(0))
    tauk(1) = trans(ga%h(1)-ga%hm, zopt(1))
    do k=2,ga%nl
       tauk(k) = trans(ga%h(k), zopt(k))
    enddo

  end subroutine compute_layer_trans

  double precision function trans(thickness, opt_depth)
    ! Equation A.1. tau(z, Z -> trans(Z - z, z_0)

    double precision :: thickness, opt_depth
    trans = exp( -thickness/opt_depth )
    
  end function trans

  subroutine compute_mean_radiation(ga, tabsat, gamma, zopt, uprad, dnrad)
      
    type(box_type), intent(in) :: ga
    double precision, intent(in) :: tabsat(ga%nl)
    double precision, intent(in) :: gamma
    double precision, intent(in) :: zopt(0:ga%nl)
    double precision, intent(out) :: uprad(ga%nl), dnrad(ga%nl)

    integer :: nz 
    parameter ( nz=10001 )

    double precision :: delz, hbot, htop, upint, dnint, fdn(nz)
    double precision :: fup(nz), zz
    double precision :: T(nz), tau_up(nz), tau_down(nz)
    integer :: i, k

    double precision :: stefan,sigov2
    parameter ( stefan=5.67040D-8, sigov2=0.5d0*stefan )

    ! Computes sig*I_k/2*zk up and down as per A.26 and A.41
    
    ! Layer 1 up- and downgoing mean radiation
    hbot = ga%hm
    htop = ga%h(1)
    delz = ( htop - hbot )/dble(nz-1)
    do i=1,nz
       zz = hbot + dble(i-1)*delz
       T(i) = tabsat(1) - gamma*zz
       tau_up(i)   = trans(htop - zz, zopt(1))
       tau_down(i) = trans(zz - hbot, zopt(1))
    enddo

    fup(:) = sigov2*T(:)**4 * tau_up(:)
    fdn(:) = sigov2*T(:)**4 * tau_down(:)
    upint = trapin(fup, nz, delz)
    dnint = trapin(fdn, nz, delz)
    uprad(1) = upint/zopt(1)
    dnrad(1) = dnint/zopt(1)

    ! Upper layers up- and downgoing mean radiation
    do k=2,ga%nl
       hbot = htop
       htop = hbot + ga%h(k)
       delz = ga%h(k)/dble(nz-1)
       do i=1,nz
          zz = hbot + dble(i-1)*delz
          T(i) = tabsat(k) - gamma*zz
          tau_up(i) = trans(htop - zz, zopt(k))
          tau_down(i) = trans(zz - hbot, zopt(k))
       enddo

       fup(:) = sigov2*T(:)**4 * tau_up(:)
       fdn(:) = sigov2*T(:)**4 * tau_down(:)
       upint = trapin(fup, nz, delz)
       dnint = trapin(fdn, nz, delz)
       uprad(k) = upint/zopt(k)
       dnrad(k) = dnint/zopt(k)
    enddo

  end subroutine compute_mean_radiation

  subroutine compute_ml_mean_temp(ga, rad, uprad, tauk, uprad_m, atm_tmbar)

    type(box_type), intent(in) :: ga
    type(radiate_type), intent(in) :: rad
    double precision, intent(in) :: uprad(ga%nl)
    double precision, intent(in) :: tauk(0:ga%nl)
    double precision, intent(inout) :: uprad_m
    double precision, intent(out) :: atm_tmbar
      
    integer :: nz, nitmax
    parameter ( nz=10001, nitmax=200)
    double precision :: stefan, sigov2, tmbtol
    parameter ( stefan=5.67040D-8, sigov2=0.5d0*stefan,tmbtol=1.0d-13)

    double precision :: deltm, delz
    integer :: iter, i, k
    double precision :: fup(nz), T(nz), tau_up(nz), zz, rhstat

    rhstat = 0.0d0
    do k=1,ga%nl
       rhstat = rhstat*tauk(k) + uprad(k)
    enddo
    rhstat = (-rhstat - rad%fsbar)/product(tauk(1:))

    ! Compute atmosphere m. l. mean temperature
    ! rhstat now contains required value of
    ! integral for mixed layer upgoing radiation
    ! Adjust tmbara to get required value of integral
    atm_tmbar = 300.0d0
    deltm = 0.0d0
    delz = ga%hm/dble(nz-1)
    iter = 0
100 continue

    do i=1,nz
       zz = dble(i-1)*delz
       T(i) = atm_tmbar - rad%gamma*zz
       tau_up(i) = trans(ga%hm - zz, rad%zopt(0))
    enddo
    fup(:) = T(:)**4 * tau_up(:)
    uprad_m = trapin(fup, nz, delz)
    uprad_m = uprad_m*sigov2/rad%zopt(0)
    ! Convert error to (approximate) temperature change
    deltm = 0.25d0*(rhstat - uprad_m)*atm_tmbar/uprad_m
    atm_tmbar = atm_tmbar + 0.75d0*deltm
    iter = iter + 1

    if ( iter.gt.nitmax ) then
       print *,' iteration for tmbara not converged'
       print *,' max. no. of iterations nitmax = ',nitmax
       print *,' deltm, tmbara = ',deltm,atm_tmbar
       print *,' program terminates in radiat'
       stop 1
    endif
    if ( abs( deltm ).gt.tmbtol ) goto 100

  end subroutine compute_ml_mean_temp

  double precision function compute_oml_mean_temp(rad, atm_tmbar)

    type(radiate_type), intent(in) :: rad
    double precision, intent(in) :: atm_tmbar

    double precision :: rhstoc, tocold
    integer :: iter

    integer :: nitmax
    parameter ( nitmax=200)
    double precision :: stefan, tmbtol
    parameter ( stefan=5.67040D-8, tmbtol=1.0d-13)

    ! Compute ocean m. l. mean temperature
    rhstoc = rad%xlamda*atm_tmbar +  B(atm_tmbar, 0.0d0, 0.0d0) - rad%fsbar
    ! Use atmos. temp. as initial guess for ocean temp.
    compute_oml_mean_temp = atm_tmbar

    iter = 0
150 continue
    tocold = compute_oml_mean_temp
    compute_oml_mean_temp = rhstoc/( rad%xlamda + stefan*tocold**3 )
    iter = iter + 1

    if ( iter.gt.nitmax ) then
       print *,' iteration for tmbaro not converged'
       print *,' max. no. of iterations nitmax = ',nitmax
       print *,' tocold, tmbaro = ',tocold,compute_oml_mean_temp
       print *,' program terminates in radiat'
       stop 1
    endif
    if ( abs( compute_oml_mean_temp - tocold ).gt.tmbtol ) goto 150

  end function compute_oml_mean_temp

  double precision function B(T0, z, gamma)

    double precision :: T0, z, gamma
    double precision :: sigma = 5.67040D-8
    
    B = (sigma/2.0)*(T0 - z*gamma)**4

  end function B
  
  subroutine compute_mean_flux(ga, ocn_tmbar, atm_tmbar, uprad_m, uprad, dnrad, tauk, &
       F0upbar, Fupbar, Fmdnbar, Fdnbar)

    type(box_type), intent(in) :: ga
    double precision, intent(in) :: ocn_tmbar
    double precision, intent(in) :: atm_tmbar
    double precision, intent(in) :: uprad_m
    double precision, intent(in) :: uprad(ga%nl), dnrad(ga%nl)
    double precision, intent(in) :: tauk(0:ga%nl)
    double precision, intent(out) :: F0upbar
    double precision, intent(out) :: Fupbar(0:ga%nl)
    double precision, intent(out) :: Fmdnbar
    double precision, intent(out) :: Fdnbar(ga%nl)
    
    integer :: k

    double precision :: stefan, sigov2
    parameter ( stefan=5.67040D-8, sigov2=0.5d0*stefan )

    ! Upgoing mean state fluxes
    ! Ocean
    F0upbar = stefan*ocn_tmbar**4 ! A.10
    ! Atmos. mixed layer
    Fupbar(0) = uprad_m       ! A.16
    ! Upper layers
    do k=1,ga%nl
       Fupbar(k) = Fupbar(k-1)*tauk(k) + uprad(k) ! A.28
    enddo

    ! Downgoing mean state fluxes
    ! Top layer
    Fdnbar(ga%nl) = -dnrad(ga%nl) ! A.43
    ! Other QG layers
    do k=ga%nl-1,1,-1
       Fdnbar(k) = Fdnbar(k+1)*tauk(k) - dnrad(k) ! A.43
    enddo
    ! Atmos. mixed layer
    Fmdnbar = -sigov2*atm_tmbar**4 ! A.22

  end subroutine compute_mean_flux

  subroutine write_mean_rad_params(ga, rad, atm_tmbar, &
       F0upbar, Fmdnbar, tauk, uprad, dnrad, Fupbar, Fdnbar)

    type(box_type), intent(in) :: ga
    type(radiate_type), intent(in) :: rad
    double precision, intent(in) :: atm_tmbar
    double precision, intent(in) :: F0upbar, Fmdnbar
    double precision, intent(in) :: tauk(0:ga%nl), uprad(ga%nl)
    double precision, intent(in) :: dnrad(ga%nl), Fupbar(0:ga%nl), Fdnbar(ga%nl)

    integer :: k

    double precision :: tmbtol
    parameter ( tmbtol=1.0d-13 )

    print *,' '
    print *,' Mean state radiation parameters:'
    print *,' --------------------------------'
    print 205, '  Mean forcing fsbar (W m^-2) = ',rad%fsbar
    print 205, '  Pert. ampl. fspamp (W m^-2) = ',rad%fspamp
    print 205, '  Pert. coefft fspco (W m^-2) = ',rad%fspco
    print 204, '  Optic. depth in aml. zopt(0) (m) = ',rad%zopt(0)
    do k=1,ga%nl
       print 224, '  Optic. depth (m) in layer',k,' = ',rad%zopt(k)
    enddo
    print 214, '  Lapse rate gamma   (K m^-1) = ',rad%gamma
    print *,' '
    print 209, '  A.m.l. transmissivity  taum = ',tauk(0)
    do k=1,ga%nl
       print 229, '  Layer', k ,'  transmissivity tau = ', tauk(k)
    enddo
    print 209, '  Transmissivity prod. tupmul = ',product(tauk(1:))
    do k=1,ga%nl
       print 228, '  Layer', k ,' integs uprad, dnrad = ', uprad(k),dnrad(k)
    enddo
    print *,' '
    print 214, '  Tolerance for m.l. temp (K) = ',tmbtol
    print 207, '  Atmos. mixed layer T (K, C) = ', atm_tmbar,atm_tmbar-273.15d0
    print *,' '
    print 206, '  Ocean m.l. radiat.  F0upbar = ',F0upbar
    print 206, '  Mixed layer  Fmbar up, down = ',Fupbar(0),Fmdnbar
    do k=1,ga%nl
       print 226, '  Layer', k ,'       Fbar up, down = ', Fupbar(k),Fdnbar(k)
    enddo
    print 214, '  Fractional error in OLR     = ', abs( Fupbar(ga%nl) + rad%fsbar )/abs(rad%fsbar)

204 format(a,9f13.4)
205 format(a,9f13.5)
206 format(a,9f13.6)
207 format(a,9f13.7)
209 format(a,9f13.9)
214 format(a,1p,9d13.4)
224 format(a,i2,a,9f13.4)
226 format(a,i2,a,9f13.6)
228 format(a,i2,a,9f13.8)
229 format(a,i2,a,9f13.9)
      
  end subroutine write_mean_rad_params

  subroutine compute_perturb_rad_params(ga, Fupbar, Fdnbar, uprad, dnrad, tauk, ocn_tmbar, atm_tmbar, tabsat, rad)

    type(box_type), intent(in) :: ga
    double precision, intent(in) :: Fupbar(0:ga%nl)
    double precision, intent(in) :: Fdnbar(ga%nl)
    double precision, intent(in) :: uprad(ga%nl)
    double precision, intent(in) :: dnrad(ga%nl)
    double precision, intent(in) :: tauk(0:ga%nl)
    double precision, intent(in) :: ocn_tmbar
    double precision, intent(in) :: atm_tmbar
    double precision, intent(in) :: tabsat(ga%nl)
    type(radiate_type), intent(inout) :: rad

    integer :: nz
    parameter ( nz=10001 )

    double precision :: stefan
    ! stefan is the Stefan-Boltzmann constant (usually denoted by sigma)
    parameter ( stefan=5.67040D-8 )

    double precision :: delz, upint, zz
    double precision :: fup(nz)
    integer :: i,k,l
    double precision :: Z(ga%nl), Btop(0:ga%nl), Bbot(0:ga%nl) 

    ! Fixed heights of top of each layer
    Z(1) = ga%h(1)
    do k=2,ga%nl
       Z(k) = Z(k-1) + ga%h(k)
    enddo

    ! Boltzmann radiation at the top of each layer (where 0 = mixed layer)...
    Btop(0) = B(atm_tmbar, ga%hm, rad%gamma)
    do k=1,ga%nl
       Btop(k) = B(tabsat(k), Z(k), rad%gamma)
    enddo
    ! ... and at the bottom.
    Bbot(0) = B(atm_tmbar, 0.0d0, rad%gamma)
    Bbot(1) = B(tabsat(1), ga%hm, rad%gamma)
    do k=2, ga%nl
       Bbot(k) = B(tabsat(k), Z(k-1), rad%gamma)
    enddo
    
    ! Perturbation (linearised) radiation parameters
    ! ==============================================
    ! Initialise Aup and Adown arrays
    rad%Aup(:,:) = 0.0d0
    rad%Adown(:,:) = 0.0d0

    ! Upgoing
    rad%Aup(0,0) = (Btop(0) - Fupbar(0) )/rad%zopt(0)
    do k=1,ga%nl-1
       ! Top layer has no eta coefft at top
       rad%Aup(k,k) = ( -tauk(k)*Fupbar(k-1) - uprad(k) + Btop(k) )/rad%zopt(k)
    enddo
    do k=1,ga%nl
       rad%Aup(k,k-1) = tauk(k)*( rad%Aup(k-1,k-1) + Fupbar(k-1)/rad%zopt(k) - Bbot(k)/rad%zopt(k))
    enddo
    do k=2,ga%nl
       do l=0,k-2
          rad%Aup(k,l) = rad%Aup(k-1,l)*tauk(k)
       enddo
    enddo

    ! Upgoing
    ! Ocean mixed layer
    rad%D0up = 4.0d0*stefan*ocn_tmbar**3
    ! Atmospheric mixed layer
    delz = ga%hm/dble(nz-1)
    do i=1,nz
       zz = dble(i-1)*delz
       fup(i) = (atm_tmbar-rad%gamma*zz)**3*exp(-(ga%hm - zz)/rad%zopt(0))
    enddo
    upint = trapin(fup, nz, delz)
    rad%Dup(0) = 2.0d0*stefan*upint/rad%zopt(0)
    ! Layer 1
    rad%Dup(1) = rad%Dup(0)*tauk(1)
    ! Upper layers
    do k=2,ga%nl
       rad%Dup(k) = rad%Dup(k-1)*tauk(k)
    enddo

    ! Downgoing
    do k=1,ga%nl-1
       rad%Adown(k,k-1) = ( Fdnbar(k+1)*tauk(k) - dnrad(k) + Bbot(k) )/rad%zopt(k)
    enddo
    rad%Adown(ga%nl,ga%nl-1)=(Bbot(ga%nl) - dnrad(ga%nl) )/rad%zopt(ga%nl)

    ! Intermediate layers
    do k=ga%nl-1,1,-1
       do l=k+1,ga%nl-1
          rad%Adown(k,l) = rad%Adown(k+1,l)*tauk(k)
       enddo
       rad%Adown(k, k ) = tauk(k)*( rad%Adown(k+1,k) - Fdnbar(k+1)/rad%zopt(k) - Btop(k)/rad%zopt(k))
    enddo

    ! Atmospheric mixed layer
    rad%Dmdown = - 2.0d0*stefan*atm_tmbar**3

  end subroutine compute_perturb_rad_params

  subroutine write_rad_coeffs(ga, rad)

    type(box_type), intent(in) :: ga
    type(radiate_type), intent(in) :: rad

    integer :: k,l

    print *,' '
    print *,' Linearised radiation coeffts:'
    print *,' -----------------------------'
    print *,' QG internal interface eta coeffts:'
    do k=1,ga%nl
       print 235, '  Layer', k ,' coeffts    Aup(k,l) = ', (rad%Aup(k,l),l=1,ga%nl-1)
    enddo
    do k=1,ga%nl
       print 235, '  Layer', k ,' coeffts  Adown(k,l) = ', (rad%Adown(k,l),l=1,ga%nl-1)
    enddo
    print *,' Mixed layer eta coeffts:'
    print 209, '  Mixed layer coefft     Bmup = ',rad%Aup(0,0)
    print 209, '  Layer 1 coefft       B1down = ',rad%Adown(1,0)
    do k=1,ga%nl
       print 229, '  Layer', k ,' coefft          Bup = ', rad%Aup(k,0)
    enddo
    print *,' Temperature perturbation coeffts:'
    print 209, '  Oceanic coefft         D0up = ',rad%D0up
    print 209, '  Mixed layer coefft     Dmup = ',rad%Dup(0)
    print 209, '  Mixed layer coefft   Dmdown = ',rad%Dmdown
    do k=1,ga%nl
       print 229, '  Layer', k ,' coefft          Dup = ',rad%Dup(k)
    enddo

209 format(a,9f13.9)
229 format(a,i2,a,9f13.9)
235 format(a,i2,a,1p,9d13.5)

  end subroutine write_rad_coeffs

  subroutine compute_rad_bal_coeffs(ga, rbtmoc, rbtmat, rad)

    type(box_type), intent(in) :: ga
    type(radiate_type), intent(inout) :: rad
    double precision, intent(inout) :: rbtmoc, rbtmat

    integer :: ipivrb(ga%nl),iwork(ga%nl)
    double precision :: rbalar(ga%nl,ga%nl),rballu(ga%nl,ga%nl)
    double precision :: balrhs(ga%nl)
    double precision :: berr, ferr
    integer :: info
    double precision :: rbafac(ga%nl),work(3*ga%nl)
    
    integer :: i, k

    ! Radiation balance initialisation coefficients
    ! =============================================
    ! We assume the mixed layer interface displacement etam
    ! and atmospheric topography are both zero, and
    ! that the entrainments between QG layers are all zero
    ! The first ga%nl-1 entries in the solution vector are the
    ! coefficients of the internal interface displacements;
    ! the last entry is the coefficient of Tm'
    ! Zero the matrices
    rbalar(:,:) = 0.0d0
    rballu(:,:) = 0.0d0
    ! Compute the matrices from the radiation coefficients
    ! 1st equation is balance for layer 1
    do i=1,ga%nl-1
       rbalar(1,i) = rad%Adown(1,i)
    enddo
    rbalar(1,ga%nl) = rad%Dup(0)
    ! Option 1: only e(1) is nonzero
    do k=2,ga%nl-1
       do i=1,ga%nl-1
          rbalar(k,i) = rad%Adown(k+1,i) + rad%Aup(k,i)
       enddo
       rbalar(k,ga%nl) = rad%Dup(k)
    enddo
    ! Final equation is OLR at top of atmosphere
    do i=1,ga%nl-1
       rbalar(ga%nl,i) = rad%Aup(ga%nl,i)
    enddo
    rbalar(ga%nl,ga%nl) = rad%Dup(ga%nl)
    ! Setup RHS & make copy for LU factorisation
    balrhs(:) = -1.0d0
    rbafac(:) = balrhs(:)
    rballu(:,:) = rbalar(:,:)

    ! Compute the LU factorization of rbalar
    ! DGETRF = NAG routine F07ADF
    call DGETRF (ga%nl, ga%nl, rballu, ga%nl, ipivrb, info)
    if ( info.ne.0 ) then
       print *,'  DGETRF in radiat returns info = ',info
       print *,'  program terminates in radiat'
       stop 1
    endif
    print *,' '

    ! Solve matrix equation for radiative balance coeffts using LAPACK
    ! Solve the linear system using the LU factorised matrix rballu
    ! DGETRS = NAG routine F07AEF
    call DGETRS ('Norm', ga%nl, 1, rballu, ga%nl, ipivrb, rbafac, ga%nl, info)
    if ( info.ne.0 ) then
       print *,'  DGETRS in radiat returns info = ',info
       print *,'  program terminates in radiat'
       stop 1
    endif
    ! Improve the solution in rbafac by iterative refinement
    ! DGERFS = NAG routine F07AHF
    call DGERFS ('Norm', ga%nl, 1, rbalar, ga%nl, rballu, ga%nl, &
         ipivrb, balrhs, ga%nl, rbafac, ga%nl, ferr, berr, &
         work, iwork, info)
    if ( info.ne.0 ) then
       print *,'  DGERFS in radiat returns info = ',info
       print *,'  program terminates in radiat'
       stop 1
    endif

    do k=1,ga%nl-1
       rad%rbetat(k) = rbafac(k)
    enddo
    rbtmat = rbafac(ga%nl)
    rbtmoc = ((rad%xlamda - rad%Dmdown)*rbtmat - 1.0d0)/(rad%xlamda + rad%D0up)

    print *, ' Radiation balance matrices:'
    print *
    print *, ' rbalar:'
    do k=1,ga%nl
       print '(2x,i2,1x,1p,5d17.9)', k,(rbalar(k,i),i=1,ga%nl)
    enddo
    print *
    print *, ' rballu:'
    do k=1,ga%nl
       print '(2x,i2,1x,1p,5d17.9)', k,(rballu(k,i),i=1,ga%nl)
    enddo
    print *
    print *,' Radiative balance initialisation coeffts:'
    print *,' -----------------------------------------'
    print *,' QG internal interface eta coeffts:'
    do k=1,ga%nl-1
       print 228, '  Layer', k ,' coefft       rbetat = ',rad%rbetat(k)
    enddo
    print 208, '  Atmos. m.l. T coefft rbtmat = ',rbtmat
    print 208, '  Ocean  m.l. T coefft rbtmoc = ',rbtmoc

208 format(a,9f13.8)
228 format(a,i2,a,9f13.8)

  end subroutine compute_rad_bal_coeffs

end module radsubs
