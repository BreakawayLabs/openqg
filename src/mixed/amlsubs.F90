module amlsubs

  use ncutils, only: nc_open, nc_get_double, nc_close
  use box, only: box_type
  use topog, only: topog_type
  use mixed, only: mixed_type, temp_type, zonal_temp
  use radsubs, only: radiate_type
  use numerics, only: dPdx, dPdy, avg_T, map_T_to_P
  use numerics, only: c_grid_advection
  use mixed_layer, only: adv_diff
  use mixed, only: init_mixed
  use grid, only: grid_type

  implicit none

  private

  type entrain_coeff_type

     double precision, allocatable :: aface(:)
     double precision :: bface, dface
     !     aface, bface, cface, dface = coefficients for calculating atmos.
     !     entrainment (now assumed to be entirely across interface 1)
     !     from eta, etam, topography and aTm' respectively

  end type entrain_coeff_type


  type atmos_mixed_type

     logical :: active = .false.
     type(box_type) :: b
     type(topog_type) :: topat
     type(mixed_type) :: ast
     type(mixed_type) :: hmixa
     double precision, allocatable :: xc1ast(:,:)
     double precision, allocatable :: fnetat(:,:)
     type(entrain_coeff_type) :: ent_coeff

     double precision, allocatable :: fnetoc(:,:)
     type(box_type):: go
     type(grid_type) :: g

     double precision :: hmamin
     double precision :: hmadmp
     double precision :: xcexp

     type(temp_type) :: temp

     type(radiate_type) :: rad

     double precision :: slhfav,oradav,arocav,arlaav
     
     ! slhfav, oradav are heat fluxes averaged over ocean (W m^-2)
     ! arocav, arlaav are atmospheric radiative fluxes (upwards +ve)
     ! averaged over ocean and land respectively (W m^-2)
     ! All computed in xforc, so no duplication.
     
  end type atmos_mixed_type

  public entrain_coeff_type
  public atmos_mixed_type
  public init_atmos_ml
  public step_aml

contains

  type(atmos_mixed_type) function init_atmos_ml(b, topat, hmamin, hmadmp, xcexp, temp, rad, ahmd, at2d, at4d)

    type(box_type), intent(in) :: b
    type(topog_type), intent(in) :: topat
    double precision, intent(in) :: hmamin
    double precision, intent(in) :: hmadmp
    double precision, intent(in) :: xcexp
    type(temp_type), intent(in) :: temp
    type(radiate_type), intent(in) :: rad
    double precision, intent(in) :: ahmd
    double precision, intent(in) :: at2d
    double precision, intent(in) :: at4d

    init_atmos_ml%b = b
    init_atmos_ml%topat = topat

    init_atmos_ml%hmamin = hmamin
    init_atmos_ml%hmadmp = hmadmp
    init_atmos_ml%xcexp = xcexp

    init_atmos_ml%temp = temp
    init_atmos_ml%rad = rad

    allocate(init_atmos_ml%xc1ast(b%nxt,b%nyt))
    call set_xc1ast(b, temp, init_atmos_ml%xc1ast)

    allocate(init_atmos_ml%fnetat(b%nxt,b%nyt))

    init_atmos_ml%ast   = init_mixed(b, 0.0d0, .false., .false., 0.0d0, 0.0d0, &
               at2d, at4d,  1.0d0)
    init_atmos_ml%hmixa = init_mixed(b, b%hm,  .true.,  .true.,  b%hm,  b%hm, &
               ahmd, 0.0d0, 1.0d0)

    init_atmos_ml%ent_coeff = init_entrain_coeff(b%nl, temp, rad)

    init_atmos_ml%active = .true.
       
  end function init_atmos_ml

  type(entrain_coeff_type) function init_entrain_coeff(nl, attemp, rad)

    integer, intent(in) :: nl
    type(temp_type), intent(in) :: attemp
    type(radiate_type), intent(in) :: rad

    double precision :: rrcpdt
    integer :: l

    allocate(init_entrain_coeff%aface(0:nl-1))

    ! Entrainment factors: aface(l), bface, cface and dface such that
    ! e(1) = Sum(l)[ aface(l)*eta(l) ] + bface*etam + cface*aD + dface*aTm'
    ! We are assuming all entrainment is across interface 1.
    rrcpdt = -1.0d0/(attemp%rho_cp*attemp%dT_12)
    do l=0,nl-1
       init_entrain_coeff%aface(l) = rrcpdt*( 0.0d0        + rad%Adown(1,l) - rad%Aup(nl,l) )
    enddo
    init_entrain_coeff%bface       = rrcpdt*( rad%Aup(0,0) + rad%Adown(1,0) - rad%Aup(nl,0) )
    init_entrain_coeff%dface       = rrcpdt*( rad%Dup(0)   + 0.0d0          - rad%Dup(nl) )
    print *,' '
    print *,' Entrainment coefficients:'
    print *,' -------------------------'
    print *,' QG internal interface 1 coeffts:'
    print 215, '  eta  coefficients  aface(l) = ', &
         (init_entrain_coeff%aface(l),l=1,nl-1)
    print 215, '  etam+D coefficient    bface = ',init_entrain_coeff%bface
    print 215, '  aTm    coefficient    dface = ',init_entrain_coeff%dface

215 format(a,1p,9d13.5)

  end function init_entrain_coeff

  !! Next section is for the xcexp experiments
  !! There are a few options which require adjustments to
  !! the comments to make them current
  
  subroutine set_xc1ast(ga, attemp, xc1ast)

    type(box_type), intent(in) :: ga
    type(temp_type), intent(in) :: attemp
    double precision, intent(out) :: xc1ast(ga%nxt,ga%nyt)

    integer :: i, j, tempid

    character subnam*(*)
    parameter ( subnam = 'set_xc1ast' )

    !! Option 1: use astbar as determined by radiative equilibrium:
    if (.true.) then
       do j=1,ga%nyt
          xc1ast(:,j) = zonal_temp(attemp, j)
       enddo
    endif

    !! Options 2 & 3: both require reading in a file from the
    !! local directory which holds averaged ast information:
    if (.false.) then
       print *,' Reading in mean AST for xcexp experiment'
       tempid = nc_open('avges.nc', subnam)
       xc1ast = nc_get_double(tempid, 'ast', subnam)
       call nc_close(tempid)
    endif

    !! Option 2: Directly use the mean AST fields:
    if (.false.) then
    endif

    !! Option 3: Use the AST from data, but average
    !! it in the x-direction, so that there are no 'bumps':
    if (.false.) then
       do j=1,ga%nyt
          do i=2,ga%nxt
             xc1ast(1,j) = xc1ast(1,j) + xc1ast(i,j)
          enddo
          xc1ast(1,j) = xc1ast(1,j)/dble(ga%nxt)
          do i=2,ga%nxt
             xc1ast(i,j) = xc1ast(1,j)
          enddo
       enddo
    endif

  end subroutine set_xc1ast

  subroutine step_aml(p, uek, vek, wekt, eta, compute_ent, tdta, &
       ent, aml)

    ! Timestep atmospheric mixed layer height - equation (7.16)
    ! and mixed layer temperature - equation (7.17).
    ! Also compute entrainment between atmospheric layers
    ! - equation (7.18), plus convective correction (7.19).

    type(atmos_mixed_type), intent(inout) :: aml
    double precision, intent(in) :: p(aml%b%nxp,aml%b%nyp)
    double precision, intent(in) :: uek(aml%b%nxp,aml%b%nyp-1)
    double precision, intent(in) :: vek(aml%b%nxp-1,aml%b%nyp)
    double precision, intent(in) :: wekt(aml%b%nxt,aml%b%nyt)
    double precision, intent(in) :: eta(aml%b%nxp,aml%b%nyp,aml%b%nl-1)
    logical, intent(in) :: compute_ent
    double precision, intent(in) :: tdta
    double precision, intent(inout) :: ent(aml%b%nxp,aml%b%nyp,aml%b%nl-1)

    double precision :: tmrhs(aml%b%nxt,aml%b%nyt),dhdt_adv_vel(aml%b%nxt,aml%b%nyt), &
         hnew(aml%b%nxt,aml%b%nyt), astnew(aml%b%nxt,aml%b%nyt),  dT_conv(aml%b%nxt,aml%b%nyt)
    double precision :: dp_dy(aml%b%nxp,aml%b%nyp-1)
    double precision :: dp_dx(aml%b%nxp-1,aml%b%nyp)

    if (.not. aml%active) stop 1

    ! Initialise rh sides with advective, diffusive
    ! and Del-4th terms - first 3 terms in equation
    ! (7.16) and first four terms in equation (7.17)
    call dPdx(p(:,:), aml%b, dp_dx)
    call dPdy(p(:,:), aml%b, dp_dy)
    tmrhs(:,:) = adv_diff(aml%ast, aml%b, dp_dx, dp_dy, uek, vek)
    dhdt_adv_vel(:,:) = adv_diff(aml%hmixa, aml%b, dp_dx, dp_dy, uek, vek)

    call compute_new_aml(aml, wekt, tdta, tmrhs, dhdt_adv_vel, &
         astnew, hnew, dT_conv)

    if (compute_ent) then
       ent(:,:,:) = compute_atmos_entrainment(aml, aml%b, tdta, dT_conv, eta)
    endif

    aml%ast%datam(:,:) = aml%ast%data(:,:)
    aml%ast%data(:,:) = astnew(:,:)
    aml%hmixa%datam(:,:) = aml%hmixa%data(:,:)
    aml%hmixa%data(:,:) = hnew(:,:)

  end subroutine step_aml

  subroutine compute_new_aml(aml, wekt, tdta, tmrhs, dhdt_adv_vel, &
       astnew, hnew, dT_conv)

    type(atmos_mixed_type), intent(in) :: aml
    double precision, intent(in) :: wekt(aml%b%nxt,aml%b%nyt)
    double precision, intent(in) :: tdta
    double precision, intent(inout) :: tmrhs(aml%b%nxt,aml%b%nyt),dhdt_adv_vel(aml%b%nxt,aml%b%nyt)
    double precision, intent(out) :: astnew(aml%b%nxt,aml%b%nyt)
    double precision, intent(out) :: hnew(aml%b%nxt,aml%b%nyt)
    double precision, intent(out) :: dT_conv(aml%b%nxt,aml%b%nyt)

    integer :: i,j
    double precision :: hmainv,hdrcdt, &
         dhdt_diab,dh_fix,dT_fix(aml%b%nxt,aml%b%nyt), &
         rrcpat
    double precision :: dT(aml%b%nxt,aml%b%nyt)

    !! H-mixed
    dT(:,:) = aml%temp%T1_rel - aml%ast%datam(:,:)
    ! Add forcing term and diabatic effect - the 5th and
    ! 6th terms in equation (7.17), then timestep ast
    ! Also do last term in equation (7.16) and step that
    hdrcdt = aml%hmadmp/aml%temp%rho_cp
    do j=1,aml%b%nyt
       do i=1,aml%b%nxt
          ! Predict new value of hmixa - equation (7.16)
          if ( dT(i,j) >= 2.0d0*hdrcdt*tdta ) then
             dhdt_diab = hdrcdt*( aml%hmixa%datam(i,j) - aml%b%hm )/dT(i,j)
             hnew(i,j) = aml%hmixa%datam(i,j) + tdta*(dhdt_adv_vel(i,j) - dhdt_diab)
             dh_fix = max( aml%hmamin - hnew(i,j), 0.0d0 )
             hnew(i,j) = hnew(i,j) + dh_fix
             dT_fix(i,j) = dh_fix*dT(i,j)/aml%hmixa%datam(i,j)
          else
             hnew(i,j) = aml%b%hm
             dT_fix(i,j) = 0.0d0
          endif
       enddo
    enddo

    !! AST

    ! Predict new ast - equation (7.17)
    hmainv = 1.0d0/aml%b%hm
    rrcpat = 1.0d0/aml%temp%rho_cp
    tmrhs(:,:) = tmrhs(:,:) + rrcpat*aml%fnetat(:,:)/aml%hmixa%datam(:,:) - hmainv*wekt(:,:)*aml%ast%datam(:,:)
    astnew(:,:) = aml%ast%datam(:,:) + tdta*tmrhs(:,:) + dT_fix(:,:)

    ! Check for convection & if necessary correct layer 1/2
    ! entrainment and mixed layer temperature - equation (7.19)
    ! dtanew should be >= 0 (stable case)
    ! Correction is nonzero only if dtanew < 0
    dT_conv(:,:) = min(0.0d0, aml%temp%T1_rel - astnew(:,:))
    astnew(:,:) = astnew(:,:) + dT_conv(:,:)

  end subroutine compute_new_aml

  function compute_atmos_entrainment(aml, ga, tdta, dT_conv, eta)

    type(atmos_mixed_type), intent(in) :: aml
    type(box_type), intent(in) :: ga
    double precision, intent(in) :: tdta
    double precision, intent(in) :: dT_conv(ga%nxt,ga%nyt)
    double precision, intent(in) :: eta(ga%nxp,ga%nyp,ga%nl-1)
    double precision :: compute_atmos_entrainment(ga%nxp,ga%nyp,ga%nl-1)

    double precision :: xfa(ga%nxt,ga%nyt)
    integer :: k

    ! Check for convection & if necessary correct layer 1/2
    ! entrainment and mixed layer temperature - equation (7.19)
    ! dT_conv should be >= 0 (stable case)
    ! Correction is nonzero only if dT_conv < 0

    ! Find layer 1/2 entrainment at T points - equation (7.18)
    ! Just terms in eta-m and Tm' (naturally at T points) here;
    ! other terms (naturally at p points) are added later
    xfa(:,:) = aml%xcexp*aml%ent_coeff%bface*( aml%hmixa%datam(:,:) - ga%hm ) &
                         +  aml%ent_coeff%dface*( aml%xcexp*aml%ast%datam(:,:) + (1.0d0 - aml%xcexp)*aml%xc1ast(:,:) ) &
                         + (aml%xcexp/( tdta*aml%temp%dT_12 ))*aml%hmixa%data(:,:)*dT_conv(:,:)

    compute_atmos_entrainment(:,:,:) = 0.0d0
    call map_T_to_P(xfa, ga, compute_atmos_entrainment(:,:,1))

    ! Add eta and topography contributions
    ! which evaluate naturally at p points
    ! ------------------------------------
    ! Entrainment factors: aface(l), xbfac, cface and dface such that
    ! e(1) = Sum(l)[ aface(l)*eta(l) ] + xbfac*etam + cface*aD + dface*aTm'
    ! = Sum(l)[ afacdp(l)*dp(l) ] + xbfac*etam + cface*aD + dface*aTm'
    ! We are assuming all entrainment is across interface 1.
    do k=1,ga%nl-1
       compute_atmos_entrainment(:,:,1) = compute_atmos_entrainment(:,:,1) + aml%ent_coeff%aface(k)*eta(:,:,k)
    enddo
    compute_atmos_entrainment(:,:,1) = compute_atmos_entrainment(:,:,1) + aml%ent_coeff%bface*aml%topat%dtop(:,:)

  end function compute_atmos_entrainment

end module amlsubs
