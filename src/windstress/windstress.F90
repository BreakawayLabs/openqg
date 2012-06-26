module windstress

  use intsubs, only: trapin
  use box, only: box_type
  use grid, only: grid_type, new_load_grid
  use bicubic, only: bcudata, bcuini, new_bicubic
  use ekman, only: ekman_type
  use qg, only: qg_type
  use numerics, only: map_P_to_x_face, map_P_to_y_face, map_T_to_P, dpdx_2bc, dpdy_2bc, dp2dx2_bc, dp2dy2_bc
  use ncutils, only: nc_open, nc_close, nc_get_double, nc_get_int, nc_get_dim, nc_get_text

  implicit none

  private

  type windstress_type

     logical :: active = .false.

     type(box_type) :: go
     type(box_type) :: ga
     type(grid_type) :: g
     type(bcudata) :: bcu
     double precision :: cdat
     logical :: tau_udiff

  end type windstress_type

  type windstress_config_type

     logical :: coupled = .false.

     double precision :: cdat
     ! cdat = Quadratic drag coefft. (nondim.)
     logical :: tau_udiff

     character, allocatable :: ocn_stress(:)
     character, allocatable :: atm_stress(:)

  end type windstress_config_type

  public windstress_config_type
  public load_windstress
  public windstress_type
  public init_windstress
  public xforc

contains

  type(windstress_config_Type) function load_windstress(filename)

    character (len=*), intent(in) :: filename

    integer :: windstress_id, strlen

    character :: subnam*(*)
    parameter ( subnam = 'load_windstress' )

    ! Namelist specified values
    windstress_id = nc_open(filename, subnam)

    load_windstress%coupled = nc_get_int(windstress_id, 'coupled', subnam) == 1

    if (load_windstress%coupled) then
       load_windstress%cdat = nc_get_double(windstress_id, 'cdat', subnam)
       load_windstress%tau_udiff = nc_get_int(windstress_id, 'tau_udiff', subnam) /= 0
       write(*,213) '  Quad. drag coefft Cd (nond) = ',load_windstress%cdat
213 format(a,1p,9d13.3)
    else
       strlen = nc_get_dim(windstress_id, 'strlen', subnam)
       allocate(load_windstress%ocn_stress(strlen))
       allocate(load_windstress%atm_stress(strlen))
       call nc_get_text(windstress_id, 'ocn_stress', load_windstress%ocn_stress, subnam)
       call nc_get_text(windstress_id, 'atm_stress', load_windstress%atm_stress, subnam)
    endif

  end function load_windstress

  type(windstress_type) function init_windstress(go, ga, cdat, tau_udiff)

    type(box_type), intent(in) :: go
    type(box_type), intent(in) :: ga
    double precision, intent(in) :: cdat
    logical, intent(in) :: tau_udiff

    init_windstress%go = go
    init_windstress%ga = ga
    init_windstress%g = new_load_grid(go, ga)
    init_windstress%bcu = bcuini(init_windstress%g)
    init_windstress%cdat = cdat
    init_windstress%tau_udiff = tau_udiff

    init_windstress%active = .true.

    if (tau_udiff) then
       print *,' Windstress depends on atmos.-ocean vel. diff.'
    else
       print *,' Windstress depends on atmos. velocity only'
    endif

  end function init_windstress

  subroutine xforc(stress, qgo, qga, &
    eko, eka)

    ! Compute forcings: stresses (7.1-7.5), Ekman velocities
    ! on both T and p grids (7.6-7.7), and net diabatic
    ! forcings of oceanic and atmos. mixed layers (7.8-7.10).

    ! Modified at v1.4.0 to compute tau at oceanic resolution, to
    ! enable tau optionally to depend on ocean mixed layer velocity
    ! (depending on setting of preprocessor option tau_udiff)

    type(windstress_type), intent(in) :: stress
    type(qg_type), intent(inout) :: qgo
    type(qg_type), intent(inout) :: qga

    type(ekman_type), intent(inout) :: eko
    type(ekman_type), intent(inout) :: eka

    double precision :: u1at(qga%b%nxp,qga%b%nyp),v1at(qga%b%nxp,qga%b%nyp)
    double precision :: dpdx(qga%b%nxp,qga%b%nyp),dpdy(qga%b%nxp,qga%b%nyp)
    double precision :: u2(qga%b%nxp,qga%b%nyp), cuohf(qga%b%nxp,qga%b%nyp), &
         shrmag(qga%b%nxp,qga%b%nyp), cdochi(qga%b%nxp,qga%b%nyp)
    double precision :: cdhfaa, qpreaa

    type(box_type) :: ga
    ga = qga%b

    call dPdx_2bc(qga%pm(:,:,1), ga, qga%bcco, dpdx)
    call dPdy_2bc(qga%pm(:,:,1), ga, qga%bcco, dpdy)
    u1at(:,:) = -dpdy(:,:)/ga%fnot
    v1at(:,:) =  dpdx(:,:)/ga%fnot

    cdhfaa = (stress%cdat/ga%fnot)/ga%hm

    qpreaa = 1.0d0/( abs(cdhfaa)*sqrt(2.0d0) )
    u2(:,:) = 4.0d0*cdhfaa*cdhfaa*(u1at(:,:)*u1at(:,:) + v1at(:,:)*v1at(:,:))
    shrmag(:,:) = qpreaa*sqrt( -1.0d0 + sqrt( 1.0d0 + u2(:,:) ) )
    cuohf(:,:) = cdhfaa*shrmag(:,:)
    cdochi(:,:) = stress%cdat*shrmag(:,:)/( 1.0d0 + cuohf(:,:)*cuohf(:,:) )
    eka%taux(:,:) = cdochi(:,:)*( u1at(:,:) - cuohf(:,:)*v1at(:,:) )
    eka%tauy(:,:) = cdochi(:,:)*( v1at(:,:) + cuohf(:,:)*u1at(:,:) )

    call compute_windstress(qgo, qga, stress%go, stress%ga, stress%g, &
         stress%tau_udiff, stress%bcu, stress%cdat, &
         eko, eka)

  end subroutine xforc

  subroutine compute_windstress(qgo, qga, go, ga, g, tau_udiff, bcu, cdat, &
       eko, eka)

    type(qg_type), intent(inout) :: qgo, qga
    type(box_type), intent(in) :: go, ga
    type(grid_type), intent(in) :: g
    logical, intent(in) :: tau_udiff
    type(bcudata), intent(in) :: bcu
    double precision, intent(in) :: cdat
    type(ekman_type), intent(inout) :: eko, eka

    double precision :: raoro
    integer :: ia, ja
    ! These are the offsets of ocean p points in the ocean-resolution
    ! atmospheric arrays, and their start and end indices therein; also
    ! the boundaries of the cells surrounding the timestepped q points
    double precision :: tauxaor(go%nxp,go%nyp)
    double precision :: tauyaor(go%nxp,go%nyp)
    double precision :: u1ator(go%nxp,go%nyp),v1ator(go%nxp,go%nyp)

    call atmos_geos_vel(g, ga, go, tau_udiff, bcu, qgo, qga, &
         u1ator, v1ator)

    raoro = qga%rho/qgo%rho
    call atmos_windstress(ga, go, tau_udiff, cdat, u1ator, v1ator, raoro, &
         tauxaor, tauyaor)

    ! Sample tau onto standard resolution atmos. grid
    ! ===============================================
    ! Copy across values at those points common to the two grids
    if (tau_udiff) then
       do ja=g%ny1,g%ny1+g%nypoar-1
          do ia=g%nx1,g%nx1+g%nxpoar-1
             eka%taux(ia,ja) = tauxaor(1+(ia-g%nx1)*g%ndxr, 1+(ja-g%ny1)*g%ndxr)
             eka%tauy(ia,ja) = tauyaor(1+(ia-g%nx1)*g%ndxr, 1+(ja-g%ny1)*g%ndxr)
          enddo
       enddo
    endif

    ! Rescale atmospheric stresses tauxaor, tauyaor
    ! to derive oceanic stresses tauxo and tauyo
    ! =============================================
    ! taux, tauy are actually dynamic stresses, so need to multiply
    ! by density ratio across interface - equation (7.5)
    eko%taux(:,:) = raoro*tauxaor(:,:)
    eko%tauy(:,:) = raoro*tauyaor(:,:)

  end subroutine compute_windstress

  subroutine atmos_geos_vel(g, ga, go, tau_udiff, bcu, qgo, qga, &
       u1ator, v1ator)

    type(grid_type), intent(in) :: g
    type(box_type), intent(in) :: ga
    type(box_type), intent(in) :: go
    logical, intent(in) :: tau_udiff
    type(bcudata), intent(in) :: bcu
    type(qg_type), intent(in) :: qgo
    type(qg_type), intent(in) :: qga
    double precision, intent(out) :: u1ator(go%nxp,go%nyp),v1ator(go%nxp,go%nyp)

    double precision :: dpdx(ga%nxp,ga%nyp),dpdy(ga%nxp,ga%nyp)
    double precision :: dp2_dx2(ga%nxp,ga%nyp),dp2_dy2(ga%nxp,ga%nyp)
    double precision :: dp2_dxy(ga%nxp,ga%nyp),dp3_dxxy(ga%nxp,ga%nyp)
    double precision :: dp3_dxyy(ga%nxp,ga%nyp)

    double precision :: u1at(g%nxpoar,g%nypoar),v1at(g%nxpoar,g%nypoar)
    double precision :: ux(g%nxpoar,g%nypoar),uy(g%nxpoar,g%nypoar),uxy(g%nxpoar,g%nypoar)
    double precision :: vx(g%nxpoar,g%nypoar),vy(g%nxpoar,g%nypoar),vxy(g%nxpoar,g%nypoar)

    ! Compute atmospheric geostrophic velocity at p points
    ! ====================================================
    ! Bilinearly interpolated fields aren't suitable for further
    ! differentiation, so we must derive au1 and av1 from ap1 at
    ! the normal atmospheric resolution, and then interpolate
    ! them onto the higher (oceanic) resolution grid.

    call dPdx_2bc(qga%pm(:,:,1), ga, qga%bcco, dpdx)
    call dPdy_2bc(qga%pm(:,:,1), ga, qga%bcco, dpdy)

    dp2_dx2 = dP2dx2_bc(qga%pm(:,:,1), ga, qga%bcco)
    dp2_dy2 = dP2dy2_bc(qga%pm(:,:,1), ga, qga%bcco)

    call dPdx_2bc(dpdy, ga, qga%bcco, dp2_dxy)
    call dPdx_2bc(dp2_dxy, ga, qga%bcco, dp3_dxxy)
    call dPdx_2bc(dp2_dy2, ga, qga%bcco, dp3_dxyy)

    ! Interpolate atmospheric geostrophic velocity to ocean resolution
    ! ================================================================
    u1at(:,:) = -dpdy(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)/ga%fnot
    v1at(:,:) =  dpdx(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)/ga%fnot

    ux(:,:) =  -(ga%dx/ga%fnot)*dp2_dxy(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)
    uy(:,:) =  -(ga%dy/ga%fnot)*dp2_dy2(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)
    uxy(:,:) = -(ga%dx*ga%dy/ga%fnot)*dp3_dxyy(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)

    vx(:,:) = (ga%dx/ga%fnot)*dp2_dx2(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)
    vy(:,:) = (ga%dy/ga%fnot)*dp2_dxy(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)
    vxy(:,:) = (ga%dx*ga%dy/ga%fnot)*dp3_dxxy(g%nx1:g%nx1+g%nxpoar-1,g%ny1:g%ny1+g%nypoar-1)
    
    call new_bicubic(g, go, u1at, ux, uy, uxy, bcu, u1ator)
    call new_bicubic(g, go, v1at, vx, vy, vxy, bcu, v1ator)

    ! Logically enough, can't compute the dependence on the
    ! oceanic velocity in an atmosphere only configuration!
    if (tau_udiff) then
       call compute_udiff(go, qgo, u1ator, v1ator)
    endif

  end subroutine atmos_geos_vel

  subroutine compute_udiff(go, qgo, u1ator, v1ator)

    type(box_type), intent(in) :: go
    type(qg_type), intenT(in) :: qgo
    double precision, intent(inout) :: u1ator(go%nxp,go%nyp),v1ator(go%nxp,go%nyp)

    double precision :: dpdx(go%nxp,go%nyp),dpdy(go%nxp,go%nyp)

    ! Compute geostrophic ocean velocity at p points, and
    ! amend atmospheric velocity to Q-G layer velocity diff.
    call dPdx_2bc(qgo%pm(:,:,1), go, qgo%bcco, dpdx)
    call dPdy_2bc(qgo%pm(:,:,1), go, qgo%bcco, dpdy)
    u1ator(:,:) = u1ator(:,:) + dpdy(:,:)/go%fnot
    v1ator(:,:) = v1ator(:,:) - dpdx(:,:)/go%fnot

  end subroutine compute_udiff

  subroutine atmos_windstress(ga, go, tau_udiff, cdat, u1ator, v1ator, raoro, &
       tauxaor, tauyaor)

    type(box_type), intent(in) :: ga
    type(box_type), intent(in) :: go
    logical, intent(in) :: tau_udiff
    double precision, intent(in) :: cdat
    double precision, intent(in) :: u1ator(go%nxp,go%nyp),v1ator(go%nxp,go%nyp)
    double precision, intent(in) :: raoro

    double precision, intent(inout) :: tauxaor(go%nxp,go%nyp),tauyaor(go%nxp,go%nyp)

    double precision :: cdhfaa, qpreaa
    double precision :: cuohf(go%nxp,go%nyp), shrmag(go%nxp,go%nyp), cdchi(go%nxp,go%nyp)

    ! Work out factors in windstress equations for both
    ! atmos. vel. (aa) and velocity difference (ab) cases
    ! ---------------------------------------------------
    ! Decide which to use later in code
    if (tau_udiff) then
       cdhfaa = (cdat/ga%fnot)*( 1.0d0/ga%hm + raoro/go%hm )
    else
       cdhfaa = (cdat/ga%fnot)/ga%hm
    endif

    qpreaa = 1.0d0/( abs(cdhfaa)*sqrt(2.0d0) )

    !$OMP PARALLEL WORKSHARE SHARED(cdhfaa, u1ator, v1ator, shrmag, qpreaa, cuohf,  &
    !$OMP      cdchi, cdat, tauxaor, tauyaor)
    shrmag(:,:) = sqrt( -1.0d0 + sqrt( 1.0d0 + 4.0d0*cdhfaa*cdhfaa*(u1ator(:,:)*u1ator(:,:) + v1ator(:,:)*v1ator(:,:)) ) )
    cuohf(:,:) = (qpreaa*cdhfaa)*shrmag(:,:)
    cdchi(:,:) = (qpreaa*cdat)*shrmag(:,:)/( 1.0d0 + cuohf(:,:)*cuohf(:,:) )
    tauxaor(:,:) = cdchi(:,:)*( u1ator(:,:) - cuohf(:,:)*v1ator(:,:) )
    tauyaor(:,:) = cdchi(:,:)*( v1ator(:,:) + cuohf(:,:)*u1ator(:,:) )
    !$OMP END PARALLEL WORKSHARE

  end subroutine atmos_windstress

end module windstress
