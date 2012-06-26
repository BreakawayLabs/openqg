module omlsubs

  use box, only: box_type
  use mixed, only: mixed_type, init_mixed, temp_type, zonal_temp
  use numerics, only: dPdx, dPdy, avg_T, map_T_to_P
  use mixed_layer, only: adv_diff

  implicit none

  private

  type ocean_mixed_type

     logical :: active = .false.
     type(box_type) :: b
     type(temp_type) :: temp
     type(mixed_type) :: sst

  end type ocean_mixed_type

  public ocean_mixed_type
  public init_ocean_ml
  public step_oml

contains
  
  type(ocean_mixed_type) function init_ocean_ml(b, temp, fluxn, fluxs, st2d, st4d, ycexp)

    type(box_type), intent(in) :: b
    type(temp_type), intent(in) :: temp
    logical, intent(in) :: fluxn
    logical, intent(in) :: fluxs
    double precision, intent(in) :: st2d
    double precision, intent(in) :: st4d
    double precision, intent(in) :: ycexp

    init_ocean_ml%b = b
    init_ocean_ml%temp = temp
    init_ocean_ml%sst = init_mixed(b, 0.0d0, fluxn, fluxs, &
         zonal_temp(temp, b%nyt), zonal_temp(temp, 1), &
         st2d, st4d, ycexp)

    ! Check heat-flux boundary conditions
    if (fluxs) then
       if ( b%fnot.lt.0.0d0 ) then
          print *,' '
          print *,' Southern boundary o.m.l. heat-flux activated'
          print *,' Sign of fnot -> running in southern hemisphere'
          print *,' These are inconsistent choices'
          print *,' Program terminates'
          stop
       endif
    endif
    if (fluxn) then
       if ( b%fnot.gt.0.0d0 ) then
          print *,' '
          print *,' Northern boundary o.m.l. heat-flux activated'
          print *,' Sign of fnot -> running in northern hemisphere'
          print *,' These are inconsistent choices'
          print *,' Program terminates'
          stop
       endif
    endif

    init_ocean_ml%active = .true.

  end function init_ocean_ml

  subroutine step_oml(p, uek, vek, wekt, fnetoc, compute_ent, tdto, &
       ent, oml)

    ! Timestep oceanic mixed layer temperature, eqn (7.11)
    ! Also compute entrainment between oceanic layers (7.12-7.13)

    type(ocean_mixed_type), intent(inout) :: oml
    double precision, intent(in) :: p(oml%b%nxp,oml%b%nyp)
    double precision, intent(in) :: uek(oml%b%nxp,oml%b%nyp-1)
    double precision, intent(in) :: vek(oml%b%nxp-1,oml%b%nyp)
    double precision, intent(in) :: wekt(oml%b%nxt,oml%b%nyt)
    double precision, intent(in) :: fnetoc(oml%b%nxt,oml%b%nyt)
    logical, intent(in) :: compute_ent
    double precision, intent(in) :: tdto
    double precision, intent(inout) :: ent(oml%b%nxp,oml%b%nyp,oml%b%nl)

    double precision :: sstnew(oml%b%nxt,oml%b%nyt)
    double precision :: dT_conv(oml%b%nxt,oml%b%nyt)
    double precision :: half_wekt_dTm(oml%b%nxt,oml%b%nyt)
    double precision :: dp_dy(oml%b%nxp,oml%b%nyp-1)
    double precision :: dp_dx(oml%b%nxp-1,oml%b%nyp)
    double precision :: xfo(oml%b%nxt,oml%b%nyt)

    if (.not. oml%active) stop 1

    call dPdx(p(:,:), oml%b, dp_dx)
    call dPdy(p(:,:), oml%b, dp_dy)
    half_wekt_dTm(:,:) = 0.5d0*wekt(:,:)*( oml%sst%datam(:,:) - oml%temp%T1_rel )
    ! Predict new sst - equation (7.11)
    sstnew(:,:) = oml%sst%datam(:,:) +  &
         tdto*( adv_diff(oml%sst, oml%b, dp_dx, dp_dy, uek, vek) & ! Advection/diffusion
              + (1.0d0/(oml%temp%rho_cp*oml%b%hm))*fnetoc(:,:)   & ! Thermal forcing
              + (1.0d0/oml%b%hm)*(wekt(:,:)*oml%sst%datam(:,:) - half_wekt_dTm(:,:) ) ) ! ekman forcing

    ! Check for convection & if necessary correct layer 1/2 entrainment and sst - equation (7.13)
    ! dtonew should be <= 0 (stable case). Correction is nonzero only if dtonew > 0
    dT_conv(:,:) = max(0.0d0, oml%temp%T1_rel - sstnew(:,:))
    if (compute_ent) then
       ent(:,:,:) = 0.0d0
       ! Find layer 1/2 entrainment at T points - equation (7.12)
       xfo(:,:) = (-1.0d0/oml%temp%dT_12)*half_wekt_dTm(:,:) - (oml%b%hm/(tdto*oml%temp%dT_12))*dT_conv(:,:)
       ! Correct xfo values by subtracting mean to ensure net entrainment is zero, implying
       ! no net heat flux into the deep ocean.
       xfo(:,:) = xfo(:,:) - avg_T(xfo, oml%b)
       call map_T_to_P(xfo, oml%b, ent(:,:,1))
    endif

    oml%sst%datam(:,:) = oml%sst%data(:,:)
    oml%sst%data(:,:) = sstnew(:,:) + dT_conv(:,:)

  end subroutine step_oml

end module omlsubs
