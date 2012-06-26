module mixed_layer

  use box, only: box_type
  use mixed, only: mixed_type
  use numerics, only: del2_T, c_grid_advection
  use numerics, only: map_T_to_x_face, map_T_to_y_face

  implicit none

  private

  public adv_diff

contains

  function adv_diff(mixed, b, dp_dx, dp_dy, uek, vek)

    ! Computes the advective and diffusive
    ! contributions to the RHS of the evolution equation
    ! for the mixed layer quantity.
    ! Uses time-lagged mixed layer quantity.
    ! This version includes the coupling coefficient ycexp.
    ! dp_dx, dp_dy is from the pressure in layer 1.
    ! Diffusion implemented using dummy points west & east in del2t.
    ! Apply no-flux boundary condition to Del-4th diffusion
    ! across solid boundaries, equivalent to T'_{nnn} = 0

    type(mixed_type), intent(in) :: mixed
    type(box_type), intent(in) :: b
    double precision, intent(in) :: dp_dx(b%nxp-1,b%nyp)
    double precision, intent(in) :: dp_dy(b%nxp,b%nyp-1)
    double precision, intent(in) :: uek(b%nxp,b%nyp-1)
    double precision, intent(in) :: vek(b%nxp-1,b%nyp)
    double precision :: adv_diff(b%nxt,b%nyt)

    double precision :: vel_u(b%nxp,b%nyp-1)
    double precision :: vel_v(b%nxp-1,b%nyp)
    double precision :: mixed_x(b%nxp,b%nyp-1)
    double precision :: mixed_y(b%nxp-1,b%nyp)
    double precision :: del2(b%nxt,b%nyt)
    double precision :: del4(b%nxt,b%nyt)

    if (.not. mixed%active) stop 1

    ! Compute advective velocity
    call adv_vel(mixed, b, dp_dx, dp_dy, uek, vek, vel_u, vel_v)

    ! C-grid advection scheme, second order accurate
    call map_T_to_x_face(mixed%data, b, mixed_x)
    call map_T_to_y_face(mixed%data, b, mixed%fluxs, mixed%sval, mixed%fluxn, mixed%nval, mixed_y)
    adv_diff(:,:) = c_grid_advection(vel_u, vel_v, mixed_x, mixed_y, b)

    ! Add Del-sqd and Del-4th terms to temperature evolution term
    del2 = del2_T(mixed%datam, b, mixed%fluxn, mixed%nval, mixed%fluxs, mixed%sval)
    del4 = del2_T(del2,        b, .false.,     0.0d0,      .false.,     0.0d0)
    adv_diff(:,:) = adv_diff(:,:) + mixed%d2*del2(:,:) - mixed%d4*del4(:,:)

  end function adv_diff

  subroutine adv_vel(mixed, b, dp_dx, dp_dy, uek, vek, vel_u, vel_v)

    type(mixed_type), intent(in) :: mixed
    type(box_type), intent(in) :: b
    double precision, intent(in) :: dp_dx(b%nxp-1,b%nyp)
    double precision, intent(in) :: dp_dy(b%nxp,b%nyp-1)
    double precision, intent(in) :: uek(b%nxp,b%nyp-1)
    double precision, intent(in) :: vek(b%nxp-1,b%nyp)
    double precision, intent(out) :: vel_u(b%nxp,b%nyp-1)
    double precision, intent(out) :: vel_v(b%nxp-1,b%nyp)
    
    vel_u(:,:) = (-mixed%ycexp/b%fnot)*dp_dy(:,:) + uek(:,:)
    if (.not. b%cyclic) then
       ! no normal mass flux
       vel_u(1,:) = 0.0d0
       vel_u(b%nxp,:) = 0.0d0
    endif

    vel_v(:,:) = ( mixed%ycexp/b%fnot)*dp_dx(:,:) + vek(:,:)
    if (mixed%fluxs) then
       ! Advection consistent with an outflow across the southern boundary equal to the Ekman transport,
       ! carrying fluid of a specified temperature tsbdy. p contribution to vel_v vanishes because p is uniform along bdy.
       vel_v(:,1) = vek(:,1)
    else
       ! no normal mass flux
       vel_v(:,1) = 0.0d0
    endif
    if (mixed%fluxn) then
       ! Advection consistent with an outflow across the northern boundary equal to the Ekman transport,
       ! carrying fluid of a specified temperature tnbdy. p contribution to vel_v vanishes because p is uniform along bdy.
       vel_v(:,b%nyp) = vek(:,b%nyp)
    else
       ! no normal mass flux
       vel_v(:,b%nyp) = 0.0d0
    endif

  end subroutine adv_vel

end module mixed_layer
