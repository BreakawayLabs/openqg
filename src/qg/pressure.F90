module pressure

  use intsubs, only: xintp
  use qg, only: qg_type
  use box, only: box_type
  use topog, only: topog_type

  use modes, only: p_modes_rhs, modes_to_layers
  use inhomog, only: solve_inhomog_eqn
  use homog, only: cyclic_homog, box_homog
  use constraint, only: update_mass_constr

  implicit none

  private

  public solve_pressure

contains

  subroutine solve_pressure(qg, b, topo, tdt, ent)

    type(qg_type), intent(inout) :: qg
    type(box_type), intent(in) :: b
    type(topog_type), intent(in) :: topo
    double precision, intent(in) :: tdt
    double precision, intent(in) :: ent(b%nxp,b%nyp,b%nl)

    double precision :: rhs(b%nxp,b%nyp,b%nl)
    double precision :: inhomog(b%nxp,b%nyp,b%nl)
    double precision :: homcor(b%nxp,b%nyp,b%nl)
    integer :: m
    double precision :: ent_xn(b%nl-1)

    ! Store previous pressure values
    qg%pm(:,:,:) = qg%p(:,:,:)

    ! Compute vorticity RHS for each mode - equation (8.13)
    rhs(:,:,:) = p_modes_rhs(qg%q, topo, qg%mod, qg%b)
      
    ! Solve modified Helmholtz equation to invert each mode
    ! m = 1 is the barotropic mode
    ! m=2,nlo are baroclinic modes
    do m=1,b%nl
       inhomog(:,:,m) = solve_inhomog_eqn(qg%inhom, m, rhs(:,:,m))
    enddo
    ! Have solved inhomogeneous modal problem with po = 0 on all solid boundaries


    ! Compute area integral of entrainment between
    ! layers 1 and 2 to get net diabatic effect
    ent_xn(1) = xintp(ent(:,:,1), b%nxp, b%nyp)*b%dx*b%dy
    ent_xn(2:) = 0.0d0
    call update_mass_constr(qg%con, tdt, qg%gp, ent_xn)

    ! Solve constraint equations to obtain homogenous correction
    if (b%cyclic) then
       call cyclic_homog(b, qg%tau_sign, tdt, qg%constr, inhomog, &
            qg%hom%cyc, qg%mod, &
            homcor)
    else
       call box_homog(b, inhomog, qg%hom%box, qg%con, qg%mod, homcor)
    endif

    ! Add suitable multiple of homogeneous solutions and unpack modal
    ! pressures to layer pressures - equations (7.18) and (7.19)
    qg%p = modes_to_layers(inhomog(:,:,:) + homcor(:,:,:), qg%mod, b)

  end subroutine solve_pressure

end module pressure
