module pressure

  use intsubs, only: xintp
  use qg, only: qg_type
  use box, only: box_type
  use topog, only: topog_type

  use modes, only: p_modes_rhs, modes_to_layers
  use inhomog, only: hscy, hsbx
  use homog, only: cyclic_homog, box_homog

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

    double precision :: wrk(b%nxp,b%nyp,b%nl)
    double precision :: homcor(b%nxp,b%nyp,b%nl)
    double precision :: bb_bx(b%nxt-1)
    double precision :: bb_cy(b%nxt/2+1)
    integer :: m
    double precision :: ent_xn(b%nl-1)

    ! Compute vorticity RHS for each mode - equation (8.13)
    wrk(:,:,:) = p_modes_rhs(qg%q, topo, qg%mod, qg%b)
      
    ! Solve modified Helmholtz equation to invert each mode
    ! m = 1 is the barotropic mode
    ! m=2,nlo are baroclinic modes

    do m=1,b%nl
       ! Compute Helmholtz operator for mode m
       ! Overwrite wrk with new modal pressure
       if (b%cyclic) then
          bb_cy(:) = qg%inhom%bd2(:) - qg%mod%rdm2(m)
          call hscy(qg%inhom, b, wrk(:,:,m), bb_cy)
       else
          bb_bx(:) = qg%inhom%bd2(:) - qg%mod%rdm2(m)
          call hsbx(qg%inhom, b, wrk(:,:,m), bb_bx)
       endif
    enddo
    ! Have solved inhomogeneous modal problem
    ! with po = 0 on all solid boundaries

    ! Compute area integral of entrainment between
    ! layers 1 and 2 to get net diabatic effect
    ent_xn(1) = xintp(ent(:,:,1), b%nxp, b%nyp)*b%dx*b%dy
    ent_xn(2:) = 0.0d0

    ! Solve constraint equations and add homogeneous solutions
    if (b%cyclic) then
       call cyclic_homog(b, qg%tau_sign, tdt,  &
            qg%constr, wrk,  &
            qg%hom%cyc, qg%mod, qg%con, qg%gp, &
            ent_xn, homcor)
    else
       call box_homog(b, tdt, wrk, qg%hom%box, qg%con, qg%gp, ent_xn, homcor)
    endif
    ! Add suitable multiple of homogeneous solutions and unpack modal
    ! pressures to layer pressures - equations (7.18) and (7.19)
    ! Also copy current pa to pam before overwriting with new value
    ! -------------------------------------------------------------
    ! Compute layer values
    qg%pm(:,:,:) = qg%p(:,:,:)
    qg%p = modes_to_layers(wrk(:,:,:) + homcor(:,:,:), qg%mod, b)

  end subroutine solve_pressure

end module pressure
