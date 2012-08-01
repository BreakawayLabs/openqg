module vorticity

  use qg, only: qg_type  
  use numerics, only: del2_P_bc, jacobian, int_P_dx

  implicit none

  private

  public vort_step

contains

  subroutine vort_step(qg, tdt, ent, wekp, txis, txin)

    type(qg_type), intent(inout) :: qg
    double precision, intent(in) :: tdt
    double precision, intent(in) :: ent(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision, intent(in) :: wekp(qg%b%nxp,qg%b%nyp)
    double precision, intent(in) :: txis, txin

    double precision :: qold(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision :: fohfac(qg%b%nl)
    double precision :: bdsums, bdsumn

    double precision :: dqdt(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision :: del2p(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision :: d4p(qg%b%nxp,qg%b%nyp,qg%b%nl),d6p(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision :: ah3sms,ah3smn,ah5sms,ah5smn
      
    double precision :: tmp_ajis(qg%b%nl), tmp_ajin(qg%b%nl)

    integer :: k
    ! Uses Arakawa energy and enstrophy conserving 9-point Jacobian
    if (qg%b%cyclic) then
       call jacobian(qg%p, qg%q, qg%b, dqdt, &
            qg%constr%ajis, qg%constr%ajin)
    else
       call jacobian(qg%p, qg%q, qg%b, dqdt, &
            tmp_ajis, tmp_ajin)
    endif


    ! Initialise dqdt in each layer with (Jacobian)
    ! advective term plus Del-4th and Del-6th frictional terms.

    ! Compute Del-sqd(p) at previous time level for dissipation
    del2p(:,:,:) = del2_P_bc(qg%pm(:,:,:), qg%b, qg%bcco)
    d4p(:,:,:) = del2_P_bc(del2p(:,:,:), qg%b, qg%bcco)
    d6p(:,:,:) = del2_P_bc(d4p(:,:,:),   qg%b, qg%bcco)

    do k=1,qg%b%nl
       dqdt(:,:,k) = dqdt(:,:,k) + (qg%ah2(k)/qg%b%fnot)*d4p(:,:,k) - (qg%ah4(k)/qg%b%fnot)*d6p(:,:,k)

       ! Compute third and fifth derivative contributions at
       ! zonal boundaries to the momentum constraint equations
       ! Only needed in the cyclic case
       if (qg%b%cyclic) then
          ah3sms = int_P_dx(del2p(:,2,k) - del2p(:,1,k), qg%b)/qg%b%dy
          ah3smn = int_P_dx(del2p(:,qg%b%nyp,k) - del2p(:,qg%b%nyp-1,k), qg%b)/qg%b%dy
          ah5sms = int_P_dx(d4p(:,2,k) - d4p(:,1,k), qg%b)/qg%b%dy
          ah5smn = int_P_dx(d4p(:,qg%b%nyp,k) - d4p(:,qg%b%nyp-1,k), qg%b)/qg%b%dy
          qg%constr%ap3s(k) = qg%ah2(k)*ah3sms
          qg%constr%ap3n(k) = qg%ah2(k)*ah3smn
          qg%constr%ap5s(k) = qg%ah4(k)*ah5sms
          qg%constr%ap5n(k) = qg%ah4(k)*ah5smn
       endif
    enddo


    fohfac(:) = qg%b%fnot/qg%b%h(:)
    !$OMP PARALLEL WORKSHARE SHARED(qg, fohfac, wekp, del2p, ent, dqdt, qold, tdt)

    ! Specify forcing and mass correction; then timestep
    ! These are the layer-specific terms near the top and bottom
    dqdt(:,:,1) = dqdt(:,:,1) + qg%b%dz_sign*fohfac(1)*ent(:,:,1) + qg%tau_sign*fohfac(1)*wekp(:,:) ! entrainment + wind
    dqdt(:,:,2) = dqdt(:,:,2) - qg%b%dz_sign*fohfac(2)*ent(:,:,1) ! entrainment
    dqdt(:,:,qg%topo%k_topo) = dqdt(:,:,qg%topo%k_topo) - &
         (fohfac(qg%topo%k_topo)*qg%delek/(2.0d0*abs(qg%b%fnot)))*del2p(:,:,qg%topo%k_topo) ! bottom drag

    ! Step the values of qo except at zonal boundaries
    ! The boundary values of qo will be updated later by ocqbdy,
    ! including meridional boundaries for the finite box case
    qold(:,:,:) = qg%q(:,:,:)
    qg%q(:,:,:) = qg%qm(:,:,:) + tdt*dqdt(:,:,:)
    qg%qm(:,:,:) = qold(:,:,:)
    !$OMP END PARALLEL WORKSHARE

    if (qg%b%cyclic) then
       qg%constr%txis = txis
       qg%constr%txin = txin

       qg%constr%int_Be_s(:) = 0.0d0
       qg%constr%int_Be_s(1) = qg%b%dz_sign*fohfac(1)*int_P_dx(ent(:,1,1), qg%b)
       qg%constr%int_Be_s(2) = qg%b%dz_sign*fohfac(1)*int_P_dx(ent(:,1,2), qg%b)
       qg%constr%int_Be_n(:) = 0.0d0
       qg%constr%int_Be_n(1) = -qg%b%dz_sign*fohfac(2)*int_P_dx(ent(:,qg%b%nyp,1), qg%b)
       qg%constr%int_Be_n(2) = -qg%b%dz_sign*fohfac(2)*int_P_dx(ent(:,qg%b%nyp,2), qg%b)

       ! Compute bottom drag first derivative contributions at
       ! zonal boundaries to the momentum constraint equations
       bdsums = int_P_dx(qg%pm(:,2,qg%topo%k_topo) - qg%pm(:,1,qg%topo%k_topo), qg%b)/qg%b%dy
       bdsumn = int_P_dx(qg%pm(:,qg%b%nyp,qg%topo%k_topo) - qg%pm(:,qg%b%nyp-1,qg%topo%k_topo), qg%b)/qg%b%dy
       qg%constr%bdrins = qg%b%fnot*qg%delek/(2.0d0*abs(qg%b%fnot))*bdsums
       qg%constr%bdrinn = qg%b%fnot*qg%delek/(2.0d0*abs(qg%b%fnot))*bdsumn
    endif

  end subroutine vort_step

end module vorticity
