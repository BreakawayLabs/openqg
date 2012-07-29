module vorsubs

  use box, only: box_type
  use qg, only: qg_type
  use topog, only: topog_type
  use modes, only: modes_type
  use numerics, only: del2_P_bc

  implicit none

  private

  public init_pv
  public ocqbdy

contains

  subroutine init_pv(qg, topo)

    type(qg_type), intent(inout) :: qg
    type(topog_type), intent(in) :: topo

    call qcomp (qg%q, qg%p, qg%mod%amat, qg%b, qg%bcco, topo)
    call qcomp (qg%qm,qg%pm,qg%mod%amat, qg%b, qg%bcco, topo)

  end subroutine init_pv

  subroutine qcomp (q, p, aaa, b, bcco, topo)

    ! Compute pot. vort. q from p  7.15)
    ! Works for both atmosphere and ocean

    ! Input arguments:
    ! p     : 3-D array of dynamic pressure at p points
    ! aaa   : A(nl,nl) matrix linking pressures and eta

    ! Output arguments:
    ! q    : 3-D array of vorticity at p points

    type(box_type), intent(in) :: b
    type(topog_type), intent(in) :: topo
    double precision, intent(in) :: bcco
    double precision q(b%nxp,b%nyp,b%nl),p(b%nxp,b%nyp,b%nl),aaa(b%nl,b%nl)

    integer :: i,j

    ! k = kbot layer (bottom - for topography)
    q(:,:,:) = del2_P_bc(p(:,:,:), b, bcco)/b%fnot
    q(:,:,topo%k_topo) = q(:,:,topo%k_topo) + topo%ddyn(:,:)

    do j=1,b%nyp
       q(:,j,:) = q(:,j,:) + b%beta*b%yprel(j)
       do i=1,b%nxp
          q(i,j,:) = q(i,j,:) - b%fnot*matmul(aaa, p(i,j,:))
       enddo
    enddo

  end subroutine qcomp

  subroutine ocqbdy (b, bcco, mod, p, topo, q)

     ! Derive oceanic boundary pot. vorticity qo from dynamic pressure po
     ! Equation (7.15), simplified because the tangential derivative p_tt
     ! is known to be zero on all solid boundaries, and the normal derivative
     ! p_nn is given by the mixed condition (2.27) (see also Appendix C)
     ! Does all solid boundaries, i.e. zonal boundaries, and
     ! also meridional boundaries if ocean is not cyclic.

     type(box_type), intent(in) :: b
     double precision, intent(in) :: bcco
     type(modes_type), intent(in) :: mod
     double precision, intent(in) :: p(b%nxp,b%nyp,b%nl)
     type(topog_type), intent(in) :: topo
     double precision, intent(inout) :: q(b%nxp,b%nyp,b%nl)

     integer :: i,j
     double precision :: fac

     ! Version with nondimensional bcco
     fac = bcco*b%dxm2/( 0.5d0*bcco + 1.0d0 )/b%fnot

     ! Zonal boundaries (incl. end/corner pts) - mixed BCs
     q(:,  1  ,:) = fac*( p(:,   2   ,:) - p(:,  1  ,:) )
     q(:,b%nyp,:) = fac*( p(:,b%nyp-1,:) - p(:,b%nyp,:) )
     if (.not. b%cyclic) then
        ! Meridional (internal) boundaries - mixed BCs
        q(  1  ,:,:) = fac*( p(   2  ,:,:) - p(  1 ,:,:) )
        q(b%nxp,:,:) = fac*( p(b%nxp-1,:,:) - p(b%nxp,:,:) )
     endif

     do i=1,b%nxp
        q(i,1,:) = q(i,1,:) - b%fnot*matmul(mod%amat, p(i,1,:))
        q(i,b%nyp,:) = q(i,b%nyp,:) - b%fnot*matmul(mod%amat, p(i,b%nyp,:))
     enddo
     if (.not. b%cyclic) then
        do j=1,b%nyp
           q(1,j,:) = q(1,j,:) - b%fnot*matmul(mod%amat, p(1,j,:))
           q(b%nxp,j,:) = q(b%nxp,j,:) - b%fnot*matmul(mod%amat, p(b%nxp,j,:))
        enddo
     endif

     q(:,1,:) = q(:,1,:) + b%beta*b%yprel(1)
     q(:,b%nyp,:) = q(:,b%nyp,:) + b%beta*b%yprel(b%nyp)
     if (.not. b%cyclic) then
        do j=1,b%nyp
           q(1,j,:) = q(1,j,:) + b%beta*b%yprel(j)
           q(b%nxp,j,:) = q(b%nxp,j,:) + b%beta*b%yprel(j)
        enddo
     endif

     q(:,1,topo%k_topo) = q(:,1,topo%k_topo) + topo%ddyn(:,1)
     q(:,b%nyp,topo%k_topo) = q(:,b%nyp,topo%k_topo) + topo%ddyn(:,b%nyp)
     if (.not. b%cyclic) then
        q(1,:,topo%k_topo) = q(1,:,topo%k_topo) + topo%ddyn(1,:)
        q(b%nxp,:,topo%k_topo) = q(b%nxp,:,topo%k_topo) + topo%ddyn(b%nxp,:)
     endif

   end subroutine ocqbdy

end module vorsubs
