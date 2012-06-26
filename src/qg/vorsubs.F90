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

    ! Internal points
    call qcomp (qg%q, qg%p, qg%mod%amat, qg%b, qg%bcco, topo)
    call qcomp (qg%qm,qg%pm,qg%mod%amat, qg%b, qg%bcco, topo)
    ! Zonal boundaries + meridional if not cyclic (mixed condition)
    call ocqbdy (qg%b, qg%bcco, qg%mod, qg%p, topo, qg%q)
    call ocqbdy (qg%b, qg%bcco, qg%mod, qg%pm, topo, qg%qm)
    if (qg%b%cyclic) then
       ! Meridional boundaries (periodic)
       call merqcy (qg%q, qg%p,  qg%mod%amat, qg%b, topo)         
       call merqcy (qg%qm,qg%pm, qg%mod%amat, qg%b, topo)
    endif

  end subroutine init_pv

  subroutine qcomp (q, p, aaa, b, bcco, topo)

    ! Compute pot. vort. q from p at internal points (equation 7.15)
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

    integer i,j,k
    double precision betay
    double precision :: del2p(b%nxp,b%nyp,b%nl)
    
    del2p(:,:,:) = del2_P_bc(p(:,:,:), b, bcco)

    ! k = kbot layer (bottom - for topography)
    q(:,:,:) = del2p(:,:,:)/b%fnot
    q(:,:,topo%k_topo) = q(:,:,topo%k_topo) + topo%ddyn(:,:)

    ! k = 1 layer (top(ocean) or bottom(atmos.))
    do j=2,b%nyp-1
       betay = b%beta*b%yprel(j)
       do i=2,b%nxp-1
          q(i,j,1) = q(i,j,1) + betay &
               - b%fnot*( aaa(1,1)*p(i,j,1) + aaa(1,2)*p(i,j,2) )
       enddo
    enddo
    ! Intermediate layers
    do k=2,b%nl-1
       do j=2,b%nyp-1
          betay = b%beta*b%yprel(j)
          do i=2,b%nxp-1
             q(i,j,k) = q(i,j,k) + betay &
                  - b%fnot*( aaa(k,k-1)*p(i,j,k-1) + aaa(k,k)*p(i,j,k) &
                  +aaa(k,k+1)*p(i,j,k+1) )
          enddo
       enddo
    enddo
    ! k = nl layer (bottom(ocean) or top(atmos.))
    do j=2,b%nyp-1
       betay = b%beta*b%yprel(j)
       do i=2,b%nxp-1
          q(i,j,b%nl) = q(i,j,b%nl) + betay &
               - b%fnot*(  aaa(b%nl,b%nl-1)*p(i,j,b%nl-1) &
               + aaa(b%nl, b%nl )*p(i,j, b%nl ) )
       enddo
    enddo

  end subroutine qcomp

  subroutine merqcy (q, p, aaa, b, topo)

    ! Compute pot. vort. on meridional boundaries for a cyclic domain
    ! Works for both atmosphere and ocean

    ! Input arguments:
    ! p     : 3-D array of dynamic pressure at p points
    ! aaa   : A(nl,nl) matrix linking pressures and eta
    ! ddyn  : dynamic (rescaled) topography
    ! kbot  : index of layer containing topography

    ! Output arguments:
    ! q    : 3-D array of vorticity at p points

    type(box_type), intent(in) :: b
    type(topog_type), intent(in) :: topo
    double precision q(b%nxp,b%nyp,b%nl),p(b%nxp,b%nyp,b%nl),aaa(b%nl,b%nl)

    integer j,k
    double precision dx2fac,betay

    dx2fac = b%dxm2/b%fnot

    ! k = 1 layer (top(ocean) or bottom(atmos.))
    do j=2,b%nyp-1
       betay = b%beta*b%yprel(j)
       ! Western boundary (periodic)
       q( 1 ,j, 1 ) = dx2fac* &
            (  p(1,j-1,1) + p(b%nxp-1,j,1) + p(2,j,1) &
            + p(1,j+1,1) - 4.0d0*p(1,j,1) ) + betay &
            - b%fnot*( aaa(1,1)*p(1,j,1) + aaa(1,2)*p(1,j,2) )
       ! Eastern boundary (periodic)
       q(b%nxp,j, 1 ) = q( 1 ,j, 1 )
    enddo

    ! Intermediate layers
    do k=2,b%nl-1
       do j=2,b%nyp-1
          betay = b%beta*b%yprel(j)
          ! Western boundary (periodic)
          q( 1 ,j,k) = dx2fac* &
               (  p( 1 ,j-1,k) + p(b%nxp-1,j,k) + p(2,j,k) &
               + p( 1 ,j+1,k) - 4.0d0*p( 1 ,j,k) ) + betay &
               - b%fnot*( aaa(k,k-1)*p(1,j,k-1) + aaa(k,k)*p(1,j,k) &
               +aaa(k,k+1)*p(1,j,k+1) )
          ! Eastern boundary (periodic)
          q(b%nxp,j,k) = q( 1 ,j,k)
       enddo
    enddo

    ! k = nl layer (bottom(ocean) or top(atmos.))
    do j=2,b%nyp-1
       betay = b%beta*b%yprel(j)
       ! Western boundary (periodic)
        q( 1 ,j,b%nl) = dx2fac* &
             (  p( 1 ,j-1,b%nl) + p(b%nxp-1,j,b%nl) + p(2,j,b%nl) &
             + p( 1 ,j+1,b%nl) - 4.0d0*p( 1 ,j,b%nl) ) + betay &
             - b%fnot*(  aaa(b%nl,b%nl-1)*p(1,j,b%nl-1) &
             + aaa(b%nl, b%nl )*p(1,j, b%nl ) )
        ! Eastern boundary (periodic)
        q(b%nxp,j,b%nl) = q( 1 ,j,b%nl)
     enddo

     ! k = kbot layer (bottom - for topography)
     do j=2,b%nyp-1
        q( 1 ,j,topo%k_topo) = q( 1 ,j,topo%k_topo) + topo%ddyn(1,j)
        ! Eastern boundary (periodic)
        q(b%nxp,j,topo%k_topo) = q( 1 ,j,topo%k_topo)
     enddo

   end subroutine merqcy

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

     integer j,k
     double precision fac,betays,betayn,f0Am,f0Ac,f0Ap,betay

     ! Version with nondimensional bcco
     fac = bcco*b%dxm2/( 0.5d0*bcco + 1.0d0 )/b%fnot

     ! Zonal boundaries (incl. end/corner pts) - mixed BCs
     betays = b%beta*b%yprel(  1 )
     betayn = b%beta*b%yprel(b%nyp)

     ! Top layer
     f0Ac = b%fnot*mod%amat(1,1)
     f0Ap = b%fnot*mod%amat(1,2)
     q(:,  1 ,1) =  fac*( p(:,   2  ,1) - p(:,  1 ,1) ) &
          - ( f0Ac*p(:,  1 ,1) + f0Ap*p(:,  1 ,2) ) &
          + betays
     q(:,b%nyp,1) =  fac*( p(:,b%nyp-1,1) - p(:,b%nyp,1) ) &
          - ( f0Ac*p(:,b%nyp,1) + f0Ap*p(:,b%nyp,2) ) &
          + betayn

     ! Intermediate layers
     do k=2,b%nl-1
        f0Am = b%fnot*mod%amat(k,k-1)
        f0Ac = b%fnot*mod%amat(k, k )
        f0Ap = b%fnot*mod%amat(k,k+1)
        q(:,  1 ,k) =  fac*( p(:,   2  ,k) - p(:,  1 ,k) ) &
             - (  f0Am*p(:,  1 ,k-1) + f0Ac*p(:,  1 ,k) &
             + f0Ap*p(:,  1 ,k+1) ) &
             + betays
        q(:,b%nyp,k) =  fac*( p(:,b%nyp-1,k) - p(:,b%nyp,k) ) &
             - (  f0Am*p(:,b%nyp,k-1) + f0Ac*p(:,b%nyp,k) &
             + f0Ap*p(:,b%nyp,k+1) ) &
             + betayn
     enddo

     ! Bottom layer
     f0Am = b%fnot*mod%amat(b%nl,b%nl-1)
     f0Ac = b%fnot*mod%amat(b%nl, b%nl )
     q(:,  1 ,b%nl) =  (fac*( p(:,   2  ,b%nl) - p(:,  1 ,b%nl) ) &
          - ( f0Am*p(:,  1 ,b%nl-1) + f0Ac*p(:,  1 ,b%nl) ) &
          + betays)
     q(:,b%nyp,b%nl) =  (fac*( p(:,b%nyp-1,b%nl) - p(:,b%nyp,b%nl) ) &
          - ( f0Am*p(:,b%nyp,b%nl-1) + f0Ac*p(:,b%nyp,b%nl) ) &
          + betayn)

     q(:,1,b%nl) = q(:,1,b%nl) + topo%ddyn(:,1)
     q(:,b%nyp,b%nl) = q(:,b%nyp,b%nl) + topo%ddyn(:,b%nyp)

     if (.not. b%cyclic) then
        ! Meridional (internal) boundaries - mixed BCs
        ! Top layer
        f0Ac = b%fnot*mod%amat(1,1)
        f0Ap = b%fnot*mod%amat(1,2)
        do j=2,b%nyp-1
           betay = b%beta*b%yprel(j)
           q(  1 ,j, 1 ) =  fac*( p(   2  ,j,1) - p(  1 ,j,1) ) &
                - ( f0Ac*p(  1 ,j,1) + f0Ap*p(  1 ,j,2) ) &
                + betay
           q(b%nxp,j, 1 ) =  fac*( p(b%nxp-1,j,1) - p(b%nxp,j,1) ) &
                - ( f0Ac*p(b%nxp,j,1) + f0Ap*p(b%nxp,j,2) ) &
                + betay
        enddo
        ! Intermediate layers
        do k=2,b%nl-1
           f0Am = b%fnot*mod%amat(k,k-1)
           f0Ac = b%fnot*mod%amat(k, k )
           f0Ap = b%fnot*mod%amat(k,k+1)
           do j=2,b%nyp-1
              betay = b%beta*b%yprel(j)
              q(  1 ,j,k) =  fac*( p(   2  ,j,k) - p(  1 ,j,k) ) &
                   - (  f0Am*p(  1 ,j,k-1) + f0Ac*p(  1 ,j,k) &
                   + f0Ap*p(  1 ,j,k+1) ) &
                   + betay
              q(b%nxp,j,k) =  fac*( p(b%nxp-1,j,k) - p(b%nxp,j,k) ) &
                   - (  f0Am*p(b%nxp,j,k-1) + f0Ac*p(b%nxp,j,k) &
                   + f0Ap*p(b%nxp,j,k+1) ) &
                   + betay
           enddo
        enddo

        ! Bottom layer (includes topography)
        f0Am = b%fnot*mod%amat(b%nl,b%nl-1)
        f0Ac = b%fnot*mod%amat(b%nl, b%nl )
        do j=2,b%nyp-1
           betay = b%beta*b%yprel(j)
           q(  1 ,j,b%nl) =  fac*( p(   2  ,j,b%nl) - p(  1 ,j,b%nl) ) &
                - ( f0Am*p(  1 ,j,b%nl-1) + f0Ac*p(  1 ,j,b%nl) ) &
                + betay
           q(b%nxp,j,b%nl) =  fac*( p(b%nxp-1,j,b%nl) - p(b%nxp,j,b%nl) ) &
                - ( f0Am*p(b%nxp,j,b%nl-1) + f0Ac*p(b%nxp,j,b%nl) ) &
                + betay
           q(1,j,topo%k_topo) = q(1,j,topo%k_topo) + topo%ddyn(1,j)
           q(b%nxp,j,topo%k_topo) = q(b%nxp,j,topo%k_topo) + topo%ddyn(b%nxp,j)
        enddo
     endif

   end subroutine ocqbdy

end module vorsubs
