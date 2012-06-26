module forcing

  use box, only: box_type
  use topog, only: topog_type
  use mixed, only: mixed_type
  use radsubs, only: fsprim
  use amlsubs, only: atmos_mixed_type
  use numerics, only: map_P_to_T

  implicit none
  
  private

  public compute_forcing

contains

  subroutine compute_forcing(sst_datam, go, ga, eta, &
       fnetoc, fnetat, aml)

    ! Compute forcings: stresses (7.1-7.5), Ekman velocities
    ! on both T and p grids (7.6-7.7), and net diabatic
    ! forcings of oceanic and atmos. mixed layers (7.8-7.10).

    type(box_type), intent(in) :: go
    double precision, intent(in) :: sst_datam(go%nxt,go%nyt)
    type(box_type), intent(in) :: ga
    type(atmos_mixed_type), intent(inout) :: aml
    double precision, intent(in) :: eta(aml%b%nxp,aml%b%nyp)
    double precision, intent(inout) :: fnetoc(go%nxt,go%nyt)
    double precision, intent(inout) :: fnetat(ga%nxt,ga%nyt)

    integer :: io,jo

    integer :: ia,ja
    double precision :: asto(go%nxt,go%nyt), fsp(go%nxt, go%nyt),ocfrac,fmafac,fmatop,hmafac
    double precision :: ocn_IR_up(go%nxt,go%nyt),sense_lat_flux(go%nxt,go%nyt),arlasm
    double precision :: aml_IR_down(go%nxt,go%nyt)
    integer :: natlan,natocn

    if (.not. aml%active) stop 1

    ! Compute atmospheric forcing

    ! Interpolate lagged atmospheric temperature onto oceanic grid
    call bilint (ga%xt, ga%yt, ga%nxt, ga%nyt, aml%ast%datam, &
         go%xt, go%yt, go%nxt, go%nyt, asto, 1.0d0, ga, go)

    ! Ocean/atmos infrared radiation
    ocn_IR_up(:,:) = aml%rad%D0up*sst_datam(:,:)
    aml%oradav = sum(ocn_IR_up(:,:))*go%norm

    ! Sensible and latent flux
    sense_lat_flux(:,:) = aml%rad%xlamda*( sst_datam(:,:) - asto(:,:) )
    aml%slhfav = sum(sense_lat_flux(:,:))*go%norm

    ! Atmospheric mixed layer radiation - into ocean
    aml_IR_down(:,:) = aml%rad%Dmdown*asto(:,:)
    aml%arocav = sum(aml_IR_down(:,:))*go%norm

    ! Specify atmospheric forcing everywhere
    ! Land case of (7.9) + last term of (7.8)
    ! These values will only be retained over land
    ! MONITORING - extra section for arlaav
    ! Atmospheric radiation based on astm
    !aml_IR_up - radiative forcing perturbation
    fnetat(:,:) = -aml%rad%Dup(0)*aml%ast%datam(:,:) - spread(fsprim( ga, aml%rad%fspco, ga%ytrel(:)), 1, ga%nxt)

    arlasm = sum(aml%ast%datam)
    ! Reset atmospheric forcing to zero over ocean
    natocn = 0
    do ja=aml%g%ny1,aml%g%ny1+aml%g%nyaooc-1
       do ia=aml%g%nx1,aml%g%nx1+aml%g%nxaooc-1
          fnetat(ia,ja) = fnetat(ia,ja) + (aml%rad%Dmdown - aml%rad%Dup(0))*aml%ast%datam(ia,ja)
          ! MONITORING - extra section for arlaav
          ! Count no. of atmos. cells over ocean
          arlasm = arlasm - aml%ast%datam(ia,ja)
          natocn = natocn + 1
       enddo
    enddo

    ! MONITORING - extra section for arlaav
    natlan = ga%nxt*ga%nyt - natocn
    if (natlan.eq.0) then
       aml%arlaav = 0.0d0
    else
       aml%arlaav = aml%rad%Dup(0)*arlasm/dble(natlan)
    endif

    ! Specify net ocean forcing (-ve of equation (7.10)),
    ! and modify the net atmospheric forcing over ocean
    ! points, equation (7.9) + last term of equation (7.8)
    ! fnetoc = - oF0
    ! fnetat = - aFm + aF0  (not including eta terms yet)

    ! Compute radiative forcing perturbation
    fsp(:,:) = spread(fsprim( ga, aml%rad%fspco, go%ytrel(:)), 1, go%nxt)
    ! Oceanic mixed layer diabatic forcing
    fnetoc(:,:) = -fsp(:,:) - aml_IR_down(:,:) - ocn_IR_up(:,:) - sense_lat_flux(:,:)

    ocfrac = go%dx*go%dy/(ga%dx*ga%dy)
    do jo=1,go%nyt
       ja = aml%g%ny1 + (jo-1)/aml%g%ndxr
       do io=1,go%nxt
          ! Atmospheric mixed layer radiation - in atmosphere
          ia = aml%g%nx1 + (io-1)/aml%g%ndxr
          fnetat(ia,ja) = fnetat(ia,ja) + ocfrac*( ocn_IR_up(io,jo) + sense_lat_flux(io,jo) )
       enddo
    enddo

    ! Complete atmospheric diabatic forcing at T points by
    ! adding terms involving local thickness perturbations
    ! and atmospheric topography. The eta1 and topography
    ! terms are averaged from the p grid to the T grid.
    ! These are the 1st, 2nd & 3rd terms of equation (7.8)
    fmafac = aml%rad%Adown(1,1)
    fmatop = aml%rad%Aup(0,0) + aml%rad%Adown(1,0)
    hmafac = -aml%hmadmp - aml%rad%Aup(0,0) - aml%rad%Adown(1,0)
    fnetat(:,:) = fnetat(:,:) - map_P_to_T(eta, ga, fmafac)  &
         - map_P_to_T(aml%topat%dtop, ga, fmatop) + hmafac*(aml%hmixa%datam(:,:) - ga%hm)

  end subroutine compute_forcing

  subroutine bilint (xa, ya, nxat, nyat, atmos, &
       xo, yo, nxoc, nyoc, ocean, fmult, ga, go)

    ! Performs bilinear interpolation of atmos(nxat,nyat), which
    ! is tabulated at coordinates xa(nxat), ya(nyat), to fill
    ! the array ocean(nxoc,nyoc), tabulated at xo(nxoc), yo(nyoc)
    ! Also multiplies the interpolant by the factor fmult
    ! Used to transfer data from atmospheric to oceanic grids

    integer, intent(in) ::  nxat,nyat,nxoc,nyoc
    double precision, intent(in) :: xa(nxat),ya(nyat),atmos(nxat,nyat), &
         xo(nxoc),yo(nyoc)
    double precision, intent(out) :: ocean(nxoc,nyoc)
    double precision, intent(in) :: fmult
    type(box_type), intent(in) :: ga
    type(box_type), intent(in) :: go

    integer :: io,jo,iam(go%nxp),iap(go%nxp),jam,jap
    double precision :: dxainv,dyainv,xam,wpx(go%nxp),wmx(go%nxp),wpy,wmy

    dxainv = 1.0d0/ga%dx
    dyainv = 1.0d0/ga%dy

    ! Get i-subscripts of ocean points in atmos array
    ! Compute subscripts once only; store in vector.
    ! Assumes fixed grid interval of dxa between xa values.
    do io=1,nxoc
       iam(io) = int( 1.0d0 + dxainv*( xo(io) - xa(1) ) )
       iap(io) = iam(io) + 1
       if ( iam(io).ge.1 ) then
          xam = xa(iam(io))
       else
          xam = xa(1) - ga%dx
       endif
       ! Compute x-direction weights (assumes regular grid)
       wpx(io) = dxainv*( xo(io) - xam )
       wmx(io) = 1.0d0 - wpx(io)
       ! Mend both pointers to give correct cyclic results for T points.
       ! Results will be inaccurate for p points, but this won't matter
       ! because the weight of the inaccurate value will be zero
       ! (p points never involve extrapolation, but T points can).
       iam(io) = 1 + mod( iam(io)+nxat-1, nxat )
       iap(io) = 1 + mod( iap(io)+nxat-1, nxat )
    enddo

    ! Compute y-direction weights and perform interpolation
    ! Assumes fixed grid intervals.
    do jo=1,nyoc
       jam = int( 1.0d0 + dyainv*( yo(jo) - ya(1) ) )
       jap = jam + 1
       ! Fix values for extrapolation.
       ! Boundary condition is no normal derivative.
       jam = max(jam,  1 )
       jap = min(jap,nyat)
       ! Compute y-direction weights (assumes regular grid)
       wpy = dyainv*( yo(jo) - ya(jam) )
       wmy = 1.0d0 - wpy
       do io=1,nxoc
          ocean(io,jo) = fmult*(  wmx(io)*wmy*atmos(iam(io),jam) &
               + wpx(io)*wmy*atmos(iap(io),jam) &
               + wmx(io)*wpy*atmos(iam(io),jap) &
               + wpx(io)*wpy*atmos(iap(io),jap) )
       enddo
    enddo

  end subroutine bilint

end module forcing
