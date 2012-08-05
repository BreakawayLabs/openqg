module forcing

  use box, only: box_type
  use topog, only: topog_type
  use mixed, only: mixed_type
  use radsubs, only: fsprim
  use amlsubs, only: atmos_mixed_type
  use numerics, only: map_P_to_T, avg_T, bilint

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
    asto(:,:) = bilint(aml%ast%datam, ga, go)

    ! Ocean/atmos infrared radiation
    ocn_IR_up(:,:) = aml%rad%D0up*sst_datam(:,:)
    aml%oradav = avg_T(ocn_IR_up(:,:), go)

    ! Sensible and latent flux
    sense_lat_flux(:,:) = aml%rad%xlamda*( sst_datam(:,:) - asto(:,:) )
    aml%slhfav = avg_T(sense_lat_flux(:,:), go)

    ! Atmospheric mixed layer radiation - into ocean
    aml_IR_down(:,:) = aml%rad%Dmdown*asto(:,:)
    aml%arocav = avg_T(aml_IR_down(:,:), go)

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
    if (natlan == 0) then
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


end module forcing
