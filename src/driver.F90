module driver

  use basin, only: ocn_basin_type, atm_basin_type
  use coupler, only: ocn_coupler_type, atm_coupler_type
  use qg, only: qg_type
  use core, only: step_qg
  use clock, only: clock_type
  use windstress, only: windstress_type, xforc
  use amlsubs, only: atmos_mixed_type, step_aml
  use forcing, only: compute_forcing
  use omlsubs, only: ocean_mixed_type, step_oml
  use ekman, only: ekman_from_tau

  use qg_monitor, only: qg_monitor_step, update_qg_monitor
  use tavsubs, only:  avg_step, update_avg, tavout
  use covsubs, only: covar_step, update_cov, covout
  use subsampling, only: subsamp_step, ocnc_out, atnc_out
  use print, only: print_step, print_ocn, print_atm
  use areasubs, only: update_area_avg, area_avg_step
  use valsubs, only: valids, valid_step

  use ml_monitor, only: aml_monitor_step, oml_monitor_step
  use ml_monitor, only: update_aml_monitor, update_oml_monitor

  implicit none

  private

  public dynamic_step
  public diagnostic_step
  public diagnose_final

contains

  subroutine dynamic_step(ocn, atm, stress, tdto, tdta, ocean_step, ocn_supp_step, atm_supp_step)

    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm
    type(windstress_type), intent(in) :: stress
    double precision, intenT(in) :: tdto
    double precision, intenT(in) :: tdta
    logical, intent(in) :: ocean_step
    logical, intent(in) :: ocn_supp_step
    logical, intent(in) :: atm_supp_step

    ! Atmos/Ocean QG <-> QG exchange
    if ( ocean_step ) then
       if (stress%active) then
          call xforc(stress, ocn%qg, atm%qg, ocn%cpl%ek, atm%cpl%ek)
          call ekman_from_tau(stress%go, ocn%cpl%ek, ocn%qg%tau_sign)
          call ekman_from_tau(stress%ga, atm%cpl%ek, atm%qg%tau_sign)
       endif
    endif

    ! Ocean ML -> Atmos ML exchange
    if ( ocean_step ) then
       if (atm%ml%active .and. ocn%ml%active .and. atm%cpl%sst_coupled) then
          atm%cpl%sst_datam = ocn%ml%sst%datam
       endif
    endif

    ! Atmos time step
    call step_atm(atm%ml, atm%qg, atm%cpl, tdta, ocean_step)

    ! Atmos ML -> Ocean ML exchange
    if ( ocean_step ) then
       if (ocn%cpl%fnet_coupled .and. ocn%ml%active .and. atm%ml%active) then
          ocn%cpl%fnet = atm%ml%fnetoc
       endif
    endif

    ! Ocean update
    if ( ocean_step ) then
       call step_ocn(ocn%ml, ocn%qg, ocn%cpl, tdto)
    endif

    ! Suppress computational mode of leapfrog scheme
    if ( ocn_supp_step ) then
       call suppress_ocean(ocn%qg, ocn%ml)
    endif

    if ( atm_supp_step ) then
       call suppress_atmos(atm%qg, atm%ml)
    endif

  end subroutine dynamic_step

  subroutine step_ocn(oml, qgo, ocn_cpl, tdto)

    type(ocean_mixed_type), intent(inout) :: oml
    type(qg_type), intent(inout) :: qgo
    type(ocn_coupler_type), intent(inout) :: ocn_cpl
    double precision, intent(in) :: tdto

    ! Step ocean mixed layer
    if (oml%active) then
       ! QG -> ML exchange
       if (ocn_cpl%p1_coupled .and. qgo%active) then
          ocn_cpl%p1(:,:) = qgo%p(:,:,1)
       endif
       call step_oml(ocn_cpl%p1, ocn_cpl%ek%uek, ocn_cpl%ek%vek, ocn_cpl%ek%wekt, ocn_cpl%fnet, ocn_cpl%ent_coupled, tdto, &
            ocn_cpl%ent, oml)
    endif
    if (qgo%active) then
       ! Step ocean qg model
       call step_qg(tdto, ocn_cpl%ent, ocn_cpl%ek%wekp, ocn_cpl%ek%txis, ocn_cpl%ek%txin, qgo)
    endif

  end subroutine step_ocn

  subroutine step_atm(aml, qga, atm_cpl, tdta, ocean_step)
    type(atmos_mixed_type), intent(inout) :: aml
    type(qg_type), intent(inout) :: qga
    type(atm_coupler_type), intent(inout) :: atm_cpl
    double precision, intent(in) :: tdta
    logical, intent(in) :: ocean_step

    integer :: k

    ! Step atmospheric mixed layer
    if (aml%active) then
       ! QG -> ML exchange
       if (qga%active) then
          if (atm_cpl%eta_coupled) then
             do k=1,aml%b%nl-1
                atm_cpl%eta(:,:,k) = (qga%pm(:,:,k) - qga%pm(:,:,k+1))/qga%gp(k)
             enddo
          endif
          if (atm_cpl%p1_coupled) then
             atm_cpl%p1(:,:) = qga%p(:,:,1)
          endif
       endif
       if ( ocean_step ) then
          call compute_forcing(atm_cpl%sst_datam, aml%go, aml%b, atm_cpl%eta(:,:,1), &
               aml%fnetoc, aml%fnetat, aml)
       endif
       call step_aml(atm_cpl%p1, atm_cpl%ek%uek, atm_cpl%ek%vek, atm_cpl%ek%wekt, atm_cpl%eta, atm_cpl%ent_coupled, tdta, &
            atm_cpl%ent, aml)
    endif
    if (qga%active) then
       ! Step atmospheric qg channel model
       call step_qg(tdta, atm_cpl%ent, atm_cpl%ek%wekp, atm_cpl%ek%txis, atm_cpl%ek%txin, qga)
    endif

  end subroutine step_atm

  subroutine suppress_atmos(qga, aml)
    type(qg_type), intent(inout) :: qga
    type(atmos_mixed_type), intent(inout) :: aml

    if (qga%active) then
       call suppress_comp_mode(qga)
    endif
    if (aml%active) then
       aml%ast%data(:,:) = 0.5d0*(aml%ast%data(:,:)+aml%ast%datam(:,:) )
       aml%hmixa%data(:,:) = 0.5d0*(aml%hmixa%data(:,:)+aml%hmixa%datam(:,:))
    endif
  end subroutine suppress_atmos

  subroutine suppress_ocean(qgo, oml)
    type(qg_type), intent(inout) :: qgo
    type(ocean_mixed_type), intent(inout) :: oml

    if (qgo%active) then
       call suppress_comp_mode(qgo)
    endif
    if (oml%active) then
       oml%sst%data(:,:) = 0.5d0*( oml%sst%data(:,:)+oml%sst%datam(:,:) )
    endif

  end subroutine suppress_ocean

  subroutine suppress_comp_mode(qg)

    type(qg_type), intent(inout) :: qg

    ! Average oceanic time levels
    qg%q(:,:,:) = 0.5d0*( qg%q(:,:,:)+qg%qm(:,:,:) )
    qg%p(:,:,:) = 0.5d0*( qg%p(:,:,:)+qg%pm(:,:,:) )

    ! Also average constraint variables
    qg%con%dpi(:) = 0.5d0*( qg%con%dpi(:) + qg%con%dpip(:) )
    if (qg%b%cyclic) then
       qg%constr%cs(:) = 0.5d0*( qg%constr%cs(:) + qg%constr%csp(:) )
       qg%constr%cn(:) = 0.5d0*( qg%constr%cn(:) + qg%constr%cnp(:) )
    endif

  end subroutine suppress_comp_mode

  subroutine diagnostic_step(outdir, ocn, atm, clk, nt, solnok)
    character (len=*), intent(in) :: outdir
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm
    type(clock_type), intent(in) :: clk
    integer, intent(in) :: nt
    logical, intent(out) :: solnok

    logical :: force
    integer :: ntdone
    double precision :: tsec, tday, tyrs
    double precision :: secday,daysyr,secsyr
    parameter ( secday=86400.0d0, daysyr=365.0d0, &
         secsyr=secday*daysyr )

    tsec = nt*clk%dta
    tday = tsec/secday
    tyrs = tsec/secsyr
    ntdone = nt - clk%nsteps0
    solnok = .true.

    if (valid_step(ntdone, atm%valids)) then
       call valids(solnok, atm%valids, atm%qg, atm%ml%ast, atm%cpl%ek)
    endif
    if (valid_step(ntdone, ocn%valids)) then
       call valids(solnok, ocn%valids, ocn%qg, ocn%ml%sst, ocn%cpl%ek)
    endif

    force = .not. solnok

    ! Ensure that monitoring and averaging are done synchronously
    ! as post-processing routines assume this.
    if (qg_monitor_step(ocn%qg_mon, force, ntdone)) then
       call update_qg_monitor(ocn%qg, ocn%cpl%ek, ocn%cpl%ent, tsec, ocn%qg_mon, outdir, ntdone, tyrs, clk%dto)
    endif
    if (qg_monitor_step(atm%qg_mon, force, ntdone)) then
       call update_qg_monitor(atm%qg, atm%cpl%ek, atm%cpl%ent, tsec, atm%qg_mon, outdir, ntdone, tyrs, clk%dta)
    endif

    if (ocn%ml%active .and. area_avg_step(ntdone, ocn%area_avgs)) then
       call update_area_avg(outdir, ntdone, tyrs, ocn%ml%sst, ocn%ml%b, ocn%area_avgs)
    endif
    if (atm%ml%active .and. area_avg_step(ntdone, atm%area_avgs)) then
       call update_area_avg(outdir, ntdone, tyrs, atm%ml%ast, atm%ml%b, atm%area_avgs)
    endif

    if (oml_monitor_step(ntdone, ocn%oml_mon, force)) then
       call update_oml_monitor(ocn%ml, ocn%cpl%ek, ocn%oml_mon, outdir, ntdone, tyrs)
    endif
    if (aml_monitor_step(ntdone, atm%aml_mon, force)) then
       call update_aml_monitor(atm%qg, atm%aml_mon, atm%cpl%ek, atm%ml, atm%qg_mon%etam, outdir, ntdone, tyrs)
    endif


    if (subsamp_step(ntdone, ocn%subsamp, force)) then
       call ocnc_out (outdir, ntdone, tyrs, ocn%ml%sst, ocn%qg, ocn%cpl%ek, ocn%b, ocn%subsamp)
    endif
    if (subsamp_step(ntdone, atm%subsamp, force)) then
       call atnc_out (outdir, ntdone, tyrs, atm%ml%ast, atm%ml%hmixa, atm%qg, atm%cpl%ek, atm%b, atm%subsamp)
    endif

    if ( print_step(ocn%print, ntdone, force) .or. (print_step(atm%print, ntdone, force)) ) then
       print *,' '
       write(*,'(a,i12,f11.2,f11.4)') '  Sample output at nt, tday, tyrs = ',nt,tday,tyrs
    endif
    if (print_step(ocn%print, ntdone, force)) call print_ocn(ocn%qg, ocn%ml)
    if (print_step(atm%print, ntdone, force)) call print_atm(atm%qg, atm%ml)

    ! Add contributions to running means of forcing fields
    if (avg_step(ntdone, ocn%avg)) call update_avg(ocn%qg, ocn%b, ocn%cpl%ek, ocn%ml%sst, ocn%cpl%fnet, ocn%avg)
    if (avg_step(ntdone, atm%avg)) call update_avg(atm%qg, atm%b, atm%cpl%ek, atm%ml%ast, atm%ml%fnetat, atm%avg)

    ! Accumulate contributions to covariance matrices
    if (covar_step(ntdone, ocn%cov)) call update_cov(ocn%ml%sst, ocn%qg, ocn%cov)
    if (covar_step(ntdone, atm%cov)) call update_cov(atm%ml%ast, atm%qg, atm%cov)

  end subroutine diagnostic_step

  subroutine diagnose_final(outdir, ocn, atm)

    character (len=*), intent(in) :: outdir
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm

    ! Compute and output time-averaged fields
    call tavout(outdir, ocn%b, ocn%avg)
    call tavout(outdir, atm%b, atm%avg)

    ! Write out covariance matrix
    call covout(outdir, ocn%cov)
    call covout(outdir, atm%cov)

  end subroutine diagnose_final
  
end module driver
