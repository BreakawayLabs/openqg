program openqg

  use util, only: s2s, streq
  use constants, only: PIBY2, SECDAY, SECSYR
  use units, only: m_to_km
  ! GLAM
  use box, only: box_type, init_box_from_mesh
  use grid, only: grid_type, load_grid, print_grid
  use topog, only: topog_type, load_topog
  use glam, only: load_glam
  use mesh, only: mesh_type

  ! QG
  use qg, only: qg_type, init_qg
  use coupler, only: ocn_coupler_type, init_ocn_coupler
  use coupler, only: atm_coupler_type, init_atm_coupler
  use basin, only: ocn_basin_type, atm_basin_type
  use basin, only: load_ocean_basin_config, ocn_basin_config_type
  use basin, only: load_atmos_basin_config, atm_basin_config_type


  use mixed, only: init_temp_from_rad
  use windstress, only: windstress_type
  use radsubs, only: radiate_type, fsprim
  use clock, only: clock_type

  use radsubs, only: load_rad
  use mixed, only: temp_type, zonal_temp
  use clock, only: load_clock, days_to_steps
  use ekman, only: ekman_from_tau, load_tau_from_file
  use windstress, only: xforc, init_windstress, windstress_config_type, load_windstress
  use state, only: save_qg, save_oml, save_aml, restart_qg
  use state, only: restart_oml, restart_aml
  use omlsubs, only: ocean_mixed_type, init_ocean_ml  
  use amlsubs, only: atmos_mixed_type, init_atmos_ml
  use forcing, only: compute_forcing
  use radsubs, only: radiat  
  use constraint, only: init_constraint
  use qg, only: init_foo_constr
  use driver, only: dynamic_step, diagnostic_step, diagnose_final

  ! Diagnostics
  use qg_monitor, only: qg_monitor_type, init_qg_monitor
  use tavsubs, only: tavout, init_basin_average
  use covsubs, only: basin_covar_type, init_basin_covar, covout
  use subsampling, only: subsampling_type, init_ocn_subsamp, init_atm_subsamp
  use areasubs, only: area_avg_type, init_area_avg
  use print, only: init_print
  use valsubs, only: valids_type
  use ml_monitor, only: init_oml_monitor, init_aml_monitor

  use vorsubs, only: init_pv

  implicit none

  call main()

contains

  subroutine main()

    implicit none

    type(clock_type) :: clk
    character (len=64) :: outdir  
    type(windstress_type) :: stress

    type(ocn_basin_type) :: ocn
    type(atm_basin_type) :: atm

    call init(ocn, atm, stress, clk, outdir)
    call run (ocn, atm, stress, clk, outdir)

  end subroutine main

  subroutine init(ocn, atm, stress, clk, outdir)

    type(windstress_type), intent(out) :: stress
    type(clock_type), intent(out) :: clk
    character (len=64), intent(out) :: outdir  
    type(ocn_basin_type), intent(out) :: ocn
    type(atm_basin_type), intent(out) :: atm

    type(ocn_basin_config_type) :: ocn_config
    type(atm_basin_config_type) :: atm_config
    type(mesh_type) :: mesh

    integer :: lenod
    character (len=64) :: indir
    logical :: use_ocean, use_atmos
    type(grid_type) :: g

#ifdef OPENMP
    ! Check OpenMP set up.
    call check_openmp()
#endif

    ! Parse command line arguments
    call get_command_argument(1, indir)
    call get_command_argument(2, outdir)
    lenod = index(outdir, '   ') - 1

    call load_glam(trim(indir)//"/glam.nc", use_ocean, use_atmos, mesh)
    if (use_atmos) then
       atm_config = load_atmos_basin_config(trim(indir)//"/atm_basin.nc")
    endif
    if (use_ocean) then
       ocn_config = load_ocean_basin_config(trim(indir)//"/ocn_basin.nc")
    endif

    call init_basins(indir, ocn_config, atm_config, mesh, ocn, atm)

    call init_mixed_and_qg(ocn_config, atm_config, ocn, atm)

    call init_dynamic_state(indir, ocn_config, atm_config, ocn, atm)

    call init_couplers(indir, ocn_config, atm_config, ocn, atm, stress)

    ! Diagnostic bits and pieces
    clk = load_clock(trim(indir)//"/clock.nc")

    call init_diagnostics(outdir(1:lenod), clk, ocn_config, atm_config, ocn, atm)

    ! Write a bunch of things to screen/output before beginning the actual run
    if (ocn%active .and. atm%active) then
       g = load_grid(ocn%qg%b, atm%qg%b)
       call print_grid(g, ocn%qg%b, atm%qg%b)
    endif

    call write_out_misc(outdir(1:lenod), ocn, atm, clk, g)

    ! Write Matlab-readable copy of input and derived parameters
    call write_out_param(outdir(1:lenod), ocn%b, atm%b, g, ocn%cov, atm%cov, clk, ocn%qg, atm%qg, &
         ocn%ml, atm%ml, stress, atm%cpl, ocn)

  end subroutine init

  subroutine init_basins(indir, ocn_config, atm_config, mesh, ocn, atm)
    character (len=*), intent(in) :: indir
    type(ocn_basin_config_type), intent(in) :: ocn_config
    type(atm_basin_config_type), intent(in) :: atm_config
    type(mesh_type), intent(in) :: mesh
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm

    if (atm_config%active) then
       atm%active = .true.
       allocate(atm%name(size(atm_config%name)))
       atm%name = atm_config%name
       atm%b = init_box_from_mesh(mesh, atm_config%nx, atm_config%ny, atm_config%ndx, atm_config%ndy, &
            atm_config%nx1, atm_config%ny1, atm_config%nl, atm_config%h, atm_config%hm, atm_config%dz_sign)

       call load_topog(atm%b, atm_config%topog, atm%topog, indir)
    endif
    if (ocn_config%active) then
       ocn%active = .true.
       allocate(ocn%name(size(ocn_config%name)))
       ocn%name = ocn_config%name

       ocn%b = init_box_from_mesh(mesh, ocn_config%nx, ocn_config%ny, ocn_config%ndx, ocn_config%ndy, &
            ocn_config%nx1, ocn_config%ny1, ocn_config%nl, ocn_config%h, ocn_config%hm, ocn_config%dz_sign)
       if (ocn%b%cyclic) then
          print *,' Model running in cyclic ocean configuration'
       else
          print *,' Model running in finite box ocean configuration'
       endif
       call load_topog(ocn%b, ocn_config%topog, ocn%topog, indir)
    endif

  end subroutine init_basins

  subroutine init_mixed_and_qg(ocn_config, atm_config, ocn, atm)
    type(ocn_basin_config_type), intent(in) :: ocn_config
    type(atm_basin_config_type), intent(in) :: atm_config
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm

    type(radiate_type) :: rad
    double precision :: rbtmat, rbtmoc
    double precision :: ocn_tmbar, atm_tmbar
    type(temp_type) :: octemp, attemp

    if (atm_config%use_qg) then
       atm%qg = init_qg(atm%b, atm_config%bcco, atm_config%delek, atm_config%ah2, atm_config%ah4, &
            atm_config%gp, atm_config%rho, atm_config%tau_sign, atm%topog)
    endif

    if (ocn_config%use_qg) then
       ocn%qg = init_qg(ocn%b, ocn_config%bcco, ocn_config%delek, ocn_config%ah2, ocn_config%ah4, &
            ocn_config%gp, ocn_config%rho, ocn_config%tau_sign, ocn%topog)
    endif

    if (atm_config%use_aml) then
       ! Compute mean state radiative balance and perturbation
       ! radiation coefficients A, B, C and D. Also compute atmosphere
       ! and ocean mixed layer temperatures that ensure equilibrium
       rad = load_rad(atm%b, atm_config%fsbar, atm_config%fspamp, &
            atm_config%zopt, atm_config%gamma, atm_config%xlamda, atm_config%T_abs)
       call radiat(atm%b, rad, rbtmat, rbtmoc, ocn_tmbar, atm_tmbar)
       call init_temp_from_rad(atm%b, atm%b, rad%tabsat(1), rad%tabsat(2), atm_tmbar, rbtmat, &
            rad%fspco, atm_config%rho, atm_config%cp, attemp)
       atm%ml = init_atmos_ml(atm%b, atm%topog, atm_config%hmamin, atm_config%hmadmp, atm_config%xcexp, attemp, &
            rad, atm_config%ahmd, atm_config%st2d, atm_config%st4d)
    endif

    if (ocn_config%use_oml) then
       if (ocn_config%use_rad_temp) then
          if (atm_config%use_aml) then
             call init_temp_from_rad(ocn%b, atm%b, ocn_config%T1_abs, ocn_config%T2_abs, ocn_tmbar, rbtmoc, &
                  rad%fspco, ocn_config%rho, ocn_config%cp, octemp)
             ocn%ml = init_ocean_ml(ocn%b, octemp, ocn_config%flux_n, ocn_config%flux_s, ocn_config%st2d, &
                  ocn_config%st4d, ocn_config%ycexp)
          else
             stop "use_rad_temp = .true. only valid with active atmospheric mixed layer"
          endif
       else
          stop "Currently unsupported option: use_rad_temp = .false."
       endif
    endif

  end subroutine init_mixed_and_qg

  subroutine init_dynamic_state(indir, ocn_config, atm_config, ocn, atm)
    character (len=*), intent(in) :: indir
    type(ocn_basin_config_type), intent(in) :: ocn_config
    type(atm_basin_config_type), intent(in) :: atm_config
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm

    double precision :: plfac(3),play
    integer :: j, k

    if (ocn%qg%active) then
       if (streq(ocn_config%qg_state, 'zero')) then
          ocn%qg%p = 0.0d0
          ocn%qg%pm = 0.0d0
       else
          call restart_qg(indir, s2s(ocn_config%qg_state), ocn%qg)
       endif
       call init_pv(ocn%qg, ocn%qg%topo)
       if (ocn%qg%b%cyclic) then
          ocn%qg%constr = init_constraint(ocn%qg%b, ocn%qg%mod, ocn%qg%p, ocn%qg%pm)
       endif
       call init_foo_constr(ocn%qg%b, ocn%qg%p, ocn%qg%pm, ocn%qg%con)
    endif
    if (atm%qg%active) then
       if (streq(atm_config%qg_state, 'zero')) then
          atm%qg%p = 0.0d0
          atm%qg%pm = 0.0d0
       else if (streq(atm_config%qg_state, 'rbal')) then
          ! Derive suitable multiplier of Fs' for each atmos. layer,
          ! from the eta coefficients of Fs' derived in radiat.
          ! We have nla layers but only nla-1 eta coeffts,
          ! and so need an extra constraint.
          plfac(1) = 0.0d0
          do k=2,atm%b%nl
             plfac(k) = plfac(k-1) - atm%qg%gp(k-1)*atm%ml%rad%rbetat(k-1)
          enddo

          ! Initialise atmospheric pressure
          do k=1,atm%b%nl
             do j=1,atm%b%nyp
                play = plfac(k)*fsprim( atm%b, atm%ml%rad%fspco, atm%qg%b%yprel(j) )
                atm%qg%p(:,j,k) = play
                atm%qg%pm(:,j,k) = play
             enddo
          enddo
       else
          call restart_qg(indir, s2s(atm_config%qg_state), atm%qg)
       endif
       call init_pv(atm%qg, atm%qg%topo)
       if (atm%qg%b%cyclic) then
          atm%qg%constr = init_constraint(atm%qg%b, atm%qg%mod, atm%qg%p, atm%qg%pm)
       endif
       call init_foo_constr(atm%qg%b, atm%qg%p, atm%qg%pm, atm%qg%con)
    endif

    if(ocn%ml%active) then
       if (streq(ocn_config%oml_state, 'zero')) then
          ocn%ml%sst%data = 0.0d0
          ocn%ml%sst%datam = 0.0d0
       else if (streq(ocn_config%oml_state, 'rbal')) then
          do j=1,ocn%ml%b%nyt
             ocn%ml%sst%data(:,j) = zonal_temp(ocn%ml%temp, j)
             ocn%ml%sst%datam(:,j) = ocn%ml%sst%data(:,j)
          enddo
       else
          call restart_oml(indir, s2s(ocn_config%oml_state), ocn%ml)
       endif
    endif
    if (atm%ml%active) then
       if (streq(atm_config%aml_state, 'zero')) then
          atm%ml%ast%data = 0.0d0
          atm%ml%ast%datam = 0.0d0
          atm%ml%hmixa%data = atm%ml%b%hm
          atm%ml%hmixa%datam = atm%ml%b%hm
       else if (streq(atm_config%aml_state, 'rbal')) then
          do j=1,atm%ml%b%nyt
             atm%ml%ast%data(:,j) = zonal_temp(atm%ml%temp, j)
             atm%ml%ast%datam(:,j) = atm%ml%ast%data(:,j)
          enddo
          atm%ml%hmixa%data = atm%ml%b%hm
          atm%ml%hmixa%datam = atm%ml%b%hm
       else
          call restart_aml(indir, s2s(atm_config%aml_state), atm%ml)
       endif
    endif

  end subroutine init_dynamic_state

  subroutine init_couplers(indir, ocn_config, atm_config, ocn, atm, stress)
    character (len=*), intent(in) :: indir
    type(ocn_basin_config_type), intent(in) :: ocn_config
    type(atm_basin_config_type), intent(in) :: atm_config
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm
    type(windstress_type), intent(out) :: stress

    type(windstress_config_type) :: ws

    if (ocn%qg%active .or. ocn%ml%active) then
       ocn%cpl = init_ocn_coupler(ocn%b, ocn_config%fnet_cpl, ocn_config%p1_cpl, ocn_config%ent_cpl, &
            ocn%qg%active, ocn%ml%active)
    endif

    if(atm%qg%active .or. atm%ml%active) then
       atm%cpl = init_atm_coupler(atm%b, ocn%b, atm_config%ent_cpl, atm_config%p1_cpl, &
            atm_config%eta_cpl, atm_config%sst_cpl)
    endif

    ! Initialise windstres
    ws = load_windstress(trim(indir)//"/windstress.nc")
    if (ws%coupled .and. ocn%qg%active .and. atm%qg%active) then
       stress = init_windstress(ocn%b, atm%b, ws%cdat, ws%tau_udiff)
       call xforc(stress, ocn%qg, atm%qg, ocn%cpl%ek, atm%cpl%ek)
       ! tau-udiff option only applies in coupled mode
       ! i.e. neither atmos-only nor ocean only are set
    else if (ocn%active) then
       call load_tau_from_file(indir, s2s(ws%ocn_stress), ocn%b, ocn%cpl%ek)
    else if (atm%active) then
       call load_tau_from_file(indir, s2s(ws%atm_stress), atm%b, atm%cpl%ek)
    endif
    if (ocn%active) then
       call ekman_from_tau(ocn%b, ocn%cpl%ek, ocn%qg%tau_sign)
    endif
    if (atm%active) then
       call ekman_from_tau(atm%b, atm%cpl%ek, ocn%qg%tau_sign)
    endif

    ! Initialis coupling dependent fields
    if (atm%ml%active .and. (atm%cpl%sst_coupled .or. (ocn%ml%active .and. ocn%cpl%fnet_coupled))) then
       atm%ml%go = ocn%b
    else
       atm%ml%go = atm%b
    endif
    if (atm%ml%active) then
       atm%ml%g = load_grid(atm%ml%go, atm%ml%b)
       allocate(atm%ml%fnetoc(atm%ml%go%nxt,atm%ml%go%nyt))
    endif

  end subroutine init_couplers

  subroutine init_diagnostics(outdir, clk, ocn_config, atm_config, ocn, atm)
    character (len=*), intent(in) :: outdir
    type(clock_type), intent(in) :: clk
    type(ocn_basin_config_type), intent(in) :: ocn_config
    type(atm_basin_config_type), intent(in) :: atm_config
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm

    integer :: nocmon, numoutsteps

    if (ocn%active) then
       if (ocn_config%compute_covar) then
          ocn%cov = init_basin_covar(ocn%b, ocn_config%dtcov, clk, s2s(ocn_config%cov_out), ocn_config%ns)
       endif
       if (ocn_config%compute_avges) then
          ocn%avg = init_basin_average(ocn%b, ocn_config%dtav, clk, s2s(ocn_config%avg_out))
       endif
       if (ocn_config%compute_subsample) then
          ocn%subsamp = init_ocn_subsamp(s2s(ocn_config%p_file), s2s(ocn_config%t_file), ocn_config%nsk, ocn_config%outfl, clk, &
               ocn_config%dtsubsamp, outdir, ocn%b)
       endif
       if (ocn_config%print_summary) then
          ocn%print = init_print(days_to_steps(ocn_config%dtprint, clk))
       endif
       if (ocn_config%check_valids) then
          ocn%valids%tauext = ocn_config%max_tau
          ocn%valids%wtext = ocn_config%max_wt
          ocn%valids%stext = ocn_config%max_st
          ocn%valids%pext = ocn_config%max_p
          ocn%valids%qext = ocn_config%max_q
          ocn%valids%nvalid = days_to_steps(ocn_config%dtvalid, clk)
          ocn%valids%active = .true.
       endif

       nocmon = days_to_steps(ocn_config%dtdiag, clk)
       numoutsteps = clk%ntsrun/nocmon + 1
       if (ocn_config%area_avg .and. ocn%ml%active) then
          ocn%area_avgs = init_area_avg(outdir, s2s(ocn_config%area_filename), ocn%b, ocn_config%narea, &
               ocn_config%xlo, ocn_config%xhi, ocn_config%ylo, ocn_config%yhi, nocmon, numoutsteps)
       endif
       if (ocn_config%monitor_qg .and. ocn%qg%active) then
          ocn%qg_mon = init_qg_monitor(ocn%qg%b, clk%dta*clk%nsteps0, s2s(ocn_config%qg_mon_out), &
               outdir, nocmon, numoutsteps, ocn%qg)
       endif
       if (ocn_config%monitor_ml .and. ocn%ml%active) then
          ocn%oml_mon = init_oml_monitor(nocmon, s2s(ocn_config%ml_mon_out), outdir, numoutsteps)
       endif
    endif

    if (atm%active) then
       if (atm_config%compute_covar) then
          atm%cov = init_basin_covar(atm%b, atm_config%dtcov, clk, s2s(atm_config%cov_out), atm_config%ns)
       endif
       if (atm_config%compute_avges) then
          atm%avg = init_basin_average(atm%b, atm_config%dtav, clk, s2s(atm_config%avg_out))
       endif
       if (atm_config%compute_subsample) then
          atm%subsamp = init_atm_subsamp(s2s(atm_config%p_file), s2s(atm_config%t_file), atm_config%nsk, ocn_config%outfl, clk, &
               atm_config%dtsubsamp, outdir, atm%b)
       endif
       if (atm_config%print_summary) then
          atm%print = init_print(days_to_steps(atm_config%dtprint, clk))
       endif
       if (atm_config%check_valids) then
          atm%valids%tauext = atm_config%max_tau
          atm%valids%wtext = atm_config%max_wt
          atm%valids%stext = atm_config%max_st
          atm%valids%pext = atm_config%max_p
          atm%valids%qext = atm_config%max_q
          atm%valids%nvalid = days_to_steps(atm_config%dtvalid, clk)
          atm%valids%active = .true.
       endif

       nocmon = days_to_steps(atm_config%dtdiag, clk)
       numoutsteps = clk%ntsrun/nocmon + 1
       if (atm_config%area_avg .and. atm%ml%active) then
          atm%area_avgs = init_area_avg(outdir, s2s(atm_config%area_filename), atm%b, atm_config%narea, &
               atm_config%xlo, atm_config%xhi, atm_config%ylo, atm_config%yhi, nocmon, numoutsteps)
       endif
       if (atm_config%monitor_qg .and. atm%qg%active) then
          atm%qg_mon = init_qg_monitor(atm%qg%b, clk%dta*clk%nsteps0, s2s(atm_config%qg_mon_out), &
               outdir, nocmon, numoutsteps, atm%qg)
       endif
       if (atm_config%monitor_ml .and. atm%ml%active) then
          atm%aml_mon = init_aml_monitor(nocmon, s2s(atm_config%ml_mon_out), ocn%b, atm%b, outdir, numoutsteps)
       endif
    endif

  end subroutine init_diagnostics

  subroutine run(ocn, atm, stress, clk, outdir)

    type(windstress_type), intent(in) :: stress
    type(clock_type), intent(in) :: clk
    character (len=64), intent(in) :: outdir  
    type(ocn_basin_type), intent(inout) :: ocn
    type(atm_basin_type), intent(inout) :: atm

    integer :: k,nt,lenod,ntdone
    double precision :: tday,tyrs
    logical :: solnok
    logical :: ocean_step, ocn_supp_step, atm_supp_step

    lenod = index(outdir, '   ') - 1
    if (atm%ml%active) then
       if (atm%cpl%eta_coupled) then
          do k=1,atm%ml%b%nl-1
             atm%cpl%eta(:,:,k) = (atm%qg%pm(:,:,k) - atm%qg%pm(:,:,k+1))/atm%qg%gp(k)
          enddo
       endif
       if (ocn%ml%active) then
          if (atm%cpl%sst_coupled) then
             atm%cpl%sst_datam = ocn%ml%sst%datam
          endif
          call compute_forcing(atm%cpl%sst_datam, atm%ml%go, atm%ml%b, atm%cpl%eta(:,:,1), &
               atm%ml%fnetoc, atm%ml%fnetat, atm%ml)
       endif
    endif

    nt = clk%nsteps0
    call diagnostic_step(outdir(1:lenod), ocn, atm, clk, nt, solnok)
    do nt = clk%nsteps0+1, clk%nsteps

       ocean_step    = mod(nt-1,clk%nstr) == 0
       ocn_supp_step = mod(nt-1,25*clk%nstr) == 0
       atm_supp_step = mod(nt-1,100) == 0
       call dynamic_step(ocn, atm, stress, clk%tdto, clk%tdta, ocean_step, ocn_supp_step, atm_supp_step)
       
       ! Timestep done; do checking and diagnostics as necessary
       tyrs = nt*clk%dta/SECSYR

       ! Periodically check validity of solution
       call diagnostic_step(outdir(1:lenod), ocn, atm, clk, nt, solnok)
       if ( .not.solnok ) then
          tday = nt*clk%dta/SECDAY
          print *,' valids has detected invalid values'
          print '(a,i12,f12.2,f11.3)', '  problem occurs at nt, tday, tyrs = ',nt,tday,tyrs
          call save_all(outdir(1:lenod), ocn, atm, "invalid.nc")
          print *,' program terminates'
          stop 1
       endif

       ! Occasionally dump restart file
       ntdone = nt - clk%nsteps0
       if (mod(ntdone,clk%noutre) == 0) then
          call save_all(outdir(1:lenod), ocn, atm, "restart.nc")
       endif
    enddo

    ! dump restart file
    call save_all(outdir(1:lenod), ocn, atm, "lastday.nc")
    call diagnose_final(outdir(1:lenod), ocn, atm)

    tday = nt*clk%dta/SECDAY
    print *
    print '(a,i12,f11.2,f11.4)', '  End of run at nt, tday, tyrs = ',nt,tday,tyrs

  end subroutine run

  subroutine save_all(outdir, ocn, atm, suffix)

    character (len=*), intent(in) :: outdir
    type(ocn_basin_type), intent(in) :: ocn
    type(atm_basin_type), intent(in) :: atm
    character (len=*), intent(in) :: suffix

    if (ocn%qg%active) then
       call save_qg(outdir, s2s(ocn%name)//"_qg_"//suffix, ocn%qg)
    endif
    if (atm%qg%active) then
       call save_qg(outdir, s2s(atm%name)//"_qg_"//suffix, atm%qg)
    endif
    if (ocn%ml%active) then
       call save_oml(outdir, s2s(ocn%name)//"_ml_"//suffix, ocn%ml)
    endif
    if (atm%ml%active) then
       call save_aml(outdir, s2s(atm%name)//"_ml_"//suffix, atm%ml)
    endif

  end subroutine save_all

  subroutine write_out_misc(outdir, ocn, atm, clk, g)
    character (len=*), intent(in) :: outdir
    type(ocn_basin_type), intent(in) :: ocn
    type(atm_basin_type), intent(in) :: atm
    type(clock_type), intent(in) :: clk
    type(grid_type), intent(in) :: g

    integer :: k

    print *,' '
    print *,' Control parameters:'
    print *,' ==================='
    print *,' outdir = ',outdir

    ! Print out a few interesting numbers
    print *,' '
    print *,' Oceanic parameters:'
    print *,' -------------------'
    print 201, '  No. of timesteps per day    = ',nint(SECDAY/clk%dto)
    print 203, '  Mixed layer thickness   (m) = ',ocn%b%hm
    if (ocn%ml%active .and. ocn%qg%active) then
       call diffts(2, ocn%b%nl, (/ocn%ml%sst%d2/), 1, ocn%b%dx, ocn%qg%mod%rdef)
       call diffts(4, ocn%b%nl, (/ocn%ml%sst%d4/), 1, ocn%b%dx, ocn%qg%mod%rdef)
    endif
    if (ocn%qg%active) then
       print 207, '  Reduced gravities  (m s^-2) = ', (ocn%qg%gp(k),k=1,ocn%b%nl-1)
       print 206, '  Baroclinic wavespeeds (m/s) = ', (ocn%qg%mod%c_phase(k),k=2,ocn%b%nl)
       print 206, '  Courant number(s)           = ', ( (clk%dto/ocn%b%dx)*ocn%qg%mod%c_phase(k),k=2,ocn%b%nl)
       print 204, '  Deformation radii      (km) = ', (m_to_km(ocn%qg%mod%rdef(k)),k=2,ocn%b%nl)
       print 205, '  Gridlengths per radius      = ', (ocn%qg%mod%rdef(k)/ocn%b%dx,k=2,ocn%b%nl)
       print 213, '  Del-sqd coeffts  (m^2 s^-1) = ', (ocn%qg%ah2(k),k=1,ocn%b%nl)
       call diffts(2, ocn%b%nl, ocn%qg%ah2, ocn%b%nl, ocn%b%dx, ocn%qg%mod%rdef)
       print 213, '  Del-4th coeffts  (m^4 s^-1) = ', (ocn%qg%ah4(k),k=1,ocn%b%nl)
       call diffts(4, ocn%b%nl, ocn%qg%ah4, ocn%b%nl, ocn%b%dx, ocn%qg%mod%rdef)
       print 204, '  Munk b.l. width scale  (km) = ', (m_to_km((ocn%qg%ah4(k)/ocn%b%beta)**0.2d0),k=1,ocn%b%nl)
       print 203, '  Bottom Ekm. layer thickness = ', ocn%qg%delek
       print 213, '  Bottom layer Ekman number   = ', (ocn%qg%delek/ocn%b%h(ocn%b%nl))**2
       if (ocn%qg%delek < 0.0d0) then
          print *,' Invalid -ve value of delek'         
          print *,' Program terminates'
          stop
       else if (ocn%qg%delek == 0.0d0) then
       else
          print 203, '  Spindown timescale   (days) = ',2.0d0*ocn%b%h(ocn%b%nl)/(abs(ocn%b%fnot)*ocn%qg%delek)/SECDAY
       endif
       print 213, '  Mixed BC coeff. bccooc (nd) = ',ocn%qg%bcco
    endif
    
    if (atm%active) then
       print *,' '
       print *,' Atmospheric parameters:'
       print *,' -----------------------'
       print 201, '  At ocean res., no. of cells = ',g%nxtaor,g%nytaor
       print 201, '  No. of timesteps per day    = ', nint(SECDAY/clk%dta)
       print 203, '  Mixed layer thickness   (m) = ',atm%b%hm
       call diffts(2, atm%b%nl, (/atm%ml%ast%d2/), 1, atm%b%dx, atm%qg%mod%rdef)
       call diffts(4, atm%b%nl, (/atm%ml%ast%d4/), 1, atm%b%dx, atm%qg%mod%rdef)
       call diffts(2, atm%b%nl, (/atm%ml%hmixa%d2/), 1, atm%b%dx, atm%qg%mod%rdef)
       print 207, '  Reduced gravities  (m s^-2) = ', (atm%qg%gp(k),k=1,atm%b%nl-1)
       print 206, '  Baroclinic wavespeeds (m/s) = ', (atm%qg%mod%c_phase(k),k=2,atm%b%nl)
       print 206, '  Courant number(s)           = ', ( (clk%dta/atm%b%dx)*atm%qg%mod%c_phase(k),k=2,atm%b%nl)
       print 204, '  Deformation radii      (km) = ', (m_to_km(atm%qg%mod%rdef(k)),k=2,atm%b%nl)
       print 205, '  Gridlengths per radius      = ', (atm%qg%mod%rdef(k)/atm%b%dx,k=2,atm%b%nl)
       print 213, '  Del-4th coeffts  (m^4 s^-1) = ', (atm%qg%ah4(k),k=1,atm%b%nl)
       call diffts(4, atm%b%nl, atm%qg%ah4, atm%b%nl, atm%b%dx, atm%qg%mod%rdef)
       print 213, '  Mixed BC coeff. bccoat (nd) = ',atm%qg%bcco
    endif

    print *,' '
    print *,' Coupling parameters:'
    print *,' --------------------'
    if (atm%ml%active) then
       print 205, '  Coefft. Lambda   (W m^-2/K) = ',atm%ml%rad%xlamda
    endif

201 format(a,9i13)
203 format(a,9f13.3)
204 format(a,9f13.4)
205 format(a,9f13.5)
206 format(a,9f13.6)
207 format(a,9f13.7)
213 format(a,1p,9d13.3)

  end subroutine write_out_misc

  subroutine write_out_param(outdir, go, ga, g, cov_ocn, cov_atm, clk, &
       qgo, qga, oml, aml, stress, atm_cpl, ocn)

    character (len=*), intent(in) :: outdir  
    type(box_type), intent(in) :: go
    type(box_type), intent(in) :: ga
    type(grid_type), intent(in) :: g
    type(basin_covar_type), intent(in) :: cov_ocn
    type(basin_covar_type), intent(in) :: cov_atm
    type(clock_type), intent(in) :: clk
    type(qg_type), intent(in) :: qgo
    type(qg_type), intent(in) :: qga
    type(ocean_mixed_type), intent(in) :: oml
    type(atmos_mixed_type), intent(in) :: aml
    type(windstress_type), intent(in) :: stress
    type(atm_coupler_type), intent(in) :: atm_cpl
    type(ocn_basin_type), intent(in) :: ocn

    integer :: k

    !     Write parameters out in matlab format file called
    !     input_parameters.m in the output directory

    open (10, file=outdir//'/input_parameters.m', &
         status='unknown')

    write(10,'(a)') '%%Matlab script to read in parameters'

    !     Preprocessor options set in Makefile
    write(10,123) 'oceanonly = ',(.not.atm_cpl%active), &
         ';  %% Ocean only run?'
    write(10,123) 'atmosonly = ',(.not.ocn%active), &
         ';  %% Atmos only run?'
    if (qgo%active) then
       write(10,123) 'cyclicoc = ',qgo%b%cyclic,';   %% Cyclic ocean?'
    endif
    if (oml%active) then
       write(10,123) 'hflxsb = ',oml%sst%fluxs, ';     %% S boundary heat flux?'
       write(10,123) 'hflxnb = ',oml%sst%fluxn, ';     %% N boundary heat flux?'
    endif
    write(10,123) 'tauudiff = ',stress%tau_udiff, ';   %% Use oc. vel. in tau?'

    !     Input parameters from parameter.src
    if (qgo%active) then
       write(10,122) 'nxto= ',qgo%b%nxt,';        %% Ocean x gridcells'
       write(10,122) 'nyto= ',qgo%b%nyt,';        %% Ocean y gridcells'
       write(10,122) 'nlo= ',qgo%b%nl,';         %% Ocean QG layers'
    endif
    if (qga%active) then
       write(10,122) 'nxta= ',qga%b%nxt,';        %% Atmos. x gridcells'
       write(10,122) 'nyta= ',qga%b%nyt,';        %% Atmos. y gridcells'
       write(10,122) 'nla= ',qga%b%nl,';         %% Atmos. QG layers'
    endif
    if(ocn%active .and. atm_cpl%active) then
       write(10,122) 'nxaooc= ',g%nxaooc, ';      %% Atmos. x gridcells over ocean'
       write(10,122) 'nyaooc= ',g%nyaooc, ';      %% Atmos. y gridcells over ocean'
       write(10,122) 'ndxr= ',g%ndxr, ';        %% Atmos./Ocean gridlength ratio'
       write(10,122) 'nx1= ',g%nx1,';         %% Starting index for', &
            ' ocean in atmospheric grid'
       write(10,122) 'ny1= ',g%ny1,';         %% Starting index for', &
            ' ocean in atmospheric grid'
    endif
    if (qgo%active) then
       write(10,121) 'fnot= ',qgo%b%fnot,';     %% Coriolis parameter'
       write(10,121) 'beta= ',qgo%b%beta,';     %% Beta'
    endif
    
    if (cov_ocn%active) then
       write(10,122) 'nso= ',cov_ocn%ns, &
            ';         %% Ocean x subsampling'
       write(10,122) 'nvo= ',cov_ocn%nv, &
            ';         %% Ocean subsamp vector'
       write(10,122) 'nco= ',cov_ocn%nc, &
            ';         %% Ocean covar components'
    endif
    if (cov_atm%active) then
       write(10,122) 'nsa= ',cov_atm%ns, &
            ';         %% Atmos. x subsampling'
       write(10,122) 'nva= ',cov_atm%nv, &
            ';         %% Atmos. subsamp vector'
       write(10,122) 'nca= ',cov_atm%nc, &
            ';         %% Atmos. covar components'
    endif

    !     Run parameters from input.params
    write(10,121) 'tini= ',clk%tini,';     %% Start time in years'
    write(10,121) 'trun= ',clk%trun,';     %% Run length in years'
    write(10,121) 'tend= ',clk%tini + clk%trun,';     %% Final time in years'
    write(10,121) 'dto= ',clk%dto, &
         ';      %% Ocean timestep in seconds'
    write(10,121) 'dta= ',clk%dta, &
         ';      %% Atmos. timestep in seconds'
    write(10,121) 'dxo= ',go%dx,';      %% Ocean grid spacing in km'
    write(10,121) 'dxa= ',ga%dx,';      %% Atmos. grid spacing in km'
    
    write(10,121) 'delek= ',qgo%delek, &
         ';    %% Ocean bottom Ekman thickness'
    if (qgo%active .and. qga%active) then
       write(10,121) 'cdat= ',stress%cdat,';     %% Air-Sea', &
            ' momentum drag coefficient'
    endif
    if (qga%active) then
       write(10,121) 'rhoat= ',qga%rho,';    %% Atmos. density'
    endif
    if (qgo%active) then
       write(10,121) 'rhooc= ',qgo%rho,';    %% Ocean  density'
    endif
    if (qga%active) then
       write(10,121) 'bccoat= ',qga%bcco, ';   %% Mixed BC coefficient for atmos (nondim.)'
    endif
    if (qgo%active) then
       write(10,121) 'bccooc= ',qgo%bcco, ';   %% Mixed BC coefficient for ocean (nondim.)'
    endif
    if (aml%active) then
       write(10,121) 'xcexp= ',aml%xcexp, ';    %% coupling coefficient x'
    endif
    if (oml%active) then
       write(10,121) 'ycexp= ',oml%sst%ycexp, ';    %% coupling coefficient y'
    endif

    !     Mixed layer parameters
    if (qgo%active) then
       write(10,121) 'hmoc= ',qgo%b%hm,';     %% Fixed ocean  ml depth'
    endif
    if (qga%active) then
       write(10,121) 'hmat= ',qga%b%hm,';     %% Fixed atmos. ml depth'
    endif
    if (oml%active) then
       write(10,121) 'st2d= ',oml%sst%d2, ';     %% sst lateral diffusivity'
    endif
    if (aml%active) then
       write(10,121) 'ahmd= ',aml%hmixa%d2, ';     %% hmixa lateral diffusivity'
       write(10,121) 'at2d= ',aml%ast%d2, ';     %% ast grad-2 diffusivity'
       write(10,121) 'at4d= ',aml%ast%d4, ';     %% ast grad-4 diffusivity'
    endif
    if (oml%active) then
       write(10,121) 'tsbdy= ',oml%temp%dT_NS, ';    %% o.m.l. S. bdy. temp (rel)'
    endif
    if (aml%active) then
       write(10,121) 'xlamda= ',aml%rad%xlamda, ';   %% Sensible/latent transfer'
       write(10,121) 'hmadmp= ',aml%hmadmp, ';   %% At. mixed layer h damping'
       
       ! Radiation parameters
       write(10,121) 'fsbar= ',aml%rad%fsbar, ';    %% Mean radiation forcing'
       write(10,121) 'fspamp= ',aml%rad%fspamp, ';   %% Radiation perturbation'
       write(10,121) 'zm= ',aml%rad%zopt(0), ';       %% Optical depth in a.m.l.'
       write(10,121) 'zopt= ',aml%rad%zopt(1), ';     %% Optical depth in layer 1'
       do k=2,aml%b%nl
          write(10,121) 'zopt= [zopt ',aml%rad%zopt(k),'];   %% Layers 2,n'
       enddo
       write(10,121) 'gamma= ',aml%rad%gamma,';    %% Adiabatic lapse rate'
    endif

    !     Oceanic QG layer parameters
    !     Reduced gravities, temperatures, friction coeffs and thicknesses
    if (qgo%active) then
       write(10,121) 'gpoc= ',qgo%gp(1), ';     %% Reduced gravity for ocean 1/2 interface'
       do k=2,qgo%b%nl-1
          write(10,121) 'gpoc= [gpoc ',qgo%gp(k), '];   %% Interfaces 2,n-1'
       enddo
       write(10,121) 'ah2oc= ',qgo%ah2(1),';    %% Del-sqd coefft ocean'
       do k=2,qgo%b%nl
          write(10,121) 'ah2oc= [ah2oc ',qgo%ah2(k),']; %% Layers 2,n'
       enddo
       write(10,121) 'ah4oc= ',qgo%ah4(1),';    %% Del-4th coefft ocean'
       do k=2,qgo%b%nl
          write(10,121) 'ah4oc= [ah4oc ',qgo%ah4(k),']; %% Layers 2,n'
       enddo
    endif
    if (oml%active) then
       write(10,121) 'tocc= ',oml%temp%T1_rel, ';     %% Rel. temperature for ocean layer 1'
    endif
    if (qgo%active) then
       write(10,121) 'hoc= ',qgo%b%h(1), ';      %% Thickness of ocean layer 1'
       do k=2,qgo%b%nl
          write(10,121) 'hoc= [hoc ',qgo%b%h(k),'];     %% Layers 2,n'
       enddo
    endif
    
    !     Atmospheric QG layer parameters
    !     Reduced gravities, temperatures, friction coeffs and thicknesses
    if (qga%active) then
       write(10,121) 'gpat= ',qga%gp(1), ';     %% Reduced gravity for atmos 1/2 interface'
       do k=2,qga%b%nl-1
          write(10,121) 'gpat= [gpat ',qga%gp(k), '];   %% Interfaces 2,n-1'
       enddo
       write(10,121) 'ah4at= ',qga%ah4(1),';    %% Del-4th coefft atmos'
       do k=2,qga%b%nl
          write(10,121) 'ah4at= [ah4at ',qga%ah4(k),'];   %% Layers 2,n'
       enddo
    endif
    if (aml%active) then
       write(10,121) 'tat= ',aml%temp%T1_rel, ';      %% Rel. temperature for atmos layer 1'
    endif
    if (qga%active) then
       write(10,121) 'hat= ',qga%b%h(1), ';      %% Thickness of atmos layer 1'
       do k=2,qga%b%nl
          write(10,121) 'hat= [hat ',qga%b%h(k),'];     %% Layers 2,n'
       enddo
    endif
    
    !     Derived parameters
    write(10,'(a)') '%%Derived parameters'
    !     Ocean
    if (qgo%active) then
       write(10,121) 'cphsoc= ',qgo%mod%c_phase(2), &
            ';   %% Baroclinic wavespeed for ocean mode 1'
       do k=3,qgo%b%nl
          write(10,121) 'cphsoc= [cphsoc ',qgo%mod%c_phase(k), '];   %% Higher modes'
       enddo
       write(10,121) 'rdefoc= ',qgo%mod%rdef(2), &
            ';   %% Deformation radius for ocean mode 1'
       do k=3,qgo%b%nl
          write(10,121) 'rdefoc= [rdefoc ',qgo%mod%rdef(k), '];   %% Higher modes'
       enddo
    endif
    if (oml%active) then
       write(10,121) 'dT_NS= ',oml%temp%dT_NS, ';    %% Rel. o.m.l. N/S. bndry temp. (K)'
    endif
    !     Atmosphere
    if (qga%active) then
       write(10,121) 'cphsat= ',qga%mod%c_phase(2), &
            ';   %% Baroclinic wavespeed for atmos mode 1'
       do k=3,qga%b%nl
          write(10,121) 'cphsat= [cphsat ',qga%mod%c_phase(k), '];   %% Higher modes'
       enddo
       write(10,121) 'rdefat= ',qga%mod%rdef(2), ';   %% Deformation radius for atmos mode 1'
       do k=3,qga%b%nl
          write(10,121) 'rdefat= [rdefat ',qga%mod%rdef(k), '];   %% Higher modes'
       enddo
    endif
    if (aml%active) then
       write(10,121) 'aface= ',aml%ent_coeff%aface(1), ';    %% eta    coefficient aface(1)'
       do k=2,aml%b%nl-1
          write(10,121) 'aface= [aface ',aml%ent_coeff%aface(k), &
               '];   %% Other interfaces'
       enddo
       write(10,121) 'bface= ',aml%ent_coeff%bface, ';    %% etam   coefficient bface'
       write(10,121) 'dface= ',aml%ent_coeff%dface, ';    %% aTm    coefficient dface'
    endif

    close(10)

121 format(A,1P,E13.5,A,A)
122 format(A,I10,A,A)
123 format(A,L,A,A)

  end subroutine write_out_param

  subroutine diffts (nord, nl, coeff, ncoef, dx, rdef)

    ! Computes diffusive decay timescales for circular eddies whose radii
    ! are the baroclinic Rossby radii, and for two-gridpoint noise.
    ! See section 8.6 of the Userguide for derivation of timescales.

    ! Input arguments:
    ! nord  : order of the diffusive term
    ! nl    : no. of QG layers  (=> nl-1 baroclinic modes)
    ! coeff : vector of diffusion coefficients (should be >= 0)
    ! ncoef : length of coefficient vector
    ! dx    : gridlength (m)
    ! rdef  : vector of nl modal deformation radii (m)
    !         (infinite value for barotropic mode replaced by 0.0)
    !         (all the above are unchanged on exit)


    integer, intent(in) :: nord,nl,ncoef
    double precision, intent(in) :: coeff(ncoef),dx,rdef(nl)

    integer :: nlmax
    parameter ( nlmax=9 )

    integer :: k,m
    double precision :: tdamp(nlmax),sinfac

    ! Check internal storage is sufficient
    if (nl > nlmax) then
       print *,' diffts has insufficient nlmax = ',nlmax
       print *,' called with nl = ',nl
       print *,' program terminates in diffts'
       stop
    endif

    ! Check all diffusion coefficients are non-negative
    ! (need positive coeffts for damping)
    do k=1,ncoef
       if (coeff(k) < 0.0d0) then
          print *,' diffts called with -ve diffusion coefft'
          print *,' coeff vector = ',(coeff(m),m=1,ncoef)
          print *,' program terminates in diffts'
          stop
       endif
    enddo

    ! Compute decay timescale(s) for a circular eddy
    ! at the Rossby radius for each baroclinic mode
    do m=2,nl
       sinfac = 2.0d0*sin( PIBY2*dx/rdef(m) )/dx
       ! Avoid infinities if coefft = 0
       do k=1,ncoef
          if (coeff(k) ==  0.0d0) then
             tdamp(k) = 0.0d0
          else
             tdamp(k) = 1.0d0/( sinfac**nord*coeff(k)*dble(nord)*SECDAY )
          endif
       enddo
       print 225, '  Mode',m-1,' damping time  (days) = ', &
            (tdamp(k),k=1,ncoef)
    enddo

    ! Compute decay timescale for two-gridpoint noise
    ! for each coefft, avoiding infinities if coefft = 0
    do k=1,ncoef
       if (coeff(k) == 0.0d0) then
          tdamp(k) = 0.0d0
       else
          tdamp(k) = (0.5d0*dx)**nord/coeff(k)/3600.0d0
       endif
    enddo
    print 205, '  Gridpoint timescale (hours) = ', &
         (tdamp(k),k=1,ncoef)

205 format(a,9f13.5)
225 format(a,i2,a,9f13.5)

  end subroutine diffts

#ifdef OPENMP
  subroutine check_openmp()

    integer :: nprocs, nthmax, numthr, thrnum
    integer :: OMP_GET_NUM_PROCS, OMP_GET_MAX_THREADS, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
    double precision :: dynadj
    double precision :: OMP_GET_DYNAMIC
    logical :: nested
    logical :: OMP_GET_NESTED

    integer :: j, nyp
    parameter ( nyp = 250 ) ! pick an arbitrary size example domain

    ! Examine OpenMP environment
    nprocs = OMP_GET_NUM_PROCS()
    nthmax = OMP_GET_MAX_THREADS()
    dynadj = OMP_GET_DYNAMIC()
    nested = OMP_GET_NESTED()
    print *, ' '
    print *, ' OpenMP parallelism is activated'
    print '(a,i5)', '  No. of processors available = ',nprocs
    print '(a,i3)', '  Max. no. of threads available = ',nthmax
    print *, ' Dynamic adjustment = ',dynadj
    print *, ' Nested parallelism = ',nested
    ! Test OpenMP is functioning correctly
    print *, ' Outside parallel section:'
    numthr = OMP_GET_NUM_THREADS()
    print *, ' Number of threads = ',numthr
    thrnum = OMP_GET_THREAD_NUM()
    print '(a,i4)', '  thrnum = ',thrnum

    print *, ' Test trivial parallel loop:'
    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP         PRIVATE (j,thrnum) &
    !$OMP         SHARED (numthr)
    !$OMP SINGLE
    numthr = OMP_GET_NUM_THREADS()
    print *, ' Number of threads = ',numthr
    !$OMP END SINGLE
    !$OMP DO SCHEDULE (STATIC)
    do j=1,numthr
       thrnum = OMP_GET_THREAD_NUM()
       print '(a,i6,i4)', '  j, thrnum = ',j,thrnum
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    ! The following loop is the most realistic test of OpenMP //ism
    ! Enable (if necessary) by removing the first ! from each line
    ! Comment out again once OpenMP is working correctly,
    ! to avoid lots of spurious output
    print *, ' Test typical j-range loop (ocean p-points)'
    !$OMP PARALLEL DO DEFAULT (NONE) &
    !$OMP         PRIVATE (j,thrnum) &
    !$OMP         SCHEDULE (STATIC)
    do j=1,nyp
       thrnum = OMP_GET_THREAD_NUM()
       print '(a,i6,i4)', '  j, thrnum = ',j,thrnum
    enddo
    !$OMP END PARALLEL DO

  end subroutine check_openmp
#endif

end program
