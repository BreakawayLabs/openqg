module basin

  use box, only: box_type
  use topog, only: topog_type
  use ncutils, only: nc_open, nc_get_dim, nc_get_int, nc_get_double, nc_get_text
  use omlsubs, only: ocean_mixed_type
  use amlsubs, only: atmos_mixed_type
  use qg, only: qg_type
  use qg_monitor, only: qg_monitor_type
  use subsampling, only: subsampling_type
  use coupler, only: ocn_coupler_type, atm_coupler_type
  use covsubs, only: basin_covar_type
  use tavsubs, only: basin_averages_type
  use print, only: print_type
  use areasubs, only: area_avg_type
  use valsubs, only: valids_type
  use ml_monitor, only: oml_monitor_type, aml_monitor_type

  implicit none

  private

  type ocn_basin_config_type

     logical :: active = .false.

     ! name
     character, allocatable :: name(:)

     ! box
     integer :: nx, ny, nl
     integer :: ndx, ndy
     integer :: nx1, ny1
     double precision, allocatable :: h(:)
     double precision :: hm
     integer :: dz_sign = -1

     ! Topo
     character, allocatable :: topog(:)

     ! Common
     double precision :: rho

     ! OML
     logical :: use_oml = .false.
     logical :: flux_n, flux_s = .false.
     double precision :: st2d, st4d, ycexp
     double precision :: cp
     logical :: use_rad_temp = .false.
     ! rad temp values
     double precision :: T1_abs, T2_abs
     ! not rad temp values
     double precision :: T1_rel, T2_rel, dT_NS
     character, allocatable :: oml_state(:)

     ! QG
     logical :: use_qg  = .false.
     double precision :: bcco
     double precision :: delek
     double precision, allocatable :: ah2(:), ah4(:)
     double precision, allocatable :: gp(:)
     integer :: tau_sign = 1
     character, allocatable :: qg_state(:)

     ! Coupling
     character, allocatable :: fnet_cpl(:)
     character, allocatable :: p1_cpl(:)
     character, allocatable :: ent_cpl(:)

     ! Diagnostics
     logical :: compute_covar = .false.
     double precision :: dtcov
     integer :: ns
     character, allocatable :: cov_out(:)

     logical :: compute_avges = .false.
     double precision :: dtav
     character, allocatable :: avg_out(:)

     logical :: compute_subsample = .false.
     double precision :: dtsubsamp
     character, allocatable :: p_file(:)
     character, allocatable :: t_file(:)     
     logical :: outfl(7)
     integer :: nsk

     logical :: check_valids = .false.
     double precision :: dtvalid
     double precision :: max_tau
     double precision :: max_wt
     double precision :: max_st
     double precision :: max_p
     double precision :: max_q
     
     logical :: print_summary = .false.
     double precision :: dtprint

     logical :: area_avg = .false.
     integer :: narea
     double precision, allocatable :: xlo(:)
     double precision, allocatable :: xhi(:)
     double precision, allocatable :: ylo(:)
     double precision, allocatable :: yhi(:)
     character, allocatable :: area_filename(:)

     logical :: monitor_ml = .false.
     character, allocatable :: ml_mon_out(:)

     logical :: monitor_qg = .false.
     character, allocatable :: qg_mon_out(:)

     double precision :: dtdiag

  end type ocn_basin_config_type

  type atm_basin_config_type

     logical :: active = .false.

     ! name
     character, allocatable :: name(:)

     ! box
     integer :: nx, ny, nl
     integer :: ndx, ndy
     integer :: nx1, ny1
     double precision, allocatable :: h(:)
     double precision :: hm
     integer :: dz_sign = 1

     ! Topo
     character, allocatable :: topog(:)

     ! Common
     double precision :: rho

     ! AML
     logical :: use_aml = .false.
     logical :: flux_n, flux_s = .false.
     double precision :: st2d, st4d, ahmd
     double precision :: cp
     double precision :: hmamin
     double precision :: hmadmp
     double precision :: xcexp
     ! Radiation parameters
     double precision , allocatable:: T_abs(:)
     double precision :: fsbar, fspamp
     double precision, allocatable :: zopt(:)
     double precision :: gamma, xlamda
     character, allocatable :: aml_state(:)

     ! QG
     logical :: use_qg = .false.
     double precision :: bcco
     double precision :: delek
     double precision, allocatable :: ah2(:), ah4(:)
     double precision, allocatable :: gp(:)
     integer :: tau_sign = -1
     character, allocatable :: qg_state(:)

     ! Coupling
     character, allocatable :: sst_cpl(:)
     character, allocatable :: p1_cpl(:)
     character, allocatable :: eta_cpl(:)
     character, allocatable :: ent_cpl(:)

     ! Diagnostics
     logical :: compute_covar = .false.
     double precision :: dtcov
     integer :: ns
     character, allocatable :: cov_out(:)

     logical :: compute_avges = .false.
     double precision :: dtav
     character, allocatable :: avg_out(:)

     logical :: compute_subsample = .false.
     double precision :: dtsubsamp
     character, allocatable :: p_file(:)
     character, allocatable :: t_file(:)     
     logical :: outfl(7)
     integer :: nsk

     logical :: check_valids = .false.
     double precision :: dtvalid
     double precision :: max_tau
     double precision :: max_wt
     double precision :: max_st
     double precision :: max_p
     double precision :: max_q

     logical :: print_summary = .false.
     double precision :: dtprint

     logical :: area_avg = .false.
     integer :: narea
     double precision, allocatable :: xlo(:)
     double precision, allocatable :: xhi(:)
     double precision, allocatable :: ylo(:)
     double precision, allocatable :: yhi(:)
     character, allocatable :: area_filename(:)

     logical :: monitor_ml = .false.
     character, allocatable :: ml_mon_out(:)

     logical :: monitor_qg = .false.
     character, allocatable :: qg_mon_out(:)

     double precision :: dtdiag

  end type atm_basin_config_type

  type ocn_basin_type

     logical :: active = .false.

     character, allocatable :: name(:)

     type(box_type) :: b
     type(topog_type) :: topog
     type(qg_type) :: qg
     type(ocean_mixed_type) :: ml
     type(ocn_coupler_type) :: cpl

     type(basin_covar_type) :: cov
     type(basin_averages_type) :: avg
     type(qg_monitor_type) :: qg_mon
     type(subsampling_type) :: subsamp
     type(print_type) :: print
     type(area_avg_type) :: area_avgs
     type(valids_type) :: valids
     type(oml_monitor_type) :: oml_mon

  end type ocn_basin_type

  type atm_basin_type

     logical :: active = .false.

     character, allocatable :: name(:)

     type(box_type) :: b
     type(topog_type) :: topog
     type(qg_type) :: qg
     type(atmos_mixed_type) :: ml
     type(atm_coupler_type) :: cpl

     type(basin_covar_type) :: cov
     type(basin_averages_type) :: avg
     type(qg_monitor_type) :: qg_mon
     type(subsampling_type) :: subsamp
     type(print_type) :: print
     type(area_avg_type) :: area_avgs
     type(valids_type) :: valids
     type(aml_monitor_type) :: aml_mon

  end type atm_basin_type

  public ocn_basin_type
  public atm_basin_type

  public ocn_basin_config_type
  public atm_basin_config_type
  public load_ocean_basin_config
  public load_atmos_basin_config

contains

  type(ocn_basin_config_type) function load_ocean_basin_config(filename)

    character (len=*) :: filename

    integer :: nc_id, strlen
    integer :: use_oml, flux_n, flux_s, use_rad_temp
    integer :: nl, use_qg
    character :: subnam*(*)
    parameter ( subnam = 'load_ocean_basin_config' )

    nc_id = nc_open(filename, subnam)

    strlen = nc_get_dim(nc_id, 'strlen', subnam)

    allocate(load_ocean_basin_config%name(strlen))
    call nc_get_text(nc_id, 'name', load_ocean_basin_config%name, subnam)

    load_ocean_basin_config%nx = nc_get_int(nc_id, 'nx', subnam)
    load_ocean_basin_config%ny = nc_get_int(nc_id, 'ny', subnam)
    load_ocean_basin_config%ndx = nc_get_int(nc_id, 'ndx', subnam)
    load_ocean_basin_config%ndy = nc_get_int(nc_id, 'ndy', subnam)       
    load_ocean_basin_config%nx1 = nc_get_int(nc_id, 'nx1', subnam)
    load_ocean_basin_config%ny1 = nc_get_int(nc_id, 'ny1', subnam)
    load_ocean_basin_config%nl = nc_get_dim(nc_id, 'nl', subnam)
    allocate(load_ocean_basin_config%h(load_ocean_basin_config%nl))
    load_ocean_basin_config%h(:) = nc_get_double(nc_id, 'h', load_ocean_basin_config%nl, subnam)
    load_ocean_basin_config%hm = nc_get_double(nc_id, 'hm', subnam)

    allocate(load_ocean_basin_config%topog(strlen))
    call nc_get_text(nc_id, 'topog', load_ocean_basin_config%topog, subnam)

    load_ocean_basin_config%rho = nc_get_double(nc_id, 'rho', subnam)

    use_oml = nc_get_int(nc_id, 'use_oml', subnam)
    load_ocean_basin_config%use_oml = use_oml == 1
    if (load_ocean_basin_config%use_oml) then
       flux_n = nc_get_int(nc_id, 'flux_n', subnam)
       flux_s = nc_get_int(nc_id, 'flux_s', subnam)
       load_ocean_basin_config%flux_n = flux_n == 1
       load_ocean_basin_config%flux_s = flux_s == 1

       if ( load_ocean_basin_config%flux_n .and. load_ocean_basin_config%flux_s ) then
          print *,' '
          print *,' Invalid model configuration: sb_hflux and'
          print *,' nb_hflux options cannot both be selected'
          print *,' Program terminates'
          stop 1
       else if ( load_ocean_basin_config%flux_s ) then
          print *,' Model running with modified o.m.l. southern b.c.'
       else if ( load_ocean_basin_config%flux_n ) then
          print *,' Model running with modified o.m.l. northern b.c.'
       else
          print *,' Model running with no heat flux N & S'
       endif

       load_ocean_basin_config%st2d = nc_get_double(nc_id, 'st2d', subnam)
       load_ocean_basin_config%st4d = nc_get_double(nc_id, 'st4d', subnam)
       load_ocean_basin_config%ycexp = nc_get_double(nc_id, 'ycexp', subnam)
       load_ocean_basin_config%cp = nc_get_double(nc_id, 'cp', subnam)
       
       use_rad_temp = nc_get_int(nc_id, 'use_rad_temp', subnam)
       load_ocean_basin_config%use_rad_temp = use_rad_temp == 1
       if (load_ocean_basin_config%use_rad_temp) then
          load_ocean_basin_config%T1_abs = nc_get_double(nc_id, 'T1_abs', subnam)
          load_ocean_basin_config%T2_abs = nc_get_double(nc_id, 'T2_abs', subnam)
       else
          load_ocean_basin_config%T1_rel = nc_get_double(nc_id, 'T1_rel', subnam)
          load_ocean_basin_config%T2_rel = nc_get_double(nc_id, 'T2_rel', subnam)
          load_ocean_basin_config%dT_NS = nc_get_double(nc_id, 'dT_NS', subnam)
       endif
       allocate(load_ocean_basin_config%oml_state(strlen))
       allocate(load_ocean_basin_config%fnet_cpl(strlen))
       allocate(load_ocean_basin_config%p1_cpl(strlen))
       call nc_get_text(nc_id, 'oml_state', load_ocean_basin_config%oml_state, subnam)
       call nc_get_text(nc_id, 'fnet_cpl', load_ocean_basin_config%fnet_cpl, subnam)
       call nc_get_text(nc_id, 'p1_cpl', load_ocean_basin_config%p1_cpl, subnam)
    endif

    use_qg = nc_get_int(nc_id, 'use_qg', subnam)
    load_ocean_basin_config%use_qg = use_qg == 1
    if (load_ocean_basin_config%use_qg) then
       nl = load_ocean_basin_config%nl
       load_ocean_basin_config%bcco = nc_get_double(nc_id, "bcco", subnam)
       load_ocean_basin_config%delek = nc_get_double(nc_id, "delek", subnam)
       allocate(load_ocean_basin_config%ah2(nl))
       allocate(load_ocean_basin_config%ah4(nl))
       allocate(load_ocean_basin_config%gp(nl-1))
       load_ocean_basin_config%ah2 = nc_get_double(nc_id, "ah2", nl, subnam)
       load_ocean_basin_config%ah4 = nc_get_double(nc_id, "ah4", nl, subnam)
       load_ocean_basin_config%gp = nc_get_double(nc_id, "gp", nl-1, subnam)
       allocate(load_ocean_basin_config%qg_state(strlen))
       allocate(load_ocean_basin_config%ent_cpl(strlen))
       call nc_get_text(nc_id, 'qg_state', load_ocean_basin_config%qg_state, subnam)
       call nc_get_text(nc_id, 'ent_cpl', load_ocean_basin_config%ent_cpl, subnam)
    endif

    load_ocean_basin_config%compute_covar = nc_get_int(nc_id, 'compute_covar', subnam) == 1
    if (load_ocean_basin_config%compute_covar) then
       load_ocean_basin_config%dtcov = nc_get_double(nc_id, "dtcov", subnam)
       load_ocean_basin_config%ns = nc_get_int(nc_id, "ns", subnam)
       allocate(load_ocean_basin_config%cov_out(strlen))
       call nc_get_text(nc_id, 'cov_out', load_ocean_basin_config%cov_out, subnam)
    endif

    load_ocean_basin_config%compute_avges = nc_get_int(nc_id, 'compute_avges', subnam) == 1
    if (load_ocean_basin_config%compute_avges) then
       load_ocean_basin_config%dtav = nc_get_double(nc_id, "dtav", subnam)
       allocate(load_ocean_basin_config%avg_out(strlen))
       call nc_get_text(nc_id, 'avg_out', load_ocean_basin_config%avg_out, subnam)
    endif

    load_ocean_basin_config%compute_subsample = nc_get_int(nc_id, 'compute_subsample', subnam) == 1
    if (load_ocean_basin_config%compute_subsample) then
       load_ocean_basin_config%dtsubsamp = nc_get_double(nc_id, "dtsubsamp", subnam)
       load_ocean_basin_config%nsk = nc_get_int(nc_id, "nsk", subnam)
       allocate(load_ocean_basin_config%p_file(strlen))
       allocate(load_ocean_basin_config%t_file(strlen))
       call nc_get_text(nc_id, 'p_file', load_ocean_basin_config%p_file, subnam)
       call nc_get_text(nc_id, 't_file',  load_ocean_basin_config%t_file, subnam)
       load_ocean_basin_config%outfl = nc_get_int(nc_id, 'outfl', 7, subnam) == 1
    endif

    load_ocean_basin_config%check_valids = nc_get_int(nc_id, 'check_valids', subnam) == 1
    if (load_ocean_basin_config%check_valids) then
       load_ocean_basin_config%dtvalid = nc_get_double(nc_id, "dtvalid", subnam)
       load_ocean_basin_config%max_tau = nc_get_double(nc_id, "max_tau", subnam)
       load_ocean_basin_config%max_wt = nc_get_double(nc_id, "max_wt", subnam)
       load_ocean_basin_config%max_st = nc_get_double(nc_id, "max_st", subnam)
       load_ocean_basin_config%max_p = nc_get_double(nc_id, "max_p", subnam)
       load_ocean_basin_config%max_q = nc_get_double(nc_id, "max_q", subnam)
    endif

    load_ocean_basin_config%print_summary = nc_get_int(nc_id, 'print_summary', subnam) == 1
    if (load_ocean_basin_config%print_summary) then
       load_ocean_basin_config%dtprint = nc_get_double(nc_id, "dtprint", subnam)
    endif

    load_ocean_basin_config%area_avg = nc_get_int(nc_id, 'area_avg', subnam) == 1
    if (load_ocean_basin_config%area_avg) then
       load_ocean_basin_config%narea = nc_get_dim(nc_id, 'narea', subnam)
       allocate(load_ocean_basin_config%xlo(load_ocean_basin_config%narea))
       allocate(load_ocean_basin_config%xhi(load_ocean_basin_config%narea))
       allocate(load_ocean_basin_config%ylo(load_ocean_basin_config%narea))
       allocate(load_ocean_basin_config%yhi(load_ocean_basin_config%narea))
       load_ocean_basin_config%xlo = nc_get_double(nc_id, "xlo", load_ocean_basin_config%narea, subnam)
       load_ocean_basin_config%xhi = nc_get_double(nc_id, "xhi", load_ocean_basin_config%narea, subnam)
       load_ocean_basin_config%ylo = nc_get_double(nc_id, "ylo", load_ocean_basin_config%narea, subnam)
       load_ocean_basin_config%yhi = nc_get_double(nc_id, "yhi", load_ocean_basin_config%narea, subnam)

       allocate(load_ocean_basin_config%area_filename(strlen))
       call nc_get_text(nc_id, 'area_filename', load_ocean_basin_config%area_filename, subnam)
    endif    

    load_ocean_basin_config%monitor_ml = nc_get_int(nc_id, 'monitor_ml', subnam) == 1
    if (load_ocean_basin_config%monitor_ml) then
       allocate(load_ocean_basin_config%ml_mon_out(strlen))
       call nc_get_text(nc_id, 'ml_mon_out', load_ocean_basin_config%ml_mon_out, subnam)
    endif

    load_ocean_basin_config%monitor_qg = nc_get_int(nc_id, 'monitor_qg', subnam) == 1
    if (load_ocean_basin_config%monitor_qg) then
       allocate(load_ocean_basin_config%qg_mon_out(strlen))
       call nc_get_text(nc_id, 'qg_mon_out', load_ocean_basin_config%qg_mon_out, subnam)
    endif

    load_ocean_basin_config%dtdiag = nc_get_double(nc_id, "dtdiag", subnam)

    load_ocean_basin_config%active = .true.

  end function load_ocean_basin_config

  type(atm_basin_config_type) function load_atmos_basin_config(filename)

    character (len=*) :: filename

    integer :: nc_id, strlen
    integer :: use_aml, flux_n, flux_s
    integer :: nl, use_qg
    character :: subnam*(*)
    parameter ( subnam = 'load_atmos_basin_config' )

    nc_id = nc_open(filename, subnam)

    strlen = nc_get_dim(nc_id, 'strlen', subnam)

    allocate(load_atmos_basin_config%name(strlen))
    call nc_get_text(nc_id, 'name', load_atmos_basin_config%name, subnam)

    load_atmos_basin_config%nx = nc_get_int(nc_id, 'nx', subnam)
    load_atmos_basin_config%ny = nc_get_int(nc_id, 'ny', subnam)
    load_atmos_basin_config%ndx = nc_get_int(nc_id, 'ndx', subnam)
    load_atmos_basin_config%ndy = nc_get_int(nc_id, 'ndy', subnam)       
    load_atmos_basin_config%nx1 = nc_get_int(nc_id, 'nx1', subnam)
    load_atmos_basin_config%ny1 = nc_get_int(nc_id, 'ny1', subnam)
    load_atmos_basin_config%nl = nc_get_dim(nc_id, 'nl', subnam)
    allocate(load_atmos_basin_config%h(load_atmos_basin_config%nl))
    load_atmos_basin_config%h(:) = nc_get_double(nc_id, 'h', load_atmos_basin_config%nl, subnam)
    load_atmos_basin_config%hm = nc_get_double(nc_id, 'hm', subnam)

    allocate(load_atmos_basin_config%topog(strlen))
    call nc_get_text(nc_id, 'topog', load_atmos_basin_config%topog, subnam)

    load_atmos_basin_config%rho = nc_get_double(nc_id, 'rho', subnam)

    use_aml = nc_get_int(nc_id, 'use_aml', subnam)
    load_atmos_basin_config%use_aml = use_aml == 1
    if (load_atmos_basin_config%use_aml) then
       flux_n = nc_get_int(nc_id, 'flux_n', subnam)
       flux_s = nc_get_int(nc_id, 'flux_s', subnam)
       load_atmos_basin_config%flux_n = flux_n == 1
       load_atmos_basin_config%flux_s = flux_s == 1

       if ( load_atmos_basin_config%flux_n .and. load_atmos_basin_config%flux_s ) then
          print *,' '
          print *,' Invalid model configuration: sb_hflux and'
          print *,' nb_hflux options cannot both be selected'
          print *,' Program terminates'
          stop 1
       else if ( load_atmos_basin_config%flux_s ) then
          print *,' Model running with modified a.m.l. southern b.c.'
       else if ( load_atmos_basin_config%flux_n ) then
          print *,' Model running with modified a.m.l. northern b.c.'
       else
          print *,' Model running with no heat flux N & S'
       endif

       load_atmos_basin_config%st2d = nc_get_double(nc_id, 'st2d', subnam)
       load_atmos_basin_config%st4d = nc_get_double(nc_id, 'st4d', subnam)
       load_atmos_basin_config%ahmd = nc_get_double(nc_id, 'ahmd', subnam)
       load_atmos_basin_config%cp = nc_get_double(nc_id, 'cp', subnam)
       load_atmos_basin_config%hmamin = nc_get_double(nc_id, 'hmamin', subnam)
       load_atmos_basin_config%hmadmp = nc_get_double(nc_id, 'hmadmp', subnam)
       load_atmos_basin_config%xcexp = nc_get_double(nc_id, 'xcexp', subnam)

       allocate(load_atmos_basin_config%T_abs(load_atmos_basin_config%nl))
       load_atmos_basin_config%T_abs = nc_get_double(nc_id, 'T_abs', load_atmos_basin_config%nl, subnam)
       load_atmos_basin_config%fsbar = nc_get_double(nc_id, 'fsbar', subnam)
       load_atmos_basin_config%fspamp = nc_get_double(nc_id, 'fspamp', subnam)
       load_atmos_basin_config%gamma = nc_get_double(nc_id, 'gamma', subnam)
       load_atmos_basin_config%xlamda = nc_get_double(nc_id, 'xlamda', subnam)

       allocate(load_atmos_basin_config%zopt(0:load_atmos_basin_config%nl))
       load_atmos_basin_config%zopt(0) = nc_get_double(nc_id, 'zopt_m', subnam)
       load_atmos_basin_config%zopt(1:) = nc_get_double(nc_id, 'zopt', load_atmos_basin_config%nl, 'load_rad')

       allocate(load_atmos_basin_config%aml_state(strlen))
       allocate(load_atmos_basin_config%sst_cpl(strlen))
       allocate(load_atmos_basin_config%eta_cpl(strlen))
       allocate(load_atmos_basin_config%p1_cpl(strlen))
       call nc_get_text(nc_id, 'aml_state', load_atmos_basin_config%aml_state, subnam)
       call nc_get_text(nc_id, 'sst_cpl', load_atmos_basin_config%sst_cpl, subnam)
       call nc_get_text(nc_id, 'eta_cpl', load_atmos_basin_config%eta_cpl, subnam)
       call nc_get_text(nc_id, 'p1_cpl', load_atmos_basin_config%p1_cpl, subnam)
    endif

    use_qg = nc_get_int(nc_id, 'use_qg', subnam)
    load_atmos_basin_config%use_qg = use_qg == 1
    if (load_atmos_basin_config%use_qg) then
       nl = load_atmos_basin_config%nl
       load_atmos_basin_config%bcco = nc_get_double(nc_id, "bcco", subnam)
       load_atmos_basin_config%delek = nc_get_double(nc_id, "delek", subnam)
       allocate(load_atmos_basin_config%ah2(nl))
       allocate(load_atmos_basin_config%ah4(nl))
       allocate(load_atmos_basin_config%gp(nl-1))
       load_atmos_basin_config%ah2 = nc_get_double(nc_id, "ah2", nl, subnam)
       load_atmos_basin_config%ah4 = nc_get_double(nc_id, "ah4", nl, subnam)
       load_atmos_basin_config%gp = nc_get_double(nc_id, "gp", nl-1, subnam)
       allocate(load_atmos_basin_config%qg_state(strlen))
       allocate(load_atmos_basin_config%ent_cpl(strlen))
       call nc_get_text(nc_id, 'qg_state', load_atmos_basin_config%qg_state, subnam)
       call nc_get_text(nc_id, 'ent_cpl', load_atmos_basin_config%ent_cpl, subnam)
    endif

    load_atmos_basin_config%compute_covar = nc_get_int(nc_id, 'compute_covar', subnam) == 1
    if (load_atmos_basin_config%compute_covar) then
       load_atmos_basin_config%dtcov = nc_get_double(nc_id, "dtcov", subnam)
       load_atmos_basin_config%ns = nc_get_int(nc_id, "ns", subnam)
       allocate(load_atmos_basin_config%cov_out(strlen))
       call nc_get_text(nc_id, 'cov_out', load_atmos_basin_config%cov_out, subnam)
    endif

    load_atmos_basin_config%compute_avges = nc_get_int(nc_id, 'compute_avges', subnam) == 1
    if (load_atmos_basin_config%compute_avges) then
       load_atmos_basin_config%dtav = nc_get_double(nc_id, "dtav", subnam)
       allocate(load_atmos_basin_config%avg_out(strlen))
       call nc_get_text(nc_id, 'avg_out', load_atmos_basin_config%avg_out, subnam)
    endif

    load_atmos_basin_config%compute_subsample = nc_get_int(nc_id, 'compute_subsample', subnam) == 1
    if (load_atmos_basin_config%compute_subsample) then
       load_atmos_basin_config%dtsubsamp = nc_get_double(nc_id, "dtsubsamp", subnam)
       load_atmos_basin_config%nsk = nc_get_int(nc_id, "nsk", subnam)
       allocate(load_atmos_basin_config%p_file(strlen))
       allocate(load_atmos_basin_config%t_file(strlen))
       call nc_get_text(nc_id, 'p_file', load_atmos_basin_config%p_file, subnam)
       call nc_get_text(nc_id, 't_file', load_atmos_basin_config%t_file, subnam)       
       load_atmos_basin_config%outfl = nc_get_int(nc_id, 'outfl', 7, subnam) == 1
    endif

    load_atmos_basin_config%check_valids = nc_get_int(nc_id, 'check_valids', subnam) == 1
    if (load_atmos_basin_config%check_valids) then
       load_atmos_basin_config%dtvalid = nc_get_double(nc_id, "dtvalid", subnam)
       load_atmos_basin_config%max_tau = nc_get_double(nc_id, "max_tau", subnam)
       load_atmos_basin_config%max_wt = nc_get_double(nc_id, "max_wt", subnam)
       load_atmos_basin_config%max_st = nc_get_double(nc_id, "max_st", subnam)
       load_atmos_basin_config%max_p = nc_get_double(nc_id, "max_p", subnam)
       load_atmos_basin_config%max_q = nc_get_double(nc_id, "max_q", subnam)
    endif

    load_atmos_basin_config%print_summary = nc_get_int(nc_id, 'print_summary', subnam) == 1
    if (load_atmos_basin_config%print_summary) then
       load_atmos_basin_config%dtprint = nc_get_double(nc_id, "dtprint", subnam)
    endif

    load_atmos_basin_config%area_avg = nc_get_int(nc_id, 'area_avg', subnam) == 1
    if (load_atmos_basin_config%area_avg) then
       load_atmos_basin_config%narea = nc_get_dim(nc_id, 'narea', subnam)
       allocate(load_atmos_basin_config%xlo(load_atmos_basin_config%narea))
       allocate(load_atmos_basin_config%xhi(load_atmos_basin_config%narea))
       allocate(load_atmos_basin_config%ylo(load_atmos_basin_config%narea))
       allocate(load_atmos_basin_config%yhi(load_atmos_basin_config%narea))
       load_atmos_basin_config%xlo = nc_get_double(nc_id, "xlo", load_atmos_basin_config%narea, subnam)
       load_atmos_basin_config%xhi = nc_get_double(nc_id, "xhi", load_atmos_basin_config%narea, subnam)
       load_atmos_basin_config%ylo = nc_get_double(nc_id, "ylo", load_atmos_basin_config%narea, subnam)
       load_atmos_basin_config%yhi = nc_get_double(nc_id, "yhi", load_atmos_basin_config%narea, subnam)

       allocate(load_atmos_basin_config%area_filename(strlen))
       call nc_get_text(nc_id, 'area_filename', load_atmos_basin_config%area_filename, subnam)
    endif

    load_atmos_basin_config%monitor_ml = nc_get_int(nc_id, 'monitor_ml', subnam) == 1
    if (load_atmos_basin_config%monitor_ml) then
       allocate(load_atmos_basin_config%ml_mon_out(strlen))
       call nc_get_text(nc_id, 'ml_mon_out', load_atmos_basin_config%ml_mon_out, subnam)
    endif

    load_atmos_basin_config%monitor_qg = nc_get_int(nc_id, 'monitor_qg', subnam) == 1
    if (load_atmos_basin_config%monitor_qg) then
       allocate(load_atmos_basin_config%qg_mon_out(strlen))
       call nc_get_text(nc_id, 'qg_mon_out', load_atmos_basin_config%qg_mon_out, subnam)
    endif

    load_atmos_basin_config%dtdiag = nc_get_double(nc_id, "dtdiag", subnam)

    load_atmos_basin_config%active = .true.

  end function load_atmos_basin_config

end module basin
