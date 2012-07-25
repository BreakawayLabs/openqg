module qg_monitor

  use units, only: m_to_km, m3_to_Sv
  use box, only: box_type
  use qg, only: qg_type
  use constraint, only: mass_constr_type
  use ncutils, only: nc_open_w, nc_create, nc_close, nc_def_dim, nc_def_float
  use ncutils, only: nc_def_int, nc_put_double, nc_enddef, nc_put_int
  use intsubs, only: xintp, trapin
  use ekman, only: ekman_type
  use numerics, only: map_P_to_x_face, map_P_to_y_face, dPdx, dPdy, int_P_dA

  use intsubs, only: genint

  implicit none

  private

  type qg_monitor_type

     character (len=64) :: filename
     integer :: nocmon
     logical :: active = .false.

     double precision :: wepm, wapm
     ! wepmoc, wepmat are mean Ekman velocities at p points (should be zero)
     ! wapmoc, wapmat are mean modulus of Ekman velocities at p points
     ! All the above are in (m s^-1).
     ! Computed in diagno

     double precision, allocatable :: etam(:)
     ! etamoc, etamat are mean values of eta at p points. (m)
     ! Computed in diagno

     double precision, allocatable :: ermas(:),emfr(:)

     ! ermasa, emfrat are the absolute and fractional
     ! mass errors at each atmospheric interface
     ! Computed in atinvq

     double precision, allocatable :: pavg(:)

     ! pavgoc, pavgat are average dynamic pressure of each layer (m^2 s^-2)
     ! Computed in diagno

     double precision, allocatable :: qavg(:)

     ! qavgoc, qavgat are average vorticity of each layer (s^-1)
     ! Computed in diagno

     double precision :: entm, enam

     ! entmat, enamat are the mean values of entrainment and
     ! |entrainment| for the atmosphere; entmoc, enamoc are the
     ! corresponding quantities for the ocean. All in (m s^-1).
     ! Computed in diagno

     double precision :: pken

     ! pkenoc, pkenat are part of the KE exchanges between layers (W m^-2)
     ! = rho0*gpr*Integ ( eta1*entrainment ) dA /Area
     ! Computed in diagno

     double precision, allocatable :: pket(:)
     ! pketoc, pketat are part of the KE exchanges between layers (W m^-2)
     ! = rho0*gpr*0.5d0*d/dt[ Integ ( eta(k)*eta(k) ) dA /Area ]
     ! Computed in diagno

     double precision, allocatable :: val(:)
     integer, allocatable :: pos(:)

     double precision, allocatable :: ekei(:)
     ! ekeioc, ekeiat are average KE of each layer (J m^-2)
     ! Computed in diagno


     double precision, allocatable :: ah2d(:)
     double precision, allocatable :: ah4d(:)
     ! ah2doc(k) is integrated dissipation in each layer (W m^-2)
     ! = ( rho0*Ah(k)*H(k) )
     ! * Integ ( u(k)*Del-sqd(u(k)) + v(k)*Del-sqd(v(k)) ) dA / Area
     ! Computed in diagno
     ! ah4doc(k), ah4dat(k) are integrated dissipation in each layer (W m^-2)
     ! = ( rho0*Ah(k)*H(k) )
     ! * Integ ( u(k)*Del-4th(u(k)) + v(k)*Del-4th(v(k)) ) dA / Area
     ! Computed in diagno


     double precision, allocatable :: et2m(:)
     double precision :: tsec ! We keep this to allow integrating pket from et2m
     ! et2moc, et2mat are the mean values of eta^2 at each interface
     ! Computed in diagno (usually) and constr (initially)


     double precision :: utau
     ! utauoc, utauat are the KE exchange between ocean & atmosphere (W m^-2)
     ! = rho*Integ ( v1*tauy - u1*taux ) dA / Area
     ! Computed in diagno

     double precision, allocatable :: circ(:)
     double precision :: ctot

     ! occirc is the volume transport (in Sverdrups) in each
     ! layer, and occtot is the total over all layers.
     ! Only non-zero in the cyclic ocean case
     ! Computed in diagno

     double precision :: btdg
     ! btdgoc is the mean energy dissipation due
     ! to drag in the bottom Ekman layer (W m^-2)
     ! = 0.5*rhooc*delek*|f0|*Integ ( u^2(nlo) + v^2(nlo) ) dA /Area
     ! Computed in diagno

     double precision, allocatable :: sfmin(:),sfmax(:)

     ! osfmin, osfmax are the min and max values of the ocean
     ! volume transport stream function (in Sverdrups) in each
     ! layer, where stream function psi = hoc(k)*(p - pbdy)/f0
     ! Sign convention is that in Gill 1982, opposite of usual fluid
     ! dynamics convention. Here u = -d(psi)/dy, v = d(psi)/dx
     ! Computed in diagno

  end type qg_monitor_type


  public qg_monitor_type
  public init_qg_monitor

  public qg_monnc_init
  
  public update_qg_monitor
  public qg_monitor_step

contains

  logical function qg_monitor_step(qg_mon, force, ntdone)
    type(qg_monitor_type), intent(in) :: qg_mon
    logical, intent(in) :: force
    integer, intent(in) :: ntdone

    qg_monitor_step = qg_mon%active .and. (force .or. mod(ntdone,qg_mon%nocmon) == 0)

  end function qg_monitor_step

  type(qg_monitor_type) function init_qg_monitor(b, tsec, filename, outdir, nocmon, numoutsteps, qg)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: tsec
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: outdir
    integer, intent(in) :: nocmon
    integer, intent(in) :: numoutsteps
    type(qg_type), intent(in) :: qg

    integer :: k
    double precision :: wko(b%nxp,b%nyp)

    init_qg_monitor%filename = filename

    allocate(init_qg_monitor%etam(b%nl-1))

    allocate(init_qg_monitor%ermas(b%nl-1))
    allocate(init_qg_monitor%emfr(b%nl-1))

    allocate(init_qg_monitor%pavg(b%nl))
    allocate(init_qg_monitor%qavg(b%nl))

    allocate(init_qg_monitor%val(b%nl))
    allocate(init_qg_monitor%pos(b%nl))

    allocate(init_qg_monitor%ekei(b%nl))

    allocate(init_qg_monitor%pket(b%nl-1))

    allocate(init_qg_monitor%ah2d(b%nl))
    allocate(init_qg_monitor%ah4d(b%nl))

    allocate(init_qg_monitor%et2m(b%nl-1))

    allocate(init_qg_monitor%circ(b%nl))

    allocate(init_qg_monitor%sfmin(b%nl))
    allocate(init_qg_monitor%sfmax(b%nl))

    init_qg_monitor%et2m(:) = 0.0d0
    init_qg_monitor%tsec = tsec

    do k=1,b%nl-1
       ! Compute eta^2 integral for energy diagnostics
       wko(:,:) = ( qg%p(:,:,k+1) - qg%p(:,:,k) )**2
       init_qg_monitor%et2m(k) = xintp(wko(:,:), b%nxp, b%nyp)
       init_qg_monitor%et2m(k) = init_qg_monitor%et2m(k)*b%norm/qg%gp(k)**2
    enddo
    init_qg_monitor%ermas = 0.0d0
    init_qg_monitor%emfr = 0.0d0

    call qg_monnc_init(outdir, init_qg_monitor%filename, qg, numoutsteps)
    init_qg_monitor%nocmon = nocmon
    init_qg_monitor%active = .true.

  end function init_qg_monitor

  subroutine update_qg_monitor(qg, ek, ent, tsec, qg_mon, outdir, ntdone, tyrs, dt, tdt)

    type(qg_type), intent(in) :: qg
    type(ekman_type), intent(in) :: ek
    double precision, intent(in) :: ent(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision, intent(in) :: tsec
    type(qg_monitor_type), intent(inout) :: qg_mon
    character (len=*), intent(in) :: outdir
    integer, intent(in) :: ntdone
    double precision, intent(in) :: tyrs
    double precision, intent(in) :: dt
    double precision, intent(in) :: tdt

    call qg_diagno(qg, ek, ent, tsec, tdt, qg_mon)
    call qg_monnc_out(outdir, qg, qg_mon, ntdone, tyrs)
    call scan_cfl(qg, dt)
    
  end subroutine update_qg_monitor

  subroutine qg_diagno(qg, ek, ent, tsec, tdt, mon)
    type(qg_type), intent(in) :: qg
    type(ekman_type), intent(in) :: ek
    double precision, intent(in) :: ent(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision, intent(in) :: tsec
    double precision, intent(in) :: tdt
    type(qg_monitor_type), intent(inout) :: mon

    double precision :: ptwk1(qg%b%nxp,qg%b%nyt)
    double precision :: ptwk2(qg%b%nxp,qg%b%nyt)
    double precision :: ptwk3(qg%b%nxp,qg%b%nyt)
    double precision :: ug(qg%b%nxp,qg%b%nyt)
    double precision :: tpwk1(qg%b%nxt,qg%b%nyp)
    double precision :: tpwk2(qg%b%nxt,qg%b%nyp)
    double precision :: tpwk3(qg%b%nxt,qg%b%nyp)
    double precision :: vg(qg%b%nxt,qg%b%nyp)
    double precision :: eta(qg%b%nxp,qg%b%nyp)

    integer :: j,k
    double precision :: et2now, etaint, pkeint
    double precision :: ugeos(qg%b%nxp,qg%b%nyp-1), vgeos(qg%b%nxp-1,qg%b%nyp), utaux, vtauy
    double precision :: pint, qint
    double precision :: ujet, ujeto(qg%b%nyt)
    double precision :: uke, u2diss, u4diss
    double precision :: vke, v2diss, v4diss
    double precision :: taux_t(qg%b%nxp-1,qg%b%nyp),tauy_t(qg%b%nxp,qg%b%nyp-1)
    double precision :: dp_dy(qg%b%nxp,qg%b%nyp-1)
    double precision :: dp_dx(qg%b%nxp-1,qg%b%nyp)
    double precision :: dpm_dy(qg%b%nxp,qg%b%nyp-1)
    double precision :: dpm_dx(qg%b%nxp-1,qg%b%nyp)

    type(box_type) :: b
    double precision pomin,pomax,poref(qg%b%nl),psiext

    double precision :: int_p2(qg%b%nl)

    b = qg%b

    ! Ekman velocity diagnostics
    ! Mean value of Wekman at p points
    mon%wepm = genint(ek%wekp, b%nxp, b%nyp, 0.5d0, 0.5d0)
    mon%wepm = mon%wepm*b%norm
    ! Mean value of abs( Wekman ) at p points
    mon%wapm = genint(abs(ek%wekp(:,:)), b%nxp, b%nyp, 0.5d0, 0.5d0)
    mon%wapm = mon%wapm*b%norm
    
    ! Entrainment diagnostics
    ! We are assuming all entrainment is across interface 1.
    mon%entm = genint(ent(:,:,1),b%nxp,b%nyp,0.5d0,0.5d0)
    mon%entm = mon%entm*b%norm
    mon%enam = genint(abs(ent(:,:,1)), b%nxp, b%nyp, 0.5d0, 0.5d0)
    mon%enam = mon%enam*b%norm

    ! Interface displacement eta diagnostics
    ! Infer eta at p points
    do k=1,b%nl-1
       eta(:,:) = b%dz_sign*( qg%p(:,:,k) - qg%p(:,:,k+1) )/qg%gp(k)
       ! Compute integral of eta
       etaint = genint(eta, b%nxp, b%nyp, 0.5d0, 0.5d0)
       mon%etam(k) = etaint*b%norm    
       ! Compute integral of eta^2
       et2now = genint(eta(:,:)*eta(:,:), b%nxp, b%nyp, 0.5d0, 0.5d0)
       et2now = et2now*b%norm
       if( tsec > mon%tsec) then
          mon%pket(k) = 0.5d0*qg%rho*qg%gp(k)*( et2now - mon%et2m(k) )/(tsec - mon%tsec)
       else
          mon%pket(k) = 0.0d0
       endif
       mon%et2m(k) = et2now
       mon%tsec = tsec
       ! Compute integral of eta*e
       ! We are assuming all entrainment is across interface 1.
       if (k == 1) then
          pkeint = genint(eta(:,:)*ent(:,:,1), b%nxp, b%nyp, 0.5d0, 0.5d0)
          mon%pken = qg%rho*qg%gp(1)*pkeint*b%norm
       endif
    enddo

    ! energetics

    ! KE exchange between ocean & atmosphere
    ! Infer u1 geostrophically
    ! Compute required integrand

    call map_P_to_x_face(ek%taux, b, taux_t)
    call map_P_to_y_face(ek%tauy, b, tauy_t)
    call dPdx(qg%p(:,:,1), b, dp_dx)
    call dPdy(qg%p(:,:,1), b, dp_dy)
    ptwk1(:,:) = -tauy_t(:,:)*dp_dy(:,:)/b%fnot
    utaux = genint(ptwk1, b%nxp, b%nyt, 0.5d0, 1.0d0)
    ! Infer v1 geostrophically
    ! Compute required integrand
    tpwk1(:,:) = taux_t(:,:)*dp_dx(:,:)/b%fnot
    vtauy = genint(tpwk1, b%nxt, b%nyp, 1.0d0, 0.5d0)
    mon%utau = qg%rho*( vtauy + utaux )*b%norm
       
    do k=1,b%nl
       pint = genint(qg%p(:,:,k), b%nxp, b%nyp, 0.5d0, 0.5d0)
       qint = genint(qg%q(:,:,k), b%nxp, b%nyp, 0.5d0, 0.5d0)
       mon%pavg(k) = pint*b%norm
       mon%qavg(k) = qint*b%norm
    enddo

    ! Layer kinetic energy and its time derivative, also stream function
    do k=1,b%nl
       call dPdx(qg%pm(:,:,k), b, dpm_dx)
       call dPdy(qg%pm(:,:,k), b, dpm_dy)

       ! Infer lagged u and v geostrophically; derive Del-4th
       ug(:,:) = -dpm_dy(:,:)/b%fnot
       ! Differentiate contents of ug to get Del-4th(lagged u)
       ! Return this field in array ptwk3. ptwk2 is workspace
       ! which contains Del-sqd(lagged u)
       call del4_ (ug, b%nxp, b%nyt, b%cyclic, b%dxm2, ptwk2, ptwk3)
          
       vg(:,:) = dpm_dx(:,:)/b%fnot
       ! Differentiate contents of vg to get Del-4th(lagged v)
       ! Return this field in array octwk3. octwk2 is workspace
       ! which contains Del-sqd(lagged v)
       call del4_ (vg, b%nxt, b%nyp, b%cyclic, b%dxm2, tpwk2, tpwk3)
       ! Find extrema of po for stream function
       ! Infer current u and v geostrophically
       ! Compute required integrands and jet position
       call dPdy(qg%p(:,:,k), b, dp_dy)
       ugeos(:,:) = -dp_dy(:,:)/b%fnot
       ! Compute required integrands
       call dPdx(qg%p(:,:,k), b, dp_dx)
       vgeos(:,:) = dp_dx(:,:)/b%fnot

       u2diss = genint(ugeos(:,:)*ptwk2(:,:), b%nxp, b%nyt, 0.5d0, 1.0d0)
       v2diss = genint(vgeos(:,:)*tpwk2(:,:), b%nxt, b%nyp, 1.0d0, 0.5d0)
       mon%ah2d(k) = -qg%rho*qg%ah2(k)*b%h(k)*(u2diss + v2diss)*b%norm

       u4diss = genint (ugeos(:,:)*ptwk3(:,:), b%nxp, b%nyt, 0.5d0, 1.0d0)
       v4diss = genint (vgeos(:,:)*tpwk3(:,:), b%nxt, b%nyp, 1.0d0, 0.5d0)
       mon%ah4d(k) =  qg%rho*qg%ah4(k)*b%h(k)*(u4diss + v4diss)*b%norm

       uke = genint(ugeos(:,:)*ugeos(:,:), b%nxp, b%nyt, 0.5d0, 1.0d0)
       vke = genint(vgeos(:,:)*vgeos(:,:), b%nxt, b%nyp, 1.0d0, 0.5d0)
       mon%ekei(k) = 0.5d0*qg%rho*b%h(k)*( uke + vke )*b%norm

       do j=1,b%nyt
          ujet = trapin(ugeos(:,j), b%nxp, 1.0d0)
          ujeto(j) = abs( ujet )/dble(b%nxt)
       enddo
       ! Find position and value of avged max in each layer
       mon%pos(k) = maxloc(ujeto(:), dim=1)
       mon%val(k) = maxval(ujeto(:))

       ! Compute total ocean circulation
       mon%circ(k) = m3_to_Sv(( qg%p(1,1,k) - qg%p(1,b%nyp,k) )/b%fnot)
    enddo
    mon%ctot = sum(mon%circ(:))

    ! Drag in Ekman layer at ocean bottom
    ! Infer lagged u geostrophically
    ! Compute required integrand
    call dPdx(qg%pm(:,:,qg%topo%k_topo), b, dpm_dx)
    call dPdy(qg%pm(:,:,qg%topo%k_topo), b, dpm_dy)

    ptwk1(:,:) = dpm_dy(:,:)**2/(b%fnot**2)
    tpwk1(:,:) = dpm_dx(:,:)**2/(b%fnot**2)
    u2diss = genint(ptwk1, b%nxp, b%nyt, 0.5d0, 1.0d0)
    v2diss = genint(tpwk1, b%nxt, b%nyp, 1.0d0, 0.5d0)
    mon%btdg = 0.5d0*qg%rho*qg%delek*abs(b%fnot)*( u2diss + v2diss )*b%norm

    ! Define reference pressures on equatorward side of domain
    if (b%fnot > 0.0d0) then
       do k=1,b%nl
          poref(k) = qg%p(1,1,k)
       enddo
    else if (b%fnot < 0.0d0) then
       do k=1,b%nl
          poref(k) = qg%p(1,b%nyp,k)
       enddo
    endif

    do k=1,b%nl
       ! Convert dynamic pressure extrema to volume
       ! transport stream function values (in Sverdrups)
       pomin = minval(qg%p(:,:,k))
       pomax = maxval(qg%p(:,:,k))
       
       psiext = min( pomin/b%fnot, pomax/b%fnot )
       mon%sfmin(k) = m3_to_Sv(b%h(k)*( psiext - poref(k)/b%fnot ))
       psiext = max( pomin/b%fnot, pomax/b%fnot )
       mon%sfmax(k) = m3_to_Sv(b%h(k)*( psiext - poref(k)/b%fnot ))
       ! Compute volume transport in each layer (in Sverdrups)
    enddo

    if (qg%b%cyclic) then
       int_p2(:) = int_P_dA(qg%p, qg%b)
       call check_continuity(qg%con, qg%gp, qg%b, tdt, int_p2, mon)
    endif

  end subroutine qg_diagno

  function d2dx2(arr, nx, ny, cyclic)

    double precision, intent(in) :: arr(nx,ny)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    logical, intent(in) :: cyclic
    double precision :: d2dx2(nx,ny)

    integer :: i

    do i=2,nx-1
       d2dx2(i,:) = arr(i-1,:) - 2.0d0*arr(i,:) + arr(i+1,:)
    enddo

    if (cyclic) then
       d2dx2(1,:) = arr(nx,:) - 2.0d0*arr(1,:) + arr(2,:)
       d2dx2(nx,:) = arr(nx-1,:) - 2.0d0*arr(nx,:) + arr(1,:)
    else
       d2dx2(1,:) = d2dx2(2,:)
       d2dx2(nx,:) = d2dx2(nx-1,:)
    endif

  end function d2dx2

  function d2dy2(arr, nx, ny)
    double precision, intent(in) :: arr(nx,ny)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    double precision :: d2dy2(nx,ny)

    integer :: j

    do j=2,ny-2
       d2dy2(:,j) = arr(:,j-1) - 2.0d0*arr(:,j) + arr(:,j+1)
    enddo
    d2dy2(:,1) = d2dy2(:,2)
    d2dy2(:,ny) = d2dy2(:,ny-1)

  end function d2dy2

  subroutine del4_ (arr, nx, ny, cyclic, dxm2, del2, del4)

    ! Computes Del-4th(arr) for a field arr(nx,ny) tabulated
    ! in a finite box or periodic box. No boundary condition is applied;
    ! the boundary values are computed using the appropriate
    ! forward or backward difference formulae.
    ! del2(nx,ny) is a workspace array used for Del-sqd(arr);
    ! Del-4th(arr) is returned in the array del4(nx,ny).
    ! arr is returned unchanged.
     
    integer, intent(in) :: nx,ny
    double precision, intent(in) :: arr(nx,ny)
    double precision, intent(in) ::dxm2
    logical, intent(in) :: cyclic
    double precision, intent(out) :: del2(nx,ny),del4(nx,ny)

    double precision :: d2arr_dx2(nx,ny)
    double precision :: d2arr_dy2(nx,ny)
    double precision :: d2del2_dx2(nx,ny)
    double precision :: d2del2_dy2(nx,ny)

    d2arr_dx2(:,:) = d2dx2(arr, nx, ny, cyclic)
    d2arr_dy2(:,:) = d2dy2(arr, nx, ny)

    del2(:,:) = dxm2*(d2arr_dx2(:,:) + d2arr_dy2(:,:))

    d2del2_dx2(:,:) = d2dx2(del2, nx, ny, cyclic)
    d2del2_dy2(:,:) = d2dy2(del2, nx, ny)

    del4(:,:) = dxm2*(d2del2_dx2(:,:) + d2del2_dy2(:,:))

  end subroutine del4_


  subroutine qg_monnc_init(outdir, filename, qg, numoutsteps)
    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename
    type(qg_type), intent(in) :: qg
    integer, intent(in) :: numoutsteps

    integer :: ncid

    character :: subnam*(*)
    parameter ( subnam = 'qg_monnc_init' )

    integer :: timedim, ldim, lmdim
    integer :: k
    double precision :: tmp(qg%b%nl)

    ncid = nc_create(outdir, filename, subnam)

    timedim = nc_def_dim(ncid, 'time', numoutsteps, subnam)
    ldim = nc_def_dim(ncid, 'z', qg%b%nl, subnam)
    lmdim = nc_def_dim(ncid, 'zm', qg%b%nl-1, subnam)
    
    call nc_def_float(ncid, 'time', timedim, 'years', subnam, 'Time')

    call nc_def_float(ncid, 'z', ldim, 'km', subnam, 'mid-layer depth axis')
    call nc_def_float(ncid, 'zm', lmdim, 'km', subnam, 'interface depth axis')

    call nc_def_float(ncid, 'wepm', timedim, 'm/s', subnam, 'Average Ekman velocity (p-grid)')
    call nc_def_float(ncid, 'wapm', timedim, 'm/s', subnam, 'Absolute Ekman velocity (p-grid)')
    call nc_def_float(ncid, 'etam', (/lmdim,timedim/), 'm', subnam, 'Average interface height')
    if (qg%b%cyclic) then
       call nc_def_float(ncid, 'ermas', (/lmdim,timedim/), 'm^4/s^2', subnam, 'Sum of mass error')
       call nc_def_float(ncid, 'emfr', (/lmdim,timedim/), ' ', subnam, 'Fractional mass error')
    endif
    call nc_def_float(ncid, 'pavg', (/ldim,timedim/), 'm^2 s^-2', subnam, 'Average dynamic pressure')
    call nc_def_float(ncid, 'qavg', (/ldim,timedim/), 's^-1', subnam, 'Average vorticity')
    call nc_def_float(ncid, 'entm', timedim, 'm/s', subnam, 'Average entrainment')
    call nc_def_float(ncid, 'enam', timedim, 'm/s', subnam, 'Absolute entrainment')
    call nc_def_float(ncid, 'pken', timedim, 'W/m^2', subnam, 'Average eta*e')
    call nc_def_float(ncid, 'pket', (/lmdim,timedim/), 'W/m^2', subnam, 'Average eta*eta_t')
    call nc_def_int(ncid,'pos',(/ldim,timedim/), 'gridsquare', subnam, 'Position of maximum velocity (gridsquare)')
    call nc_def_float(ncid, 'val', (/ldim,timedim/), 'm/s', subnam, 'Maximum velocity')
    call nc_def_float(ncid, 'ekei', (/ldim,timedim/), 'J/m^2', subnam, 'Average kinetic energy')
    call nc_def_float(ncid, 'ah2d', (/ldim,timedim/), 'W/m^2', subnam, 'Average dissipation (sqd)')
    call nc_def_float(ncid, 'ah4d', (/ldim,timedim/), 'W/m^2', subnam, 'Average dissipation (4th)')
    call nc_def_float(ncid, 'et2m', (/lmdim, timedim/), 'W/m^2', subnam, 'Average eta^2')
    call nc_def_float(ncid, 'utau', timedim, 'W/m^2', subnam, 'Average p1*we')
    call nc_def_float(ncid, 'circ', (/ldim,timedim/), 'Sv', subnam, 'Transport in layer')
    call nc_def_float(ncid, 'ctot', timedim, 'Sv', subnam, 'Total transport')
    call nc_def_float(ncid, 'btdg', timedim, 'W/m^2', subnam, 'Average bottom drag')
    call nc_def_float(ncid, 'sfmin', (/ldim,timedim/), 'Sv', subnam, 'Minimum streamfunction')
    call nc_def_float(ncid, 'sfmax', (/ldim,timedim/), 'Sv', subnam, 'Maximum streamfunction')


    ! Leave definition mode: entering data mode.
    call nc_enddef(ncid, subnam)

    ! Convert height into km and store in 'zo'
    tmp(1) = 0.5d0*qg%b%h(1)
    do k=2,qg%b%nl
       tmp(k) = tmp(k-1) + 0.5d0*(qg%b%h(k-1) + qg%b%h(k))
    enddo
    call nc_put_double(ncid, 'z', m_to_km(tmp), subnam)
    ! Convert height into km and store in 'zom'
    tmp(1) = qg%b%h(1)
    do k=2,qg%b%nl-1
       tmp(k) = tmp(k-1) + qg%b%h(k)
    enddo
    call nc_put_double(ncid, 'zm', m_to_km(tmp(:qg%b%nl-1)), subnam)

    call nc_close(ncid)

  end subroutine qg_monnc_init

  subroutine qg_monnc_out(outdir, qg, qg_mon, ntdone, tyrs)

    character (len=*), intent(in) :: outdir
    type(qg_type), intent(in) :: qg
    type(qg_monitor_type), intent(in) :: qg_mon
    integer, intent(in) :: ntdone
    double precision, intent(in) :: tyrs

    character :: subnam*(*)
    parameter ( subnam = 'qg_monnc_out' )

    integer :: ncid, startt

    ncid = nc_open_w(outdir, qg_mon%filename, subnam)

    ! Store current time as part of 'time' vector
    startt = ntdone/qg_mon%nocmon+1
    call nc_put_double(ncid, 'time', startt, tyrs, subnam)   
    call nc_put_double(ncid, 'wepm', startt, qg_mon%wepm, subnam)
    call nc_put_double(ncid, 'wapm', startt, qg_mon%wapm, subnam)
    call nc_put_double(ncid, 'etam', startt, qg_mon%etam, subnam)
    if (qg%b%cyclic) then
       call nc_put_double(ncid, 'ermas', startt, qg_mon%ermas, subnam)
       call nc_put_double(ncid, 'emfr', startt, qg_mon%emfr, subnam)
    endif
    call nc_put_double(ncid, 'pavg', startt, qg_mon%pavg, subnam)
    call nc_put_double(ncid, 'qavg', startt, qg_mon%qavg, subnam)
    call nc_put_double(ncid, 'entm', startt, qg_mon%entm, subnam)
    call nc_put_double(ncid, 'enam', startt, qg_mon%enam, subnam)
    call nc_put_double(ncid, 'pken', startt, qg_mon%pken, subnam)
    call nc_put_double(ncid, 'pket', startt, qg_mon%pket, subnam)
    call nc_put_int(ncid,'pos',(/1, startt/),(/qg%b%nl, 1/),qg_mon%pos,subnam)
    call nc_put_double(ncid, 'val',startt, qg_mon%val, subnam)
    call nc_put_double(ncid, 'ekei', startt, qg_mon%ekei, subnam)
    call nc_put_double(ncid, 'ah2d', startt, qg_mon%ah2d, subnam)
    call nc_put_double(ncid, 'ah4d', startt, qg_mon%ah4d, subnam)
    call nc_put_double(ncid, 'et2m', startt, qg_mon%et2m, subnam)
    call nc_put_double(ncid, 'utau', startt, qg_mon%utau, subnam)
    call nc_put_double(ncid, 'circ', startt, qg_mon%circ, subnam)
    call nc_put_double(ncid, 'ctot', startt, qg_mon%ctot, subnam)
    call nc_put_double(ncid, 'btdg', startt, qg_mon%btdg, subnam)
    call nc_put_double(ncid, 'sfmin', startt, qg_mon%sfmin, subnam)
    call nc_put_double(ncid, 'sfmax', startt, qg_mon%sfmax, subnam)

    call nc_close(ncid)

  end subroutine qg_monnc_out

  subroutine scan_cfl(qg, dt)

    ! Check how well the current state satisfies the CFL criterion.

    type(qg_type), intent(in) :: qg
    double precision, intent(in) :: dt
    
    integer :: i, j, k
    double precision :: uabs, vabs

    double precision :: dp_dy(qg%b%nxp,qg%b%nyp-1,qg%b%nl)
    double precision :: dp_dx(qg%b%nxp-1,qg%b%nyp,qg%b%nl)
    double precision :: ugeos(qg%b%nxp,qg%b%nyp-1,qg%b%nl)
    double precision :: vgeos(qg%b%nxp-1,qg%b%nyp,qg%b%nl)

    double precision :: umax(qg%b%nl),vmax(qg%b%nl)

    double precision :: cflcrit
    parameter ( cflcrit=0.8d0 )

    ! Check extreme velocity components and derive Courant numbers
    ! Infer u, v geostrophically
    do k=1,qg%b%nl
       call dPdx(qg%p(:,:,k), qg%b, dp_dx(:,:,k))
       call dPdy(qg%p(:,:,k), qg%b, dp_dy(:,:,k))
       ugeos(:,:,k) = -dp_dy(:,:,k)/qg%b%fnot
       vgeos(:,:,k) =  dp_dx(:,:,k)/qg%b%fnot
       umax(k) = maxval(abs(ugeos(:,:,k)))
       vmax(k) = maxval(abs(vgeos(:,:,k)))
    enddo

    print 217, '  max |u|(k) = ',(umax(k),k=1,qg%b%nl)
    print 217, '  CFL |u|(k) = ',(umax(k)*dt/qg%b%dx,k=1,qg%b%nl)
    print 217, '  max |v|(k) = ',(vmax(k),k=1,qg%b%nl)
    print 217, '  CFL |v|(k) = ',(vmax(k)*dt/qg%b%dx,k=1,qg%b%nl)
    ! If have bad CFL values, scan for and print locations
    do k=1,qg%b%nl
       if (umax(k)*dt/qg%b%dx >= cflcrit) then
          do j=1,qg%b%nyp-1
             do i=1,qg%b%nxp
                uabs = (dt/qg%b%dx)*abs(ugeos(i,j,k))
                if (uabs >= cflcrit) then
                   print 250, '  Bad |u|; CFL, i, j, k = ',uabs,i,j,k
                endif
             enddo
          enddo
       endif
       if (vmax(k)*dt/qg%b%dx >= cflcrit) then
          do j=1,qg%b%nyp
             do i=1,qg%b%nxp-1
                vabs = (dt/qg%b%dx)*abs(vgeos(i,j,k))
                if (vabs >= cflcrit) then
                   print 250, '  Bad |v|; CFL, i, j, k = ',vabs,i,j,k
                endif
             enddo
          enddo
       endif
    enddo

217 format(a,1p,9d15.7)
250 format(a,1p,d15.7,2i8,i6)
    
  end subroutine scan_cfl

  pure subroutine check_continuity(con, gp, b, tdt, int_p, qg_mon)
    ! Check continuity is satisfied at each interface
    ! Update continuity measures at each interface
    type(mass_constr_type), intent(in) :: con
    type(box_type), intent(in) :: b
    double precision, intent(in) :: gp(b%nl-1)
    double precision, intent(in) :: tdt
    double precision, intent(in) :: int_p(b%nl)
    type(qg_monitor_type), intent(inout) :: qg_mon

    integer :: k
    double precision :: est1, est2, edif, esum

    double precision ecrit
    parameter ( ecrit=1.0d-13 )

    do k=1,b%nl-1
       ! Choose sign of dpioc so that +ve dpioc -> +ve eta
       ! Check continuity is satisfied at each interface
       ! MONITORING - extra section for ermaso, emfroc
       ! Compute alternative estimates of new dpioc
       est1 = b%dz_sign*(int_p(k) - int_p(k+1))
       est2 = con%dpi(k)
       edif = est1 - est2
       esum = abs(est1) + abs(est2)
       qg_mon%ermas(k) = edif
       ! Compute fractional error if entrainment is significant;
       ! fraction is meaningless if est1, est2 just noisy zeros
       if (esum > (ecrit*b%xl*b%yl*tdt*gp(k))) then
          qg_mon%emfr(k) = 2.0d0*edif/esum
       else
          qg_mon%emfr(k) = 0.0d0
       endif
    enddo

  end subroutine check_continuity


end module qg_monitor
