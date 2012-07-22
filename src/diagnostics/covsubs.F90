module covsubs

  use mixed, only: mixed_type
  use qg, only: qg_type
  use ncutils, only: nc_enddef, nc_def_dim, nc_create, nc_close
  use ncutils, only: nc_def_float, nc_put_double, nc_open, nc_get_int
  use clock, only: clock_type, days_to_steps
  use box, only: box_type

  implicit none

  private

  type basin_covar_type

     logical :: active = .false.

     integer :: ns, nv, nc

     double precision, allocatable :: covp(:), avgp(:)
     double precision, allocatable :: covt(:), avgt(:)
     double precision :: swtp, swtt

     integer :: nup, nut

     ! matrices of covariance, vectors of averages
     ! and sums of weights at subsampled points, for ocean data

     integer :: ntcov

     character (len=64) :: filename

  end type basin_covar_type

  public basin_covar_type
  public init_basin_covar
  public update_cov
  public covar_step
  public covout

contains

  type(basin_covar_type) function init_basin_covar(b, dtcov, clk, filename, ns)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: dtcov
    type(clock_type), intent(in) :: clk
    character (len=*), intent(in) :: filename
    integer, intent(in) :: ns

    init_basin_covar%active = dtcov > 0.0d0
    if (init_basin_covar%active) then
       init_basin_covar%filename = filename
       init_basin_covar%ns = ns
       call finish_loading_cov(dtcov, b, clk, init_basin_covar)
    endif

  end function init_basin_covar

  subroutine finish_loading_cov(dt, b, clk, cov)

    double precision, intent(in) :: dt
    type(box_type), intent(in) :: b
    type(clock_type), intent(in) :: clk
    type(basin_covar_type), intent(inout) :: cov

    if (cov%active) then
       cov%nv = b%nxt*b%nyt/(cov%ns*cov%ns)
       cov%nc = cov%nv*(cov%nv+1)/2    
       cov%ntcov = days_to_steps(dt, clk)
       allocate(cov%covp(cov%nc))
       allocate(cov%avgp(cov%nv))
       allocate(cov%covt(cov%nc))
       allocate(cov%avgt(cov%nv))

       cov%nup = 0
       cov%nut = 0
       cov%swtp = 0.0d0
       cov%swtt = 0.0d0

       cov%covp(:) = 0.0d0
       cov%covt(:) = 0.0d0
       cov%avgp(:) = 0.0d0
       cov%avgt(:) = 0.0d0

       if (mod(b%nxt,cov%ns) /= 0 .or. mod(b%nyt,cov%ns) /= 0) then
          print *,' '
          print *,' nxt, nyt, ns = ',b%nxt,b%nyt,cov%ns
          print *,' Need nxt, nyt both to be integer multiples of'
          print *,' ns for covariance subsampling to work correctly'
          print *,' Program terminates'
          stop 1
       endif
    endif

  end subroutine finish_loading_cov

  logical pure function covar_step(nt, cov)

    integer, intent(in) :: nt
    type(basin_covar_type), intent(in) :: cov
    
    if (cov%active .and. nt > 0) then
       covar_step = mod(nt,cov%ntcov) == 0
    else
       covar_step = .false.
    endif

  end function covar_step

  subroutine update_cov(st, qg, cov)

    ! Accumulates time covariance matrices during a run

    type(mixed_type), intent(in) :: st
    type(qg_type), intent(in) :: qg
    type(basin_covar_type), intent(inout) :: cov

    integer :: ifault
    double precision :: u(cov%nv)

    if (.not. cov%active) stop 1

    ! Spatially subsample array p into vector u
    call psampl (qg%p(:,:,1), qg%b%nxp, qg%b%nyp, cov%ns, u, cov%nv)
    ! Update covariance matrix
    call dssp (u, cov%avgp, cov%covp, 1.0d0, cov%swtp, cov%nv, cov%nup, ifault)
    if (ifault /= 0) then
       print *,' dssp problem in update_cov on p; ifault = ',ifault
    endif

    if (st%active) then
       ! Spatially subsample array st into vector u
       call tsampl (st%data, qg%b%nxt, qg%b%nyt, cov%ns, u, cov%nv)
       ! Update covariance matrix
       call dssp (u, cov%avgt, cov%covt, 1.0d0, cov%swtt, cov%nv, cov%nut, ifault)
       if (ifault /= 0) then
          print *,' dssp problem in update_cov on st; ifault = ',ifault
       endif
    endif

  end subroutine update_cov

  subroutine covout(outdir, cov)

    character (len=*), intent(in) :: outdir
    type(basin_covar_type), intent(in) :: cov

    character :: subnam*(*)
    parameter ( subnam = 'covout' )

    integer :: idim
    integer :: ncdim,nvdim,covp_id,covt_id
    integer :: avgp_id,avgt_id,swtp_id,swtt_id
    integer :: covncid

    if (cov%active) then
       covncid = nc_create(outdir, cov%filename, subnam)

       ! now initialise netcdf file
       idim = nc_def_dim(covncid, 's', 1, subnam)

       ncdim = nc_def_dim(covncid, 'nc', cov%nc, subnam)
       nvdim = nc_def_dim(covncid, 'nv', cov%nv, subnam)
       covp_id = nc_def_float(covncid, 'covp', ncdim, ' ', subnam)
       covt_id = nc_def_float(covncid, 'covt', ncdim, ' ', subnam)
       avgp_id = nc_def_float(covncid, 'avgp', ncdim, ' ', subnam)
       avgt_id = nc_def_float(covncid, 'avgt', ncdim, ' ', subnam)
       swtp_id = nc_def_float(covncid, 'swtp', ncdim, ' ', subnam)
       swtt_id = nc_def_float(covncid, 'swtt', ncdim, ' ', subnam)

       !! Leave definition mode: entering data mode.
       call nc_enddef(covncid, subnam)

       !! Write data to arrays
       call nc_put_double(covncid, covp_id, cov%covp, subnam)
       call nc_put_double(covncid, covt_id, cov%covt, subnam)
       call nc_put_double(covncid, avgp_id, cov%avgp, subnam)
       call nc_put_double(covncid, avgt_id, cov%avgt, subnam)
       call nc_put_double(covncid, swtp_id, (/cov%swtp/), subnam)
       call nc_put_double(covncid, swtt_id, (/cov%swtt/), subnam)

       call nc_close(covncid, subnam)
    endif

  end subroutine covout


  subroutine tsampl (datat, nx, ny, nsi, dsamt, nsvec)

    ! Given an array datat(nx,ny) of data tabulated at T points, and
    ! a subsampling interval nsi, computes the values of dsamt(nsvec),
    ! a vector containing the required subsample of the data
    ! This will work for both channel and box configurations

    integer, intent(in) :: nx,ny,nsi,nsvec
    double precision, intent(in) :: datat(nx,ny)
    double precision, intent(out) :: dsamt(nsvec)

    integer :: is,js,ivs,id,id1,id2,jd,jd1,jd2
    double precision :: sumd,sumi
    ivs = 0
    if (nsi > 1) then
       ! Subsample the data by computing the arithmetic mean
       do js=1,(ny/nsi)
          jd1 = 1 + (js-1)*nsi
          jd2 = js*nsi
          do is=1,(nx/nsi)
             id1 = 1 + (is-1)*nsi
             id2 = is*nsi
             sumd = 0.0d0
             do jd=jd1,jd2
                sumi = 0.0d0
                do id=id1,id2
                   sumi = sumi + datat(id,jd)
                enddo
                sumd = sumd + sumi
             enddo
             ivs = (js-1)*(nx/nsi) + is
             dsamt(ivs) = sumd
          enddo
       enddo
    else if (nsi == 1) then
       ! Just copy the data into the subsample array
       do jd=1,ny
          do id=1,nx
             ivs = (jd-1)*nx + id
             dsamt(ivs) = datat(id,jd)
          enddo
       enddo
    else
       print *,' WARNING: tsampl called with invalid nsi = ',nsi
    endif
    ivs = (nx/nsi)*(ny/nsi)
    if (ivs /= nsvec) then
       print *,' WARNING: inconsistent subsample vector in tsampl'
       print *,' nsvec, final ivs = ',nsvec,ivs
    endif

  end subroutine tsampl

  subroutine psampl (datap, nx, ny, nsi, dsamp, nsvec)

    ! Given an array datap(nx,ny) of data tabulated at p points, and
    ! a subsampling interval nsi, computes the values of dsamp(nvec),
    ! a vector containing the required subsample of data

    integer, intent(in) :: nx,ny,nsi,nsvec
    double precision, intent(in) :: datap(nx,ny)
    double precision, intent(out) :: dsamp(nsvec)

    integer :: is,js,ivs,id,id1,id2,jd,jd1,jd2
    double precision :: sumd,sumi,sums,sumn
    ivs = 0
    do js=1,(ny/nsi)
       jd1 = 1 + (js-1)*nsi
       jd2 = 1 + js*nsi
       do is=1,(nx/nsi)
          id1 = 1 + (is-1)*nsi
          id2 = 1 + is*nsi
          ! Sum over inner rows of averaging area
          sumd = 0.0d0
          do jd=jd1+1,jd2-1
             sumi = 0.5d0*datap(id1,jd)
             do id=id1+1,id2-1
                sumi = sumi + datap(id,jd)
             enddo
             sumi = sumi + 0.5d0*datap(id2,jd)
             sumd = sumd + sumi
          enddo
          ! Contributions from S & N boundaries of averaging area
          sums = 0.5d0*datap(id1,jd1)
          sumn = 0.5d0*datap(id1,jd2)
          do id=id1+1,id2-1
             sums = sums + datap(id,jd1)
             sumn = sumn + datap(id,jd2)
          enddo
          sums = sums + 0.5d0*datap(id2,jd1)
          sumn = sumn + 0.5d0*datap(id2,jd2)
          ivs = (js-1)*(nx/nsi) + is
          dsamp(ivs) = sumd + 0.5d0*( sums + sumn )
       enddo
    enddo
    ivs = (nx/nsi)*(ny/nsi)
    if (ivs /= nsvec) then
       print *,' WARNING: inconsistent subsample vector in psampl'
       print *,' nsvec, final ivs = ',nsvec,ivs
    endif

  end subroutine psampl

  subroutine dssp (x, xmean, xssp, wt, sumwt, nvar, nunit, ifault)

    ! Algorithm AS 41 j.r.statist.soc.c. (1971) vol. 20 no.2
    ! This subroutine updates the mean vector xmean (length nvar)
    ! and the matrix of corrected sums of squares and products xssp
    ! (length nvar(nvar+1)/2, stored by lower triangle), when a
    ! data vector x (length nvar) with weight wt is either included
    ! (wt > 0) or excluded (wt < 0).  sumwt is the current sum of
    ! weights on entry and the updated sum on exit and nunit is
    ! the current and updated sample size.  ifault=0 indicates normal
    ! exit,  ifault=1 indicates zero or negative value of sumwt,
    ! ifault=2 indicates zero or negative nunit, ifault=3 indicates
    ! nvar < 1.  Note that x, xmean, xssp, wt and sumwt are double
    ! precision and must be declared as such in the calling program.

    double precision, intent(inout) ::  x(*)
    double precision, intent(inout) :: xmean(*)
    double precision, intent(inout) :: xssp(*)
    double precision, intent(in) ::  wt
    double precision, intent(inout) :: sumwt
    integer, intent(in) ::  nvar
    integer, intent(inout) :: nunit
    integer, intent(out) :: ifault

    double precision :: co
    parameter ( co=0.0d0 )

    integer :: i,j,k
    double precision :: b, c

    ! Check variates, weights and sample size
    ifault = 0
    if (nvar < 1) then
       ifault = 3
       return
    endif
    if (wt < 0.0d0) then
       nunit = nunit-1
    else if (wt == 0.0d0) then
       return
    else if (wt > 0.0d0) then
       nunit = nunit+1
    endif
    sumwt = sumwt+wt
    if (sumwt <= co) then
       ifault = 1
       return
    endif
    b = wt/sumwt
    if (nunit < 1) then
       ifault = 2
    else if (nunit == 1) then
       ! Initialise means and ssp for sample size  =  1
       do i = 1,nvar
          if (wt < co) then
             xmean(i) = xmean(i) + b*(x(i)-xmean(i))
          else
             xmean(i) = x(i)
          endif
          x(i) = co
          k = (i*(i-1))/2
          do j = 1,i
             k = k + 1
             xssp(k) = co
          enddo
       enddo
    else if (nunit > 1) then
       ! Update means and ssp for sample size greater than 1
       c = wt - b*wt
       do i = 1,nvar
          x(i) = x(i) - xmean(i)
          xmean(i) = xmean(i) + b*x(i)
       enddo
       do i = 1,nvar
          k = (i*(i-1))/2
          do j = 1,i
             k = k + 1
             xssp(k) = xssp(k) + c*x(i)*x(j)
          enddo
       enddo
    endif

  end subroutine dssp

end module covsubs
