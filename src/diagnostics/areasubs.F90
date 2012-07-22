module areasubs

  use units, only: m_to_km
  use ncutils, only: nc_enddef, nc_def_dim, nc_def_int, nc_def_float
  use ncutils, only: nc_put_int, nc_put_double, nc_create, nc_close, nc_open_w

  use mixed, only: mixed_type
  use box, only: box_type

  implicit none

  private

  type subarea_type

     integer :: i1t, i2t, j1t, j2t
     integer :: i1p, i2p, j1p, j2p
     double precision :: fwt, fet, fst, fnt
     double precision :: fwp, fep, fsp, fnp

  end type subarea_type

  type area_avg_type

     logical :: active = .false.
     integer :: nocmon

     type(subarea_type), allocatable :: subareas(:)
     integer :: num_areas
     character (len=64) :: filename

  end type area_avg_type

  public area_avg_type

  public area_avg_step
  public update_area_avg
  public init_area_avg

contains

  logical function area_avg_step(ntdone, area_avgs)

    integer, intent(in) :: ntdone
    type(area_avg_type), intent(in) :: area_avgs

    area_avg_step = area_avgs%active .and. mod(ntdone,area_avgs%nocmon) == 0

  end function area_avg_step

  subroutine update_area_avg (outdir, ntdone, tyrs, st, b, area_avgs)

    ! Defines several subareas of ocean and atmosphere, and
    ! computes averages of various quantities over these areas
    ! Number and limits of areas specified in file "areas.limits"

    character (len=*), intent(in) :: outdir
    double precision, intent(in) :: tyrs
    integer, intent(in) :: ntdone
    type(mixed_type), intent(in) :: st
    type(box_type), intent(in) :: b
    type(area_avg_type), intent(in) :: area_avgs

    character :: subnam*(*)
    parameter ( subnam = 'areavg' )

    integer :: m,arencid

    double precision :: tav(area_avgs%num_areas)

    integer :: startt
   
    startt = ntdone/area_avgs%nocmon + 1

    ! Compute limited-area diagnostics
    arencid = nc_open_w(outdir, area_avgs%filename, subnam)
    call nc_put_double(arencid, 'time', startt, tyrs, subnam)
    do m=1,area_avgs%num_areas
       call areint (st%data, b%nxt, b%nyt, &
            area_avgs%subareas(m)%i1t, area_avgs%subareas(m)%i2t, &
            area_avgs%subareas(m)%j1t, area_avgs%subareas(m)%j2t, &
            area_avgs%subareas(m)%fwt, area_avgs%subareas(m)%fet, &
            area_avgs%subareas(m)%fst, area_avgs%subareas(m)%fnt, &
            tav(m))
    enddo
    call nc_put_double(arencid, 'data', startt, tav, subnam)
    call nc_close(arencid, subnam)

  end subroutine update_area_avg

  type(area_avg_type) function init_area_avg(outdir, filename, b, narea, xlo, xhi, ylo, yhi, nocmon, numoutsteps)

    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename
    type(box_type), intent(in) :: b
    integer, intent(in) :: narea
    double precision, intent(in) :: xlo(narea)
    double precision, intent(in) :: xhi(narea)
    double precision, intent(in) :: ylo(narea)
    double precision, intent(in) :: yhi(narea)
    integer, intent(in) :: nocmon
    integer, intent(in) :: numoutsteps

    integer :: m, nc_id, tmp(narea)    
    integer :: timedim, nocdim, noc_id, ocdat_id, tim_id

    character :: subnam*(*)
    parameter ( subnam = 'init_area_avg' )

    init_area_avg%filename = filename

    allocate(init_area_avg%subareas(narea))
    init_area_avg%num_areas = narea

    do m=1,narea
       call init_subarea(b, xlo(m), xhi(m), ylo(m), yhi(m), m, init_area_avg%subareas(m))
    enddo
    nc_id = nc_create(outdir, init_area_avg%filename, subnam)
    timedim = nc_def_dim(nc_id, 'time', numoutsteps, subnam)
    nocdim = nc_def_dim(nc_id, 'nare', narea, subnam)
    noc_id = nc_def_int(nc_id, 'nare', (/nocdim/), ' ', subnam, 'Areas')
    tim_id = nc_def_float(nc_id, 'time', (/timedim/),'years', subnam, 'Time')
    ocdat_id = nc_def_float(nc_id, 'data', (/nocdim, timedim/),'K', subnam, 'Area averaged temperature')
    call nc_enddef(nc_id, subnam)
    do m=1,narea
       tmp(m) = m
    enddo
    call nc_put_int(nc_id, noc_id, tmp, subnam)
    call nc_close(nc_id, subnam)
    init_area_avg%nocmon = nocmon
    init_area_avg%active = .true.

  end function init_area_avg

  subroutine init_subarea(b, xlo, xhi, ylo, yhi, index, subarea)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: xlo, xhi, ylo, yhi
    integer, intent(in) :: index
    type(subarea_type), intent(inout) :: subarea

    double precision :: frlo, frhi
    double precision :: rilo, rihi, rjlo, rjhi

    ! T points
    rilo = 1.0d0 + ( xlo - 0.5d0*b%dx )/b%dx
    rihi = 1.0d0 + ( xhi - 0.5d0*b%dx )/b%dx
    frlo = mod( rilo, 1.0d0 )
    subarea%i1t = int( rilo )
    if (frlo >= 0.5d0) then
       frlo = frlo - 1.0d0
       subarea%i1t = subarea%i1t + 1
    endif
    subarea%fwt = 0.5d0 - frlo

    frhi = mod( rihi, 1.0d0 )
    subarea%i2t = int( rihi )
    if (frhi > 0.5d0) then
       frhi = frhi - 1.0d0
       subarea%i2t = subarea%i2t + 1
    endif
    subarea%fet = 0.5d0 + frhi

    rjlo = 1.0d0 + ( ylo - 0.5d0*b%dy )/b%dy
    rjhi = 1.0d0 + ( yhi - 0.5d0*b%dy )/b%dy
    frlo = mod( rjlo, 1.0d0 )
    subarea%j1t = int( rjlo )
    if (frlo >= 0.5d0) then
       frlo = frlo - 1.0d0
       subarea%j1t = subarea%j1t + 1
    endif
    subarea%fst = 0.5d0 - frlo

    frhi = mod( rjhi, 1.0d0 )
    subarea%j2t = int( rjhi )
    if (frhi > 0.5d0) then
       frhi = frhi - 1.0d0
       subarea%j2t = subarea%j2t + 1
    endif
    subarea%fnt = 0.5d0 + frhi

    ! p points
    rilo = 1.0d0 + xlo/b%dx
    rihi = 1.0d0 + xhi/b%dx
    frlo = mod( rilo, 1.0d0 )
    subarea%i1p = int( rilo )
    if (frlo >= 0.5d0) then
       frlo = frlo - 1.0d0
       subarea%i1p = subarea%i1p + 1
    endif
    subarea%fwp = 0.5d0 - frlo

    frhi = mod( rihi, 1.0d0 )
    subarea%i2p = int( rihi )
    if (frhi > 0.5d0) then
       frhi = frhi - 1.0d0
       subarea%i2p = subarea%i2p + 1
    endif
    subarea%fep = 0.5d0 + frhi

    rjlo = 1.0d0 + ylo/b%dy
    rjhi = 1.0d0 + yhi/b%dy
    frlo = mod( rjlo, 1.0d0 )
    subarea%j1p = int( rjlo )
    if (frlo >= 0.5d0) then
       frlo = frlo - 1.0d0
       subarea%j1p = subarea%j1p + 1
    endif
    subarea%fsp = 0.5d0 - frlo

    frhi = mod( rjhi, 1.0d0 )
    subarea%j2p = int( rjhi )
    if (frhi > 0.5d0) then
       frhi = frhi - 1.0d0
       subarea%j2p = subarea%j2p + 1
    endif
    subarea%fnp = 0.5d0 + frhi

    print *
    print '(a,4i9)', '  area name       = ', index
    print '(a,4f9.1)', '  xlo, xhi, ylo, yhi = ', m_to_km(xlo), m_to_km(xhi), m_to_km(ylo), m_to_km(yhi)
    print '(a,4i9)',   '  i1t, i2t, j1t, j2t = ', subarea%i1t, subarea%i2t, subarea%j1t, subarea%j2t
    print '(a,4f9.3)', '             weights = ', subarea%fwt, subarea%fet, subarea%fst, subarea%fnt
    print '(a,4i9)',   '  i1p, i2p, j1p, j2p = ', subarea%i1p, subarea%i2p, subarea%j1p, subarea%j2p
    print '(a,4f9.3)', '             weights = ', subarea%fwp, subarea%fep, subarea%fsp, subarea%fnp

    if (xlo < 0.0d0 .or. xlo > b%xl .or. &
        xhi < 0.0d0 .or. xhi > b%xl) then
       print *,' limits outside x range for area = ',index
       print *,' xlo, xhi = ',xlo,xhi
       print *,' limits are : ',0.0d0,b%xl
       print *,' program terminates in areavg'
       stop
    endif
    if (ylo < 0.0d0 .or. ylo > b%yl .or. &
        yhi < 0.0d0 .or. yhi > b%yl) then
       print *,' limits outside y range for area = ',index
       print *,' ylo, yhi = ',ylo,yhi
       print *,' limits are : ',0.0d0,b%yl
       print *,' program terminates in areavg'
       stop
    endif
    ! Check subscript values
    if (subarea%i1t < 1 .or. subarea%i1t > b%nxt .or. &
        subarea%i2t < 1 .or. subarea%i2t > b%nxt) then
       print *,' invalid subscript for area = ',index
       print *,' i1t, i2t = ',subarea%i1t,subarea%i2t
       print *,' program terminates in areavg'
       stop
    endif
    if (subarea%j1t < 1 .or. subarea%j1t > b%nyt .or. &
        subarea%j2t < 1 .or. subarea%j2t > b%nyt) then
       print *,' invalid subscript for area = ',index
       print *,' j1t, j2t = ',subarea%j1t,subarea%j2t
       print *,' program terminates in areavg'
       stop
    endif
    if (subarea%i1p < 1 .or. subarea%i1p > b%nxp .or. &
        subarea%i2p < 1 .or. subarea%i2p > b%nxp) then
       print *,' invalid subscript for area = ',index
       print *,' i1p, i2p = ',subarea%i1p,subarea%i2p
       print *,' program terminates in areavg'
       stop
    endif
    if (subarea%j1p < 1 .or. subarea%j1p > b%nyp .or. &
        subarea%j2p < 1 .or. subarea%j2p > b%nyp) then
       print *,' invalid subscript for area = ',index
       print *,' j1p, j2p = ',subarea%j1p,subarea%j2p
       print *,' program terminates in areavg'
       stop
    endif

  end subroutine init_subarea

  subroutine areint (array, nx, ny, ilo, ihi, jlo, jhi, facw, face, facs, facn, answer)

    ! Computes area integral of subsection of array(nx,ny) from
    ! i=ilo to i=ihi and j=jlo to j=jhi. Weighting factors facw,
    ! face, facs and facn are applied to the boundary points.
    ! Modified version with reduced accumulator error

    integer, intent(in) :: nx,ny,ilo,ihi,jlo,jhi
    double precision, intent(in) :: array(nx,ny),facw,face,facs,facn
    double precision, intent(out) :: answer

    integer :: i,j
    double precision :: asums,wsums,asumi,wsumi,arsum,wtsum,asumn,wsumn

    arsum = 0.0d0
    wtsum = 0.0d0

    ! Southern boundary
    asums = facw*array(ilo,jlo)
    wsums = facw
    ! Northern boundary
    asumn = facw*array(ilo,jhi)
    wsumn = facw

    ! Inner latitudes
    do j=jlo+1,jhi-1
       asumi = facw*array(ilo,j)
       wsumi = facw
       do i=ilo+1,ihi-1
          asumi = asumi + array(i,j)
          wsumi = wsumi + 1.0d0
       enddo
       asumi = asumi + face*array(ihi,j)
       wsumi = wsumi + face
       arsum = arsum + asumi
       wtsum = wtsum + wsumi
    enddo
    do i=ilo+1,ihi-1
       ! Southern boundary
       asums = asums + array(i,jlo)
       wsums = wsums + 1.0d0
       ! Northern boundary
       asumn = asumn + array(i,jhi)
       wsumn = wsumn + 1.0d0
    enddo

    ! Southern boundary
    asums = asums + face*array(ihi,jlo)
    wsums = wsums + face
    ! Northern boundary
    asumn = asumn + face*array(ihi,jhi)
    wsumn = wsumn + face

    answer = ( facs*asums + facn*asumn + arsum ) &
         /( facs*wsums + facn*wsumn + wtsum )

  end subroutine areint

end module areasubs
