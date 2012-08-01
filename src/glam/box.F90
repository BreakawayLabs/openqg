module box

  use units, only: m_to_km
  use mesh, only: mesh_type  
  use ncutils, only: nc_open, nc_get_double, nc_get_int, nc_get_dim

  implicit none

  private

  type box_type
     ! Specified values
     integer :: nxt = 0
     integer :: nyt = 0
     integer :: nl = 0
     double precision :: dx = 0.0d0
     double precision :: dy = 0.0d0
     double precision :: x0, y0
     double precision, allocatable :: h(:)
     double precision :: hm
     integer :: dz_sign

     logical :: cyclic

     ! Derived values
     double precision, allocatable :: yp(:)
     double precision, allocatable :: yprel(:)
     double precision, allocatable :: yt(:)
     double precision, allocatable :: ytrel(:)

     double precision, allocatable :: xp(:)
     double precision, allocatable :: xt(:)

     double precision :: fnot, beta
     integer :: nxp = 0
     integer :: nyp = 0
     double precision :: xl, yl
     double precision :: dxm2
     double precision :: rdxf0

  end type box_type

  public box_type
  public print_box_summary
  public load_box_from_filename
  public init_box_from_mesh

contains

  type(box_type) function init_box_from_mesh(mesh, nx, ny, ndx, ndy, nx1, ny1, nl, h, hm, dz_sign)

    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: nx, ny, nl
    integer, intent(in) :: ndx, ndy
    integer, intent(in) :: nx1, ny1
    double precision, intent(in) :: h(nl), hm
    integer, intent(in) :: dz_sign

    init_box_from_mesh%nxt = nx*ndx
    init_box_from_mesh%nyt = ny*ndy

    init_box_from_mesh%x0 = (nx1-1)*mesh%dx
    init_box_from_mesh%y0 = (ny1-1)*mesh%dy

    init_box_from_mesh%dx = mesh%dx/ndx
    init_box_from_mesh%dy = mesh%dy/ndy

    init_box_from_mesh%nl = nl
    allocate(init_box_from_mesh%h(nl))
    init_box_from_mesh%h(:) = h(:)
    init_box_from_mesh%hm = hm

    init_box_from_mesh%dz_sign = dz_sign

    init_box_from_mesh%cyclic = nx == mesh%nxt
    init_box_from_mesh%fnot = mesh%fnot
    init_box_from_mesh%beta = mesh%beta

    call compute_derived_values(init_box_from_mesh)
    call print_box_summary(init_box_from_mesh, " ")

  end function init_box_from_mesh

  type(box_type) function load_box_from_filename(mesh, filename)

    type(mesh_type), intent(in) :: mesh
    character (len=*), intent(in) :: filename

    integer :: box_id
    integer :: nx, ny, ndx, ndy, nx1, ny1

    box_id = nc_open(filename, 'Load grid')
    
    nx = nc_get_int(box_id, 'nx', 'load_box')
    ny = nc_get_int(box_id, 'ny', 'load_box')
    ndx = nc_get_int(box_id, 'ndx', 'load_box')
    ndy = nc_get_int(box_id, 'ndy', 'load_box')

    load_box_from_filename%nxt = nx*ndx
    load_box_from_filename%nyt = ny*ndy
       
    nx1 = nc_get_int(box_id, 'nx1', 'load_box')
    ny1 = nc_get_int(box_id, 'ny1', 'load_box')
    load_box_from_filename%x0 = (nx1-1)*mesh%dx
    load_box_from_filename%y0 = (ny1-1)*mesh%dy

    load_box_from_filename%dx = mesh%dx/ndx
    load_box_from_filename%dy = mesh%dy/ndy

    load_box_from_filename%nl = nc_get_dim(box_id, 'nl', 'load_box')
    allocate(load_box_from_filename%h(load_box_from_filename%nl))
    load_box_from_filename%h(:) = nc_get_double(box_id, 'h', load_box_from_filename%nl, 'load_box')
    load_box_from_filename%hm = nc_get_double(box_id, 'hm', 'load_box')

    load_box_from_filename%dz_sign = nc_get_int(box_id, 'dz_sign', 'load_box')

    load_box_from_filename%cyclic = nx == mesh%nxt
    load_box_from_filename%fnot = mesh%fnot
    load_box_from_filename%beta = mesh%beta

    call compute_derived_values(load_box_from_filename)
    call print_box_summary(load_box_from_filename, filename )

  end function load_box_from_filename

  pure subroutine compute_derived_values(b)

    type(box_type), intent(inout) :: b

    integer :: i, j

    b%nxp = b%nxt + 1
    b%nyp = b%nyt + 1

    b%xl = b%nxt*b%dx
    b%yl = b%nyt*b%dy

    b%dxm2 = 1.0d0/(b%dx*b%dx)

    allocate(b%yp(b%nyp))
    allocate(b%yprel(b%nyp))
    do j=1,b%nyp
       b%yp(j) = b%y0 + (j-1)*b%dy
       b%yprel(j) = -0.5d0*b%yl + (j-1)*b%dy 
    enddo
    allocate(b%yt(b%nyt))
    allocate(b%ytrel(b%nyt))
    do j=1,b%nyt
       b%yt(j) = b%yp(j) + 0.5d0*b%dy
       b%ytrel(j) = b%yprel(j) + 0.5d0*b%dy
    enddo
    allocate(b%xp(b%nxp))
    do i=1,b%nxp
       b%xp(i) = b%x0 + (i-1)*b%dx
    enddo
    allocate(b%xt(b%nxt))
    do i=1,b%nxt
       b%xt(i) = b%xp(i) + 0.5d0*b%dx
    enddo

    b%rdxf0 = 1.0d0/(b%dx*b%fnot) 

  end subroutine compute_derived_values


  subroutine print_box_summary(b, title)
    type(box_type), intent(in) :: b
    character (len=*), intent(in) :: title

    print *, "Box Summary: ", title
    print *, "Indexing"
    print *, "--------"
    print *, "T-grid (nxt,nyt,nl): ", b%nxt, b%nyt, b%nl
    print *, "P-grid (nxp,nyp,nl): ", b%nxp, b%nyp, b%nl
    print *, ""
    print *, "Physical size"
    print *, "-------------"
    print *, "T-grid size (dx, dy) (km): ", m_to_km(b%dx), m_to_km(b%dy)
    print *, "S-W coordinate (x0, y0) (km): ", m_to_km(b%x0), m_to_km(b%y0)
    print *, "Total size (x, y) (km): ", m_to_km(b%xl), m_to_km(b%yl)
    print *, "QG layer depths (m): ", b%h
    print *, "Mixed layer depth (m): ", b%hm
    print *, "Cyclic: ", b%cyclic
    print *, "Southwest P-cell (x, y) (km): ", b%xp(1), b%yp(1)
    print *, "Northeast P-cell (x, y) (km): ", b%xp(b%nxp), b%yp(b%nyp)
    print *, "Southwest T-cell (x, y) (km): ", b%xt(1), b%yt(1)
    print *, "Northeast T-cell (x, y) (km): ", b%xt(b%nxt), b%yt(b%nyt)
    print *, ""
    print *, "Coriolis"
    print *, "--------"
    print *, "f0: ", b%fnot
    print *, "beta: ", b%beta
    print *, ""

  end subroutine print_box_summary


end module box
