module mesh
  
  use ncutils, only: nc_open, nc_get_int, nc_get_double

  implicit none

  private

  type mesh_type

     integer :: nxt, nyt
     double precision :: dx, dy
     double precision :: x0, y0

     logical :: cyclic
     double precision :: x1, lat

     double precision :: fnot, beta

  end type mesh_type

  public mesh_type
  public init_mesh
  public load_mesh
  public print_mesh_summary

contains

  type(mesh_type) pure function init_mesh(nxt, nyt, dx, dy, x0, y0, cyclic, x1, lat)

    integer, intent(in) :: nxt, nyt
    double precision, intent(in) :: dx, dy
    double precision, intent(in) :: x0, y0
    logical, intent(in) :: cyclic
    double precision, intent(in) :: x1
    double precision, intent(in) :: lat

    init_mesh%nxt = nxt
    init_mesh%nyt = nyt
    init_mesh%dx = dx
    init_mesh%dy = dy
    init_mesh%x0 = x0
    init_mesh%y0 = y0
    init_mesh%cyclic = cyclic
    init_mesh%x1 = x1
    init_mesh%lat = lat

    call compute_derived_values(init_mesh)

  end function init_mesh

  type(mesh_type) function load_mesh(filename)

    character (len=*), intent(in) :: filename
    integer :: cyclic
    integer :: mesh_id

    mesh_id = nc_open(filename, 'Load mesh') 
    load_mesh%nxt = nc_get_int(mesh_id, 'nxt', 'load_mesh')
    load_mesh%nyt = nc_get_int(mesh_id, 'nyt', 'load_mesh')

    load_mesh%x0 = nc_get_double(mesh_id, 'x0', 'load_mesh')
    load_mesh%y0 = nc_get_double(mesh_id, 'y0', 'load_mesh')

    load_mesh%dx = nc_get_double(mesh_id, 'dx', 'load_mesh')
    load_mesh%dy = nc_get_double(mesh_id, 'dy', 'load_mesh')
 
    cyclic = nc_get_int(mesh_id, 'cyclic', 'load_mesh')
    load_mesh%cyclic = (cyclic == 1)
    load_mesh%lat = nc_get_double(mesh_id, 'lat', 'load_mesh')
    load_mesh%x1 = nc_get_double(mesh_id, 'x1', 'load_mesh')

    call compute_derived_values(load_mesh)

  end function load_mesh

  pure subroutine compute_derived_values(mesh)
    type(mesh_type), intent(inout) :: mesh

    double precision :: pi, T, omega, phi, R0, xl, R

    pi = 3.14159265358979323d0
    T = (23*60 + 56)*60 + 4.098903691d0 ! Mean Solar Day
    omega = 2*pi/T                    ! rotation of earth
    phi = pi*mesh%lat/180.0d0      ! latitude in radians
    R0 = 6378.1e3 ! Estimate of the radius of the earth from Google.
    if (mesh%cyclic) then
       xl = mesh%nxt*mesh%dx
       R = xl/(2*pi*cos(phi))                  ! Effective radius of earth
    else if (mesh%x1 == 0.0d0) then
       ! Actual earth size
       R = R0*cos(phi)
    else
       R = mesh%x1/(2*pi*cos(phi))                  ! Effective radius of earth
    endif

    mesh%fnot = 2*omega*sin(phi)
    mesh%beta = 2*omega*cos(phi)/R

  end subroutine compute_derived_values

  subroutine print_mesh_summary(mesh, title)
    type(mesh_type), intent(in) :: mesh
    character (len=*), intent(in) :: title

    print *, "Mesh Summary: ", title
    print *, "Indexing"
    print *, "--------"
    print *, "T-grid (nx,ny): ", mesh%nxt, mesh%nyt
    print *, ""
    print *, "Physical size"
    print *, "-------------"
    print *, "T-grid size (dx, dy) (km): ", mesh%dx/1.0d3, mesh%dy/1.0d3
    print *, "S-W coordinate (x0, y0) (km): ", mesh%x0/1.0d3, mesh%y0/1.0d3
    print *, "Total size (x, y) (km): ", (mesh%nxt*mesh%dx)/1.0d3, (mesh%nyt*mesh%dy)/1.0d3
    print *, "Cyclic: ", mesh%cyclic
    print *, ""
    print *, "Coriolis"
    print *, "--------"
    print *, "Latitude (deg): ", mesh%lat
    print *, "f0: ", mesh%fnot
    print *, "beta: ", mesh%beta
    print *, ""

  end subroutine print_mesh_summary


end module mesh
