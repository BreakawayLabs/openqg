module mesh

  use units, only: m_to_km
  use ncutils, only: nc_open, nc_get_int, nc_get_double
  use constants, only: PI, R0, OMEGA

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

  pure subroutine compute_derived_values(mesh)
    type(mesh_type), intent(inout) :: mesh

    double precision :: phi, xl, R

    phi = PI*mesh%lat/180.0d0      ! latitude in radians
    if (mesh%cyclic) then
       xl = mesh%nxt*mesh%dx
       R = xl/(2*PI*cos(phi))                  ! Effective radius of earth
    else if (mesh%x1 == 0.0d0) then
       ! Actual earth size
       R = R0*cos(phi)
    else
       R = mesh%x1/(2*PI*cos(phi))                  ! Effective radius of earth
    endif

    mesh%fnot = 2*OMEGA*sin(phi)
    mesh%beta = 2*OMEGA*cos(phi)/R

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
    print *, "T-grid size (dx, dy) (km): ", m_to_km(mesh%dx), m_to_km(mesh%dy)
    print *, "S-W coordinate (x0, y0) (km): ", m_to_km(mesh%x0), m_to_km(mesh%y0)
    print *, "Total size (x, y) (km): ", m_to_km(mesh%nxt*mesh%dx), m_to_km(mesh%nyt*mesh%dy)
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
