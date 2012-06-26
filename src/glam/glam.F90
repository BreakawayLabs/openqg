module glam

  use ncutils, only: nc_open, nc_get_int, nc_get_double
  use mesh, only: mesh_type, init_mesh

  implicit none

  private

  public load_glam

contains

  subroutine load_glam(filename, use_ocean, use_atmos, m)
    character (len=*), intent(in) :: filename
    logical, intent(out) :: use_ocean
    logical, intent(out) :: use_atmos
    type(mesh_type), intent(out) :: m

    integer :: glam_id

    integer :: nxt, nyt
    double precision :: dx, dy, x0, y0
    logical :: cyclic
    double precision :: x1, lat

    glam_id = nc_open(filename, 'Load glam') 
    nxt = nc_get_int(glam_id, 'nxt', 'load_glam')
    nyt = nc_get_int(glam_id, 'nyt', 'load_glam')

    x0 = nc_get_double(glam_id, 'x0', 'load_glam')
    y0 = nc_get_double(glam_id, 'y0', 'load_glam')

    dx = nc_get_double(glam_id, 'dx', 'load_glam')
    dy = nc_get_double(glam_id, 'dy', 'load_glam')
 
    cyclic = nc_get_int(glam_id, 'cyclic', 'load_glam') == 1
    lat = nc_get_double(glam_id, 'lat', 'load_glam')
    x1 = nc_get_double(glam_id, 'x1', 'load_glam')

    m = init_mesh(nxt, nyt, dx, dy, x0, y0, cyclic, x1, lat)

    use_ocean = nc_get_int(glam_id, 'use_ocean', 'load_glam') == 1
    use_atmos = nc_get_int(glam_id, 'use_atmos', 'load_glam') == 1

    if ( .not. use_ocean .and. .not. use_atmos) then
       print *,' '
       print *,' Invalid model configuration: ocean_only and'
       print *,' atmos_only options cannot both be selected'
       print *,' Program terminates'
       stop 1
    else if ( .not. use_atmos) then
       print *,' Model running in ocean_only configuration'
    else if ( .not. use_ocean) then
       print *,' Model running in atmos_only configuration'
    else
       print *,' Model running in full coupled configuration'
    endif

  end subroutine load_glam

end module glam
