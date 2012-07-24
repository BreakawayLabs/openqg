module ncutils

  use netcdf

  implicit none

  private

  public handle_err
  public nc_close
  public nc_create
  public nc_open
  public nc_open_w

  public nc_def_dim
  public nc_get_dim

  public nc_get_int


  public nc_get_double

  public nc_get_text

  public nc_enddef

  interface nc_get_int
     module procedure nc_get_int
     module procedure nc_get_int_1d
  end interface nc_get_int

  interface nc_get_double
     module procedure nc_get_double
     module procedure nc_get_double_1d
     module procedure nc_get_double_2d
     module procedure nc_get_double_3d
  end interface nc_get_double

  public nc_def_int
  interface nc_def_int
     module procedure nc_def_int
     module procedure nc_def_int_nd
  end interface nc_def_int
  
  public nc_def_double
  interface nc_def_double
     module procedure nc_def_double
     module procedure nc_def_double_nd
  end interface nc_def_double

  public nc_def_float
  interface nc_def_float
     module procedure nc_def_float
     module procedure nc_def_float_nd
  end interface nc_def_float

  public nc_put_int
  interface nc_put_int
     module procedure nc_put_int
     module procedure nc_put_int_simple
  end interface

  public nc_put_double

  interface nc_put_double
     module procedure nc_put_double_block_1d_id
     module procedure nc_put_double_block_2d_id
     module procedure nc_put_double_block_3d_id
     module procedure nc_put_double_block_1d_name
     module procedure nc_put_double_block_2d_name
     module procedure nc_put_double_block_3d_name
     module procedure nc_put_double_line_array_1d_id
     module procedure nc_put_double_line_array_1d_name
     module procedure nc_put_double_line_array_2d_id
     module procedure nc_put_double_line_array_2d_name
     module procedure nc_put_double_line_array_3d_id
     module procedure nc_put_double_line_array_3d_name
     module procedure nc_put_double_line_scalar_id
     module procedure nc_put_double_line_scalar_name
  end interface nc_put_double

contains

  subroutine handle_err(ncstat, fromst, msg)

    integer, intent(in) ::  ncstat
    character (len=*), intent(in), optional :: fromst
    character (len=*), intent(in), optional :: msg

    ! fromst is an optional string indicating where the call came
    ! from that caused the netCDF problem (e.g. subroutine name)

    ! Routine which interprets errors from netCDF output functions,
    ! prints them to standard output and then kills the whole run.
    if (ncstat /= NF90_NOERR) then
       if ( present(fromst) .and. present(msg) ) then
          print *, trim(fromst)//': '//trim( nf90_strerror(ncstat))//': '//msg
       else if ( present(fromst) ) then
          print *, trim(fromst)//': '//trim( nf90_strerror(ncstat) )
       else
          print *, trim( nf90_strerror(ncstat) )
       endif
       stop 1
    endif

  end subroutine handle_err

  subroutine nc_close(ncid, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in), optional :: subnam

    integer :: status

    status = nf90_close(ncid)
    if (status /= NF90_NOERR) call handle_err(status, subnam)

  end subroutine nc_close

  integer function nc_create(outdir, filename, subnam)
    
    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename
    character (len=*), intent(in), optional :: subnam

    integer :: status

    status = nf90_create(trim(outdir)//'/'//filename, NF90_CLOBBER, nc_create)
    if (status /= NF90_NOERR) call handle_err(status, subnam, trim(outdir)//'/'//filename)

  end function nc_create

  integer function nc_open(filename, subnam)

    character (len=*), intent(in) :: filename
    character (len=*), intent(in), optional :: subnam

    integer :: status

    status = nf90_open(filename, NF90_NOWRITE, nc_open)
    if (status /= NF90_NOERR) call handle_err(status, subnam, filename)

  end function nc_open

  integer function nc_open_w(outdir, filename, subnam)

    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename
    character (len=*), intent(in), optional :: subnam

    integer :: status

    status = nf90_open(trim(outdir)//'/'//filename, NF90_WRITE, nc_open_w)
    if (status /= NF90_NOERR) call handle_err(status, subnam, filename)

  end function nc_open_w

  double precision function nc_get_double(ncid, varname, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_get_var(ncid, varid, nc_get_double)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end function nc_get_double

  function nc_get_double_1d(ncid, varname, n, subnam)

    integer, intent(in) :: n
    double precision :: nc_get_double_1d(n)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_get_var(ncid, varid, nc_get_double_1d)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end function nc_get_double_1d

  function nc_get_double_2d(ncid, varname, n, m, subnam)

    integer, intent(in) :: n,m
    double precision :: nc_get_double_2d(n,m)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_get_var(ncid, varid, nc_get_double_2d)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end function nc_get_double_2d

  function nc_get_double_3d(ncid, varname, n, m, zzzz, subnam)

    integer, intent(in) :: n,m, zzzz
    double precision :: nc_get_double_3d(n,m,zzzz)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_get_var(ncid, varid, nc_get_double_3d)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end function nc_get_double_3d

  integer function nc_def_dim(ncid, dimname, n, subnam)
    
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: dimname
    integer, intent(in) :: n
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat

    ncstat = nf90_def_dim(ncid, dimname, n, nc_def_dim)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end function nc_def_dim

  integer function nc_get_dim(ncid, dimname, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: dimname
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    ncstat = nf90_inq_dimid(ncid, dimname, varid)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)
    ncstat = nf90_inquire_dimension(ncid, varid, len=nc_get_dim)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end function nc_get_dim

  subroutine nc_enddef(ncid, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat

    ncstat = nf90_enddef(ncid)
    if (ncstat /= NF90_NOERR ) call handle_err(ncstat, subnam)
    
  end subroutine nc_enddef

  
  subroutine nc_def_var(ncid, varname, dim, type, units, subnam, longname)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dim
    integer, intent(in) :: type
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    integer :: ncstat, var_id

    ncstat = nf90_def_var(ncid, varname, type, dim, var_id)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)
    ncstat = nf90_put_att(ncid, var_id, 'units', units)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)
    if (present(longname)) then
       ncstat = nf90_put_att(ncid, var_id, 'long_name', longname)
       if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)
    endif

  end subroutine nc_def_var

  subroutine nc_def_var_nd(ncid, varname, dims, type, units, subnam, longname)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dims(:)
    integer, intent(in) :: type
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    integer :: ncstat, var_id

    ncstat = nf90_def_var(ncid, varname, type, dims, var_id)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)
    ncstat = nf90_put_att(ncid, var_id, 'units', units)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)
    if (present(longname)) then
       ncstat = nf90_put_att(ncid, var_id, 'long_name', longname)
       if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)
    endif

  end subroutine nc_def_var_nd

  subroutine nc_def_float(ncid, varname, dim, units, subnam, longname)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dim
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    call nc_def_var(ncid, varname, dim, NF90_FLOAT, units, subnam, longname)

  end subroutine nc_def_float

  subroutine nc_def_float_nd(ncid, varname, dims, units, subnam, longname)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dims(:)
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    call nc_def_var_nd(ncid, varname, dims, NF90_FLOAT, units, subnam, longname)

  end subroutine nc_def_float_nd

  subroutine nc_def_int(ncid, varname, dim, units, subnam, longname)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dim
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    call nc_def_var(ncid, varname, dim, NF90_INT, units, subnam, longname)

  end subroutine nc_def_int

  subroutine nc_def_int_nd(ncid, varname, dims, units, subnam, longname)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dims(:)
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    call nc_def_var_nd(ncid, varname, dims, NF90_INT, units, subnam, longname)

  end subroutine nc_def_int_nd

  subroutine nc_def_double(ncid, varname, dim, units, subnam, longname)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dim
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    call nc_def_var(ncid, varname, dim, NF90_DOUBLE, units, subnam, longname)

  end subroutine nc_def_double

  subroutine nc_def_double_nd(ncid, varname, dims, units, subnam, longname)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: dims(:)
    character (len=*), intent(in) :: units
    character (len=*), intent(in), optional :: subnam
    character (len=*), intent(in), optional :: longname

    call nc_def_var_nd(ncid, varname, dims, NF90_DOUBLE, units, subnam, longname)

  end subroutine nc_def_double_nd

  integer function nc_get_int(ncid, varname, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_get_var(ncid, varid, nc_get_int)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end function nc_get_int

  function nc_get_int_1d(ncid, varname, n, subnam)

    integer, intent(in) :: n
    integer :: nc_get_int_1d(n)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_get_var(ncid, varid, nc_get_int_1d)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end function nc_get_int_1d

  subroutine nc_get_text(ncid, varname, data, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character, allocatable, intent(inout) :: data(:)
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_get_var(ncid, varid, data, count=shape(data))
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end subroutine nc_get_text

  subroutine nc_put_int(nc_id, varname, start, count, data, subnam)
    integer, intent(in) :: nc_id
    character (len=*), intent(in) :: varname
    integer, intent(in) :: start(:)
    integer, intent(in) :: count(:)
    integer, intent(in) :: data(:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(nc_id, varname, subnam)
    ncstat = nf90_put_var(nc_id, varid, data, start, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_int

  subroutine nc_put_int_simple(nc_id, varname, data, subnam)
    integer, intent(in) :: nc_id
    character (len=*), intent(in) :: varname
    integer, intent(in) :: data(:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(nc_id, varname, subnam)
    ncstat = nf90_put_var(nc_id, varid, data)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_int_simple


  subroutine nc_put_double_block_1d_id(ncid, varid, data, subnam)

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    double precision, intent(in) :: data(:)
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat

    ncstat = nf90_put_var(ncid, varid, data)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_block_1d_id

  subroutine nc_put_double_block_2d_id(ncid, varid, data, subnam)

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    double precision, intent(in) :: data(:,:)
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat

    ncstat = nf90_put_var(ncid, varid, data)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_block_2d_id

  subroutine nc_put_double_block_3d_id(ncid, varid, data, subnam)

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    double precision, intent(in) :: data(:,:,:)
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat

    ncstat = nf90_put_var(ncid, varid, data)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_block_3d_id

  subroutine nc_put_double_block_1d_name(ncid, varname, data, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    double precision, intent(in) :: data(:)
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_put_var(ncid, varid, data)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_block_1d_name

  subroutine nc_put_double_block_2d_name(ncid, varname, data, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    double precision, intent(in) :: data(:,:)
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_put_var(ncid, varid, data)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_block_2d_name

  subroutine nc_put_double_block_3d_name(ncid, varname, data, subnam)

    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    double precision, intent(in) :: data(:,:,:)
    character (len=*), intent(in), optional :: subnam

    integer :: ncstat, varid

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_put_var(ncid, varid, data)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_block_3d_name

! nc_put_double_line_scalar_id (ncid, varid, start, data, subname) (nc_put_float are a bit like this...)
! nc_put_double_line_array_id 
! nc_put_double_line_scalar_name (ncid, varname, start, data, subname)
! nc_put_double_line_array_name

  subroutine nc_put_double_line_array_1d_name(ncid, varname, start, data, subnam)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: start
    double precision, intent(in) :: data(:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat, varid
    integer :: startt(2), count(2)

    startt = (/1, start/)
    count = (/size(data), 1/)

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_put_var(ncid, varid, data, startt, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_array_1d_name

  subroutine nc_put_double_line_array_1d_id(ncid, varid, start, data, subnam)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start
    double precision, intent(in) :: data(:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat
    integer :: startt(2), count(2)

    startt = (/1, start/)
    count = (/size(data), 1/)

    ncstat = nf90_put_var(ncid, varid, data, startt, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_array_1d_id

  subroutine nc_put_double_line_array_2d_name(ncid, varname, start, data, subnam)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: start
    double precision, intent(in) :: data(:,:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat, varid
    integer :: startt(3), count(3)

    startt = (/1, 1, start/)
    count = (/shape(data), 1/)

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_put_var(ncid, varid, data, startt, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_array_2d_name

  subroutine nc_put_double_line_array_2d_id(ncid, varid, start, data, subnam)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start
    double precision, intent(in) :: data(:,:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat
    integer :: startt(3), count(3)

    startt = (/1, 1, start/)
    count = (/shape(data), 1/)

    ncstat = nf90_put_var(ncid, varid, data, startt, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_array_2d_id

  subroutine nc_put_double_line_array_3d_name(ncid, varname, start, data, subnam)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: start
    double precision, intent(in) :: data(:,:,:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat, varid
    integer :: startt(4), count(4)

    startt = (/1, 1, 1, start/)
    count = (/shape(data), 1/)

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_put_var(ncid, varid, data, startt, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_array_3d_name

  subroutine nc_put_double_line_array_3d_id(ncid, varid, start, data, subnam)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start
    double precision, intent(in) :: data(:,:,:)
    character (len=*), intent(in) :: subnam

    integer :: ncstat
    integer :: startt(4), count(4)

    startt = (/1, 1, 1, start/)
    count = (/shape(data), 1/)

    ncstat = nf90_put_var(ncid, varid, data, startt, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_array_3d_id

  subroutine nc_put_double_line_scalar_name(ncid, varname, start, data, subnam)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    integer, intent(in) :: start
    double precision, intent(in) :: data
    character (len=*), intent(in) :: subnam

    integer :: ncstat, varid
    integer :: startt(2), count(2)

    startt = (/start, 1/)
    count = (/1, 1/)

    varid = nc_inq_varid(ncid, varname, subnam)
    ncstat = nf90_put_var(ncid, varid, (/data/), startt, count)
    if (ncstat /= NF90_NOERR ) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_scalar_name

  subroutine nc_put_double_line_scalar_id(ncid, varid, start, data, subnam)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start
    double precision, intent(in) :: data
    character (len=*), intent(in) :: subnam

    integer :: ncstat
    integer :: startt(2), count(2)

    startt = (/1, start/)
    count = (/1, 1/)

    ncstat = nf90_put_var(ncid, varid, (/data/), startt, count)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam)

  end subroutine nc_put_double_line_scalar_id

  integer function nc_inq_varid(ncid, varname, subnam)
    integer, intent(in) :: ncid
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: subnam

    integer :: ncstat

    ncstat = nf90_inq_varid(ncid, varname, nc_inq_varid)
    if (ncstat /= NF90_NOERR) call handle_err(ncstat, subnam, varname)

  end function nc_inq_varid
  
end module ncutils
