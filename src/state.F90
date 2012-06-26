module state

  use util, only: streq, s2s
  use box, only: box_type
  use qg, only: qg_type
  use mixed, only: zonal_temp
  use radsubs, only: fsprim
  use ncutils, only: nc_open, nc_get_dim, nc_get_double, nc_def_dim, nc_create
  use ncutils, only: nc_def_double, nc_put_double, nc_enddef, nc_close
  use vorsubs, only: init_pv
  use omlsubs, only: ocean_mixed_type
  use amlsubs, only: atmos_mixed_type

  implicit none

  private

  public save_qg
  public save_aml
  public save_oml

  public restart_qg
  public restart_aml
  public restart_oml

contains

  subroutine restart_qg(indir, filename, qg)

    character (len=*), intent(in) :: indir
    character (len=*), intent(in) :: filename
    type(qg_type), intent(inout) :: qg

    integer :: nc_id
    character subnam*(*)
    parameter ( subnam = 'restart_qg' )

    nc_id = nc_open(trim(indir)//"/"//trim(filename), subnam)    
    
    qg%p = nc_get_double(nc_id, 'p', qg%b%nxp, qg%b%nyp, qg%b%nl, subnam)
    qg%pm = nc_get_double(nc_id, 'pm', qg%b%nxp, qg%b%nyp, qg%b%nl, subnam)

  end subroutine restart_qg

  subroutine restart_oml(indir, filename, oml)

    character (len=*), intent(in) :: indir
    character (len=*), intent(in) :: filename
    type(ocean_mixed_type), intent(inout) :: oml

    integer :: nc_id
    character subnam*(*)
    parameter ( subnam = 'restart_oml' )

    nc_id = nc_open(trim(indir)//"/"//trim(filename), subnam)    
    
    oml%sst%data = nc_get_double(nc_id, 'sst', oml%b%nxt, oml%b%nyt, subnam)
    oml%sst%datam = nc_get_double(nc_id, 'sstm', oml%b%nxt, oml%b%nyt, subnam)

  end subroutine restart_oml

  subroutine restart_aml(indir, filename, aml)

    character (len=*), intent(in) :: indir
    character (len=*), intent(in) :: filename
    type(atmos_mixed_type), intent(inout) :: aml

    integer :: nc_id
    character subnam*(*)
    parameter ( subnam = 'restart_aml' )

    nc_id = nc_open(trim(indir)//"/"//trim(filename), subnam)    
    
    aml%ast%data = nc_get_double(nc_id, 'ast', aml%b%nxt, aml%b%nyt, subnam)
    aml%ast%datam = nc_get_double(nc_id, 'astm', aml%b%nxt, aml%b%nyt, subnam)
    aml%hmixa%data = nc_get_double(nc_id, 'hmixa', aml%b%nxt, aml%b%nyt, subnam)
    aml%hmixa%datam = nc_get_double(nc_id, 'hmixam', aml%b%nxt, aml%b%nyt, subnam)

  end subroutine restart_aml

  subroutine save_qg(outdir, filename, qg)

    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename
    type(qg_type), intent(in) :: qg

    integer :: nc_id, xdim, ydim, zdim, tmp_id
    character subnam*(*)
    parameter ( subnam = 'save_qg' )

    nc_id = nc_create(outdir, filename)

    xdim = nc_def_dim(nc_id, "nxp", qg%b%nxp, subnam)
    ydim = nc_def_dim(nc_id, "nyp", qg%b%nyp, subnam)
    zdim = nc_def_dim(nc_id, "nl", qg%b%nl, subnam)

    tmp_id = nc_def_double(nc_id, 'x', xdim, 'km', subnam)
    tmp_id = nc_def_double(nc_id, 'y', ydim, 'km', subnam)
    tmp_id = nc_def_double(nc_id, 'p', (/xdim,ydim,zdim/), 'Current pressure', subnam)
    tmp_id = nc_def_double(nc_id, 'pm', (/xdim,ydim,zdim/), 'Previous pressure', subnam)

    call nc_enddef(nc_id)

    call nc_put_double(nc_id, 'x', qg%b%xp/1.0d3, subnam)
    call nc_put_double(nc_id, 'y', qg%b%yp/1.0d3, subnam)
    call nc_put_double(nc_id, 'p', qg%p, subnam)
    call nc_put_double(nc_id, 'pm', qg%pm, subnam)

    call nc_close(nc_id)

  end subroutine save_qg

  subroutine save_oml(outdir, filename, oml)

    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename
    type(ocean_mixed_type), intent(in) :: oml

    integer :: nc_id, xdim, ydim, tmp_id
    character subnam*(*)
    parameter ( subnam = 'save_oml' )

    nc_id = nc_create(outdir, filename)    

    xdim = nc_def_dim(nc_id, "nxt", oml%b%nxt, subnam)
    ydim = nc_def_dim(nc_id, "nyt", oml%b%nyt, subnam)

    tmp_id = nc_def_double(nc_id, 'x', xdim, 'km', subnam)
    tmp_id = nc_def_double(nc_id, 'y', ydim, 'km', subnam)
    tmp_id = nc_def_double(nc_id, 'sst', (/xdim,ydim/), 'Current mixed layer temp', subnam)
    tmp_id = nc_def_double(nc_id, 'sstm', (/xdim,ydim/), 'Previous mixed layer temp', subnam)

    call nc_enddef(nc_id)

    call nc_put_double(nc_id, 'x', oml%b%xt/1.0d3, subnam)
    call nc_put_double(nc_id, 'y', oml%b%yt/1.0d3, subnam)
    call nc_put_double(nc_id, 'sst', oml%sst%data, subnam)
    call nc_put_double(nc_id, 'sstm', oml%sst%datam, subnam)

    call nc_close(nc_id)

  end subroutine save_oml

  subroutine save_aml(outdir, filename, aml)

    character (len=*), intent(in) :: outdir
    character (len=*), intent(in) :: filename
    type(atmos_mixed_type), intent(in) :: aml

    integer :: nc_id, xdim, ydim, tmp_id
    character subnam*(*)
    parameter ( subnam = 'save_aml' )

    nc_id = nc_create(outdir, filename)    

    xdim = nc_def_dim(nc_id, "nxt", aml%b%nxt, subnam)
    ydim = nc_def_dim(nc_id, "nyt", aml%b%nyt, subnam)

    tmp_id = nc_def_double(nc_id, 'x', xdim, 'km', subnam)
    tmp_id = nc_def_double(nc_id, 'y', ydim, 'km', subnam)
    tmp_id = nc_def_double(nc_id, 'ast', (/xdim,ydim/), 'Current mixed layer temp', subnam)
    tmp_id = nc_def_double(nc_id, 'astm', (/xdim,ydim/), 'Previous mixed layer temp', subnam)
    tmp_id = nc_def_double(nc_id, 'hmixa', (/xdim,ydim/), 'Current mixed layer height', subnam)
    tmp_id = nc_def_double(nc_id, 'hmixam', (/xdim,ydim/), 'Previous mixed layer height', subnam)

    call nc_enddef(nc_id)

    call nc_put_double(nc_id, 'x', aml%b%xt/1.0d3, subnam)
    call nc_put_double(nc_id, 'y', aml%b%yt/1.0d3, subnam)
    call nc_put_double(nc_id, 'ast', aml%ast%data, subnam)
    call nc_put_double(nc_id, 'astm', aml%ast%datam, subnam)
    call nc_put_double(nc_id, 'hmixa', aml%hmixa%data, subnam)
    call nc_put_double(nc_id, 'hmixam', aml%hmixa%datam, subnam)

    call nc_close(nc_id)

  end subroutine save_aml

end module state
