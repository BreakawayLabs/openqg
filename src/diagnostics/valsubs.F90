module valsubs

  use mixed, only: mixed_type
  use qg, only: qg_type
  use ekman, only: ekman_type
  use box, only: box_type
  use omlsubs, only: ocean_mixed_type
  use amlsubs, only: atmos_mixed_type

  implicit none

  type valids_type

     logical :: active = .false.
     integer :: nvalid

     double precision :: tauext
     double precision :: wtext
     double precision :: stext
     double precision :: pext
     double precision :: qext

  end type valids_type
  
  public valids_type
  public valids
  public valid_step

contains

  logical function valid_step(ntdone, valids) 
    integer, intent(in) :: ntdone
    type(valids_type), intent(in) :: valids

    valid_step = valids%active .and. mod(ntdone,valids%nvalid).eq.0

  end function valid_step

  ! Scan prognostic variable arrays for extreme values,
  ! and set argument to .false. if any are located.
  ! Diagnosis and termination of program can
  ! then be controlled by the calling program.
  ! Intended to be cheap to call, to enable
  ! frequent checking of solution validity.

  subroutine valids(solnok, val, qg, st, ek)

    logical, intent(inout) :: solnok
    type(valids_type), intent(in) :: val
    type(qg_type), intent(in) :: qg
    type(mixed_type), intent(in) :: st
    type(ekman_type), intent(in) :: ek

    double precision :: pmin, pmax
    double precision :: qmin, qmax
    double precision :: stmin, stmax
    double precision :: wekmin, wekmax
    double precision :: txmin, txmax
    double precision :: tymin, tymax

    pmin = minval(qg%p)
    pmax = maxval(qg%p)
    qmin = minval(qg%q)
    qmax = maxval(qg%q)
    if (st%active) then
       stmin = minval(st%data)
       stmax = maxval(st%data)
    endif
    wekmin = minval(ek%wekt)
    wekmax = maxval(ek%wekt)
    txmin = minval(ek%taux)
    txmax = maxval(ek%taux)
    tymin = minval(ek%tauy)
    tymax = maxval(ek%tauy)

    if ( abs(pmin) >= val%pext .or. abs(pmax) >= val%pext ) then
       print *,' '
       print *,' Invalid values of p found by valids'
       print *,' p: min, max = ', pmin, pmax
       solnok = .false.
       ! Scan for location of relevant extremum
       if ( abs(pmin) >= val%pext ) then
          call scan3D (qg%p, qg%b%nxp, qg%b%nyp, qg%b%nl, pmin, 'p minimum')
       endif
       if ( abs(pmax) >= val%pext ) then
          call scan3D (qg%p, qg%b%nxp, qg%b%nyp, qg%b%nl, pmax, 'p maximum')
       endif
    endif
    if ( abs(qmin) >= val%qext .or. abs(qmax) >= val%qext ) then
       print *,' '
       print *,' Invalid values of q found by valids'
       print *,' q: min, max = ', qmin, qmax
       solnok = .false.
       ! Scan for location of relevant extremum
       if ( abs(qmin) >= val%qext ) then
          call scan3D (qg%q, qg%b%nxp, qg%b%nyp, qg%b%nl, qmin, 'q minimum')
       endif
       if ( abs(qmax) >= val%qext ) then
          call scan3D (qg%q, qg%b%nxp, qg%b%nyp, qg%b%nl, qmax, 'q maximum')
       endif
    endif
    if (st%active) then
       if ( abs(stmin) >= val%stext .or. abs(stmax) >= val%stext ) then
          print *,' '
          print *,' Invalid values of st found by valids'
          print *,' st: min, max = ', stmin, stmax
          solnok = .false.
          ! Scan for location of relevant extremum
          if ( abs(stmin) >= val%stext ) then
             call scan2D (st%data, qg%b%nxt, qg%b%nyt, stmin, 'ST minimum')
          endif
          if ( abs(stmax) >= val%stext ) then
             call scan2D (st%data, qg%b%nxt, qg%b%nyt, stmax, 'ST maximum')
          endif
       endif
    endif
    if ( abs(wekmin) >= val%wtext .or. abs(wekmax) >= val%wtext ) then
       print *,' '
       print *,' Invalid values of wekt by valids'
       print *,' wekt: min, max = ', wekmin, wekmax
       solnok = .false.
       ! Scan for location of relevant extremum
       if ( abs(wekmin) >= val%wtext ) then
          call scan2D (ek%wekt, qg%b%nxt, qg%b%nyt, wekmin, 'wekt minimum')
       endif
       if ( abs(wekmax) >= val%wtext ) then
          call scan2D (ek%wekt, qg%b%nxt, qg%b%nyt, wekmax, 'wekt maximum')
       endif
    endif
    if ( abs(txmin).ge.val%tauext .or. abs(txmax).ge.val%tauext .or. &
         abs(tymin).ge.val%tauext .or. abs(tymax).ge.val%tauext ) then
       print *,' '
       print *,' Invalid values of taux or tauy found by valids'
       print *,' taux: min, max = ',txmin,txmax
       print *,' tauy: min, max = ',tymin,tymax
       solnok = .false.
    endif

  end subroutine valids

  subroutine scan2D (array, nx, ny, extrem, string)

    ! Scan a 2D array of dimensions (nx,ny) for the location of the
    ! value extrem, and print the location and surrounding values.
    ! 'string' identifies the extremum being located

    integer, intent(in) :: nx,ny
    double precision, intent(in) ::  array(nx,ny),extrem
    character, intent(in) ::  string*(*)

    integer :: i,j,im,jm,ilo,ihi

    im = 1
    jm = 1

    do j=1,ny
       do i=1,nx
          if ( array(i,j).eq.extrem ) then
             write(*,'(2x,a,a,2i7)') string,' located at i, j = ',i,j
             im = i
             jm = j
          endif
       enddo
    enddo
    ! Print neighbourhood of extremum
    ilo = max( 1,im-2)
    ihi = min(nx,im+2)
    write(*,'(3x,5i14)') (i,i=ilo,ihi)
    do j=min(ny,jm+3),max(1,jm-3),-1
       write(*,'(i7,1p,5d14.5)') j,(array(i,j),i=ilo,ihi)
    enddo

  end subroutine scan2D

  subroutine scan3D (array, nx, ny, nl, extrem, string)

    ! Scan a 3D array of dimensions (nx,ny,nl) for the location of the
    ! value extrem, and print the location and surrounding values.
    ! 'string' identifies the extremum being located

    integer, intent(in) :: nx,ny,nl
    double precision, intent(in) :: array(nx,ny,nl),extrem
    character, intent(in) :: string*(*)

    integer :: i,j,k,im,jm,km,ilo,ihi
    integer :: bad_count

    im = 1
    jm = 1
    km = 1
    bad_count = 0

    do k=1,nl
       do j=1,ny
          do i=1,nx
             if ( bad_count < 20 .and. array(i,j,k).eq.extrem ) then
                write(*,'(2x,a,a,3i7)') &
                     string,' located at i, j, k = ',i,j,k
                im = i
                jm = j
                km = k
                bad_count = bad_count + 1
             endif
          enddo
       enddo
    enddo
    ! Print neighbourhood of extremum
    ilo = max( 1,im-2)
    ihi = min(nx,im+2)
    write(*,'(3x,5i14)') (i,i=ilo,ihi)
    do j=min(ny,jm+3),max(1,jm-3),-1
       write(*,'(i7,1p,5d14.5)') j,(array(i,j,km),i=ilo,ihi)
    enddo

  end subroutine scan3D

end module valsubs
