module ekman

  use box, only: box_type
  use numerics, only: map_P_to_x_face, map_P_to_y_face, map_T_to_P, int_P_dx
  use ncutils, only: nc_open, nc_close, nc_get_double

  implicit none

  private

  type ekman_type

     double precision, allocatable :: taux(:,:), tauy(:,:)

     double precision, allocatable :: wekt(:,:)
     double precision, allocatable :: wekp(:,:)
     ! wekpo is the ocean Ekman velocity at p gridpoints (m s^-1)
     ! wekto is the ocean Ekman velocity at T gridpoints (m s^-1)

     double precision, allocatable :: uek(:,:)
     double precision, allocatable :: vek(:,:)

     ! txis, txin are the integrals of the windstress
     ! component taux along the southern and northern boundaries of the
     ! internal q domain (i.e. 1/2 gridlength in from physical boundaries)
     double precision :: txis, txin

  end type ekman_type


  public ekman_type
  public init_ekman
  public ekman_from_tau
  public load_tau_from_file

contains

  type(ekman_type) function init_ekman(b)
    
    type(box_type), intent(in) :: b

    allocate(init_ekman%taux(b%nxp,b%nyp))
    allocate(init_ekman%tauy(b%nxp,b%nyp))

    allocate(init_ekman%wekt(b%nxt,b%nyt))
    allocate(init_ekman%wekp(b%nxp,b%nyp))

    allocate(init_ekman%uek(b%nxp,b%nyt))
    allocate(init_ekman%vek(b%nxt,b%nyp))

  end function init_ekman

  subroutine ekman_from_tau(b, ek, tau_sign)

    type(box_type), intent(in) :: b
    type(ekman_type), intent(inout) :: ek
    integer, intent(in) :: tau_sign

    double precision :: taux_t(b%nxp-1,b%nyp),tauy_t(b%nxp,b%nyp-1)

    if (b%cyclic) then
       ! Calculate contributions of integrated Ekman velocity
       ! at internal p points to the momentum constraints.
       ! ----------------------------------------------------
       ! Compute by applying Stokes' theorem 1/2 (atmos. res)
       ! gridlength in from the Southern & Northern boundaries,
       ! i.e. surrounding those p-cells which are timestepped.
       ek%txis = 0.5d0*int_P_dx(ek%taux(:,1)       + ek%taux(:,2),     b)
       ek%txin = 0.5d0*int_P_dx(ek%taux(:,b%nyp-1) + ek%taux(:,b%nyp), b)
    endif

    call map_P_to_x_face(ek%taux, b, taux_t)
    call map_P_to_y_face(ek%tauy, b, tauy_t)

    ek%uek(:,:) = ( tau_sign/(b%fnot*b%hm))*tauy_t(:,:)
    ek%vek(:,:) = (-tau_sign/(b%fnot*b%hm))*taux_t(:,:)

    ! Compute oceanic Ekman velocity at T points - equation (7.7)
    ek%wekt(:,:) = ((1.0d0/(b%dx*b%fnot))*(tauy_t(2:,:) - tauy_t(:b%nxt,:))  - &
                    (1.0d0/(b%dy*b%fnot))*(taux_t(:,2:) - taux_t(:,:b%nyt)) )
    call map_T_to_P(ek%wekt, b, ek%wekp)

  end subroutine ekman_from_tau

  subroutine load_tau_from_file(indir, filename, b, ek)
    character (len=64), intent(in) :: indir
    character (len=*), intent(in) :: filename
    type(box_type), intent(in) :: b
    type(ekman_type), intent(inout) :: ek

    integer :: tempid

    character :: subnam*(*)
    parameter ( subnam = 'load_tau_from_file' )

    print *,' Mean forcing for ocean only case read from netCDF'
    tempid = nc_open(trim(indir)//'/'//filename, subnam)
    ek%taux = nc_get_double(tempid, 'tauxo', b%nxp, b%nyp, subnam)
    ek%tauy = nc_get_double(tempid, 'tauyo', b%nxp, b%nyp, subnam)
    call nc_close(tempid, subnam)

  end subroutine load_tau_from_file

end module ekman
