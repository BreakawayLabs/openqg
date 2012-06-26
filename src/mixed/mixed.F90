module mixed

  use box, only: box_type
  use ncutils, only: nc_open, nc_get_double
  use numerics, only: sin_lat

  implicit none

  private

  type mixed_type
     
     logical :: active = .false.

     double precision, allocatable :: data(:,:)
     double precision, allocatable :: datam(:,:)

     ! Advection
     double precision :: ycexp
     logical :: fluxn, fluxs
     double precision :: nval, sval

     ! Diffusion
     double precision :: d2, d4

  end type mixed_type

  type temp_type

     type(box_type) :: b

     ! Vertical temperature profile
     double precision :: T1_rel ! T_bar - T_(1)
     double precision :: dT_12  ! T(1) - T(2)

     ! Zonal temperature profile
     double precision :: dT_NS  ! (T(N) - T(S))/2
     double precision :: y0
     double precision :: ymax

     double precision :: rho_cp
     
  end type temp_type

     ! tabsoc is the absolute (potential) temp. of each ocean layer (K)
     ! toc/tat is the temperature anomaly of each oceanic/atmos layer,
     !         relative to the mean state radiative equilibrium (K)
     ! tmbaro/tmbara is the mean ocean/atmos mixed layer absolute temperature (K)

     ! tsbdy is the (relative) southern boundary temp-
     ! erature specified for the oceanic mixed later (K),
     ! only active if the sb_hflux option is activated.
     ! tnbdy is the (relative) northern boundary temp-
     ! erature specified for the oceanic mixed later (K),
     ! only active if the nb_hflux option is activated.
     !
     ! sstbar is the radiative equilibrium SST anomaly (K)
     ! (a function of latitude only)
     ! astbar is the radiative equilibrium a.m.l. temp. anomaly (K)
     ! (a function of latitude only)


  public mixed_type
  public temp_type

  public init_mixed

  public init_temp
  public init_temp_from_rad
  public zonal_temp

contains

  type(mixed_type) function init_mixed(b, init_val, fluxn, fluxs, nval, sval, d2, d4, ycexp)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: init_val
    logical, intent(in) :: fluxn, fluxs
    double precision, intent(in) :: nval, sval
    double precision, intent(in) :: d2, d4
    double precision, intent(in) :: ycexp

    allocate(init_mixed%data(b%nxt,b%nyt))
    allocate(init_mixed%datam(b%nxt,b%nyt))
    init_mixed%data(:,:) = init_val
    init_mixed%datam(:,:) = init_val

    init_mixed%fluxn = fluxn
    init_mixed%fluxs = fluxs

    init_mixed%nval = nval
    init_mixed%sval = sval

    init_mixed%d2 = d2
    init_mixed%d4 = d4
    init_mixed%ycexp = ycexp

    init_mixed%active = .true.

  end function init_mixed

  subroutine init_temp_from_rad(b, ga, T1, T2, Tbar, rbtm, fspco, rho, cp, temp)

    type(box_type), intent(in) :: b
    type(box_type), intent(in) :: ga
    double precision, intent(in) :: T1
    double precision, intent(in) :: T2
    double precision, intent(in) :: Tbar
    double precision, intent(in) :: rbtm
    double precision, intent(in) :: fspco
    double precision, intent(in) :: rho
    double precision, intent(in) :: cp
    type(temp_type), intent(out) :: temp

    temp%b = b
    temp%T1_rel = T1 - Tbar
    temp%dT_12 = T1 - T2
    temp%dT_NS = abs(rbtm*fspco)
    temp%y0 = ga%yp(1) + 0.5d0*ga%yl
    temp%ymax = ga%yp(1) + ga%yl
    temp%rho_cp = rho*cp

  end subroutine init_temp_from_rad

  subroutine init_temp(b, T1_rel, T2_rel, dT_NS, rho, cp, temp)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: T1_rel
    double precision, intent(in) :: T2_rel
    double precision, intent(in) :: dT_NS
    double precision, intent(in) :: rho
    double precision, intent(in) :: cp
    type(temp_type), intent(out) :: temp

    temp%b = b
    temp%T1_rel = T1_rel
    temp%dT_12 = T1_rel - T2_rel
    temp%dT_NS = dT_NS
    temp%y0 = b%yp(1) + 0.5d0*b%yl
    temp%ymax = b%yp(1) + b%yl
    temp%rho_cp = rho*cp

  end subroutine init_temp

  double precision function zonal_temp(temp, j)

    type(temp_type), intent(in) :: temp
    integer, intent(in) :: j

    zonal_temp = sin_lat(temp%b, temp%b%yt(j), temp%y0, temp%ymax, temp%dT_NS)

  end function zonal_temp

end module mixed
