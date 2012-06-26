module qg

  use util, only: strlen
  use ncutils, only: nc_open, nc_get_double, nc_get_dim, nc_close
  use box, only: box_type
  use numerics, only: int_P_dA
  use intsubs, only: trapin
  use topog, only: topog_type
  use constraint, only: constraint_type, core_constr_type, init_core_constr

  use modes, only: modes_type, init_modes, eigmod
  use homog, only: homog_type, init_homog, homsol
  use inhomog, only: inhomog_type, init_inhomog

  implicit none

  private

  type qg_type
     logical :: active = .false.

     ! bcco = mixed boundary coefft. for fluid (nondim.)
     ! bcco = 0  =>  free slip/no stress
     ! bcco -> Inf  =>  no slip

     ! gpoc = reduced gravity across oceanic interfaces (m s^-2)
     ! ah2oc = Del-sqd damping coefft for ocean (m^2 s^-1)
     ! ah4oc = Del-4th damping coefft for ocean (m^4 s^-1)
     ! delek = bottom Ekman layer thickness (m)

     ! Specified physical parameters.
     double precision :: bcco
     double precision :: delek
     double precision, allocatable :: ah2(:), ah4(:)
     double precision, allocatable :: gp(:)
     double precision :: rho

     integer :: tau_sign

     ! Specified GLAM objects
     type(box_type) :: b
     type(topog_type) :: topo

     ! Derivied static parameters
     type(modes_type) :: mod
     type(homog_type) :: hom
     type(inhomog_type) :: inhom

     ! Dynamic state
     double precision, allocatable :: q(:,:,:), qm(:,:,:)
     double precision, allocatable :: p(:,:,:), pm(:,:,:)
     type(constraint_type) :: constr
     type(core_constr_type) :: con

  end type qg_type

  public qg_type
  public init_foo_constr
  public init_qg

contains

  type(qg_type) function init_qg(go, bcco, delek, ah2, ah4, gp, rho, tau_sign, topo)
    type(box_type), intent(in) :: go
    double precision, intent(in) :: bcco
    double precision, intent(in) :: delek
    double precision, intent(in) :: ah2(go%nl)
    double precision, intent(in) :: ah4(go%nl)
    double precision, intent(in) :: gp(go%nl-1)
    double precision, intent(in) :: rho
    integer, intent(in) :: tau_sign
    type(topog_type), intent(in) :: topo

    init_qg%b = go
    init_qg%bcco = bcco
    init_qg%delek = delek
    allocate(init_qg%ah2(go%nl))
    allocate(init_qg%ah4(go%nl))
    allocate(init_qg%gp(go%nl))
    init_qg%ah2 = ah2
    init_qg%ah4 = ah4
    init_qg%gp = gp
    init_qg%rho = rho
    init_qg%tau_sign = tau_sign

    init_qg%topo = topo
    allocate(init_qg%q(init_qg%b%nxp,init_qg%b%nyp,init_qg%b%nl))
    allocate(init_qg%qm(init_qg%b%nxp,init_qg%b%nyp,init_qg%b%nl))
    allocate(init_qg%p(init_qg%b%nxp,init_qg%b%nyp,init_qg%b%nl))
    allocate(init_qg%pm(init_qg%b%nxp,init_qg%b%nyp,init_qg%b%nl))

    init_qg%q(:,:,:) = 0.0d0
    init_qg%qm(:,:,:) = 0.0d0
    init_qg%p(:,:,:) = 0.0d0
    init_qg%pm(:,:,:) = 0.0d0

    init_qg%mod = init_modes(init_qg%b%nl)
    init_qg%hom = init_homog(init_qg%b)
    init_qg%inhom = init_inhomog(init_qg%b)

    call eigmod(init_qg%b, init_qg%gp, init_qg%mod)
    call homsol(init_qg%b, init_qg%mod, init_qg%hom, init_qg%inhom)

    init_qg%con = init_core_constr(init_qg%b%nl)

    init_qg%active = .true.

  end function init_qg

  subroutine init_foo_constr(b, p, pm, con)
    type(box_type), intent(in) :: b
    double precision, intent(in) :: p(b%nxp,b%nyp,b%nl)
    double precision, intent(in) :: pm(b%nxp,b%nyp,b%nl)
    type(core_constr_type), intent(inout) :: con

    double precision wk(b%nxp,b%nyp,b%nl)

    integer :: k

    do k=1,b%nl-1
       ! Choose sign of dpioc so that +ve dpioc -> +ve eta
       wk(:,:,1) = b%dz_sign*(pm(:,:,k) - pm(:,:,k+1))
       wk(:,:,2) = b%dz_sign*(p(:,:,k)  - p(:,:,k+1))
       con%dpip(k) = int_P_dA(wk(:,:,1), b)
       con%dpi(k) = int_P_dA(wk(:,:,2), b)
    enddo
    
  end subroutine init_foo_constr

end module qg

