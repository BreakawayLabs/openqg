module coupler

  use util, only: streq, s2s
  use box, only: box_type
  use ekman, only:  ekman_type, init_ekman

  implicit none

  private

  type ocn_coupler_type

     logical :: active = .false.

     ! External forcing
     double precision, allocatable :: fnet(:,:)
     type(ekman_type) :: ek

     ! Internal exchange
     double precision, allocatable :: p1(:,:)
     double precision, allocatable :: ent(:,:,:)

     logical :: fnet_coupled = .false.
     logical :: p1_coupled = .false.
     logical :: ent_coupled = .false.

  end type ocn_coupler_type

  type atm_coupler_type

     logical :: active = .false.

     ! Internal exchange
     double precision, allocatable :: eta(:,:,:)
     double precision, allocatable :: p1(:,:)
     double precision, allocatable :: ent(:,:,:)

     ! External forcing
     type(ekman_type) :: ek
     double precision, allocatable :: sst_datam(:,:)

     logical :: sst_coupled = .false.
     logical :: p1_coupled = .false.
     logical :: eta_coupled = .false.
     logical :: ent_coupled = .false.

  end type atm_coupler_type

  public ocn_coupler_type
  public atm_coupler_type

  public init_ocn_coupler
  public init_atm_coupler

contains  

  type(ocn_coupler_type) function init_ocn_coupler(b, fnet_cpl, p1_cpl, ent_cpl, qgo_active, oml_active)

    type(box_type), intent(in) :: b
    character, allocatable, intent(in) :: fnet_cpl(:)
    character, allocatable, intent(in) :: p1_cpl(:)
    character, allocatable, intent(in) :: ent_cpl(:)
    logical, intent(in) :: qgo_active
    logical, intent(in) :: oml_active

    init_ocn_coupler%ek = init_ekman(b)

    allocate(init_ocn_coupler%fnet(b%nxt,b%nyt))
    allocate(init_ocn_coupler%p1(b%nxp,b%nyp))
    allocate(init_ocn_coupler%ent(b%nxp,b%nyp,b%nl))

    init_ocn_coupler%fnet(:,:) = 0.0d0
    init_ocn_coupler%p1(:,:) = 0.0d0
    init_ocn_coupler%ent(:,:,:) = 0.0d0

    if (oml_active) then
       if (streq(fnet_cpl, "coupled")) then
          init_ocn_coupler%fnet_coupled = .true.
          init_ocn_coupler%fnet(:,:) = 0.0d0
       else if (streq(fnet_cpl, "zero")) then
          init_ocn_coupler%fnet_coupled = .false.
          init_ocn_coupler%fnet(:,:) = 0.0d0
       else
          init_ocn_coupler%fnet_coupled = .false.
          ! FIXME
       endif

       if (streq(p1_cpl, "coupled")) then
          init_ocn_coupler%p1_coupled = .true.
          init_ocn_coupler%p1(:,:) = 0.0d0
       else if (streq(p1_cpl, "zero")) then
          init_ocn_coupler%p1_coupled = .false.
          init_ocn_coupler%p1(:,:) = 0.0d0
       else
          init_ocn_coupler%p1_coupled = .false.
          ! FIXME
       endif
    endif

    if (qgo_active) then
       if (streq(ent_cpl, "coupled")) then
          init_ocn_coupler%ent_coupled = .true.
          init_ocn_coupler%ent(:,:,:) = 0.0d0
       else if (streq(ent_cpl, "zero")) then
          init_ocn_coupler%ent_coupled = .false.
          init_ocn_coupler%ent(:,:,:) = 0.0d0
       else
          init_ocn_coupler%ent_coupled = .false.
          ! FIXME
       endif
    endif

    init_ocn_coupler%active = .true.

  end function init_ocn_coupler

  type(atm_coupler_type) function init_atm_coupler(b, go, ent_cpl, p1_cpl, eta_cpl, sst_cpl)

    type(box_type), intent(in) :: b
    type(box_type), intenT(in) :: go
    character, allocatable, intent(in) :: ent_cpl(:)
    character, allocatable, intent(in) :: p1_cpl(:)
    character, allocatable, intent(in) :: eta_cpl(:)
    character, allocatable, intent(in) :: sst_cpl(:)

    init_atm_coupler%ek = init_ekman(b)

    if (streq(ent_cpl, "coupled")) then
       init_atm_coupler%ent_coupled = .true.
       allocate(init_atm_coupler%ent(b%nxp,b%nyp,b%nl))
       init_atm_coupler%ent(:,:,:) = 0.0d0
    else if (streq(ent_cpl, "zero")) then
       init_atm_coupler%ent_coupled = .false.
       allocate(init_atm_coupler%ent(b%nxp,b%nyp,b%nl))
       init_atm_coupler%ent(:,:,:) = 0.0d0
    else
       init_atm_coupler%ent_coupled = .false.
       ! FIXME
    endif

    if (streq(p1_cpl, "coupled")) then
       init_atm_coupler%p1_coupled = .true.
       allocate(init_atm_coupler%p1(b%nxp,b%nyp))
       init_atm_coupler%p1(:,:) = 0.0d0
    else if (streq(p1_cpl, "zero")) then
       init_atm_coupler%p1_coupled = .false.
       allocate(init_atm_coupler%p1(b%nxp,b%nyp))
       init_atm_coupler%p1(:,:) = 0.0d0
    else
       init_atm_coupler%p1_coupled = .false.
       ! FIXME
    endif

    if (streq(eta_cpl, "coupled")) then
       init_atm_coupler%eta_coupled = .true.
       allocate(init_atm_coupler%eta(b%nxp,b%nyp,b%nl-1))
       init_atm_coupler%eta(:,:,:) = 0.0d0
    else if (streq(eta_cpl, "zero")) then
       init_atm_coupler%eta_coupled = .false.
       allocate(init_atm_coupler%eta(b%nxp,b%nyp,b%nl-1))
       init_atm_coupler%eta(:,:,:) = 0.0d0
    else
       init_atm_coupler%eta_coupled = .false.
       ! FIXME
    endif

    if (streq(sst_cpl, "coupled")) then
       init_atm_coupler%sst_coupled = .true.
       allocate(init_atm_coupler%sst_datam(go%nxt,go%nyt))
       init_atm_coupler%sst_datam(:,:) = 0.0d0
    else if (streq(sst_cpl, "zero")) then
       init_atm_coupler%sst_coupled = .false.
       allocate(init_atm_coupler%sst_datam(go%nxt,go%nyt))
       init_atm_coupler%sst_datam(:,:) = 0.0d0
    else
       init_atm_coupler%sst_coupled = .false.
       ! FIXME
    endif

    init_atm_coupler%active = .true.

  end function init_atm_coupler

end module coupler
