module numerics

  use intsubs, only: trapin
  use box, only: box_type
  use constants, only: PI

  implicit none

  private

  ! Linear mapping 
  public map_T_to_P
  public map_P_to_T
  public map_T_to_y_face
  public map_T_to_x_face
  public map_P_to_y_face
  public map_P_to_x_face

  public avg_T
  public avg_P
  public int_P_dA
  public int_P_x_face_dA
  public int_P_y_face_dA

  public jacobian ! Second order (no boundary conditions)
  public c_grid_advection

  ! 1st deriv
  public dPdx ! First order
  public dPdy 
  public dPdx_2 ! Second order
  public dPdy_2
  public dPdx_2bc ! Second order w/ boundary conditions
  public dPdy_2bc

  ! 2nd deriv
  public dP2dx2_bc ! Second order 2nd deriv w/ boundary conditions
  public dP2dy2_bc
  public del2_P_bc

  public del2_T ! Place del^2 onto T points from u/v points with boundary values.

  public sin_lat

contains

  pure function map_P_to_T(p_data, b, fac)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp, b%nyp)
    double precision, intent(in) :: fac
    double precision :: map_P_to_T(b%nxt,b%nyt)

    map_P_to_T(:,:) = (fac*0.25d0)*(p_data(:b%nxp-1,:b%nyp-1) + p_data(2:,:b%nyp-1) + p_data(:b%nxp-1,2:) + p_data(2:,2:))

  end function map_P_to_T


  subroutine map_T_to_P(t_data, b, p_data)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: t_data(b%nxp-1,b%nyp-1)
    double precision, intent(out) :: p_data(b%nxp,b%nyp)

    integer :: i, j

    ! Intenal points
    p_data(2:b%nxp-1,2:b%nyp-1) = 0.25*(t_data(2:b%nxp-1,2:b%nyp-1) + t_data(:b%nxp-2,2:b%nyp-1) + &
                                        t_data(2:b%nxp-1,:b%nyp-2) + t_data(:b%nxp-2,:b%nyp-2))

    ! N/S boundaries
    do i=2,b%nxp-1
       p_data(i,1)     = 0.5d0*(t_data(i-1,1)       + t_data(i,1))
       p_data(i,b%nyp) = 0.5d0*(t_data(i-1,b%nyp-1) + t_data(i,b%nyp-1))
    enddo

    ! E/W boundaries
    if (b%cyclic) then
       do j=2,b%nyp-1
          p_data(1,j) = 0.25d0*(t_data(1,j) + t_data(b%nxp-1,j) + t_data(1,j-1) + t_data(b%nxp-1,j-1))
          p_data(b%nxp,j) = p_data(1,j)
       enddo
    else
       do j=2,b%nyp-1
          p_data(1,  j) = 0.5d0*(t_data(1,    j-1) + t_data(1,    j))
          p_data(b%nxp,j) = 0.5d0*(t_data(b%nxp-1,j-1) + t_data(b%nxp-1,j))
       enddo
    endif

    ! Corners
    if (b%cyclic) then
       p_data(1,    1)     = 0.5d0*(t_data(1,1)       + t_data(b%nxp-1,1))
       p_data(1,    b%nyp) = 0.5d0*(t_data(1,b%nyp-1) + t_data(b%nxp-1,b%nyp-1))
       p_data(b%nxp,1)     = p_data(1,1)
       p_data(b%nxp,b%nyp) = p_data(1,b%nyp)
    else
       p_data(1,    1)     = t_data(1,      1)
       p_data(b%nxp,1)     = t_data(b%nxp-1,1)
       p_data(1,    b%nyp) = t_data(1,      b%nyp-1)
       p_data(b%nxp,b%nyp) = t_data(b%nxp-1,b%nyp-1)
    endif

  end subroutine map_T_to_P

  pure function del2_T(t_data, b, north_flux, north_flux_val, south_flux, south_flux_val)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: t_data(b%nxt,b%nyt)
    logical, intent(in) :: north_flux
    double precision, intent(in) :: north_flux_val
    logical, intent(in) :: south_flux
    double precision, intent(in) :: south_flux_val
    double precision :: del2_T(b%nxt,b%nyt)

    double precision :: dt_x(b%nxt+1,b%nyt)
    double precision :: dt_y(b%nxt,b%nyt+1)

    dt_x(2:b%nxt,:) = t_data(2:b%nxt,:) - t_data(:b%nxt-1,:)
    if (b%cyclic) then
       dt_x(1,    :) = t_data(1,:) - t_data(b%nxt,:)
       dt_x(b%nxt+1,:) = dt_x(1,:)
    else
       dt_x(1,    :) = 0.0d0
       dt_x(b%nxt+1,:) = 0.0d0
    endif
    
    if (south_flux) then
       dt_y(:,1) = t_data(:,1) - south_flux_val
    else
       dt_y(:,1) = 0.0d0
    endif
    dt_y(:,2:b%nyt) = t_data(:,2:b%nyt) - t_data(:,:b%nyt-1)
    if (north_flux) then
       dt_y(:,b%nyt+1) = north_flux_val - t_data(:,b%nyt)
    else
       dt_y(:,b%nyt+1) = 0.0d0
    endif

    del2_T(:,:) = (1.0d0/(b%dx*b%dx))*(dt_x(2:,:) - dt_x(:b%nxt,:) + dt_y(:,2:) - dt_y(:,:b%nyt) )

  end function del2_T

  double precision pure function int_T_dA(t_data, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: t_data(b%nxt, b%nyt)

    int_T_dA = sum(t_data(:,:))*b%dx*b%dy

  end function int_T_dA

  double precision function avg_T(t_data, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: t_data(b%nxt, b%nyt)

    avg_T = sum(t_data(:,:))/(b%nxt*b%nyt)

  end function avg_T

  double precision pure function int_P_dA(p_data, b)
    
    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp, b%nyp)

    ! All points
    int_P_dA = sum(p_data(:,:))
    ! Remove N/S overhang
    int_P_dA = int_P_dA - 0.5d0*(sum(p_data(:,1)) + sum(p_data(:,b%nyp)))
    ! Remove E/W overhang
    int_P_dA = int_P_dA - 0.5d0*(sum(p_data(1,:)) + sum(p_data(b%nxp,:)))
    ! Add back in corners
    int_P_dA = int_P_dA + 0.25d0*(p_data(1,1) + p_data(1,b%nyp) + p_data(b%nxp,1) + p_data(b%nxp,b%nyp))
    ! Multiply by grid size
    int_P_dA = int_P_dA*b%dx*b%dy

  end function int_P_dA

  double precision function int_P_y_face_dA(y_face_data, b)
    
    type(box_type), intent(in) :: b
    double precision, intent(in) :: y_face_data(b%nxp,b%nyt)

    double precision :: int_dx(b%nyt)
    integer :: j
    if (b%cyclic) then
       do j=1,b%nyt
          int_dx(j) = sum(y_face_data(:b%nxp-1,j))*b%dx
       enddo
    else
       do j=1,b%nyt
          int_dx(j) = (sum(y_face_data(2:b%nxp-1,j)) &
                       + 0.5d0*y_face_data(1,j)      &
                       + 0.5d0*y_face_data(b%nxp,j))*b%dx
       enddo
    endif

    int_P_y_face_dA = sum(int_dx(:))*b%dy

  end function int_P_y_face_dA

  double precision function int_P_x_face_dA(x_face_data, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: x_face_data(b%nxt,b%nxp)
    
    double precision :: int_dy(b%nxt)
    integer :: i
    do i=1,b%nxt
       int_dy(i) = (sum(x_face_data(i,2:b%nyp-1)) &
                    + 0.5d0*x_face_data(i,1)      &
                    + 0.5d0*x_face_data(i,b%nyp))*b%dy
    enddo

    int_P_x_face_dA = sum(int_dy(:))*b%dx

  end function int_P_x_face_dA

  double precision function avg_P(p_data, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp, b%nyp)

    avg_P = sum(p_data(:,:))/(b%nxp*b%nyp)

  end function avg_P

  subroutine map_T_to_y_face(t_data, b, zero_south, south_val, zero_north, north_val, y_face)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: t_data(b%nxt,b%nyt)
    logical, intent(in) :: zero_south
    double precision, intent(in) :: south_val
    logical, intent(in) :: zero_north
    double precision, intent(in) :: north_val
    double precision, intent(out) :: y_face(b%nxt,b%nyt+1)

    integer :: j

    if (zero_south) then
       y_face(:,1) = t_data(:,1)
    else
       y_face(:,1) = 0.5d0*(t_data(:,1) + south_val)
    endif
    do j=2,b%nyt
       y_face(:,j) = 0.5d0*(t_data(:,j-1) + t_data(:,j))
    enddo
    if (zero_north) then
       y_face(:,b%nyt+1) = t_data(:,b%nyt)
    else
       y_face(:,b%nyt+1) = 0.5d0*(t_data(:,b%nyt) + north_val)
    endif

  end subroutine map_T_to_y_face

  subroutine map_T_to_x_face(t_data, b, x_face)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: t_data(b%nxt, b%nyt)
    double precision, intent(out) :: x_face(b%nxt+1,b%nyt)

    if (b%cyclic) then
       x_face(1,:) = 0.5d0*(t_data(1,:) + t_data(b%nxt,:))
    else
       x_face(1,:) = t_data(1,:)
    endif
    x_face(2:b%nxt,:) = 0.5d0*(t_data(:b%nxt-1,:) + t_data(2:,:))
    if (b%cyclic) then
       x_face(b%nxt+1,:) = x_face(1,:)
    else
       x_face(b%nxt+1,:) = t_data(b%nxt,:)
    endif

  end subroutine map_T_to_x_face

  subroutine map_P_to_x_face(p_data, b, x_face)
    
    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(out) :: x_face(b%nxp-1,b%nyp)

    x_face(:,:) = 0.5d0*(p_data(2:,:) + p_data(:b%nxp-1,:))

  end subroutine map_P_to_x_face

  subroutine map_P_to_y_face(p_data, b, y_face)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(out) :: y_face(b%nxp,b%nyp-1)

    y_face(:,:) = 0.5d0*(p_data(:,2:) + p_data(:,:b%nyp-1))

  end subroutine map_P_to_y_face

  subroutine dPdx(p_data, b, dp_dx)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(out) :: dp_dx(b%nxp-1,b%nyp)

    dp_dx(:,:) = (1.0d0/b%dx)*(p_data(2:,:) - p_data(:b%nxp-1,:))
    
  end subroutine dPdx

  subroutine dPdy(p_data, b, dp_dy)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(out) :: dp_dy(b%nxp,b%nyp-1)

    dp_dy(:,:) = (1.0d0/b%dy)*(p_data(:,2:) - p_data(:,:b%nyp-1))

  end subroutine dPdy

  pure function dPdx_2(p_data, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision :: dPdx_2(b%nxp,b%nyp)

    dPdx_2(2:b%nxp-1,:) = p_data(3:,:) - p_data(:b%nxp-2,:)
    if (b%cyclic) then
       dPdx_2(1,:) = p_data(2,:) - p_data(b%nxp-1,:)
    else
       dPdx_2(1,:) = 0.0d0
    endif
    dPdx_2(b%nxp,:) = dPdx_2(1,:)

  end function dPdx_2

  pure function dPdy_2(p_data, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)

    double precision :: dPdy_2(b%nxp,b%nyp)

    dPdy_2(:,1) = 0.0d0
    dPdy_2(:,2:b%nyp-1) = p_data(:,3:) - p_data(:,:b%nyp-2)
    dPdy_2(:,b%nyp) = 0.0d0

  end function dPdy_2

  subroutine jacobian(p, q, b, jacob, ajis, ajin)

    ! Jacobian advection term J(q,p) plus
    ! ah2fac times Del-4th(p) minus ah4fac times Del-6th(p).
    ! Uses Arakawa energy and enstrophy conserving 9-point Jacobian

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p(b%nxp,b%nyp,b%nl)
    double precision, intent(in) :: q(b%nxp,b%nyp,b%nl)
    double precision, intent(out) :: jacob(b%nxp,b%nyp,b%nl)
    double precision, intent(out) :: ajis(b%nl)
    double precision, intent(out) :: ajin(b%nl)

    double precision :: aj5sms,aj5smn,aj9sms,aj9smn

    double precision :: dp_dx(b%nxp,b%nyp,b%nl), dp_dy(b%nxp,b%nyp,b%nl)
    double precision :: dq_dx(b%nxp,b%nyp,b%nl), dq_dy(b%nxp,b%nyp,b%nl)
    double precision :: fac

    integer :: k

    fac = 1.0d0/(3.0d0*(2.0d0*b%dx)*(2.0d0*b%dy)*b%fnot)
    do k=1,b%nl
       dp_dx(:,:,k) = dPdx_2(p(:,:,k), b)
       dp_dy(:,:,k) = dPdy_2(p(:,:,k), b)
       dq_dx(:,:,k) = dPdx_2(q(:,:,k), b)
       dq_dy(:,:,k) = dPdy_2(q(:,:,k), b)
     
       ! Compute advective (Jacobian) and diffusive contributions to dq/dt
       ! Initialise dq/dt with Jacobian advection and diffusive terms
       jacob(:,:,k) = fac*( dq_dx(:,:,k)*dp_dy(:,:,k) - dq_dy(:,:,k)*dp_dx(:,:,k) + &
            dPdx_2(q(:,:,k)*dp_dy(:,:,k), b) - dPdy_2(q(:,:,k)*dp_dx(:,:,k), b) + &
            dPdy_2(p(:,:,k)*dq_dx(:,:,k), b) - dPdx_2(p(:,:,k)*dq_dy(:,:,k), b) )
    enddo
    do k=1,b%nl
       if (b%cyclic) then
          ! Southern boundary Jacobian term
          aj5sms = trapin(q(:,1,k)*dp_dx(:,2,k), b%nxp, b%dx)/(2.0d0*b%dx)
          aj9sms = trapin(q(:,2,k)*dp_dx(:,2,k), b%nxp, b%dx)/(2.0d0*b%dx)
          ajis(k) = ( aj5sms + 2.0d0*aj9sms )/6.0d0
          
          ! Northern boundary Jacobian term
          aj5smn = trapin(q(:,b%nyp,k)*dp_dx(:,b%nyp-1,k), b%nxp, b%dx)/(2.0d0*b%dx)
          aj9smn = trapin(q(:,b%nyp-1,k)*dp_dx(:,b%nyp-1,k), b%nxp, b%dx)/(2.0d0*b%dx)
          ajin(k) = -( aj5smn + 2.0d0*aj9smn )/6.0d0
       endif   
    enddo

  end subroutine jacobian

 
  ! C-grid advection scheme, second order accurate
  pure function c_grid_advection(vel_u, vel_v, t_data_x_face, t_data_y_face, b)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: vel_u(b%nxt+1,b%nyt)
    double precision, intent(in) :: vel_v(b%nxt,b%nyt+1)
    double precision, intent(in) :: t_data_x_face(b%nxt+1,b%nyt)
    double precision, intent(in) :: t_data_y_face(b%nxt,b%nyt+1)
    double precision :: c_grid_advection(b%nxt,b%nyt)

    double precision :: uT(b%nxt+1,b%nyt),vT(b%nxt,b%nyt+1)

    uT(:,:) = vel_u(:,:)*t_data_x_face(:,:)
    vT(:,:) = vel_v(:,:)*t_data_y_face(:,:)

    c_grid_advection(:,:) = (-1.0d0/b%dx)*( uT(2:b%nxt+1,:) - uT(1:b%nxt,:) ) + (-1.0d0/b%dy)*( vT(:,2:b%nyt+1) - vT(:,1:b%nyt) )

  end function c_grid_advection

  function del2_P_bc(p_data, b, alpha)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp,b%nl)
    double precision, intent(in) :: alpha

    double precision :: del2_P_bc(b%nxp,b%nyp,b%nl)

    integer :: k

    ! Compute Del-sqd(p) at previous time level for dissipation
    do k=1,b%nl
       del2_P_bc(:,:,k) = dP2dx2_bc(p_data(:,:,k), b, alpha) + dP2dy2_bc(p_data(:,:,k), b, alpha)
    enddo

  end function del2_P_bc

  function dP2dx2_bc(p_data, b, alpha)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(in) :: alpha

    double precision :: dP2dx2_bc(b%nxp,b%nyp)

    double precision :: fac
    fac = 1.0d0/(b%dx*b%dx)
    dP2dx2_bc(2:b%nxp-1,:) = fac*(p_data(1:b%nxp-2,:) + p_data(3:b%nxp,:) - 2.0d0*p_data(2:b%nxp-1,:))
    if (.not. b%cyclic) then
       dP2dx2_bc(1,:)     =  fac*alpha/( 0.5d0*alpha + 1.0d0 )*(p_data(2,:) - p_data(1,:))
       dP2dx2_bc(b%nxp,:) = -fac*alpha/( 0.5d0*alpha + 1.0d0 )*(p_data(b%nxp,:) - p_data(b%nxp-1,:))
    else
       dP2dx2_bc(1,:)     = fac*(p_data(b%nxp-1,:) + p_data(2,:) - 2.0d0*p_data(1,:))
       dP2dx2_bc(b%nxp,:) = dP2dx2_bc(1,:)
    endif
    dP2dx2_bc(:,1) = 0.0d0
    dP2dx2_bc(:,b%nyp) = 0.0d0
    
  end function dP2dx2_bc

  function dP2dy2_bc(p_data, b, alpha)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(in) :: alpha

    double precision :: dP2dy2_bc(b%nxp, b%nyp)

    double precision :: fac
    fac = 1.0d0/(b%dy*b%dy)
    dP2dy2_bc(:,2:b%nyp-1) = fac*(p_data(:,1:b%nyp-2) + p_data(:,3:b%nyp) - 2.0d0*p_data(:,2:b%nyp-1))
    dP2dy2_bc(:,1)     =  fac*alpha/( 0.5d0*alpha + 1.0d0 )*(p_data(:,2) - p_data(:,1))
    dP2dy2_bc(:,b%nyp) = -fac*alpha/( 0.5d0*alpha + 1.0d0 )*(p_data(:,b%nyp) - p_data(:,b%nyp-1))
    if (.not. b%cyclic) then
       dP2dy2_bc(1,:) = 0.0d0
       dP2dy2_bc(b%nxp,:) = 0.0d0
    endif

  end function dP2dy2_bc

  subroutine dPdx_2bc(p_data, b, alpha, dpdx)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(in) :: alpha
    double precision, intent(inout) :: dpdx(b%nxp,b%nyp)

    integer :: i

    do i=2,b%nxp-1
       dpdx(i,:) = (p_data(i+1,:) - p_data(i-1,:))/(2.0d0*b%dx)
    enddo
    if (.not. b%cyclic) then
       dpdx(1,:) = (p_data(2,:) - p_data(1,:))/(b%dx*(0.5d0*alpha + 1.0d0))
       dpdx(b%nxp,:) = (p_data(b%nxp,:) - p_data(b%nxp-1,:))/(b%dx*(0.5d0*alpha + 1.0d0))
    else
       dpdx(1,:) = (p_data(2,:) - p_data(b%nxp-1,:))/(2.0d0*b%dx)
       dpdx(b%nxp,:) = dpdx(1,:)
    endif
    dpdx(:,1) = 0.0d0
    dpdx(:,b%nyp) = 0.0d0
    
  end subroutine dPdx_2bc

  subroutine dPdy_2bc(p_data, b, alpha, dpdy)

    type(box_type), intent(in) :: b
    double precision, intent(in) :: p_data(b%nxp,b%nyp)
    double precision, intent(in) :: alpha
    double precision, intent(inout) :: dpdy(b%nxp,b%nyp)

    integer :: j

    do j=2,b%nyp-1
       dpdy(:,j) = (p_data(:,j+1) - p_data(:,j-1))/(2.0d0*b%dy)
    enddo
    dpdy(:,1)   = (p_data(:,2)   - p_data(:,1))/(b%dy*(0.5d0*alpha + 1.0d0))
    dpdy(:,b%nyp) = (p_data(:,b%nyp) - p_data(:,b%nyp-1))/(b%dy*(0.5d0*alpha + 1.0d0))
    if (.not. b%cyclic) then
       dpdy(1,:) = 0.0d0
       dpdy(b%nxp,:) = 0.0d0
    endif
    
  end subroutine dPdy_2bc


  double precision pure function sin_lat(b, y, y0, ymax, bar)
            
    type(box_type), intent(in) :: b
    double precision, intent(in) ::  y
    double precision, intent(in) ::  y0
    double precision, intent(in) ::  ymax
    double precision, intent(in) :: bar
   
    sin_lat = -0.50d0*sign(bar, b%fnot)*sin( (PI/2.0d0)*(y - y0)/(ymax - y0) )

  end function sin_lat

end module numerics
