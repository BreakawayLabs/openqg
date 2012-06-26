module bicubic

  use box, only: box_type
  use grid, only: grid_type

  implicit none
  
  private

  type bcudata

     double precision, allocatable :: stfn(:,:,:)
     double precision, allocatable :: stfn_temp(:,:)
     double precision, allocatable :: st(:,:)
     double precision :: stinv(16,16)

  end type bcudata

  public bcudata
  public bcuini
  public new_bicubic

contains

  subroutine new_bicubic(g, go, u, ux, uy, uxy, bcu, output)

    type(grid_type), intent(in) :: g
    type(box_type), intent(in) :: go
    double precision, intent(in) :: u(g%nxpoar,g%nypoar)
    double precision, intent(in) :: ux(g%nxpoar,g%nypoar)
    double precision, intent(in) :: uy(g%nxpoar,g%nypoar)
    double precision, intent(in) :: uxy(g%nxpoar,g%nypoar)
    type(bcudata), intent(in) :: bcu
    double precision, intent(out) :: output(go%nxp,go%nyp)

    double precision :: d(16, g%nxtoar,g%nytoar)
    double precision :: out_temp2((g%ndxr+1), (g%ndxr+1), g%nxtoar, g%nytoar)
 
    d = bc1(u, ux, uy, uxy, g)

    out_temp2 = bc2(d, bcu, g)

    output = bc3(g, go, out_temp2)

  end subroutine new_bicubic

  function bc3(g, go, out_temp2)
    type(grid_type), intent(in) :: g
    type(box_type), intent(in) :: go
    double precision, intent(in) :: out_temp2((g%ndxr+1), (g%ndxr+1), g%nxtoar, g%nytoar)

    double precision :: bc3(go%nxp,go%nyp)

    integer :: i, j, iout, jout

    do j=1,g%nytoar
       jout = 1 + (j-1)*g%ndxr
       do i=1,g%nxtoar
          iout = 1 + (i-1)*g%ndxr
          bc3(iout:iout+g%ndxr,jout:jout+g%ndxr) = out_temp2(:,:,i,j)
       enddo
    enddo

  end function bc3

  function bc2(d, bcu, g)
    type(grid_type), intent(in) :: g
    type(bcudata), intent(in) :: bcu
    double precision, intent(in) :: d(16,g%nxtoar,g%nytoar)
    double precision :: bc2((g%ndxr+1), (g%ndxr+1), g%nxtoar, g%nytoar)

    bc2 = reshape(matmul(bcu%st(:,:), reshape(d(:,:,:), (/16, g%nxtoar*g%nytoar/))), shape(bc2))

  end function bc2

  function bc1(u, ux, uy, uxy, g)
    type(grid_type), intent(in) :: g
    double precision, intent(in) :: u(g%nxpoar,g%nypoar)
    double precision, intent(in) :: ux(g%nxpoar,g%nypoar)
    double precision, intent(in) :: uy(g%nxpoar,g%nypoar)
    double precision, intent(in) :: uxy(g%nxpoar,g%nypoar)

    double precision :: bc1(16, g%nxtoar,g%nytoar)
    integer :: k, ip, jq

    k=1
    do jq=0,1
       do ip=0,1
          bc1(k,   :,:) = u  (1+ip:g%nxtoar+ip,1+jq:g%nytoar+jq)
          bc1(4+k, :,:) = ux (1+ip:g%nxtoar+ip,1+jq:g%nytoar+jq)
          bc1(8+k, :,:) = uy (1+ip:g%nxtoar+ip,1+jq:g%nytoar+jq)
          bc1(12+k,:,:) = uxy(1+ip:g%nxtoar+ip,1+jq:g%nytoar+jq)
          k = k + 1
       enddo
    enddo

  end function bc1


  type(bcudata) pure function bcuini(g)

    ! Computes transformation matrices between data values and bicubic
    ! fit coefficients, and the coefficient multipliers at the finer
    ! resolution gridpoints, and makes them available via a common block

    type(grid_type), intent(in) :: g

    double precision :: ss,tt
    integer :: ip,jp,i,j,k

    allocate(bcuini%stfn(0:g%ndxr,0:g%ndxr,16))
    allocate(bcuini%stfn_temp((g%ndxr+1)*(g%ndxr+1),16))
    allocate(bcuini%st((g%ndxr+1)*(g%ndxr+1), 16))

    ! Setup the coefficient multipliers on the refined grid
    ! =====================================================
    ! ss and tt are normalised coordinates within the gridcell
    ! Loop through the fine mesh gridpoints (ip, jp)
    do jp=0,g%ndxr
       tt = dble(jp)/dble(g%ndxr)
       do ip=0,g%ndxr
          ss = dble(ip)/dble(g%ndxr)
          ! Loop through the bicubic terms, saving
          ! the results as a 16-element vector
          k = 0
          do j=0,3
             do i=0,3
                k = k + 1
                bcuini%stfn(ip,jp,k) = ss**(i)*tt**(j)
                bcuini%stfn_temp(ip + jp*(1+g%ndxr) + 1,k) = ss**(i)*tt**(j)
             enddo
          enddo
       enddo
    enddo

    bcuini%stinv(:,1) =  (/1, 0, -3,  2, 0, 0,  0,  0, -3,  0,  9, -6,  2,  0, -6,  4/)
    bcuini%stinv(:,2) =  (/0, 0,  3, -2, 0, 0,  0,  0,  0,  0, -9,  6,  0,  0,  6, -4/)
    bcuini%stinv(:,3) =  (/0, 0,  0,  0, 0, 0,  0,  0,  3,  0, -9,  6, -2,  0,  6, -4/)
    bcuini%stinv(:,4) =  (/0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  9, -6,  0,  0, -6,  4/)
    bcuini%stinv(:,5) =  (/0, 1, -2,  1, 0, 0,  0,  0,  0, -3,  6, -3,  0,  2, -4,  2/)
    bcuini%stinv(:,6) =  (/0, 0, -1,  1, 0, 0,  0,  0,  0,  0,  3, -3,  0,  0, -2,  2/)
    bcuini%stinv(:,7) =  (/0, 0,  0,  0, 0, 0,  0,  0,  0,  3, -6,  3,  0, -2,  4, -2/)
    bcuini%stinv(:,8) =  (/0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  3,  0,  0,  2, -2/)
    bcuini%stinv(:,9) =  (/0, 0,  0,  0, 1, 0, -3,  2, -2,  0,  6, -4,  1,  0, -3,  2/)
    bcuini%stinv(:,10) = (/0, 0,  0,  0, 0, 0,  3, -2,  0,  0, -6,  4,  0,  0,  3, -2/)
    bcuini%stinv(:,11) = (/0, 0,  0,  0, 0, 0,  0,  0, -1,  0,  3, -2,  1,  0, -3,  2/)
    bcuini%stinv(:,12) = (/0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  2,  0,  0,  3, -2/)
    bcuini%stinv(:,13) = (/0, 0,  0,  0, 0, 1, -2,  1,  0, -2,  4, -2,  0,  1, -2,  1/)
    bcuini%stinv(:,14) = (/0, 0,  0,  0, 0, 0, -1,  1,  0,  0,  2, -2,  0,  0, -1,  1/)
    bcuini%stinv(:,15) = (/0, 0,  0,  0, 0, 0,  0,  0,  0, -1,  2, -1,  0,  1, -2,  1/)
    bcuini%stinv(:,16) = (/0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  1, -1,  0,  0, -1,  1/)

    bcuini%st(:,:) = matmul(bcuini%stfn_temp(:,:), bcuini%stinv(:,:))

  end function bcuini

end module bicubic
