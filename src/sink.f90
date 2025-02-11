module sink
  use Types, only: wp=>dp
  use soot_params
  use collision
  use evap
  implicit none

contains

  subroutine calcSink(coag_kernel)

    real(wp), dimension(:), allocatable :: wall_coeffs, l0_coeffs, e0_coeffs
    procedure(cvg_kernel), pointer, intent(in) :: coag_kernel
    !! Which collision function to use
    integer :: i

    allocate(wall_coeffs(NT))
    allocate(l0_coeffs(NT))
    allocate(e0_coeffs(NT))

    if (ND .gt. 0) then
       do i=1,ND
          wall_coeffs(i) = 1.0_wp/((i)**(1.0_wp/3.0_wp))
       end do
    end if
    if ((ND+1) .le. NT) then
       wall_coeffs(ND+1) = 1.5_wp/s_len(ND+1) * (r_limits(ND+1)**(2.0_wp/3.0_wp) - l_limit**(2.0_wp/3.0_wp))
    end if
    if ((ND+2) .le. NT) then
       do i=(ND+2),NT
          wall_coeffs(i) = 1.5_wp/s_len(i)*(r_limits(i)**(2.0_wp/3.0_wp) - r_limits(i-1)**(2.0_wp/3.0_wp))
       end do
    end if

    if (ND .gt. 0) then
       do i=1,ND
          l0_coeffs(i) = 1.0_wp/((i)**(1.0_wp/2.0_wp))
       end do
    end if
    if ((ND+1) .le. NT) then
       l0_coeffs(ND+1) = 2/s_len(ND+1) * (r_limits(ND+1)**(1.0_wp/2.0_wp) - l_limit**(1.0_wp/2.0_wp))
    end if
    if ((ND+1) .le. NT) then
       do i=(ND+2),NT
          l0_coeffs(i) = 2.0_wp/s_len(i)*(r_limits(i)**(1.0_wp/2.0_wp) - r_limits(i-1)**(1.0_wp/2.0_wp))
       end do
    end if

    e0_coeffs(1) = 0

    ! Discrete
    if (ND .gt. 1) then
       e0_coeffs(2) = (1.0_wp/2.0_wp)*calcEvap(coag_kernel, A, 2)*2.0_wp*0.5_wp
    end if
    if (ND .gt. 2) then
       do i=3,ND
          e0_coeffs(i) = 1.0_wp/i*calcEvap(coag_kernel, A, i)
       end do
    end if

    ! Sectional


  end subroutine calcSink


end module sink
