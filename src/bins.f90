module bins
  use Types, only: wp=>dp
  use constants
  private
  real(wp) :: l_limit
  !! Left/lower limit of sectional bins
  real(wp), dimension(:), allocatable :: r_limits
  !! Right/upper limits of sectional bins

  real(wp), dimension(:), allocatable :: s_len
  !! Length of each bin, in number of molecules
  real(wp), dimension(:), allocatable :: s_avg
  !! Average particle mass in each bin, in kg
  real(wp), dimension(:), allocatable :: convert_num_to_log
  !! Conversion factor from number per bin to dN/dlog10(Dp)
  real(wp), dimension(:), allocatable :: dpBins

  public :: find_bin, allocate_bins, initialize_bins

contains

  subroutine allocate_bins (NT)
    integer, intent(in) :: NT

    allocate(r_limits(NT))
    r_limits = 0.0_wp
    allocate(s_len(NT))
    s_len = 0.0_wp
    allocate(s_avg(NT))
    s_avg = 0.0_wp
    allocate(convert_num_to_log(NT))
    convert_num_to_log = 0.0_wp
    allocate(dpBins(NT))
    dpBins = 0.0_wp
  end subroutine allocate_bins

  subroutine initialize_bins (DIS_SECT, ND, NS, NT, geo_factor, v_mono)
    integer, intent(in) :: DIS_SECT, ND, NS, NT
    real(wp), intent(in) :: geo_factor, v_mono
    integer i

    l_limit = 0

    if ((DIS_SECT .eq. 1) .or. (DIS_SECT .eq. 3)) then
       do i=1,ND
          r_limits(i) = i
          s_avg(i) = i
          convert_num_to_log(i) = 3.0_wp/(log10((i*0.5_wp)/(i-0.5_wp)))
       end do
       s_len = 1.0
       l_limit = ND+0.5 
    end if

    if (NS .gt. 0) then
       r_limits(ND+1) = l_limit+geo_factor
       s_len(ND+1) = r_limits(ND+1)-l_limit
       s_avg(ND+1) = s_len(ND+1)/(log(geo_factor))
       convert_num_to_log(ND+1) = 3.0_wp/log10(geo_factor)
    end if
    if (NS .gt. 1) then
       do i = ND + 2, NT
          r_limits(i) = r_limits(i-1)*geo_factor
          s_len(i) = r_limits(i)-r_limits(i-1)
          s_avg = s_len(i)/log(geo_factor)
          convert_num_to_log(i) = 3.0_wp/log10(geo_factor)
       end do
    end if

    dpBins = (s_avg*v_mono)**(1.0_wp/3.0_wp)*(6.0_wp/PI)**(1.0_wp/3.0_wp)
  end subroutine initialize_bins



  function find_bin (qty)
    !! Given a quantity, figure out which bin it goes.
    !! If this returns 0, then it is bigger than the largest bin
    integer :: find_bin
    integer, intent(in) :: qty
    integer i

    find_bin = 0
    do i=1,size(r_limits)
       if (qty .le. r_limits(i)) then
          find_bin = i
          exit
       end if
    end do
    return
  end function find_bin

  subroutine set_l_limit (limit)
    real(wp), intent(in) :: limit
    l_limit = limit
  end subroutine set_l_limit


end module bins
