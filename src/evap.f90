module evap
  use Types, only: wp=>dp
  use collision


  public :: calcEvap
contains

  function calcEvap (coag_kernel, A, k) result(val)
    procedure(cvg_kernel), pointer, intent(in) :: coag_kernel
    real(wp), intent(in) :: A
    integer, intent(in) :: k
    real(wp) :: val
    val = coag_kernel(1,k-1)*exp(1.5_wp*A*(k**(2.0_wp/3.0_wp)-(k-1)**(2.0_wp/3.0_wp)));
    return
  end function calcEvap
end module evap
