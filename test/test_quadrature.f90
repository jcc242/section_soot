module test_quadrature
  use Types, only: wp=>dp
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use quadrature, only: integral
  use dqag, only: quad_fun
  implicit none
  private

  public :: collect_quadrature

  contains

    subroutine collect_quadrature (testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
           new_unittest("integral", test_integral), &
           new_unittest("integral2", test_integral2) &
           ]
      
    end subroutine collect_quadrature

    subroutine test_integral (error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: val
      real(wp), parameter :: xmin = -1.0, xmax = 1.0, rel_tol=1.e-15, abs_tol=1.e-15
      procedure(quad_fun), pointer :: g=>x2
      val = integral(g, xmin, xmax, rel_tol, abs_tol)

      call check(error, val, 0.6666666_wp, thr=5.e-6_wp)
      if (allocated(error)) return

      
    end subroutine test_integral

    subroutine test_integral2 (error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: val
      real(wp), parameter :: xmin = -1.0, xmax = 1.0, rel_tol=1.e-15, abs_tol=1.e-15
      procedure(quad_fun), pointer :: g=>fexp
      val = integral(g, xmin, xmax, rel_tol, abs_tol)

      call check(error, val, 2.3504_wp, thr=5.e-5_wp)
      if (allocated(error)) return
      
    end subroutine test_integral2

    function x2 (x)
      real(wp) :: x2
      real(wp), intent(in) :: x

      x2 = x*x
      return
    end function x2

    function fexp (x)
      real(wp) :: fexp
      real(wp), intent(in) :: x

      fexp = exp(x)
      return
    end function fexp


end module test_quadrature
