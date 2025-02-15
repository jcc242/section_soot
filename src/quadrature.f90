module quadrature
  use Types, wp=>dp
!!$  use dqag, only: qag, quad_fun

  public :: integral

  interface
     function quad_fun(x) result(beta)
       use Types, wp=>dp
       real(wp), intent(in) :: x
       real(wp) :: beta
     end function quad_fun
  end interface

  type quad_fun_ptr
     procedure(quad_fun), nopass, pointer :: integrand
  end type quad_fun_ptr

contains

  function integral (fun, xmin, xmax, rel_tol, abs_tol) result(q)
    real(wp), intent(in) :: xmin, xmax, rel_tol, abs_tol
!!$    procedure(quad_fun), pointer :: fun
!!$    type(quad_fun_ptr), target :: fun
    real( kind = 4), external :: fun
    real(wp) :: q, abs_error
    integer :: num_eval, err
    integer, parameter :: key = 6
    external qag

    call qag(fun, xmin, xmax, abs_tol, rel_tol, key, q,&
         &abs_error, num_eval, err)

  end function integral


end module quadrature
