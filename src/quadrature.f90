module quadrature
  use Types, wp=>dp
  use dqag, only: qag, quad_fun

  public :: integral


contains

  function integral (fun, xmin, xmax, rel_tol, abs_tol) result(q)
    real(wp), intent(in) :: xmin, xmax, rel_tol, abs_tol
    procedure(quad_fun), pointer :: fun
    real(wp) :: q, abs_error
    integer :: num_eval, err
    integer, parameter :: key = 6

    call qag(fun, xmin, xmax, abs_tol, rel_tol, key, q,&
         &abs_error, num_eval, err)

  end function integral


end module quadrature
