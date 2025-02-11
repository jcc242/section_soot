module quadrature

  public :: integral

  abstract interface
     function quad_fun(i, j) result(beta)
       use Types, only: wp=>dp
       integer, intent(in) :: i,j
       real(wp) :: beta
     end function cvg_kernel
  end interface

contains

 function integral (fun, xmin, xmax, rel_tol, abs_tol) result(q)
   real(wp), intent(in) :: xmin, xmax, rel_error, abs_error
   procedure(quad_fun), pointer :: fun
   real(wp) :: q, abs_error
   external dqag
   integer :: num_eval, err, last
   integer, parameter :: limit = 100, key = 6, len_w=limit*4
   integer, dimension(limit) :: iwork
   real(wp), dimension(len_w) :: work

   call dqag(fun, xmin, xmax, abs_tol, rel_tol, key, q&
        &abs_error, num_eval, err, limit, len_w, last, iwork, work)

 end function integral
 

end module quadrature
