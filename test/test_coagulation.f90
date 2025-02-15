module test_coagulation
  use Types, only: wp=>dp
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use collision
  use coagulation
  use quadrature
  use bins, only: coagSnkCoefs, coagSrcIds, coagSrcCoefs
  use input, only: read_input
  use soot_params
  implicit none
  private

  public :: collect_coagulation

  contains

    subroutine collect_coagulation (testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
           new_unittest("test_coag_discrete", test_coag_discrete) &
           ]
      
    end subroutine collect_coagulation

    subroutine test_coag_discrete (error)
      type(error_type), allocatable, intent(out) :: error
      character(len=100) :: filename
      procedure(cvg_kernel), pointer :: kernel
      integer :: i,j
      real(wp), dimension(10,10) :: truth_coagSrcIds
      real(wp), dimension(10,10) :: truth_coagSrcCoefs
      real(wp), dimension(10,10) :: truth_coagSnkCoefs

      filename = "params/test_coag_discrete.toml"
      truth_coagSrcIds = transpose(reshape((/2, 3, 4, 5, 6, 7, 8, 9,10,10&
           &, 3, 4, 5, 6, 7, 8, 9,10,10,10&
           &, 4, 5, 6, 7, 8, 9,10,10,10,10&
           &, 5, 6, 7, 8, 9,10,10,10,10,10&
           &, 6, 7, 8, 9,10,10,10,10,10,10&
           &, 7, 8, 9,10,10,10,10,10,10,10&
           &, 8, 9,10,10,10,10,10,10,10,10&
           &, 9,10,10,10,10,10,10,10,10,10&
           &,10,10,10,10,10,10,10,10,10,10&
           &,10,10,10,10,10,10,10,10,10,10/),shape(truth_coagSrcIds)))

      truth_coagSrcCoefs = transpose(reshape((/&
           &8.4652e-10,7.0203e-10,6.8710e-10,7.0004e-10,7.2232e-10,7.4827e-10&
           &,7.7567e-10,8.0352e-10,8.3135e-10,8.5893e-10&
           &,7.0203e-10,4.7509e-10,4.1561e-10,3.9400e-10,3.8650e-10,3.8561e-10&
           &,3.8823e-10,3.9288e-10,3.9875e-10,4.0538e-10&
           &,6.8710e-10,4.1561e-10,3.3887e-10,3.0597e-10,2.8957e-10,2.8102e-10&
           &,2.7677e-10,2.7509e-10,2.7505e-10,2.7609e-10&
           &,7.0004e-10,3.9400e-10,3.0597e-10,2.6663e-10,2.4557e-10,2.3325e-10&
           &,2.2573e-10,2.2112e-10,2.1838e-10,2.1691e-10&
           &,7.2232e-10,3.8650e-10,2.8957e-10,2.4557e-10,2.2139e-10,2.0666e-10&
           &,1.9715e-10,1.9080e-10,1.8650e-10,1.8359e-10&
           &,7.4827e-10,3.8561e-10,2.8102e-10,2.3325e-10,2.0666e-10,1.9018e-10&
           &,1.7926e-10,1.7172e-10,1.6637e-10,1.6251e-10&
           &,7.7567e-10,3.8823e-10,2.7677e-10,2.2573e-10,1.9715e-10,1.7926e-10&
           &,1.6725e-10,1.5881e-10,1.5268e-10,1.4814e-10&
           &,8.0352e-10,3.9288e-10,2.7509e-10,2.2112e-10,1.9080e-10,1.7172e-10&
           &,1.5881e-10,1.4964e-10,1.4290e-10,1.3782e-10&
           &,8.3135e-10,3.9875e-10,2.7505e-10,2.1838e-10,1.8650e-10,1.6637e-10&
           &,1.5268e-10,1.4290e-10,1.3565e-10,1.3013e-10&
           &,8.5893e-10,4.0538e-10,2.7609e-10,2.1691e-10,1.8359e-10,1.6251e-10&
           &,1.4814e-10,1.3782e-10,1.3013e-10,1.2424e-10&
      &/),shape(truth_coagSrcCoefs)))

      truth_coagSnkCoefs = transpose(reshape((/&
           &4.2326e-10,2.3401e-10,1.7177e-10,1.4001e-10,1.2039e-10,1.0690e-10&
           &,9.6958e-11,8.9280e-11,8.3135e-11,7.8084e-11&
           &,4.6802e-10,2.3754e-10,1.6624e-10,1.3133e-10,1.1043e-10,9.6404e-11&
           &,8.6274e-11,7.8575e-11,7.2499e-11,6.7564e-11&
           &,5.1532e-10,2.4936e-10,1.6943e-10,1.3113e-10,1.0859e-10,9.3675e-11&
           &,8.3032e-11,7.5025e-11,6.8761e-11,6.3714e-11&
           &,5.6003e-10,2.6266e-10,1.7484e-10,1.3332e-10,1.0914e-10,9.3299e-11&
           &,8.2084e-11,7.3706e-11,6.7194e-11,6.1975e-11&
           &,6.0193e-10,2.7607e-10,1.8098e-10,1.3643e-10,1.1069e-10,9.3937e-11&
           &,8.2147e-11,7.3385e-11,6.6607e-11,6.1198e-11&
           &,6.4137e-10,2.8921e-10,1.8735e-10,1.3995e-10,1.1272e-10,9.5089e-11&
           &,8.2736e-11,7.3593e-11,6.6546e-11,6.0942e-11&
           &,6.7871e-10,3.0196e-10,1.9374e-10,1.4365e-10,1.1501e-10,9.6525e-11&
           &,8.3626e-11,7.4110e-11,6.6797e-11,6.0997e-11&
           &,7.1424e-10,3.1430e-10,2.0007e-10,1.4741e-10,1.1742e-10,9.8124e-11&
           &,8.4697e-11,7.4819e-11,6.7245e-11,6.1252e-11&
           &,7.4822e-10,3.2625e-10,2.0628e-10,1.5119e-10,1.1989e-10,9.9819e-11&
           &,8.5882e-11,7.5651e-11,6.7823e-11,6.1641e-11&
           &,7.8084e-10,3.3782e-10,2.1238e-10,1.5494e-10,1.2240e-10,1.0157e-10&
           &,8.7139e-11,7.6565e-11,6.8490e-11,6.2122e-11&
      &/),shape(truth_coagSnkCoefs)))



      call read_input(trim(filename))
      kernel => getKernel(2, temperature, pressure, rho, m_mono)
      call calcCoagulation(kernel, coagSnkCoefs, coagSrcIds, coagSrcCoefs)

      do i = 1,NT
         do j = 1,NT
            call check(error, coagSrcIds(i,j), truth_coagSrcIds(i,j), thr=5.e-14_wp)
            call check(error, coagSrcCoefs(i,j), truth_coagSrcCoefs(i,j), thr=5.e-14_wp)
            call check(error, coagSnkCoefs(i,j), truth_coagSnkCoefs(i,j), thr=5.e-14_wp)
         end do
      end do
      
      if (allocated(error)) return
    end subroutine test_coag_discrete

end module test_coagulation
