module test_collision
  use Types, only: wp=>dp
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use collision
  use input, only: read_input
  use soot_params, only: temperature, pressure, rho, m_mono
  implicit none
  private

  public :: collect_collision

  contains

    subroutine collect_collision (testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
           new_unittest("test_fuchs", test_fuchs), &
           new_unittest("test_free", test_free) &
           ]
      
    end subroutine collect_collision

    subroutine test_free (error)
      type(error_type), allocatable, intent(out) :: error
      integer :: i,j,idx
      real(wp) :: val
      character(len=100) :: filename
      procedure(cvg_kernel), pointer :: kernel
      real(wp), dimension(10), parameter :: truth = (/&
           &4.23263e-10,&
           &4.68024e-10,&
           &4.75097e-10,&
           &5.1533e-10,&
           &4.98736e-10,&
           &5.08313e-10,&
           &5.60039e-10,&
           &5.25339e-10,&
           &5.24541e-10,&
           &5.33278e-10/)


      filename = "params/test_free.toml"

      call read_input(trim(filename))

      kernel => getKernel(2, temperature, pressure, rho, m_mono)

      idx = 1
      do i=1,ND
         do j=1,i
            val = kernel(i,j)
            call check(error, val, truth(idx), thr=5.e-6_wp)
            if (allocated(error)) return
            idx = idx + 1
         end do
      end do
      
    end subroutine test_free

    subroutine test_fuchs (error)
      type(error_type), allocatable, intent(out) :: error
      integer :: i,j, idx
      real(wp) :: val
      character(len=100) :: filename
      procedure(cvg_kernel), pointer :: kernel
      real(wp), dimension(10), parameter :: truth = (/&
           &4.2326e-10,&
           &4.68019e-10,&
           &4.7509e-10,&
           &5.15323e-10,&
           &4.98727e-10,&
           &5.08301e-10,&
           &5.6003e-10,&
           &5.25327e-10,&
           &5.24527e-10,&
           &5.33262e-10/)
      filename = "params/test_fuchs.toml"

      call read_input(trim(filename))

      kernel => getKernel(1, temperature, pressure, rho, m_mono)

      idx = 1
      do i=1,ND
         do j=1,i
            val = kernel(i,j)
            call check(error, val, truth(idx), thr=5.e-12_wp)
            if (allocated(error)) return
            idx = idx + 1
         end do
      end do
      

      
    end subroutine test_fuchs

end module test_collision
