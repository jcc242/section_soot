!> Driver for unit testing
program tester
  use Types, only: stderr
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type, &
       & select_suite, run_selected, get_argument
  use test_quadrature, only : collect_quadrature
  implicit none
  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
       new_testsuite("quadrature", collect_quadrature) &
       ]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then
     is = select_suite(testsuites, suite_name)
     if (is > 0 .and. is <= size(testsuites)) then
        if (allocated(test_name)) then
           write(stderr, fmt) "Suite:", testsuites(is)%name
           call run_selected(testsuites(is)%collect, test_name, stderr, stat)
           if (stat < 0) then
              error stop 1
           end if
        else
           write(stderr, fmt) "Testing:", testsuites(is)%name
           call run_testsuite(testsuites(is)%collect, stderr, stat)
        end if
     else
        write(stderr, fmt) "Available testsuites"
        do is = 1, size(testsuites)
           write(stderr, fmt) "-", testsuites(is)%name
        end do
        error stop 1
     end if
  else
     do is = 1, size(testsuites)
        write(stderr, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, stderr, stat)
     end do
  end if

  if (stat > 0) then
     write(stderr, '(i0, 1x, a)') stat, "test(s) failed!"
     error stop 1
  end if

end program tester
