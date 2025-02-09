program main
  use soot_params, only: read_input
  use section_soot, only: calcCoagulation
  implicit none

  call read_input
  call calcCoagulation()
end program main
