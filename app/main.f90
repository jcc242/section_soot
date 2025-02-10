program main
  use Types, only: wp=>dp
  use soot_params, only: read_input
  use section_soot, only: calcCoagulation
  use collision, only: getKernel, cvg_kernel
  implicit none

  real(wp) :: T, P, rho, m_mono

  procedure(cvg_kernel), pointer :: kernel
  call read_input
  kernel => getKernel(T, P, rho, m_mono) 
  call calcCoagulation(kernel)
end program main
