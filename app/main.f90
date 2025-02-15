program main
  use Types, only: wp=>dp
  use input, only: read_args, read_input
  use coagulation, only: calcCoagulation
  use collision, only: getKernel, cvg_kernel
  use sink, only: calcSink
  use bins, only: coagSnkCoefs, coagSrcIds, coagSrcCoefs
  use soot_params, only: temperature, pressure, rho, m_mono, kernel_type=>kernel
  implicit none
  procedure(cvg_kernel), pointer :: kernel
  character(len=100) :: filename

  call read_args(filename)
  call read_input(filename)

  kernel => getKernel(kernel_type, temperature, pressure, rho, m_mono) 
  call calcCoagulation(kernel, coagSnkCoefs, coagSrcIds, coagSrcCoefs)

  call calcSink(kernel)
end program main


