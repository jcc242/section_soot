module input
  use soot_params, only: read_file

contains

subroutine read_input(filename)
  character(*) :: filename
  logical :: res
  inquire(file=trim(filename), exist=res)
  if (res) then
     write(*,*) "Reading input from: ",trim(filename)
     call read_file(trim(filename))
  else
     write(*,*) "Could not find file: ",trim(filename)
     error stop
  end if
end subroutine read_input

subroutine read_args (filename)
  character(len=100), intent(out) :: filename
  integer :: num_args

  num_args = command_argument_count()
  if (num_args .eq. 0) then
     write(*,*) "No input file given, using default params.toml"
     filename = "params/params.toml"
  else if (num_args .gt. 1) then
     write(*,*) "More than one argument given, only checking first for filename"
     call get_command_argument(1, filename)
  else
     call get_command_argument(1, filename)
  end if
end subroutine read_args


end module input
