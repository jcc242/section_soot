module section_soot
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, section_soot!"
  end subroutine say_hello
end module section_soot
