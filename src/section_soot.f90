module section_soot
  use Types, only: wp=>dp
  use soot_params
  use collision
  implicit none
  private

  public :: calcCoagulation
contains
  subroutine calcCoagulation(coag_kernel)
    !! Calculates the mass based coagulation coeffcient for every two size bins.

    !! Geometric factor between neighboring sectional bins
    procedure(cvg_kernel), pointer, intent(in) :: coag_kernel
    !! Which collision function to use
    real(wp), dimension(NT, NT) :: coagSnkCoefs
    real(wp), dimension(NT, NT) :: coagSrcIds
    real(wp), dimension(NT, NT) :: coagSrcCoefs
    real(wp) :: beta


    integer :: i, j, nk

    if ((DIS_SECT .eq. 1) .or. (DIS_SECT .eq. 2)) then
       do i=1,ND
          do j=1,i
             ! Step one, find mass of collisions
             beta = coag_kernel(i,j)
             coagSnkCoefs(i,j) = beta/j
             coagSnkCoefs(j,i) = beta/i

             nk = i + j

             ! Step two, find bin new particle goes into
          end do
       end do
    end if


    

  end subroutine calcCoagulation
end module section_soot
