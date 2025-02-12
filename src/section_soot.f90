module section_soot
  use Types, only: wp=>dp
  use soot_params
  implicit none

  public :: GDE

contains

  function GDE (t, conc, R, kappa, coagSnkCoefs, coagSrcIds, coagSrcCoefs, wall_coefs, l0_coefs, e0_coefs, w, l0, e0, m0, CONST MONOMER, COAG_OFF)
    real(wp), dimension(:), allocatable :: GDE
    real(wp), intent(in) :: t, R, kappa, coagSnkCoefs, coagSrcIds, coagSrcCoefs, wall_coefs, l0_coefs, e0_coefs, w, l0, e0, m0, CONST MONOMER, COAG_OFF
    real(wp), dimension(:), intent(in) :: conc

    integer :: i,j
    real(wp) :: pre_factor

    allocate(GDE(NT))
    GDE = 0.0

    if (COAG_OFF .eqv. .true.) then
       coagSnkCoefs(2:end,2:end) = 0.0
       coagSrcCoefs(2:end,2:end) = 0.0
    end if

    do i=1,NT
       do j=1,i
          if (i .eq. j) then
             pre_factor = 0.5
          else
             pre_factor = 1.0
          end if
          if (coagSrcIds(i,j) .ne. 0) then
             GDE(coagSrcIds(i,j)) = GDE(coagSrcIds(i,j)) + pre_factor*coagSrcCoefs(i,j)*conc(i)*conc(j)
             if (coagSrcIds(i,j) .lt. NT) then
                GDE(coagSrcIds(i,j)+1) = GDE(coagSrcIds(i,j)+1) + pre_factor*(coagSnkCoefs(i,j)+coagSnkCoefs(j,i)-coagSrcCoefs(i,j))*conc(i)*conc(j)
             end if
          end if
       end do
    end do

    
          
    

  end function GDE

end module section_soot

