module coagulation
  use Types, only: wp=>dp
  use soot_params
  use collision
!!$  use bins
  use quadrature, only: integral, quad_fun
  implicit none
  private

  public :: calcCoagulation

  procedure(cvg_kernel), pointer :: m_coag_kernel

contains

  subroutine calcCoagulation(coag_kernel, coagSnkCoefs, coagSrcIds, coagSrcCoefs)
    !! Calculates the mass based coagulation coeffcient for every two size bins.

    !! Geometric factor between neighboring sectional bins
    procedure(cvg_kernel), pointer, intent(in) :: coag_kernel
    !! Which collision function to use
    real(wp), dimension(NT, NT), intent(out) :: coagSnkCoefs
    real(wp), dimension(NT, NT), intent(out) :: coagSrcIds
    real(wp), dimension(NT, NT), intent(out) :: coagSrcCoefs
    real(wp) :: beta

    integer :: i, j, nk, idx

    m_coag_kernel=>coag_kernel

    if ((DIS_SECT .eq. 1) .or. (DIS_SECT .eq. 3)) then
       do i=1,ND
          do j=1,i
             ! Step one, find mass of collisions
             beta = coag_kernel(real(i, kind = 8),real(j, kind = 8))
             coagSnkCoefs(i,j) = beta/j
             coagSnkCoefs(j,i) = beta/i

             nk = i + j

             ! Step two, find bin new particle goes into
             idx = find_bin(nk)
             if (idx .gt. 0) then
                ! We found a bin
                coagSrcIds(i,j) = idx
                coagSrcCoefs(i,j) = coagSnkCoefs(i,j) + coagSnkCoefs(j,i)
             else if (idx .eq. 0) then
                ! We are larger than the largest bin
                if (grow_bey_bound .eqv. .false.) then
                   ! Dump in the largest bin
                   coagSrcIds(i,j) = NT
                   coagSrcCoefs(i,j) = coagSnkCoefs(i,j) + coagSnkCoefs(j,i)
                end if
             end if
             coagSrcIds(j,i) = coagSrcIds(i,j)
             coagSrcCoefs(j,i) = coagSrcCoefs(i,j)
          end do
       end do
    end if ! Discrete bins

!!$    if ((DIS_SECT .eq. 1) .or. (DIS_SECT .eq. 2)) then
!!$      ! Sectional bins 
!!$       do i=ND+1:NT
!!$          if (i-1 .eq. ND) then
!!$             coagSnkCoefs(1,j) = integral(coag_kernel(1,))

  end subroutine calcCoagulation

  function igral1 (x) result(val)
    real(wp) :: val
    real(wp), intent(in) :: x
    val = m_coag_kernel(1.0_wp,x)/x
  end function igral1

end module coagulation
