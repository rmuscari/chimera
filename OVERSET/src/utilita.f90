!===============================================================================
!  subroutines di utilità generale
!===============================================================================

!===============================================================================
! restituisce una unità logica libera utilizzabile per open
subroutine getunit(Free_Unit)
use prec
implicit none
integer(kind=I4P) :: Free_Unit ! Free logic unit.
integer(kind=I4P) :: n1        ! Counter.
integer(kind=I4P) :: ios       ! Inquiring flag.
logical(kind=I4P) :: lopen     ! Inquiring flag.

Free_Unit = -1      ! initializing free logic unit
n1=30               ! initializing counter
do
   inquire(unit=n1,opened=lopen,iostat=ios) ! verify logic units
   if (ios==0) then
      if (.not.lopen) then
         Free_Unit = n1                      ! assignment of free unit
         exit
      endif
   endif
   n1=n1+1                                      ! updating counter
enddo
end subroutine
