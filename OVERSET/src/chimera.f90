!===============================================================================
!  definisce una nuova cella di tipo chimera
!===============================================================================
subroutine newchi(igr,ibl,i,j,k,td,nd,md)
use prec
use moddef, only : donor
use modpar, only : donmax,blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k,td,nd
integer (kind=I4P) :: n
type (donor) :: md(donmax)

! se non c'è almeno un donatore si ferma tutto
if (nd.lt.1) stop 'nd <= 0 in newchi'

if (allocated(blochi(ibl,igr)%cella(i,j,k)%q)) &
   deallocate(blochi(ibl,igr)%cella(i,j,k)%q)

allocate(blochi(ibl,igr)%cella(i,j,k)%q(nd))

blochi(ibl,igr)%cella(i,j,k)%t = td
blochi(ibl,igr)%cella(i,j,k)%n = nd
do n=1,nd
   blochi(ibl,igr)%cella(i,j,k)%q(n) = md(n)
end do

end subroutine

!===============================================================================
!  cancella una cella di tipo chimera (per ora non serve, ma non si sa mai)
!===============================================================================
subroutine delchi(igr,ibl,i,j,k)
use prec
use modpar, only : stdflag,blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k

blochi(ibl,igr)%cella(i,j,k)%t = stdflag
blochi(ibl,igr)%cella(i,j,k)%n = 0
if (allocated(blochi(ibl,igr)%cella(i,j,k)%q)) &
   deallocate(blochi(ibl,igr)%cella(i,j,k)%q)

end subroutine

!===============================================================================
!  copia una cella chimera
!===============================================================================
subroutine copychi(igr,ibl,ia,ja,ka,ib,jb,kb)
use prec
use modpar, only : stdflag,blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,ia,ja,ka,ib,jb,kb
integer (kind=I4P) :: n,nd

! se la cella di origine non è chimera ferma tutto
if (blochi(ibl,igr)%cella(ia,ja,ka)%t .eq. stdflag) then
   stop 'NON chimera in copychi'
end if

blochi(ibl,igr)%cella(ib,jb,kb)%t = blochi(ibl,igr)%cella(ia,ja,ka)%t
blochi(ibl,igr)%cella(ib,jb,kb)%n = blochi(ibl,igr)%cella(ia,ja,ka)%n

if (allocated(blochi(ibl,igr)%cella(ib,jb,kb)%q)) &
   deallocate(blochi(ibl,igr)%cella(ib,jb,kb)%q)

nd = blochi(ibl,igr)%cella(ia,ja,ka)%n 
allocate(blochi(ibl,igr)%cella(ib,jb,kb)%q(nd))
do n=1,nd
   blochi(ibl,igr)%cella(ib,jb,kb)%q(n) = blochi(ibl,igr)%cella(ia,ja,ka)%q(n)
end do

end subroutine
