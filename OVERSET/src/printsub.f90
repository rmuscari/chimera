!===============================================================================
!  preprocessore CHIMERA per XSHIP - scrittura matrici CC
!===============================================================================
subroutine printcc(igr)
use prec
use moddef
use modpar, only : nccnat,blochi,nbl,fout,nout
implicit none
integer (kind=I4P), intent(in) :: igr
character (len=50) :: ftmp
integer (kind=I4P) :: ibl,i,j,k,m,n
integer (kind=I4P) :: ii0,iin,jj0,jjn,kk0,kkn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  apertura file di output
write(ftmp,'(a,i2.2)') trim(fout),igr
open(nout,file=ftmp,form='unformatted',status='unknown')
rewind(nout)

!  scrittura numero di blocchi
write(nout) nbl

!  scrittura dimensioni dei blocchi
do ibl=1,nbl
   write(nout) blochi(ibl,igr)%ni, blochi(ibl,igr)%nj, blochi(ibl,igr)%nk, &
               blochi(ibl,igr)%gc(1:6)
end do

!  scrittura delle matrici con il tipo
do ibl=1,nbl
   ii0 = - blochi(ibl,igr)%gc(1) + 1
   iin = blochi(ibl,igr)%ni + blochi(ibl,igr)%gc(2)
   jj0 = - blochi(ibl,igr)%gc(3) + 1
   jjn = blochi(ibl,igr)%nj + blochi(ibl,igr)%gc(4)
   kk0 = - blochi(ibl,igr)%gc(5) + 1
   kkn = blochi(ibl,igr)%nk + blochi(ibl,igr)%gc(6)
   write(nout) (((blochi(ibl,igr)%cella(i,j,k)%t,&
                  i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
end do

!  per le celle chimera (tipo>nccnat) scrittura dei dati
do ibl=1,nbl
   ii0 = - blochi(ibl,igr)%gc(1) + 1
   iin = blochi(ibl,igr)%ni + blochi(ibl,igr)%gc(2)
   jj0 = - blochi(ibl,igr)%gc(3) + 1
   jjn = blochi(ibl,igr)%nj + blochi(ibl,igr)%gc(4)
   kk0 = - blochi(ibl,igr)%gc(5) + 1
   kkn = blochi(ibl,igr)%nk + blochi(ibl,igr)%gc(6)
   do k=kk0,kkn
      do j=jj0,jjn
         do i=ii0,iin
            if (blochi(ibl,igr)%cella(i,j,k)%t.gt.nccnat) then
               m = blochi(ibl,igr)%cella(i,j,k)%n
               write(nout)  blochi(ibl,igr)%cella(i,j,k)%n
               write(nout) (blochi(ibl,igr)%cella(i,j,k)%q(n)%b,&
                            blochi(ibl,igr)%cella(i,j,k)%q(n)%i,&
                            blochi(ibl,igr)%cella(i,j,k)%q(n)%j,&
                            blochi(ibl,igr)%cella(i,j,k)%q(n)%k,&
                            blochi(ibl,igr)%cella(i,j,k)%q(n)%w,n=1,m)
            end if
         end do
      end do
   end do
end do

!  chiusura file output
close(nout)

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
! stampa debug per visualizzazione scatole
!===============================================================================
subroutine printbox
use prec
use moddef
use modpar
implicit none
integer (kind=I4P) :: p,i,j,k
type (point) :: q(8),rotate

if (nbox.lt.1) return

open(197,file="scatole.dat",status='unknown',form='formatted')
rewind(197)
write(197,'(a)') 'variables="x" "y" "z"'
do p=1,nbox

   write(197,'(a,i2.2,a)') 'zone T="box #',p,'"'
   if (tbox(p).eq.1) then
      q(1) = obox(p) - ilex(p)*ibox(p) - jlex(p)*jbox(p) - klex(p)*kbox(p)
      q(2) = obox(p) + ilex(p)*ibox(p) - jlex(p)*jbox(p) - klex(p)*kbox(p)
      q(3) = obox(p) - ilex(p)*ibox(p) + jlex(p)*jbox(p) - klex(p)*kbox(p)
      q(4) = obox(p) + ilex(p)*ibox(p) + jlex(p)*jbox(p) - klex(p)*kbox(p)
      q(5) = obox(p) - ilex(p)*ibox(p) - jlex(p)*jbox(p) + klex(p)*kbox(p)
      q(6) = obox(p) + ilex(p)*ibox(p) - jlex(p)*jbox(p) + klex(p)*kbox(p)
      q(7) = obox(p) - ilex(p)*ibox(p) + jlex(p)*jbox(p) + klex(p)*kbox(p)
      q(8) = obox(p) + ilex(p)*ibox(p) + jlex(p)*jbox(p) + klex(p)*kbox(p)
      write(197,'(a)') 'I=2, J=2, K=2'
      write(197,'(a)') 'DATAPACKING=POINT'
      write(197,*) q(1:8)
   else if (tbox(p).eq.2) then
      write(197,'(a)') 'I=2, J=37, K=2'
      write(197,'(a)') 'DATAPACKING=POINT'
      do k=1,2
         q(1) = obox(p) + real(3-2*k,8)*ilex(p)*ibox(p)
         do j=1,37
            q(2) = rotate(real(j-1,8)*10d0,jbox(p),kbox(p))
            do i=1,2
               q(3) = q(1)+real(i-1,8)*jlex(p)*q(2)
               write(197,*) q(3)
            end do
         end do
      end do
   else if (tbox(p).eq.3) then
      write(197,'(a)') 'I=2, J=37, K=2'
      write(197,'(a)') 'DATAPACKING=POINT'
      do k=1,2
         q(1) = obox(p) + real(3-2*k,8)*ilex(p)*ibox(p)
         do j=1,37
            q(2) = rotate(real(j-1,8)*10d0,jbox(p),kbox(p))
            do i=1,2
               q(3) = q(1)+(real(2-i,8)*jlex(p)+real(i-1,8)*klex(p))*q(2)
               write(197,*) q(3)
            end do
         end do
      end do
   else if (tbox(p).eq.9) then
      q(1:8) = vbox(1:8,p)
      write(197,'(a)') 'I=2, J=2, K=2'
      write(197,'(a)') 'DATAPACKING=POINT'
      write(197,*) q(1:8)
   end if

end do
close(197)

return
end

!===============================================================================
! ruota un vettore di un angolo assegnato
function rotate(alfa,v,w)
use prec
use moddef
implicit none
type (point) :: rotate
real (kind=R8P), intent(in)  :: alfa
type (point), intent(in) :: v,w
real (kind=R8P) :: alra
alra = alfa*atan(1d0)/45d0
rotate = cos(alra)*v + sin(alra)*w
end function

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  preprocessore CHIMERA per XSHIP - scrittura matrice BUFDIS
!===============================================================================
subroutine printdis(igr)
use prec
use modpar, only : nout,nbl,blochi
implicit none

character (len=50) :: ftmp
integer (kind=I4P) :: igr,ibl, nii,njj,nkk, i,j,k
real (kind=R8P) :: wrk = 0d0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  apertura file di output
write(ftmp,'(a,i2.2)') 'dist.',igr
open(nout,file=ftmp,form='unformatted',status='unknown',position='rewind')

! scrive 7 reali (alla stregua di Xnavis)
write(nout) (wrk,ibl=1,7)

!  scrittura dimensioni dei blocchi
do ibl=1,nbl
   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk
   write(nout) nii,njj,nkk
end do

!  scrittura delle matrici ICC
do ibl=1,nbl
   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk
   write(nout) (((blochi(ibl,igr)%cella(i,j,k)%dist,&
                  i=0,nii+1),j=0,njj+1),k=0,nkk+1)
end do

!  chiusura file output
close(nout)

end subroutine
