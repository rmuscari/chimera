!===============================================================================
! Determina se un punto è dentro o fuori un blocco contando le intersezioni dei
! raggi che partono dal punto con i confini del blocco stesso.
! Prima conta le intersezioni del raggio z>zc
! Poi conta le intersezioni del raggio y>yc
! Se il numero di intersezioni coincide sono abbastanza sicuro di non aver fatto
! cazzate. Altrimenti conto le intersezioni del raggio x>xc.
! In un mondo perfetto basterebbe contare le intersezioni incontrate da un solo
! raggio, ma in qualunque direzione mi muova c'è sempre la possibilità che la
! proiezione di una faccetta sul piano normale si intrecci mandando tutto a
! puttane
!===============================================================================
subroutine intersec(m,igr,ibl,nrmx,nrmy,nrmz,ing,c)
use prec
use moddef, only : point
use modpar, only : blogrd
implicit none
integer (kind=I4P) :: m,igr,ibl
real (kind=R8P), intent(in), dimension(6) :: nrmx,nrmy,nrmz
type (point), intent(in) :: ing(2,6),c

integer (kind=I4P) :: mx,my,mz

call intersez(m,igr,ibl,nrmz,ing,c)
m = mod(m,2)

return ! non usa le altre direzioni

call intersez(mz,igr,ibl,nrmz,ing,c)
mz = mod(mz,2)

call intersey(my,igr,ibl,nrmy,ing,c)
my = mod(my,2)

if (my.eq.mz) then
   m = mz
   return
end if

call intersex(mx,igr,ibl,nrmx,ing,c)
mx = mod(mx,2)

if (mx.eq.my) then
   m = mx
else if (mx.eq.mz) then
   m = mz
end if

print*,'intersec!',mx,my,mz

end subroutine

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

!===============================================================================
! conta il numero di intersezioni del raggio x=xc,y=yc,z>zc con il confine
! del blocco
!===============================================================================
subroutine intersez(m,igr,ibl,nrm,ing,c)
use prec
use moddef, only : point
use modpar, only : blogrd
implicit none
real (kind=R8P), parameter :: eps = 1d-4
integer (kind=I4P) :: m
integer (kind=I4P), intent(in) :: ibl,igr
real (kind=R8P), intent(in) :: nrm(6)
type (point), intent(in) :: ing(2,6),c
type (point) :: face(2,2)
integer (kind=I4P) :: i,j,k,nf,nii,njj,nkk

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk

m = 0
do nf=1,6
   ! se il punto esce dall'ingombro della faccia, salta a quella successiva
   if (c%x .lt. ing(1,nf)%x .or.&
       c%y .lt. ing(1,nf)%y .or.&
       c%x .gt. ing(2,nf)%x .or.&
       c%y .gt. ing(2,nf)%y .or.&
       c%z .gt. ing(2,nf)%z) cycle
   ! se la componente massima lungo z, fra tutte le faccette della faccia, del
   ! versore normale è minore di eps, vuol dire che la faccia è parallela
   ! all'asse z e la procedura è mal condizionata
   if (nrm(nf).lt.eps) cycle
   if (nf.lt.3) then
      i = (nf-1)*nii
      do k=1,nkk
         do j=1,njj
            face = blogrd(ibl,igr)%nodo(i,j-1:j,k-1:k)
            call inter3z(m,face,c)
         end do
      end do
   else if (nf.gt.4) then
      k = (nf-5)*nkk
      do j=1,njj
         do i=1,nii
            face = blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k)
            call inter3z(m,face,c)
         end do
      end do
   else
      j = (nf-3)*njj
      do k=1,nkk
         do i=1,nii
            face = blogrd(ibl,igr)%nodo(i-1:i,j,k-1:k)
            call inter3z(m,face,c)
         end do
      end do
   end if
end do

end subroutine

!===============================================================================
! controlla se il raggio verso le "z" crescenti interseca la faccia "f"
!===============================================================================
subroutine inter3z(m,f,c)
use prec
use moddef, only : point
implicit none
integer (kind=I4P) :: m,n
type (point) :: f(2,2),c
real (kind=R8P), dimension(5) :: xpi,ypi

! non può esserci intersezione
if (maxval(f(:,:)%x).lt.c%x) return
if (minval(f(:,:)%x).gt.c%x) return
if (maxval(f(:,:)%y).lt.c%y) return
if (minval(f(:,:)%y).gt.c%y) return
if (maxval(f(:,:)%z).lt.c%z) return

! se zmin(faccia) < z(c) < zmax(faccia) deve prendere qualche precauzione
! cerca l'intersezione del raggio con la faccia e controlla se z(inter) > z(c)
if (minval(f(:,:)%z).lt.c%z) then
   call interz(n,f,c)
   if (n.eq.0) return
end if

! controlla che la proiezione sul piano z=0 del punto "c" sia interna alla
! proiezione della faccia "f"
xpi(1) = f(1,1)%x - c%x
ypi(1) = f(1,1)%y - c%y
xpi(2) = f(1,2)%x - c%x
ypi(2) = f(1,2)%y - c%y
xpi(3) = f(2,2)%x - c%x
ypi(3) = f(2,2)%y - c%y
xpi(4) = f(2,1)%x - c%x
ypi(4) = f(2,1)%y - c%y
call inter2d(n,xpi,ypi)
if (n.eq.1)  m = m+1

return
end

subroutine interz(n,f,c)
use prec
use moddef
implicit none
real (kind=R8P), parameter :: eps=1.0d-50
integer (kind=I4P) :: n
type (point) :: f(2,2),c
type (point) :: dia,dib,cen
real (kind=R8P) :: r,s,den

n = 0
dia = f(2,2)-f(1,1)
dib = f(1,2)-f(2,1)

den = dia%y*dib%x - dia%x*dib%y
if (abs(den).lt.eps) return

cen = c-mean(f)
r = (cen%y*dib%x - cen%x*dib%y)/den
s = (cen%x*dia%y - cen%y*dia%x)/den

if ((r*dia%z+s*dib%z).gt.cen%z) n = 1

return
end

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

!===============================================================================
! conta il numero di intersezioni del raggio x=xc,z=zc,y>yc con il confine
! del blocco
!===============================================================================
subroutine intersey(m,igr,ibl,nrm,ing,c)
use prec
use moddef, only : point
use modpar, only : blogrd
implicit none
real (kind=R8P), parameter :: eps = 1d-4
integer (kind=I4P) :: m
integer (kind=I4P), intent(in) :: ibl,igr
real (kind=R8P), intent(in) :: nrm(6)
type (point), intent(in) :: ing(2,6),c
type (point) :: face(2,2)
integer (kind=I4P) :: i,j,k,nf,nii,njj,nkk

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk

m = 0
do nf=1,6
   ! se il punto esce dall'ingombro della faccia, salta a quella successiva
   if (c%x .lt. ing(1,nf)%x .or.&
       c%z .lt. ing(1,nf)%z .or.&
       c%x .gt. ing(2,nf)%x .or.&
       c%y .gt. ing(2,nf)%y .or.&
       c%z .gt. ing(2,nf)%z) cycle
   ! se la componente massima lungo y, fra tutte le faccette della faccia, del
   ! versore normale è minore di eps, vuol dire che la faccia è parallela
   ! all'asse y e la procedura è mal condizionata
   if (nrm(nf).lt.eps) cycle
   if (nf.lt.3) then
      i = (nf-1)*nii
      do k=1,nkk
         do j=1,njj
            face = blogrd(ibl,igr)%nodo(i,j-1:j,k-1:k)
            call inter3y(m,face,c)
         end do
      end do
   else if (nf.gt.4) then
      k = (nf-5)*nkk
      do j=1,njj
         do i=1,nii
            face = blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k)
            call inter3y(m,face,c)
         end do
      end do
   else
      j = (nf-3)*njj
      do k=1,nkk
         do i=1,nii
            face = blogrd(ibl,igr)%nodo(i-1:i,j,k-1:k)
            call inter3y(m,face,c)
         end do
      end do
   end if
end do

end subroutine

!===============================================================================
! controlla se il raggio verso le "y" crescenti interseca la faccia "f"
!===============================================================================
subroutine inter3y(m,f,c)
use prec
use moddef, only : point
implicit none
integer (kind=I4P) :: m,n
type (point) :: f(2,2),c
real (kind=R8P), dimension(5) :: xpi,zpi

! non può esserci intersezione
if (maxval(f(:,:)%x).lt.c%x) return
if (minval(f(:,:)%x).gt.c%x) return
if (maxval(f(:,:)%y).lt.c%y) return
if (maxval(f(:,:)%z).lt.c%z) return
if (minval(f(:,:)%z).gt.c%z) return

! se ymin(faccia) < y(c) < ymax(faccia) deve prendere qualche precauzione
! cerca l'intersezione del raggio con la faccia e controlla se y(inter) > y(c)
if (minval(f(:,:)%y).lt.c%y) then
   call intery(n,f,c)
   if (n.eq.0) return
end if

! controlla che la proiezione sul piano y=0 del punto "c" sia interna alla
! proiezione della faccia "f"
xpi(1) = f(1,1)%x - c%x
zpi(1) = f(1,1)%z - c%z
xpi(2) = f(1,2)%x - c%x
zpi(2) = f(1,2)%z - c%z
xpi(3) = f(2,2)%x - c%x
zpi(3) = f(2,2)%z - c%z
xpi(4) = f(2,1)%x - c%x
zpi(4) = f(2,1)%z - c%z
call inter2d(n,zpi,xpi)
if (n.eq.1)  m = m+1

return
end

subroutine intery(n,f,c)
use prec
use moddef
implicit none
real (kind=R8P), parameter :: eps=1.0d-50
integer (kind=I4P) :: n
type (point) :: f(2,2),c
type (point) :: dia,dib,cen
real (kind=R8P) :: r,s,den

n = 0
dia = f(2,2)-f(1,1)
dib = f(1,2)-f(2,1)

den = dia%x*dib%z - dia%z*dib%x
if (abs(den).lt.eps) return

cen = c-mean(f)
r = (cen%x*dib%z - cen%z*dib%x)/den
s = (cen%z*dia%x - cen%x*dia%z)/den

if ((r*dia%y+s*dib%y).gt.cen%y) n = 1

return
end

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

!===============================================================================
! conta il numero di intersezioni del raggio x=xc,z=zc,y>yc con il confine
! del blocco
!===============================================================================
subroutine intersex(m,igr,ibl,nrm,ing,c)
use prec
use moddef, only : point
use modpar, only : blogrd
implicit none
real (kind=R8P), parameter :: eps = 1d-4
integer (kind=I4P) :: m
integer (kind=I4P), intent(in) :: ibl,igr
real (kind=R8P), intent(in) :: nrm(6)
type (point), intent(in) :: ing(2,6),c
type (point) :: face(2,2)
integer (kind=I4P) :: i,j,k,nf,nii,njj,nkk

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk

m = 0
do nf=1,6
   ! se il punto esce dall'ingombro della faccia, salta a quella successiva
   if (c%y .lt. ing(1,nf)%y .or.&
       c%z .lt. ing(1,nf)%z .or.&
       c%x .gt. ing(2,nf)%x .or.&
       c%y .gt. ing(2,nf)%y .or.&
       c%z .gt. ing(2,nf)%z) cycle
   ! se la componente massima lungo x, fra tutte le faccette della faccia, del
   ! versore normale è minore di eps, vuol dire che la faccia è parallela
   ! all'asse x e la procedura è mal condizionata
   if (nrm(nf).lt.eps) cycle
   if (nf.lt.3) then
      i = (nf-1)*nii
      do k=1,nkk
         do j=1,njj
            face = blogrd(ibl,igr)%nodo(i,j-1:j,k-1:k)
            call inter3x(m,face,c)
         end do
      end do
   else if (nf.gt.4) then
      k = (nf-5)*nkk
      do j=1,njj
         do i=1,nii
            face = blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k)
            call inter3x(m,face,c)
         end do
      end do
   else
      j = (nf-3)*njj
      do k=1,nkk
         do i=1,nii
            face = blogrd(ibl,igr)%nodo(i-1:i,j,k-1:k)
            call inter3x(m,face,c)
         end do
      end do
   end if
end do

end subroutine

!===============================================================================
! controlla se il raggio verso le "x" crescenti interseca la faccia "f"
!===============================================================================
subroutine inter3x(m,f,c)
use prec
use moddef, only : point
implicit none
integer (kind=I4P) :: m,n
type (point) :: f(2,2),c
real (kind=R8P), dimension(5) :: ypi,zpi

! non può esserci intersezione
if (maxval(f(:,:)%x).lt.c%x) return
if (maxval(f(:,:)%y).lt.c%y) return
if (minval(f(:,:)%y).gt.c%y) return
if (maxval(f(:,:)%z).lt.c%z) return
if (minval(f(:,:)%z).gt.c%z) return

! se xmin(faccia) < x(c) < xmax(faccia) deve prendere qualche precauzione
! cerca l'intersezione del raggio con la faccia e controlla se x(inter) > x(c)
if (minval(f(:,:)%x).lt.c%x) then
   call interx(n,f,c)
   if (n.eq.0) return
end if

! controlla che la proiezione sul piano x=0 del punto "c" sia interna alla
! proiezione della faccia "f"
ypi(1) = f(1,1)%y - c%y
zpi(1) = f(1,1)%z - c%z
ypi(2) = f(1,2)%y - c%y
zpi(2) = f(1,2)%z - c%z
ypi(3) = f(2,2)%y - c%y
zpi(3) = f(2,2)%z - c%z
ypi(4) = f(2,1)%y - c%y
zpi(4) = f(2,1)%z - c%z
call inter2d(n,ypi,zpi)
if (n.eq.1)  m = m+1

return
end

subroutine interx(n,f,c)
use prec
use moddef
implicit none
real (kind=R8P), parameter :: eps=1.0d-50
integer (kind=I4P) :: n
type (point) :: f(2,2),c
type (point) :: dia,dib,cen
real (kind=R8P) :: r,s,den

n = 0
dia = f(2,2)-f(1,1)
dib = f(1,2)-f(2,1)

den = dia%z*dib%y - dia%y*dib%z
if (abs(den).lt.eps) return

cen = c-mean(f)
r = (cen%z*dib%y - cen%y*dib%z)/den
s = (cen%y*dia%z - cen%z*dia%y)/den

if ((r*dia%x+s*dib%x).gt.cen%x) n = 1

return
end

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

!===============================================================================
! controlla che la proiezione del punto "c" sia interna alla proiezione della
! faccia "f"
!===============================================================================
subroutine inter2d(n,xpi,ypi)
use prec
implicit none

integer (kind=I4P) :: n
real (kind=R8P)    :: xpi(5),ypi(5),y
integer (kind=I4P) :: i,ip1

!...............................................................................

! numero intersezioni
n = 0

! "c" è già nell'origine, per comodità definisco 5==1
xpi(5) = xpi(1)
ypi(5) = ypi(1)

! verifico se sta sul primo vertice
if ( (xpi(1).eq.0.0d0).and.(ypi(1).eq.0.0d0) ) then
   n = 1
   return
end if

! verifico se sta su uno degli altri vertici
do i=2,4
   if ( (xpi(i).eq.0.0d0).and.(ypi(i).eq.0.0d0) ) return
end do

! verifico se sta sul primo o sul quarto lato
do i=1,4,3
   ip1 = i+1
   if (xpi(i)*xpi(ip1).lt.0.0d0) then
      y = ypi(i) - xpi(i)*(ypi(ip1)-ypi(i))/(xpi(ip1)-xpi(i))
      if (y.eq.0.0d0) then
         n = 1
         return
      end if
   else if (xpi(i).eq.0.0d0 .and. xpi(ip1).eq.0.0d0) then
      if (ypi(i)*ypi(ip1).lt.0.0d0) then
         n = 1
         return
      end if
   end if
end do

! verifico se sta sul secondo o sul terzo lato
do i=2,3
   ip1 = i+1
   if (xpi(i)*xpi(ip1).lt.0.0d0) then
      y = ypi(i) - xpi(i)*(ypi(ip1)-ypi(i))/(xpi(ip1)-xpi(i))
      if (y.eq.0.0d0) return
   else if (xpi(i).eq.0.0d0 .and. xpi(ip1).eq.0.0d0) then
      if (ypi(i)*ypi(ip1).lt.0.0d0) return
   end if
end do

!  ! verifico se è all'interno del quadrilatero
!  ! questo metodo non funziona se la proiezione del quadrangolo sul piano z=0
!  ! si intreccia
!  pvr = -xpi(1)*(ypi(2)-ypi(1)) + ypi(1)*(xpi(2)-xpi(1))
!  do i=2,4
!     ip1 = i+1
!     pvn = -xpi(i)*(ypi(ip1)-ypi(i)) + ypi(i)*(xpi(ip1)-xpi(i))
!     if (pvr*pvn.lt.0d0) return
!  end do
!  n = 1

! verifico se è all'interno del quadrilatero
do i=1,4
   ip1 = i+1
   if (xpi(i).eq.0d0) then
      ! il raggio passa per il primo estremo del segmento: aumento n e continuo
      n = n+1
   else if (xpi(i)*xpi(ip1).lt.0d0) then
      ! il raggio interseca il segmento "aperto" (senza gli estremi):
      ! aumento n e continuo
      y = ypi(i) - xpi(i)*(ypi(ip1)-ypi(i))/(xpi(ip1)-xpi(i))
      if (y.gt.0d0) n = n+1
   end if
end do
n = mod(n,2)

return
end
