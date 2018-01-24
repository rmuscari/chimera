!===============================================================================
!  lettura reticolo
!===============================================================================
subroutine iogrid
use prec
use modpar
implicit none
integer (kind=I4P) :: igr,ibl

!  reticolo senza cornice e senza blocchi perforanti
!  va letto se fgr>1 altrimenti sbaglia il calcolo dei limiti
if (lgrid .or. fgr.gt.1) call readgridnocc

if (lgrid) then

   !  generazione cornici per reticolo fitto
   igr = 1
   call calframe(igr,2)

   !  generazione reticoli radi
   do igr=2,ngr

      !  campionatura reticolo livello precedente
      do ibl=1,nbl
         call collectblock(igr,ibl)
      end do

      !  ridefinizione facce, spigoli e vertici per c.c. di adiacenza
      call calframe(igr,1)

   end do

else

   ! i reticoli con cornici sono già scritti ...
   call readgridcc

end if

!  ricerca dei limiti inferiore e superiore per i blocchi e i patch
call limits

! calcola una volta per tutte le coordinate dei centri cella
call cencel

!  scrittura reticolo con cornici
if (lgrid)  call writegridcc

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  lettura reticolo senza cornici
!===============================================================================
subroutine readgridnocc
use prec
use moddef
use modpar, only : blogrd,nbl,ngrd,fgrd
implicit none
integer (kind=I4P) :: ibl,jbl,i,j,k,nii,njj,nkk

!  apre file reticolo
open(ngrd,file=fgrd,form='unformatted',status='old')
rewind(ngrd)

!  lettura numero e dimensione blocchi (già fatta in readpar)
read(ngrd) ibl
do jbl=1,ibl
   read(ngrd) i,j,k
end do

!  lettura coordinate
do ibl=1,nbl
   nii = blogrd(ibl,1)%ni
   njj = blogrd(ibl,1)%nj
   nkk = blogrd(ibl,1)%nk
   read(ngrd) (((blogrd(ibl,1)%nodo(i,j,k)%x,i=0,nii),j=0,njj),k=0,nkk)
   read(ngrd) (((blogrd(ibl,1)%nodo(i,j,k)%y,i=0,nii),j=0,njj),k=0,nkk)
   read(ngrd) (((blogrd(ibl,1)%nodo(i,j,k)%z,i=0,nii),j=0,njj),k=0,nkk)
end do

!  chiusura file reticolo
close(ngrd)

end subroutine

!===============================================================================
!  lettura reticolo con cornici
!===============================================================================
subroutine readgridcc
use prec
use modpar, only : nout,fgr,ngr,nbl, blogrd, fout
implicit none
integer (kind=I4P) :: igr,ibl
integer (kind=I4P) :: i,j,k, ii0,iin,jj0,jjn,kk0,kkn
character (len=50) :: ftmp

do igr=fgr,ngr

   !  apertura file reticolo
   write(ftmp,'(a,i2.2,a)') trim(fout),igr,'.grd'
   open(nout,file=ftmp,form='unformatted',status='old')
   rewind(nout)

   !  lettura numero di blocchi
   read(nout) ibl

   !  lettura dimensioni dei blocchi
   do ibl=1,nbl
      read(nout) i,j,k
   end do

   !  lettura reticolo
   do ibl=1,nbl
      ii0 =                    - blogrd(ibl,igr)%gc(1)
      iin = blogrd(ibl,igr)%ni + blogrd(ibl,igr)%gc(2)
      jj0 =                    - blogrd(ibl,igr)%gc(3)
      jjn = blogrd(ibl,igr)%nj + blogrd(ibl,igr)%gc(4)
      kk0 =                    - blogrd(ibl,igr)%gc(5)
      kkn = blogrd(ibl,igr)%nk + blogrd(ibl,igr)%gc(6)
      read(nout) (((blogrd(ibl,1)%nodo(i,j,k)%x,i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
      read(nout) (((blogrd(ibl,1)%nodo(i,j,k)%y,i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
      read(nout) (((blogrd(ibl,1)%nodo(i,j,k)%z,i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
   end do

   !  chiusura file reticolo
   close(nout)

end do

end subroutine

!===============================================================================
!  scrittura reticolo con cornici
!===============================================================================
subroutine writegridcc
use prec
use modpar, only : nout,fgr,ngr,nbl, blogrd, fout
implicit none
integer (kind=I4P) :: igr,ibl
integer (kind=I4P) :: i,j,k, ii0,iin,jj0,jjn,kk0,kkn
character (len=50) :: ftmp

do igr=fgr,ngr

   !  apertura file reticolo
   write(ftmp,'(a,i2.2,a)') trim(fout),igr,'.grd'
   open(nout,file=ftmp,form='unformatted',status='unknown')
   rewind(nout)

   !  scrittura numero di blocchi
   write(nout) nbl

   !  scrittura dimensioni dei blocchi
   do ibl=1,nbl
      write(nout) blogrd(ibl,igr)%ni,&
                  blogrd(ibl,igr)%nj,&
                  blogrd(ibl,igr)%nk,&
                  blogrd(ibl,igr)%gc(1:6)
   end do

   !  scrittura reticolo
   do ibl=1,nbl
      ii0 =                    - blogrd(ibl,igr)%gc(1)
      iin = blogrd(ibl,igr)%ni + blogrd(ibl,igr)%gc(2)
      jj0 =                    - blogrd(ibl,igr)%gc(3)
      jjn = blogrd(ibl,igr)%nj + blogrd(ibl,igr)%gc(4)
      kk0 =                    - blogrd(ibl,igr)%gc(5)
      kkn = blogrd(ibl,igr)%nk + blogrd(ibl,igr)%gc(6)
      write(nout) (((blogrd(ibl,igr)%nodo(i,j,k)%x,&
                     i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
      write(nout) (((blogrd(ibl,igr)%nodo(i,j,k)%y,&
                     i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
      write(nout) (((blogrd(ibl,igr)%nodo(i,j,k)%z,&
                     i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
   end do

   !  chiusura file output
   close(nout)

end do

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  genera i reticoli radi a partire da quelli più fitti
!===============================================================================
subroutine collectblock(igr,ibl)
use prec
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl
integer (kind=I4P) :: nii,njj,nkk, i,j,k, ii,jj,kk
nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk
do k=-1,nkk+1
   kk = 2*k
   do j=-1,njj+1
      jj = 2*j
      do i=-1,nii+1
         ii = 2*i
         blogrd(ibl,igr)%nodo(i,j,k) = blogrd(ibl,igr-1)%nodo(ii,jj,kk)
      end do
   end do
end do
end subroutine

!===============================================================================
!  ricerca dei punti estremali dei blocchi
!===============================================================================
subroutine limits
use prec
use moddef
use modpar, only : ngr,nbl,npa,bp,fgr,boulay,nr_sottoblocchi,&
                   i_centro,j_centro,k_centro,ingombro,ingofacc,ingopatch,&
                   blogrd
implicit none
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: mgr,ibl,i,j,k,p,nl(6)
integer (kind=I4P) :: subbl,subni,subnj,subnk,macroi,macroj,macrok
type (point)     :: c,inf,sup,ing,suq

mgr = 2**(ngr-1)

!  limiti dei blocchi SENZA cornice
!
do ibl=1,nbl

   nii = blogrd(ibl,1)%ni
   njj = blogrd(ibl,1)%nj
   nkk = blogrd(ibl,1)%nk

   subni = nii/mgr
   subnj = njj/mgr
   subnk = nkk/mgr
   subbl = 1

   do macrok=1,mgr
   do macroj=1,mgr
   do macroi=1,mgr

      inf =  1.0d20
      sup = -1.0d20

      i_centro(subbl,ibl) = (subni*(2*macroi-1))/2
      j_centro(subbl,ibl) = (subnj*(2*macroj-1))/2
      k_centro(subbl,ibl) = (subnk*(2*macrok-1))/2

      do k=(macrok-1)*subnk,macrok*subnk
      do j=(macroj-1)*subnj,macroj*subnj
      do i=(macroi-1)*subni,macroi*subni

         c = blogrd(ibl,1)%nodo(i,j,k)

         if (c%x.lt.inf%x)  inf%x = c%x
         if (c%y.lt.inf%y)  inf%y = c%y
         if (c%z.lt.inf%z)  inf%z = c%z

         if (c%x.gt.sup%x)  sup%x = c%x
         if (c%y.gt.sup%y)  sup%y = c%y
         if (c%z.gt.sup%z)  sup%z = c%z

      end do
      end do
      end do

      ingombro(1,subbl,ibl) = inf
      ingombro(2,subbl,ibl) = sup
      subbl = subbl+1

   end do
   end do
   end do

   ingombro(1,0,ibl)%x = minval(ingombro(1,1:nr_sottoblocchi,ibl)%x)
   ingombro(1,0,ibl)%y = minval(ingombro(1,1:nr_sottoblocchi,ibl)%y)
   ingombro(1,0,ibl)%z = minval(ingombro(1,1:nr_sottoblocchi,ibl)%z)
   ingombro(2,0,ibl)%x = maxval(ingombro(2,1:nr_sottoblocchi,ibl)%x)
   ingombro(2,0,ibl)%y = maxval(ingombro(2,1:nr_sottoblocchi,ibl)%y)
   ingombro(2,0,ibl)%z = maxval(ingombro(2,1:nr_sottoblocchi,ibl)%z)

   ! ingombro facce
   !
   inf =  1.0d20
   sup = -1.0d20
   ing =  1.0d20
   suq = -1.0d20
   do k=0,nkk
   do j=0,njj
      i=0
      c = blogrd(ibl,1)%nodo(i,j,k)

      if (c%x.lt.inf%x)  inf%x = c%x
      if (c%y.lt.inf%y)  inf%y = c%y
      if (c%z.lt.inf%z)  inf%z = c%z

      if (c%x.gt.sup%x)  sup%x = c%x
      if (c%y.gt.sup%y)  sup%y = c%y
      if (c%z.gt.sup%z)  sup%z = c%z

      i=nii
      c = blogrd(ibl,1)%nodo(i,j,k)

      if (c%x.lt.ing%x)  ing%x = c%x
      if (c%y.lt.ing%y)  ing%y = c%y
      if (c%z.lt.ing%z)  ing%z = c%z

      if (c%x.gt.suq%x)  suq%x = c%x
      if (c%y.gt.suq%y)  suq%y = c%y
      if (c%z.gt.suq%z)  suq%z = c%z
   end do
   end do
   ingofacc(1,1,ibl) = inf
   ingofacc(2,1,ibl) = sup
   ingofacc(1,2,ibl) = ing
   ingofacc(2,2,ibl) = suq

   inf =  1.0d20
   sup = -1.0d20
   ing =  1.0d20
   suq = -1.0d20
   do k=0,nkk
   do i=0,nii
      j=0
      c = blogrd(ibl,1)%nodo(i,j,k)

      if (c%x.lt.inf%x)  inf%x = c%x
      if (c%y.lt.inf%y)  inf%y = c%y
      if (c%z.lt.inf%z)  inf%z = c%z

      if (c%x.gt.sup%x)  sup%x = c%x
      if (c%y.gt.sup%y)  sup%y = c%y
      if (c%z.gt.sup%z)  sup%z = c%z

      j=njj
      c = blogrd(ibl,1)%nodo(i,j,k)

      if (c%x.lt.ing%x)  ing%x = c%x
      if (c%y.lt.ing%y)  ing%y = c%y
      if (c%z.lt.ing%z)  ing%z = c%z

      if (c%x.gt.suq%x)  suq%x = c%x
      if (c%y.gt.suq%y)  suq%y = c%y
      if (c%z.gt.suq%z)  suq%z = c%z
   end do
   end do
   ingofacc(1,3,ibl) = inf
   ingofacc(2,3,ibl) = sup
   ingofacc(1,4,ibl) = ing
   ingofacc(2,4,ibl) = suq

   inf =  1.0d20
   sup = -1.0d20
   ing =  1.0d20
   suq = -1.0d20
   do j=0,njj
   do i=0,nii
      k=0
      c = blogrd(ibl,1)%nodo(i,j,k)

      if (c%x.lt.inf%x)  inf%x = c%x
      if (c%y.lt.inf%y)  inf%y = c%y
      if (c%z.lt.inf%z)  inf%z = c%z

      if (c%x.gt.sup%x)  sup%x = c%x
      if (c%y.gt.sup%y)  sup%y = c%y
      if (c%z.gt.sup%z)  sup%z = c%z

      k=nkk
      c = blogrd(ibl,1)%nodo(i,j,k)

      if (c%x.lt.ing%x)  ing%x = c%x
      if (c%y.lt.ing%y)  ing%y = c%y
      if (c%z.lt.ing%z)  ing%z = c%z

      if (c%x.gt.suq%x)  suq%x = c%x
      if (c%y.gt.suq%y)  suq%y = c%y
      if (c%z.gt.suq%z)  suq%z = c%z
   end do
   end do
   ingofacc(1,5,ibl) = inf
   ingofacc(2,5,ibl) = sup
   ingofacc(1,6,ibl) = ing
   ingofacc(2,6,ibl) = suq

end do

! ingombro patches
!
do p=1,npa
   ibl = bp(1,p)
   nl  = bp(5:10,p)
   inf =  1.0d20
   sup = -1.0d20
   ! ancora non ho calcolato boulaf
   ing = 2d0*sqrt(boulay)
   do k=nl(5),nl(6)
   do j=nl(3),nl(4)
   do i=nl(1),nl(2)
      c = blogrd(ibl,1)%nodo(i,j,k)

      if (c%x.lt.inf%x)  inf%x = c%x
      if (c%y.lt.inf%y)  inf%y = c%y
      if (c%z.lt.inf%z)  inf%z = c%z

      if (c%x.gt.sup%x)  sup%x = c%x
      if (c%y.gt.sup%y)  sup%y = c%y
      if (c%z.gt.sup%z)  sup%z = c%z
   end do
   end do
   end do
   ingopatch(1,p) = inf-ing
   ingopatch(2,p) = sup+ing
end do

end subroutine

!===============================================================================
! calcola le coordinate dei centri cella
!===============================================================================
subroutine cencel
use prec
use moddef, only : operator(+), operator(*)
use modpar, only : fgr,ngr,nbl,blomet,blogrd
implicit none
integer (kind=I4P) :: igr,ibl,i,j,k
integer (kind=I4P) :: ii0,iin,jj0,jjn,kk0,kkn

do igr=fgr,ngr
   do ibl=1,nbl
      ii0 = - blogrd(ibl,igr)%gc(1) + 1
      iin = blogrd(ibl,igr)%ni + blogrd(ibl,igr)%gc(2)
      jj0 = - blogrd(ibl,igr)%gc(3) + 1
      jjn = blogrd(ibl,igr)%nj + blogrd(ibl,igr)%gc(4)
      kk0 = - blogrd(ibl,igr)%gc(5) + 1
      kkn = blogrd(ibl,igr)%nk + blogrd(ibl,igr)%gc(6)
      do k=kk0,kkn
         do j=jj0,jjn
            do i=ii0,iin
               blomet(ibl,igr)%cella(i,j,k)%cen = 0.125d0 * &
                                       ( blogrd(ibl,igr)%nodo(i  ,j  ,k  ) &
                                       + blogrd(ibl,igr)%nodo(i-1,j  ,k  ) &
                                       + blogrd(ibl,igr)%nodo(i  ,j-1,k  ) &
                                       + blogrd(ibl,igr)%nodo(i-1,j-1,k  ) &
                                       + blogrd(ibl,igr)%nodo(i  ,j  ,k-1) &
                                       + blogrd(ibl,igr)%nodo(i-1,j  ,k-1) &
                                       + blogrd(ibl,igr)%nodo(i  ,j-1,k-1) &
                                       + blogrd(ibl,igr)%nodo(i-1,j-1,k-1) )
            end do
         end do
      end do
   end do
end do

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  calcolo metrica
!===============================================================================
subroutine metrix(igr)
use prec
use moddef
use modpar
implicit none
integer (kind=I4P) :: igr,ibl
real (kind=R8P) :: volch

!  stampa di controllo
write(*,'(a,i3,a)') 'inizio metrix  (igr =',igr,')'

!  per il "controllo di qualità" sul calcolo dei volumi ...
volch = 0.0d0

do ibl=1,nbl
   !  calcola le normali alle facce
   call calnorm(igr,ibl)

   ! controlla se le facce dei blocchi sono parallele ad un asse cartesiano
   call parallelismo_facce(igr,ibl)

   !  calcolo volumi
   call calvol(igr,ibl,volch)
end do

!  incongruenza massima nel calcolo dei volumi per il livello attuale
write(ndeb,'(a,i3.3,a)') 'reticolo #',igr,'# {{{1'
write(ndeb,'(a,i3.3,a,$)') 'reticolo #',igr,'#  incongruenza massima = '
write(ndeb,'(f5.1,"%")') 1.0d2*volch

end subroutine

!===============================================================================
!  calcola le normali alle facce della cella
!===============================================================================
subroutine calnorm(igr,ibl)
use prec
use moddef, only : point
use modpar, only : blogrd,blomet
implicit none
integer (kind=I4P), intent(in) :: igr,ibl
integer (kind=I4P) :: i,j,k
integer (kind=I4P) :: nii,njj,nkk,gcc(6)
real    (kind=R8P) :: si,sj,sk,snormale
type (point) :: anormale
type (point), dimension(2,2) :: tmpp,tmpq

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk
gcc = blogrd(ibl,igr)%gc

i = max(nii/2,1)
j = max(njj/2,1)
k = max(nkk/2,1)

tmpp = blogrd(ibl,igr)%nodo(i-1,j-1:j,k-1:k)
tmpq = blogrd(ibl,igr)%nodo(i,j-1:j,k-1:k)
si = snormale(tmpp,tmpq)

tmpp = blogrd(ibl,igr)%nodo(i-1:i,j-1,k-1:k)
tmpq = blogrd(ibl,igr)%nodo(i-1:i,j,k-1:k)
sj = snormale(tmpp,tmpq)

tmpp = blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k-1)
tmpq = blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k)
sk = snormale(tmpp,tmpq)

!-------------------------------------------------------------------------------
!  direzione "I"
!-------------------------------------------------------------------------------

do k=-gcc(5)+1,nkk+gcc(6)
   do j=-gcc(3)+1,njj+gcc(4)
      do i=-gcc(1)+1,nii+gcc(2)
         tmpp = blogrd(ibl,igr)%nodo(i,j-1:j,k-1:k)
         blomet(ibl,igr)%cella(i,j,k)%sni = anormale(tmpp,si)
      end do
   end do
end do

!-------------------------------------------------------------------------------
!  direzione "J"
!-------------------------------------------------------------------------------

do k=-gcc(5)+1,nkk+gcc(6)
   do i=-gcc(1)+1,nii+gcc(2)
      do j=-gcc(3)+1,njj+gcc(4)
         tmpp = blogrd(ibl,igr)%nodo(i-1:i,j,k-1:k)
         blomet(ibl,igr)%cella(i,j,k)%snj = anormale(tmpp,sj)
      end do
   end do
end do

!-------------------------------------------------------------------------------
!  direzione "K"
!-------------------------------------------------------------------------------

do j=-gcc(3)+1,njj+gcc(4)
   do i=-gcc(1)+1,nii+gcc(2)
      do k=-gcc(5)+1,nkk+gcc(6)
         tmpp = blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k)
         blomet(ibl,igr)%cella(i,j,k)%snk = anormale(tmpp,sk)
      end do
   end do
end do

end subroutine

!...............................................................................

function snormale(v,w)
use prec
use moddef
implicit none
real (kind=R8P)                :: snormale
type (point), dimension(2,2) :: v,w
type (point)                 :: r,s,t,u
r = w(2,2) - w(1,1)
s = w(2,1) - w(1,2)
t = mean(w(:,:)) - mean(v(:,:))
u = r.cross.s
snormale = 0.5d0
if ((u.dot.t).lt.0.0d0)  snormale = -0.5d0
end function

function anormale(w,segno)
use prec
use moddef
implicit none
type (point)  :: anormale
type (point)  :: w(2,2),r,s
real (kind=R8P) :: segno
r = w(2,2) - w(1,1)
s = w(2,1) - w(1,2)
anormale = segno*(r.cross.s)
end function

!===============================================================================
!  Per ciascuna faccia calcola il massimo di abs(n{x,y,z})/norma(n)
!  cerca di valutare se un faccia è interamente parallela ad un asse cartesiano
!  perché in questo caso l'algoritmo dei raggi per valutare se un punto è
!  interno o esterno ad un blocco fallisce miseramente.
!===============================================================================
subroutine parallelismo_facce(igr,ibl)
use prec
use moddef, only : point,vers0
use modpar, only : ingofax,ingofay,ingofaz,blogrd,blomet
implicit none
integer (kind=I4P), intent(in) :: igr,ibl
integer (kind=I4P) :: i,j,k,nf
integer (kind=I4P) :: nii,njj,nkk
real (kind=R8P) :: ingx,ingy,ingz
type (point) :: nrm

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk
do nf=1,6
   ingx =  0d0
   ingy =  0d0
   ingz =  0d0
   if (nf.lt.3) then
      i = (nf-1)*nii
      do k=1,nkk
         do j=1,njj
            nrm = vers0(blomet(ibl,igr)%cella(i,j,k)%sni)
            ingx = max(ingx,abs(nrm%x))
            ingy = max(ingy,abs(nrm%y))
            ingz = max(ingz,abs(nrm%z))
         end do
      end do
   else if (nf.gt.4) then
      k = (nf-5)*nkk
      do j=1,njj
         do i=1,nii
            nrm = vers0(blomet(ibl,igr)%cella(i,j,k)%snk)
            ingx = max(ingx,abs(nrm%x))
            ingy = max(ingy,abs(nrm%y))
            ingz = max(ingz,abs(nrm%z))
         end do
      end do
   else
      j = (nf-3)*njj
      do k=1,nkk
         do i=1,nii
            nrm = vers0(blomet(ibl,igr)%cella(i,j,k)%snj)
            ingx = max(ingx,abs(nrm%x))
            ingy = max(ingy,abs(nrm%y))
            ingz = max(ingz,abs(nrm%z))
         end do
      end do
   end if
   ingofax(nf,ibl) = ingx
   ingofay(nf,ibl) = ingy
   ingofaz(nf,ibl) = ingz
end do

end subroutine

!===============================================================================
!  calcola volumi
!===============================================================================
subroutine calvol(igr,ibl,volc)
use prec
use moddef
use modpar, only : blomet,blogrd,blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl
character (len=80) :: str
integer (kind=I4P) :: i,j,k,l
integer (kind=I4P) :: nii,njj,nkk
real    (kind=R8P) :: volx,voly,volz,volume
real    (kind=R8P) :: vdif,volc
type (point)     :: c(6),n(6), volm

!...............................................................................

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk

!  inizio del loop sulle celle
do k=1,nkk
   do j=1,njj
      do i=1,nii

         c(1) = 0.25d0 * ( blogrd(ibl,igr)%nodo(i-1,j-1,k-1) &
                         + blogrd(ibl,igr)%nodo(i-1,j  ,k-1) &
                         + blogrd(ibl,igr)%nodo(i-1,j-1,k  ) &
                         + blogrd(ibl,igr)%nodo(i-1,j  ,k  ) )
         c(2) = 0.25d0 * ( blogrd(ibl,igr)%nodo(i  ,j-1,k-1) &
                         + blogrd(ibl,igr)%nodo(i  ,j  ,k-1) &
                         + blogrd(ibl,igr)%nodo(i  ,j-1,k  ) &
                         + blogrd(ibl,igr)%nodo(i  ,j  ,k  ) )
         c(3) = 0.25d0 * ( blogrd(ibl,igr)%nodo(i-1,j-1,k-1) &
                         + blogrd(ibl,igr)%nodo(i  ,j-1,k-1) &
                         + blogrd(ibl,igr)%nodo(i-1,j-1,k  ) &
                         + blogrd(ibl,igr)%nodo(i  ,j-1,k  ) )
         c(4) = 0.25d0 * ( blogrd(ibl,igr)%nodo(i-1,j  ,k-1) &
                         + blogrd(ibl,igr)%nodo(i  ,j  ,k-1) &
                         + blogrd(ibl,igr)%nodo(i-1,j  ,k  ) &
                         + blogrd(ibl,igr)%nodo(i  ,j  ,k  ) )
         c(5) = 0.25d0 * ( blogrd(ibl,igr)%nodo(i-1,j-1,k-1) &
                         + blogrd(ibl,igr)%nodo(i  ,j-1,k-1) &
                         + blogrd(ibl,igr)%nodo(i-1,j  ,k-1) &
                         + blogrd(ibl,igr)%nodo(i  ,j  ,k-1) )
         c(6) = 0.25d0 * ( blogrd(ibl,igr)%nodo(i-1,j-1,k  ) &
                         + blogrd(ibl,igr)%nodo(i  ,j-1,k  ) &
                         + blogrd(ibl,igr)%nodo(i-1,j  ,k  ) &
                         + blogrd(ibl,igr)%nodo(i  ,j  ,k  ) )

         n(1) = -blomet(ibl,igr)%cella(i-1,j,k)%sni
         n(2) =  blomet(ibl,igr)%cella(i,j,k)%sni
         n(3) = -blomet(ibl,igr)%cella(i,j-1,k)%snj
         n(4) =  blomet(ibl,igr)%cella(i,j,k)%snj
         n(5) = -blomet(ibl,igr)%cella(i,j,k-1)%snk
         n(6) =  blomet(ibl,igr)%cella(i,j,k)%snk

         volm = 0.0d0
         do l=1,6
            volm = volm + c(l)*n(l)
         end do

         volx = volm%x
         voly = volm%y
         volz = volm%z

         volume = (volx+voly+volz)/3.0d0

         if (volm.le.0.0d0) then
            str = 'Blocco       ' // trim(castr(ibl))
            write(*,'(/,a)') trim(str)
            str = 'ni,nj,nk     ' // trim(castr(nii)) // '  ' &
                                  // trim(castr(njj)) // '  ' &
                                  // trim(castr(nkk))
            write(*,'(a)') trim(str)
            str = ' i, j, k     ' // trim(castr(i)) // '  ' &
                                  // trim(castr(j)) // '  ' &
                                  // trim(castr(k))
            write(*,'(a)') trim(str)
            write(*,'(a,e11.4)')      '     volx = ',volx
            write(*,'(a,e11.4)')      '     voly = ',voly
            write(*,'(a,e11.4)')      '     volz = ',volz

            do l=1,6
               write(*,'(3(e15.7))') c(l)
            end do
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i-1,j-1,k-1)
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i  ,j-1,k-1)
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i-1,j  ,k-1)
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i  ,j  ,k-1)
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i-1,j-1,k  )
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i  ,j-1,k  )
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i-1,j  ,k  )
            write(*,'(3(e16.9))') blogrd(ibl,igr)%nodo(i  ,j  ,k  )
            stop
         else
            vdif = max(abs(volx-voly),abs(volx-volz),abs(voly-volz))
            volc = max(volc,vdif/volume)
         end if

         blochi(ibl,igr)%cella(i,j,k)%size = volume

      end do
   end do
end do

end subroutine

!===============================================================================
!  calcola le cornici dei blocchi
!===============================================================================
subroutine calframe(igr,c)
use prec
use moddef
use modpar
implicit none
integer (kind=I4P), intent(in) :: igr,c
integer (kind=I4P) :: ibl,tgr,nf,fcc,abl,p,q
integer (kind=I4P) :: nlp(6),nlq(6)

!...............................................................................
!  per tener conto delle dimensioni dei patch nel reticolo igr
tgr = 2**(igr-1)

!...............................................................................
!  TOPPE SULLE FACCE - determina le celle fittizie ...

do p=1,npa

   ibl = bp(1,p)
   nf  = bp(2,p)
   fcc = bp(3,p)
   nlp = bp(5:10,p)/tgr

   ! c.c. naturali
   !
   if (fcc.le.nccnat) then
      call estrface(igr,ibl,nf,nlp,c)

   ! c.c. adiacenza uno a uno
   !
   else if (fcc.gt.100) then
      q   = bp(4,p)
      abl = bp(1,q)
      nlq = bp(5:10,q)/tgr

      call adiaface(igr,ibl,abl,nf,fcc,nlp,nlq)

   ! c.c. chimera
   !
   else
      call estrface(igr,ibl,nf,nlp,c)

   end if
end do

!...............................................................................
!  SPIGOLI

! calcola gli spigoli fittizi come medie dei nodi adiacenti
do ibl=1,nbl
   call meanspig(igr,ibl)
end do

!...............................................................................
!  VERTICI

!  calcola i vertici fittizi come medie dei nodi adiacenti
do ibl=1,nbl
   call meanvert(igr,ibl)
end do

end subroutine

!===============================================================================
!  calcola le posizioni dei punti fittizi esterni al blocco tramite
!  assegnazione diretta dal blocco adiacente
!===============================================================================
subroutine adiaface(igr,ibl,abl,nf,fcc,nl,nla)
use prec
use moddef
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,abl
integer (kind=I4P) :: nf,fcc
integer (kind=I4P) :: nl(6),nla(6)
integer (kind=I4P) :: nii,njj,nkk,gcc(6)
integer (kind=I4P) :: li,lj,lk, m, l0,l1,ln
integer (kind=I4P) :: i,j,k, ia,ja,ka, ib,jb,kb

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk
gcc = blogrd(ibl,igr)%gc

li =     fcc      / 100
lj = mod(fcc,100) / 10
lk = mod(fcc,10)
m  = mod(nf,2)

!-------------------------------------------------------------------------------
if (nf.lt.3) then      ! direzione "I"
!-------------------------------------------------------------------------------
   l0 = m*(-gcc(1))   +   (1-m)*(nii+1)
   l1 = m*( -1 )   +   (1-m)*(nii+gcc(2))
   ln =                (1-m)*nii
   do k=nl(5),nl(6)
      call tangind(lk,k-nl(5),ia,ja,ka,nla)
      do j=nl(3),nl(4)
         call tangind(lj,j-nl(3),ia,ja,ka,nla)
         do i=l0,l1
            call normind(li,m,i-ln,ia,ja,ka,ib,jb,kb,nla)
            blogrd(ibl,igr)%nodo(i,j,k) = blogrd(ibl,igr)%nodo(ln,j,k) + &
                                        ( blogrd(abl,igr)%nodo(ia,ja,ka) &
                                        - blogrd(abl,igr)%nodo(ib,jb,kb) )
         end do
      end do
   end do

!-------------------------------------------------------------------------------
else if (nf.gt.4) then ! direzione "K"
!-------------------------------------------------------------------------------
   l0 = m*(-gcc(5))   +   (1-m)*(nkk+1)
   l1 = m*( -1 )   +   (1-m)*(nkk+gcc(6))
   ln =                (1-m)*nkk
   do j=nl(3),nl(4)
      call tangind(lj,j-nl(3),ia,ja,ka,nla)
      do i=nl(1),nl(2)
         call tangind(li,i-nl(1),ia,ja,ka,nla)
         do k=l0,l1
            call normind(lk,m,k-ln,ia,ja,ka,ib,jb,kb,nla)
            blogrd(ibl,igr)%nodo(i,j,k) = blogrd(ibl,igr)%nodo(i,j,ln) + &
                                        ( blogrd(abl,igr)%nodo(ia,ja,ka) &
                                        - blogrd(abl,igr)%nodo(ib,jb,kb) )
         end do
      end do
   end do

!-------------------------------------------------------------------------------
else                   ! direzione "J"
!-------------------------------------------------------------------------------
   l0 = m*(-gcc(3))   +   (1-m)*(njj+1)
   l1 = m*( -1 )   +   (1-m)*(njj+gcc(4))
   ln =                (1-m)*njj
   do k=nl(5),nl(6)
      call tangind(lk,k-nl(5),ia,ja,ka,nla)
      do i=nl(1),nl(2)
         call tangind(li,i-nl(1),ia,ja,ka,nla)
         do j=l0,l1
            call normind(lj,m,j-ln,ia,ja,ka,ib,jb,kb,nla)
            blogrd(ibl,igr)%nodo(i,j,k) = blogrd(ibl,igr)%nodo(i,ln,k) + &
                                        ( blogrd(abl,igr)%nodo(ia,ja,ka) &
                                        - blogrd(abl,igr)%nodo(ib,jb,kb) )
         end do
      end do
   end do

end if

end subroutine

!===============================================================================
subroutine normind(l,m,r,ia,ja,ka,ib,jb,kb,nla)
use prec
implicit none
integer (kind=I4P), intent(in)  :: l,m,r,nla(6)
integer (kind=I4P), intent(out) :: ia,ja,ka, ib,jb,kb

ib = ia
jb = ja
kb = ka

if      (l.eq.1) then
   ib = m*nla(2)
   ia = ib + r
else if (l.eq.2) then
   ib = (1-m)*nla(2)
   ia = ib - r
else if (l.eq.3) then
   jb = m*nla(4)
   ja = jb + r
else if (l.eq.4) then
   jb = (1-m)*nla(4)
   ja = jb - r
else if (l.eq.5) then
   kb = m*nla(6)
   ka = kb + r
else if (l.eq.6) then
   kb = (1-m)*nla(6)
   ka = kb - r
end if
end subroutine

!===============================================================================
subroutine tangind(l,r,ia,ja,ka,nla)
use prec
implicit none
integer (kind=I4P) :: l,r,ia,ja,ka,nla(6)
if      (l.eq.1) then
   ia = nla(1)+r
else if (l.eq.2) then
   ia = nla(2)-r
else if (l.eq.3) then
   ja = nla(3)+r
else if (l.eq.4) then
   ja = nla(4)-r
else if (l.eq.5) then
   ka = nla(5)+r
else if (l.eq.6) then
   ka = nla(6)-r
end if
end subroutine

!===============================================================================
!  calcola le posizioni dei punti fittizi esterni al blocco
!  estrapolando lungo le linee cordinate (c'è anche la possibilità, commentata,
!  di estrapolare lungo le normali alle facce)
!===============================================================================
subroutine estrface(igr,ibl,nf,nl,c)
use prec
use moddef
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl
integer (kind=I4P) :: nf,nl(6),c
integer (kind=I4P) :: i,j,k, m,s,ss
integer (kind=I4P) :: nii,njj,nkk
type (point), allocatable :: win(:,:)

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk

!  allocazione memoria di lavoro
m = max(nii,njj,nkk)
allocate(win(0:4,(m+1)*(m+1)))

s = 2*mod(nf,2)-1
ss = 2*s

!-------------------------------------------------------------------------------
if (nf.lt.3) then      ! direzione "I"
!-------------------------------------------------------------------------------
   m = 0
   if (nii.eq.1) ss = s
   do k=nl(5),nl(6)
      do j=nl(3),nl(4)
         m = m+1
         win(0,m) = blogrd(ibl,igr)%nodo(nl(1)+ss,j,k)
         win(1,m) = blogrd(ibl,igr)%nodo(nl(1)+s ,j,k)
         win(2,m) = blogrd(ibl,igr)%nodo(nl(1)   ,j,k)
      end do
   end do

   ! Per Xnavis va usato "estrlin" !
   ! call estrort(nl(4)-nl(3),nl(6)-nl(5),win)
   call estrlin(nl(4)-nl(3),nl(6)-nl(5),win)

   m = 0
   if (c.eq.1) then
      ! sui reticoli igr>1 estrapolo solo i nodi sulla seconda cornice
      ! i nodi sulla prima cornice sono "collezionati" dal reticolo fitto
      do k=nl(5),nl(6)
         do j=nl(3),nl(4)
            m = m+1
            blogrd(ibl,igr)%nodo(nl(1)-2*s,j,k) = win(4,m)
         end do
      end do
   else
      ! sul reticolo fitto devo estrapolare tutte le cornici
      do k=nl(5),nl(6)
         do j=nl(3),nl(4)
            m = m+1
            blogrd(ibl,igr)%nodo(nl(1)-s  ,j,k) = win(3,m)
            blogrd(ibl,igr)%nodo(nl(1)-2*s,j,k) = win(4,m)
         end do
      end do
   end if

!-------------------------------------------------------------------------------
else if (nf.gt.4) then ! direzione "K"
!-------------------------------------------------------------------------------
   m = 0
   if (nkk.eq.1) ss = s
   do j=nl(3),nl(4)
      do i=nl(1),nl(2)
         m = m+1
         win(0,m) = blogrd(ibl,igr)%nodo(i,j,nl(5)+ss)
         win(1,m) = blogrd(ibl,igr)%nodo(i,j,nl(5)+s )
         win(2,m) = blogrd(ibl,igr)%nodo(i,j,nl(5)   )
      end do
   end do

   ! Per Xnavis va usato "estrlin" !
   ! call estrort(nl(2)-nl(1),nl(4)-nl(3),win)
   call estrlin(nl(2)-nl(1),nl(4)-nl(3),win)

   m = 0
   if (c.eq.1) then
      ! sui reticoli igr>1 estrapolo solo i nodi sulla seconda cornice
      ! i nodi sulla prima cornice sono "collezionati" dal reticolo fitto
      do j=nl(3),nl(4)
         do i=nl(1),nl(2)
            m = m+1
            blogrd(ibl,igr)%nodo(i,j,nl(5)-2*s) = win(4,m)
         end do
      end do
   else
      ! sul reticolo fitto devo estrapolare tutte le cornici
      do j=nl(3),nl(4)
         do i=nl(1),nl(2)
            m = m+1
            blogrd(ibl,igr)%nodo(i,j,nl(5)-s  ) = win(3,m)
            blogrd(ibl,igr)%nodo(i,j,nl(5)-2*s) = win(4,m)
         end do
      end do
   end if

!-------------------------------------------------------------------------------
else                   ! direzione "J"
!-------------------------------------------------------------------------------
   m = 0
   if (njj.eq.1) ss = s
   do k=nl(5),nl(6)
      do i=nl(1),nl(2)
         m = m+1
         win(0,m) = blogrd(ibl,igr)%nodo(i,nl(3)+ss,k)
         win(1,m) = blogrd(ibl,igr)%nodo(i,nl(3)+s ,k)
         win(2,m) = blogrd(ibl,igr)%nodo(i,nl(3)   ,k)
      end do
   end do

   ! Per Xnavis va usato "estrlin" !
   ! call estrort(nl(2)-nl(1),nl(6)-nl(5),win)
   call estrlin(nl(2)-nl(1),nl(6)-nl(5),win)

   m = 0
   if (c.eq.1) then
      ! sui reticoli igr>1 estrapolo solo i nodi sulla seconda cornice
      ! i nodi sulla prima cornice sono "collezionati" dal reticolo fitto
      do k=nl(5),nl(6)
         do i=nl(1),nl(2)
            m = m+1
            blogrd(ibl,igr)%nodo(i,nl(3)-2*s,k) = win(4,m)
         end do
      end do
   else
      ! sul reticolo fitto devo estrapolare tutte le cornici
      do k=nl(5),nl(6)
         do i=nl(1),nl(2)
            m = m+1
            blogrd(ibl,igr)%nodo(i,nl(3)-s  ,k) = win(3,m)
            blogrd(ibl,igr)%nodo(i,nl(3)-2*s,k) = win(4,m)
         end do
      end do
   end if

!-------------------------------------------------------------------------------
end if
!-------------------------------------------------------------------------------

!  libera memoria di lavoro
deallocate(win)

end subroutine

!===============================================================================
subroutine estrlin(ni,nj,w)
use prec
use moddef
implicit none
integer (kind=I4P) :: ni,nj, i,j
real (kind=R8P)             :: r,d
type (point) :: w(0:4,0:ni,0:nj),s
do j=0,nj
   do i=0,ni
      d = dist(w(1,i,j),w(0,i,j))
      r = 1.0d0
      if (d.gt.0.0d0) r = dist(w(2,i,j),w(1,i,j))/d
      r = sqrt(r)
      s = r*(w(2,i,j)-w(1,i,j))
      w(3,i,j) = w(2,i,j) + s
      w(4,i,j) = w(3,i,j) + s
   end do
end do
end subroutine

!===============================================================================
subroutine estrort(ni,nj,w)
use prec
use moddef
implicit none
integer (kind=I4P) :: ni,nj, i,j
type (point)     :: w(0:4,0:ni,0:nj),r,s
real (kind=R8P)    :: t

!  calcola le normali nei centri faccia
do j=1,nj
   do i=1,ni
      r = (w(2,i,j)-w(2,i-1,j-1)).cross.(w(2,i-1,j)-w(2,i,j-1))
      s = mean(w(2,i-1:i,j-1:j)) - mean(w(1,i-1:i,j-1:j))
      t = r.dot.r
      ! se uno spigolo è degenere, |r| è nullo :-)
      if (t.gt.1.0d-20) then
         w(3,i,j) = r * ((r.dot.s)/t)
      else
         w(3,i,j) = s * (1.0d-6/(s.dot.s))
      end if
   end do
end do

!  le trasferisce sui vertici delle celle
do j=1,nj-1
   do i=1,ni-1
      w(4,i,j) = mean(w(3,i:i+1,j:j+1))
   end do
end do

do j=1,nj-1
   w(4,0 ,j) = 0.5d0*(w(3,1 ,j)+w(3,1 ,j+1))
   w(4,ni,j) = 0.5d0*(w(3,ni,j)+w(3,ni,j+1))
end do

do i=1,ni-1
   w(4,i,0 ) = 0.5d0*(w(3,i,1 )+w(3,i+1,1 ))
   w(4,i,nj) = 0.5d0*(w(3,i,nj)+w(3,i+1,nj))
end do

w(4,0 ,0 ) = w(3,1 ,1 )
w(4,ni,0 ) = w(3,ni,1 )
w(4,0 ,nj) = w(3,1 ,nj)
w(4,ni,nj) = w(3,ni,nj)

!  calcola i punti estrapolati
do j=0,nj
   do i=0,ni
      r = w(4,i,j)
      w(3,i,j) = w(2,i,j) + r
      w(4,i,j) = w(3,i,j) + r
   end do
end do

end subroutine

!===============================================================================
!  calcola i punti negli spigoli fittizi mediando i nodi adiacenti
!===============================================================================
subroutine meanspig(igr,ibl)
use prec
use moddef
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl
integer (kind=I4P) :: nii,njj,nkk,gc(6)
integer (kind=I4P) :: i,j,k,m
integer (kind=I4P) :: i0,i1,is
integer (kind=I4P) :: j0,j1,js
integer (kind=I4P) :: k0,k1,ks
type (point), allocatable, dimension(:,:) :: tmpt

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk
gc  = blogrd(ibl,igr)%gc

!-------------------------------------------------------------------------------
!   direzione "I"
!-------------------------------------------------------------------------------

do j=3,4
   m = mod(j,2)
   js = 1-2*m ; j0 = njj*(1-m) ; j1 = j0+js*gc(j)
   do k=5,6
      m = mod(k,2)
      ks = 1-2*m ; k0 = nkk*(1-m) ; k1 = k0+ks*gc(k)
      allocate(tmpt(gc(j)+1,gc(k)+1))
      do i=0,nii
         tmpt = blogrd(ibl,igr)%nodo(i,j0:j1:js,k0:k1:ks)
         call meanspig3(tmpt,gc(j),gc(k))
         blogrd(ibl,igr)%nodo(i,j0+js:j1:js,k0+ks:k1:ks) = tmpt(2:,2:)
      end do
      deallocate(tmpt)
   end do
end do

!-------------------------------------------------------------------------------
!   direzione "J"
!-------------------------------------------------------------------------------

do k=5,6
   m = mod(k,2)
   ks = 1-2*m ; k0 = nkk*(1-m) ; k1 = k0+ks*gc(k)
   do i=1,2
      m = mod(i,2)
      is = 1-2*m ; i0 = nii*(1-m) ; i1 = i0+is*gc(i)
      allocate(tmpt(gc(i)+1,gc(k)+1))
      do j=0,njj
         tmpt = blogrd(ibl,igr)%nodo(i0:i1:is,j,k0:k1:ks)
         call meanspig3(tmpt,gc(i),gc(k))
         blogrd(ibl,igr)%nodo(i0+is:i1:is,j,k0+ks:k1:ks) = tmpt(2:,2:)
      end do
      deallocate(tmpt)
   end do
end do

!-------------------------------------------------------------------------------
!   direzione "K"
!-------------------------------------------------------------------------------

do i=1,2
   m = mod(i,2)
   is = 1-2*m ; i0 = nii*(1-m) ; i1 = i0+is*gc(i)
   do j=3,4
      m = mod(j,2)
      js = 1-2*m ; j0 = njj*(1-m) ; j1 = j0+js*gc(j)
      allocate(tmpt(gc(i)+1,gc(j)+1))
      do k=0,nkk
         tmpt = blogrd(ibl,igr)%nodo(i0:i1:is,j0:j1:js,k)
         call meanspig3(tmpt,gc(i),gc(j))
         blogrd(ibl,igr)%nodo(i0+is:i1:is,j0+js:j1:js,k) = tmpt(2:,2:)
      end do
      deallocate(tmpt)
   end do
end do

end subroutine

!...............................................................................
subroutine meanspig3(pt,ni,nj)
use prec
use moddef
implicit none
integer (kind=I4P) :: ni,nj
type (point) :: pt(0:ni,0:nj)
integer (kind=I4P) :: i,j
type (point) :: da,db,dav,dbv
real (kind=R8P) :: nrm

da  = pt(ni,0) - pt(0,0)
nrm = norma(da)
dav = 0d0
if (nrm.gt.1.0d-10) then
   dav = da / nrm
end if

db = pt(0,nj) - pt(0,0)
nrm = norma(db)
dbv = 0d0
if (nrm.gt.1.0d-10) then
   dbv = db / nrm
end if

! if ((dav.dot.dbv).lt.0.05d0) then
   do j=1,nj
      do i=1,ni
         pt(i,j) = pt(0,0) + (pt(i,0)-pt(0,0)) + (pt(0,j)-pt(0,0))
      end do
   end do
!  else
!     pt(1,1) = 0.5d0*(pt(1,0)+pt(0,1))
!     pt(2,2) = 0.5d0*(pt(2,0)+pt(0,2))
!     pt(1,2) = 0.5d0*(pt(2,2)+pt(0,2))
!     pt(2,1) = 0.5d0*(pt(2,2)+pt(2,0))
!  end if

end subroutine

!===============================================================================
!  calcola i punti nei vertici fittizi mediando i nodi adiacenti
!===============================================================================
subroutine meanvert(igr,ibl)
use prec
use moddef
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl
real (kind=R8P), parameter :: r = 1.0_R8P
integer (kind=I4P) :: nii,njj,nkk,gc(6)
integer (kind=I4P) :: i,fi,is,i0,i1,ii
integer (kind=I4P) :: j,fj,js,j0,j1,jj
integer (kind=I4P) :: k,fk,ks,k0,k1,kk

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk
gc  = blogrd(ibl,igr)%gc

! loop sugli otto vertici del blocco senza cornici
do k=0,1
   fk = 5+k ; ks = 2*k-1 ; k0 = nkk*k ; k1 = k0+ks*gc(fk)
   do j=0,1
      fj = 3+j ; js = 2*j-1 ; j0 = njj*j ; j1 = j0+js*gc(fj)
      do i=0,1
         fi = 1+i ; is = 2*i-1 ; i0 = nii*i ; i1 = i0+is*gc(fi)

         ! loop sulle ghost cells del vertice considerato
         do kk=k0+ks,k1,ks
            do jj=j0+js,j1,js
               do ii=i0+is,i1,is

                  blogrd(ibl,igr)%nodo(ii,jj,kk) = &
              blogrd(ibl,igr)%nodo(i0,j0,k0)  &
         + r*(blogrd(ibl,igr)%nodo(ii,j0,k0) - blogrd(ibl,igr)%nodo(i0,j0,k0)) &
         + r*(blogrd(ibl,igr)%nodo(i0,jj,k0) - blogrd(ibl,igr)%nodo(i0,j0,k0)) &
         + r*(blogrd(ibl,igr)%nodo(i0,j0,kk) - blogrd(ibl,igr)%nodo(i0,j0,k0))

               end do
            end do
         end do

      end do
   end do
end do

end subroutine
