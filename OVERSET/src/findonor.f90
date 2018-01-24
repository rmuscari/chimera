!===============================================================================
!  ricerca della cella donatrice
!
!  gli elementi del vettore "flag" servono ad indicare
!  in INGRESSO
!      flag(1) --> ricerca nel gruppo "1" di blocchi
!      ....
!      flag(6) --> ricerca nel gruppo "6" di blocchi
!      flag(7) --> centro-cella o centro-faccia
!      flag(8) --> accetta qualsiasi donatore (anche + grande)
!      flag(9) --> ricerca del punto + vicino (se non trova donatori)
!
!  in USCITA
!     flag(10) = 0 non ha trovato donatori
!                1 non ha trovato donatori; ha cercato la cella + vicina
!                2 ha trovato un donatore (non importa se + piccolo o + grande)
!                3 ha trovato un donatore a priorita' piu' alta (e quindi non
!                  importa se + piccolo o + grande)
!                4 ha trovato un donatore + piccolo della cella chimera
!
!===============================================================================
subroutine findonor(igr,ibl,cellsize,flag,c,md,lvbad)
use prec
use moddef, only : donor,dist
use modpar
implicit none

integer (kind=I4P), intent(in) :: igr,ibl
real    (kind=R8P), intent(in) :: cellsize
type (point), intent(in) :: c
integer (kind=I4P) :: flag(10)
type (donor)     :: md(donmax)
integer (kind=I1P) :: lvbad(0:lvmax)
integer (kind=I4P) :: cou,m,n, bd,id,jd,kd
integer (kind=I4P) :: don(4,200)
real    (kind=R8P) :: sl,sd
logical (kind=I1P) :: t

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! cerca tutti i possibili donatori dando la precedenza ai blocchi dello stesso
! livello e con priorità più alta
! i blocchi a priorità più bassa vengono considerati solo per celle di cornice
do m=1,6
   cou = 0
   if (flag(m).eq.0) go to 17
   do n=nlv(m-1,ibl)+1,nlv(m,ibl)
      bd = blv(n,ibl)
      t = .false.
      call findfast(igr,bd,id,jd,kd,t,c)
      if (t) then
         cou = cou+1
         don(:,cou) = (/ bd,id,jd,kd /)
         lvbad(lv(bd)) = -1_I1P
      end if
   end do
   if (cou.gt.0) then
      call picksmaller(igr,cou,don)
      if (flag(8).eq.1) then
         ! ha trovato una donatrice, mi trovo su una cornice e non mi
         ! interessano né le dimensioni né la priorità
         flag(10) = 2
         go to 19
      else if (pri(don(1,1)).gt.pri(ibl)) then
         ! la priorità della cella donatrice è più alta; non mi interessano le
         ! dimensioni
         flag(10) = 3
         go to 19
      else
         ! verifico se la cella donatrice è più piccola della cella chimera
         sd = blochi(don(1,1),igr)%cella(don(2,1),don(3,1),don(4,1))%size
         sl = cellsize
         if (sd.lt.sl) then
            ! lo è; posso accettarla
            flag(10) = 4
            go to 19
         end if
      end if
   end if
17 continue
end do

! non ha trovato nessuna donatrice ...
if (flag(9).eq.0) then
   ! ... pazienza
   return
else
   ! ... deve cercare la più vicina
   flag(10) = 1
   sd = 1.0d30
   do bd=1,nbl
      if (bd.ne.ibl) then
         sl = 1.1d0*sd
         call nearest_on_boundary(igr,bd,id,jd,kd,sl,c)
         if (sl.lt.sd) then
            sd = sl
            don(:,1) = (/ bd,id,jd,kd /)
         end if
      end if
   end do
end if

! ha trovato una o più donatrici (e ha già scelto la più piccola)
19 continue
md(1) = donor(don(1,1),don(2,1),don(3,1),don(4,1),1.0)

! cerca il vertice più vicino della donatrice principale (attorno a cui cercare
! gli altri donatori)
call nearest_vertex(igr,md,c)

!  cerca le altre donatrici attorno alla principale
call multidonor(igr,md,c)

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  Cerca nel blocco dbl una possibile cella donatrice id,jd,kd per il punto c
!===============================================================================
subroutine findfast(igr,dbl,id,jd,kd,t,c)
use prec
use modpar
implicit none

integer (kind=I4P) :: igr,dbl,id,jd,kd
logical (kind=I1P) :: t,fuori_ingombro
type (point)     :: c
integer (kind=I4P) :: il,jl,kl,n,m,mc,step

!...............................................................................

!  se "c" esce dai limiti del blocco torna subito alla subroutine chiamante
if (fuori_ingombro(c,ingombro(1,0,dbl))) return

! ripeto il controllo considerando i sotto-blocchi in cui ho diviso dbl
do n=1,nr_sottoblocchi
   if (.not.fuori_ingombro(c,ingombro(1,n,dbl))) go to 29
end do
! "c" esce dai limiti di TUTTI i sotto-blocchi di dbl
return

! ha trovato un sotto-blocco nei cui limiti casca "c"
29 continue

!  conta numero di intersezioni raggio x=xc,y=yc,z>zc  con i confini del blocco
call intersec(m,igr,dbl,ingofax(1,dbl),ingofay(1,dbl),ingofaz(1,dbl),&
              ingofacc(1,1,dbl),c)
if (m.eq.0) return

! cerca la cella più vicina partendo dal centro del primo sotto-blocco nel cui
! ingombro cade "c"
m = 2**(igr-1)
il = i_centro(n,dbl)/m
jl = j_centro(n,dbl)/m
kl = k_centro(n,dbl)/m
do n=ngr-igr+1,1,-1
   step = 2**(n-1)
   call ricerca_per_linee_coordinate(igr,dbl,step,il,jl,kl,c)
end do

!  cerca la cella donatrice tra quelle attorno al vertice + vicino
!  questo loop fa cacare, andrebbe cambiato
mc = max(il,blogrd(dbl,igr)%ni-il,&
         jl,blogrd(dbl,igr)%nj-jl,&
         kl,blogrd(dbl,igr)%nk-kl)
do m=1,mc
   call findonor3(igr,dbl,id,jd,kd,il,jl,kl,m,t,c)
   if (t) return
end do

end subroutine

!===============================================================================
!  Trova la cella del blocco (igr,dbl) compresa in
!  [nl(1):nl(2),nl(3):nl(4),nl(5):nl(6)] che contiene "c"
!===============================================================================
subroutine findonor3(igr,ibl,id,jd,kd,il,jl,kl,m,t,c)
use prec
use moddef
use modpar, only : blogrd,blomet
implicit none

integer (kind=I4P), intent(in) :: igr,ibl
integer (kind=I4P) :: id,jd,kd,il,jl,kl,m
integer (kind=I4P) :: nii,njj,nkk
type (point) :: b,c

logical (kind=I1P) :: test,t
integer (kind=I4P) :: i,j,k

!...............................................................................

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk

do k=max(1,kl-m+1),min(nkk,kl+m)
   do j=max(1,jl-m+1),min(njj,jl+m)
      do i=max(1,il-m+1),min(nii,il+m)

         ! considera solo le facce del blocco
         test = max(i-il,j-jl,k-kl).eq.m .or. min(i-il,j-jl,k-kl).eq.(1-m)
         if (.not.test) goto 10

         ! faccia I-
         b = mean(blogrd(ibl,igr)%nodo(i-1,j-1:j,k-1:k))
         test = ((c-b).dot.(-blomet(ibl,igr)%cella(i-1,j,k)%sni)) .gt. (0.0d0)
         if (test) goto 10

         ! faccia I+
         b = mean(blogrd(ibl,igr)%nodo(i,j-1:j,k-1:k))
         test = ((c-b).dot.blomet(ibl,igr)%cella(i,j,k)%sni) .gt. (0.0d0)
         if (test) goto 10

         ! faccia J-
         b = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1,k-1:k))
         test = ((c-b).dot.(-blomet(ibl,igr)%cella(i,j-1,k)%snj)) .gt. (0.0d0)
         if (test) goto 10

         ! faccia J+
         b = mean(blogrd(ibl,igr)%nodo(i-1:i,j,k-1:k))
         test = ((c-b).dot.blomet(ibl,igr)%cella(i,j,k)%snj) .gt. (0.0d0)
         if (test) goto 10

         ! faccia K-
         b = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k-1))
         test = ((c-b).dot.(-blomet(ibl,igr)%cella(i,j,k-1)%snk)) .gt. (0.0d0)
         if (test) goto 10

         ! faccia K+
         b = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k))
         test = ((c-b).dot.blomet(ibl,igr)%cella(i,j,k)%snk) .gt. (0.0d0)
         if (test) goto 10

         t  = .true.
         id = i
         jd = j
         kd = k

         return

10       continue

      end do
   end do
end do

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  cerca il vertice della prima cella donatrice "md(1)" più vicino a "c"
!===============================================================================
subroutine nearest_vertex(igr,md,c)
use prec
use moddef, only : donor,point,dist
use modpar, only : donmax,blogrd
implicit none

integer (kind=I4P), intent(in) :: igr
integer (kind=I4P) :: nii,njj,nkk
type (donor)     :: md(donmax)
type (point)     :: c,b
integer (kind=I4P) :: bd,id,jd,kd
integer (kind=I4P) :: i0,i1,iv,jv,kv
integer (kind=I4P) :: i,j,k
real    (kind=R8P) :: dn,dl

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! per comodità
bd = md(1)%b
id = md(1)%i
jd = md(1)%j
kd = md(1)%k

nii = blogrd(bd,igr)%ni
njj = blogrd(bd,igr)%nj
nkk = blogrd(bd,igr)%nk

! della cella donatrice cerca il vertice più vicino
dn = 1.0d30
iv = id
jv = jd
kv = kd
do k=kd-1,kd
   do j=jd-1,jd
      do i=id-1,id
         b = blogrd(bd,igr)%nodo(i,j,k)
         dl = dist(b,c)
         if (dl.lt.dn) then
            dn = dl
            iv = i
            jv = j
            kv = k
         else if (dl.eq.dn) then
            ! a parità di distanza privilegia il vertice che NON si trova sulle
            ! cornici
            i0 = iv*(nii-iv)*jv*(njj-jv)*kv*(nkk-kv)
            i1 = i *(nii-i )*j *(njj-j )*k *(nkk-k )
            if (i0.eq.0 .and. i1.gt.0) then
               iv = i
               jv = j
               kv = k
            end if
         end if
      end do
   end do
end do

md(2)%b = bd
md(2)%i = iv
md(2)%j = jv
md(2)%k = kv

return
end subroutine nearest_vertex

!===============================================================================
!  assegna come celle donatrici le 8 celle adiacenti al vertice in input
!===============================================================================
subroutine multidonor(igr,md,c)
use prec
use moddef, only : donor,point
use modpar, only : donmax,blogrd
implicit none
integer (kind=I4P), parameter :: nd = 8

integer (kind=I4P), intent(in) :: igr
type (donor)     :: md(donmax)
type (point)     :: c,tri(donmax)
integer (kind=I4P) :: iv,jv,kv
integer (kind=I4P) :: bd,m,i,j,k,s
real    (kind=R8P) :: u(donmax),dmax,ddon

type (donor) :: a(donmax),q

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! per comodità
bd = md(2)%b
iv = md(2)%i
jv = md(2)%j
kv = md(2)%k

! loop sulle 8 celle che insistono sul vertice passato in input
m = 0
do k=kv,kv+1
   do j=jv,jv+1
      do i=iv,iv+1
         m = m+1
         call getcenter(igr,bd,i,j,k,tri(m))
         a(m)%b = bd
         ! DONATORI ANCHE SULLE CORNICI
!        a(m)%i = i
!        a(m)%j = j
!        a(m)%k = k
         ! NIENTE DONATORI SULLE CORNICI
         a(m)%i = min(blogrd(bd,igr)%ni,max(i,1))
         a(m)%j = min(blogrd(bd,igr)%nj,max(j,1))
         a(m)%k = min(blogrd(bd,igr)%nk,max(k,1))
      end do
   end do
end do

! interpolazione trilineare
call trilinear(tri,c,u)

! memorizza i donatori con i pesi
! mette la cella più grande nella prima posizione
dmax = 0d0
s = 1
do m=1,nd
   md(m)%b = a(m)%b
   md(m)%i = a(m)%i
   md(m)%j = a(m)%j
   md(m)%k = a(m)%k
   md(m)%w = real(u(m),R8P)
   call getsize(igr,a(m)%b,a(m)%i,a(m)%j,a(m)%k,ddon)
   if (ddon.gt.dmax) then
      dmax = ddon
      s = m
   end if
end do
q     = md(s)
md(s) = md(1)
md(1) = q

return
end

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  Trova il vertice del blocco (igr,ibl) più vicino a "bc"
!===============================================================================
subroutine ricerca_per_linee_coordinate(igr,ibl,step,il,jl,kl,bc)
use prec
use moddef, only : point,dist
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,step
integer (kind=I4P) :: il,jl,kl,nl(6)
integer (kind=I4P) :: it,jt,kt
logical (kind=I1P) :: t1,t2
real    (kind=R8P) :: dl,dt
type    (point)  :: bc,pt

pt = blogrd(ibl,igr)%nodo(il,jl,kl)
nl = (/ 0,blogrd(ibl,igr)%ni, 0,blogrd(ibl,igr)%nj, 0,blogrd(ibl,igr)%nk /)
dl = dist(bc,pt)
t1 = .true.
do while (t1)
   t1 = .false.

! cambia un indice alla volta
!
   t2 = .true.
   do while (t2)
      t2 = .false.
      it = min(il+step,nl(2))
      pt = blogrd(ibl,igr)%nodo(it,jl,kl)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = max(il-step,nl(1))
      pt = blogrd(ibl,igr)%nodo(it,jl,kl)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      jt = min(jl+step,nl(4))
      pt = blogrd(ibl,igr)%nodo(il,jt,kl)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; jl = jt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      jt = max(jl-step,nl(3))
      pt = blogrd(ibl,igr)%nodo(il,jt,kl)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; jl = jt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      kt = min(kl+step,nl(6))
      pt = blogrd(ibl,igr)%nodo(il,jl,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      kt = max(kl-step,nl(5))
      pt = blogrd(ibl,igr)%nodo(il,jl,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

! cambia due indici alla volta
!
   t2 = .true.
   do while (t2)
      t2 = .false.
      it = min(il+step,nl(2))
      jt = min(jl+step,nl(4))
      kt = kl
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = min(il+step,nl(2))
      jt = max(jl-step,nl(3))
      kt = kl
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = max(il-step,nl(1))
      jt = min(jl+step,nl(4))
      kt = kl
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = max(il-step,nl(1))
      jt = max(jl-step,nl(3))
      kt = kl
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = min(il+step,nl(2))
      jt = jl
      kt = min(kl+step,nl(6))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = min(il+step,nl(2))
      jt = jl
      kt = max(kl-step,nl(5))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = max(il-step,nl(1))
      jt = jl
      kt = min(kl+step,nl(6))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = max(il-step,nl(1))
      jt = jl
      kt = max(kl-step,nl(5))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = il
      jt = min(jl+step,nl(4))
      kt = min(kl+step,nl(6))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = il
      jt = min(jl+step,nl(4))
      kt = max(kl-step,nl(5))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = il
      jt = max(jl-step,nl(3))
      kt = min(kl+step,nl(6))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

   t2 = .true.
   do while (t2)
      t2 = .false.
      it = il
      jt = max(jl-step,nl(3))
      kt = max(kl-step,nl(5))
      pt = blogrd(ibl,igr)%nodo(it,jt,kt)
      dt = dist(bc,pt)
      if (dt.lt.dl) then
         dl = dt ; il = it ; jl = jt ; kl = kt ; t1 = .true. ; t2 = .true.
      end if
   end do

end do

end subroutine

!===============================================================================
!  Trova la cella della cornice del blocco (igr,bd) più vicina a "c"
!===============================================================================
subroutine nearest_on_boundary(igr,bd,id,jd,kd,sl,c)
use prec
use moddef, only : point,dist
use modpar, only : blomet
implicit none
integer (kind=I4P) :: igr,bd,id,jd,kd
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: i,j,k,step
real (kind=R8P)    :: sl,ss
type (point)     :: c,b

nii = blomet(bd,igr)%ni
njj = blomet(bd,igr)%nj
nkk = blomet(bd,igr)%nk

step = max(nkk-1,1)
do k=1,nkk,step
   do j=1,njj
      do i=1,nii
         b = blomet(bd,igr)%cella(i,j,k)%cen
         ss = dist(c,b)
         if (ss.lt.sl) then
            id = i
            jd = j
            kd = k
            sl = ss
         end if
      end do
   end do
end do

step = max(njj-1,1)
do k=2,nkk-1
   do j=1,njj,step
      do i=1,nii
         b = blomet(bd,igr)%cella(i,j,k)%cen
         ss = dist(c,b)
         if (ss.lt.sl) then
            id = i
            jd = j
            kd = k
            sl = ss
         end if
      end do
   end do
end do

step = max(nii-1,1)
do k=2,nkk-1
   do j=2,njj-1
      do i=1,nii,step
         b = blomet(bd,igr)%cella(i,j,k)%cen
         ss = dist(c,b)
         if (ss.lt.sl) then
            id = i
            jd = j
            kd = k
            sl = ss
         end if
      end do
   end do
end do

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

function fuori_ingombro(c,ing)
use prec
use moddef, only : point
implicit none
logical (kind=I1P) :: fuori_ingombro
type (point), intent(in) :: c,ing(2)
fuori_ingombro = c%x.lt.ing(1)%x .or. c%y.lt.ing(1)%y .or. c%z.lt.ing(1)%z &
            .or. c%x.gt.ing(2)%x .or. c%y.gt.ing(2)%y .or. c%z.gt.ing(2)%z
return
end function fuori_ingombro

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

subroutine picksmaller(igr,cou,don)
use prec
implicit none
integer (kind=I4P) :: igr,cou,don(4,*)
integer (kind=I4P) :: n
real    (kind=R8P) :: sl,sd

if (cou.eq.1) return

call getsize(igr,don(1,1),don(2,1),don(3,1),don(4,1),sd)
do n=2,cou
   call getsize(igr,don(1,n),don(2,n),don(3,n),don(4,n),sl)
   if (sl.lt.sd) then
      sd = sl
      don(1,1) = don(1,n)
      don(2,1) = don(2,n)
      don(3,1) = don(3,n)
      don(4,1) = don(4,n)
   end if
end do

return
end subroutine picksmaller
