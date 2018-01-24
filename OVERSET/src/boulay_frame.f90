!===============================================================================
! ridefinisce le celle di cornice che si trovano nello strato limite
!===============================================================================
subroutine strato_limite_cornice(igr,ibl,nf,nl,td,fd)
!{{{3
use prec
use modpar
implicit none
integer (kind=I4P) :: igr,ibl,nf,nl(6),td,fd
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: i,j,k,m
integer (kind=I4P), dimension(4) :: imin,imax,jmin,jmax,kmin,kmax
integer (kind=I4P), dimension(4) :: niii,njjj,nkkk,iski,jski,kski

nii = blochi(ibl,igr)%ni
njj = blochi(ibl,igr)%nj
nkk = blochi(ibl,igr)%nk

if (fd.gt.0) go to 47

! donatori per il centro cella
select case (nf)

case (1:2)

   ! facce i=0 o i=ni
   jmin = (/ nl(3),nl(3),nl(3),nl(4) /)
   jmax = (/ nl(4),nl(4),nl(3),nl(4) /)
   kmin = (/ nl(5),nl(6),nl(5),nl(5) /)
   kmax = (/ nl(5),nl(6),nl(6),nl(6) /)
   iski = (/ 0,0,0,0 /)
   jski = (/ 0,0,1,-1 /)
   kski = (/ 1,-1,0,0 /)
   niii = (/ 0,0,0,0 /)
   njjj = (/ 0,0,njj,njj /)
   nkkk = (/ nkk,nkk,0,0 /)
   ! loop sui 4 lati della faccia
   do m=1,4
      do k=kmin(m),kmax(m)
         do j=jmin(m),jmax(m)
            do i=nl(1),nl(2)
               call sub1_framebl(igr,ibl,niii(m),njjj(m),nkkk(m),&
                                 iski(m),jski(m),kski(m),i,j,k,td)
            end do
         end do
      end do
   end do

case (3:4)

   ! facce j=0 o j=nj
   imin = (/ nl(1),nl(1),nl(1),nl(2) /)
   imax = (/ nl(2),nl(2),nl(1),nl(2) /)
   kmin = (/ nl(5),nl(6),nl(5),nl(5) /)
   kmax = (/ nl(5),nl(6),nl(6),nl(6) /)
   iski = (/ 0,0,1,-1 /)
   jski = (/ 0,0,0,0 /)
   kski = (/ 1,-1,0,0 /)
   niii = (/ 0,0,nii,nii /)
   njjj = (/ 0,0,0,0 /)
   nkkk = (/ nkk,nkk,0,0 /)
   ! loop sui 4 lati della faccia
   do m=1,4
      do k=kmin(m),kmax(m)
         do i=imin(m),imax(m)
            do j=nl(3),nl(4)
               call sub1_framebl(igr,ibl,niii(m),njjj(m),nkkk(m),&
                                 iski(m),jski(m),kski(m),i,j,k,td)
            end do
         end do
      end do
   end do

case (5:6)

   ! facce k=0 o k=nk
   imin = (/ nl(1),nl(1),nl(1),nl(2) /)
   imax = (/ nl(2),nl(2),nl(1),nl(2) /)
   jmin = (/ nl(3),nl(4),nl(3),nl(3) /)
   jmax = (/ nl(3),nl(4),nl(4),nl(4) /)
   iski = (/ 0,0,1,-1 /)
   jski = (/ 1,-1,0,0 /)
   kski = (/ 0,0,0,0 /)
   niii = (/ 0,0,nii,nii /)
   njjj = (/ njj,njj,0,0 /)
   nkkk = (/ 0,0,0,0 /)
   ! loop sui 4 lati della faccia
   do m=1,4
      do j=jmin(m),jmax(m)
         do i=imin(m),imax(m)
            do k=nl(5),nl(6)
               call sub1_framebl(igr,ibl,niii(m),njjj(m),nkkk(m),&
                                 iski(m),jski(m),kski(m),i,j,k,td)
            end do
         end do
      end do
   end do

end select

return
!...............................................................................
47 continue

! donatori per il centro faccia
select case (nf)

case (1:2)

   ! facce i=0 o i=ni
   i  = nl(nf)
   jmin = (/ nl(3),nl(3),nl(3),nl(4) /)
   jmax = (/ nl(4),nl(4),nl(3),nl(4) /)
   kmin = (/ nl(5),nl(6),nl(5),nl(5) /)
   kmax = (/ nl(5),nl(6),nl(6),nl(6) /)
   iski = (/ 0,0,0,0 /)
   jski = (/ 0,0,1,-1 /)
   kski = (/ 1,-1,0,0 /)
   niii = (/ 0,0,0,0 /)
   njjj = (/ 0,0,njj,njj /)
   nkkk = (/ nkk,nkk,0,0 /)
   ! loop sui 4 lati della faccia
   do m=1,4
      do k=kmin(m),kmax(m)
         do j=jmin(m),jmax(m)
            call sub2_framebl(igr,ibl,niii(m),njjj(m),nkkk(m),&
                              iski(m),jski(m),kski(m),i,j,k,nf,td)
         end do
      end do
   end do

case (3:4)

   ! facce j=0 o j=nj
   j  = nl(nf)
   imin = (/ nl(1),nl(1),nl(1),nl(2) /)
   imax = (/ nl(2),nl(2),nl(1),nl(2) /)
   kmin = (/ nl(5),nl(6),nl(5),nl(5) /)
   kmax = (/ nl(5),nl(6),nl(6),nl(6) /)
   iski = (/ 0,0,1,-1 /)
   jski = (/ 0,0,0,0 /)
   kski = (/ 1,-1,0,0 /)
   niii = (/ 0,0,nii,nii /)
   njjj = (/ 0,0,0,0 /)
   nkkk = (/ nkk,nkk,0,0 /)
   ! loop sui 4 lati della faccia
   do m=1,4
      do k=kmin(m),kmax(m)
         do i=imin(m),imax(m)
            call sub2_framebl(igr,ibl,niii(m),njjj(m),nkkk(m),&
                              iski(m),jski(m),kski(m),i,j,k,nf,td)
         end do
      end do
   end do

case (5:6)

   ! facce k=0 o k=nk
   k  = nl(nf)
   imin = (/ nl(1),nl(1),nl(1),nl(2) /)
   imax = (/ nl(2),nl(2),nl(1),nl(2) /)
   jmin = (/ nl(3),nl(4),nl(3),nl(3) /)
   jmax = (/ nl(3),nl(4),nl(4),nl(4) /)
   iski = (/ 0,0,1,-1 /)
   jski = (/ 1,-1,0,0 /)
   kski = (/ 0,0,0,0 /)
   niii = (/ 0,0,nii,nii /)
   njjj = (/ njj,njj,0,0 /)
   nkkk = (/ 0,0,0,0 /)
   ! loop sui 4 lati della faccia
   do m=1,4
      do j=jmin(m),jmax(m)
         do i=imin(m),imax(m)
            call sub2_framebl(igr,ibl,niii(m),njjj(m),nkkk(m),&
                              iski(m),jski(m),kski(m),i,j,k,nf,td)
         end do
      end do
   end do

end select

return
!3}}}
end subroutine strato_limite_cornice

!===============================================================================
!===============================================================================
subroutine sub1_framebl(igr,ibl,nii,njj,nkk,is,js,ks,ii,jj,kk,td)
!{{{4
use prec
use moddef, only : donor,operator(-)
use modpar
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,nii,njj,nkk,is,js,ks,ii,jj,kk,td

integer (kind=I4P) :: test_indici
integer (kind=I4P) :: i,j,k,nd,p,pp,famn,nwv(9)
integer (kind=I4P) :: ip,jp,kp, iq,jq,kq,dq,qq
real (kind=R8P)    :: dw
type (donor)     :: md(donmax)
type (point)     :: e1,e2,c,n

! distanza della prima cella dalla parete
dw = blochi(ibl,igr)%cella(ii,jj,kk)%dist

! determina il patch di parete che genera lo strato limite ...
p = blochi(ibl,igr)%cella(ii,jj,kk)%wall
if (p.lt.1)  return

! ... ne prende la famiglia ...
i = mod(bp(4,p),100)
if (i.gt.0) then
   ! ... se la famiglia esiste (i>0), blt è lo spessore di s.l. della famiglia
   blt = boulaf(i)
else
   ! ... altrimenti è lo spessore di s.l. del patch
   blt = boulap(p)
end if
if (dw.gt.blt) return

! i,j,k servono come indici di appoggio, perché se modifico ii,jj,kk si incazza
! come una vespa
i = ii
j = jj
k = kk

! se la direzione normale al patch non coincide con quella in input, ritorna
select case (bp(2,p))
   case (1) ; if (is.ne.1) return
   case (2) ; if (is.ne.-1) return
   case (3) ; if (js.ne.1) return
   case (4) ; if (js.ne.-1) return
   case (5) ; if (ks.ne.1) return
   case (6) ; if (ks.ne.-1) return
end select

! individua, se esistono, patch di parete nella stessa famiglia di "p"
call patch_relatives(p,famn)
if (famn.eq.0) return

! centro cella
call getcenter(igr,ibl,i,j,k,c)

! normale approssimativa
call getcenter(igr,ibl,i+is,j+js,k+ks,e1)
n = e1-c

! vertice più vicino scelto tra i patch selezionati
call strato_limite_4frame(igr,famn,nwv,e1,e2,c,n)

! se ha trovato un opportuno "nwv", cerca i donatori per le celle dello
! s.l. a partire da quel vertice, altrimenti ritorna
if (nwv(1).gt.0) then

   dw = blochi(ibl,igr)%cella(i,j,k)%dist
   pp = blochi(ibl,igr)%cella(i,j,k)%wall
   test_indici = (i-1)*(nii-i)*abs(is) &  ! se lo s.l. è spesso rispetto alle
               + (j-1)*(njj-j)*abs(js) &  ! dimensioni del blocco può succedere
               + (k-1)*(nkk-k)*abs(ks)    ! di sforare un indice
   do while ((dw.le.blt) .and. (pp.eq.p) .and. (test_indici.ge.0))

      ! se pp.ne.p vuol dire che la parete più vicina alla cella non è più
      ! quella per la quale sto calcolando lo strato limite (ad es. mi trovo in
      ! corrispondenza di un angolo)

      ! cerca iterativamente il vertice più vicino a "c" muovendosi prima lungo
      ! la coordinata uscente dalla parete e poi nel piano parallelo
      c = blomet(ibl,igr)%cella(i,j,k)%cen
      dq = 1
      qq = 1
      do while (dq.ne.0 .and. qq.lt.4)
         ! cerca il vertice con distanza da parete più prossima a partire da nwv
         ! lo memorizza in ip,jp,kp lasciando invariato nwv
         call ricerca_lungo_normale(igr,nwv,ip,jp,kp,dw)
         ! sul livello determinato da ip o jp o kp (dipende dalla faccia)
         ! determina la faccetta che contiene la proiezione di "c"
         ! memorizza il vertice della faccetta più vicino a "c" in ip,jp,kp
         ! aggiorna gli indici di nwv "paralleli" al piano (ad es. se la parete
         ! è sul piano J=0 aggiorna nwv(3) e nwv(5))
         iq = ip
         jq = jp
         kq = kp
         call ricerca_sul_piano(igr,nwv,ip,jp,kp,e1,e2,c,n)
         if (nwv(1).lt.1) return
         dq = abs(ip-iq)+abs(jp-jq)+abs(kp-kq)
         qq = qq+1
      end do

      ! individua donatori (attorno al vertice più vicino) e calcola pesi
      call trilinear_strato_limite(igr,nwv,ip,jp,kp,nd,md,dw,e1,e2,c)

      ! crea la cella chimera
      call newchi(igr,ibl,i,j,k,td,nd,md)

      ! passa al punto successivo di s.l.
      i = i+is
      j = j+js
      k = k+ks
      dw = blochi(ibl,igr)%cella(i,j,k)%dist
      pp = blochi(ibl,igr)%cella(i,j,k)%wall
      test_indici = (i-1)*(nii-i)*abs(is) &
                  + (j-1)*(njj-j)*abs(js) &
                  + (k-1)*(nkk-k)*abs(ks)
   end do

end if

return
!4}}}
end subroutine sub1_framebl

!===============================================================================
!===============================================================================
subroutine sub2_framebl(igr,ibl,nii,njj,nkk,is,js,ks,ii,jj,kk,nf,td)
!{{{5
use prec
use moddef
use modpar
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,nii,njj,nkk,is,js,ks,ii,jj,kk,nf,td

integer (kind=I4P) :: test_indici
integer (kind=I4P) :: i,j,k,ix,jx,kx,nd,p,pp,famn,nwv(9)
integer (kind=I4P) :: ip,jp,kp, iq,jq,kq,dq,qq
real (kind=R8P)    :: dw
type (donor)     :: md(donmax)
type (point)     :: e1,e2,c,n

! distanza della prima cella dalla parete
dw = blochi(ibl,igr)%cella(ii,jj,kk)%dist
if (dw.gt.(1d3*boulay)) return

! determina il patch di parete che genera lo strato limite ...
p = blochi(ibl,igr)%cella(ii,jj,kk)%wall
! ... ne prende la famiglia ...
i = mod(bp(4,p),100)
if (i.gt.0) then
   ! ... se la famiglia esiste (i>0), blt è lo spessore di s.l. della famiglia
   blt = boulaf(i)
else
   ! ... altrimenti è lo spessore di s.l. del patch
   blt = boulap(p)
end if
if (dw.gt.blt) return

! i,j,k servono come indici di appoggio, perché se modifico ii,jj,kk si incazza
! come una vespa
i = ii
j = jj
k = kk

! se la direzione normale al patch non coincide con quella in input, ritorna
! en passant, assegna gli indici per determinare la cella di cornice più esterna
ix = 0
jx = 0
kx = 0
select case (bp(2,p))
   case (1)
      if (is.ne.1) return
      ix = -1
   case (2)
      if (is.ne.-1) return
      ix = +1
   case (3)
      if (js.ne.1) return
      jx = -1
   case (4)
      if (js.ne.-1) return
      jx = +1
   case (5)
      if (ks.ne.1) return
      kx = -1
   case (6)
      if (ks.ne.-1) return
      kx = +1
end select

! individua, se esistono, patch di parete nella stessa famiglia di "p"
call patch_relatives(p,famn)
if (famn.eq.0) return

! centro faccia
call centro_faccia_cornice(igr,ibl,i,j,k,nf,c)

! normale approssimativa
call centro_faccia_cornice(igr,ibl,i+is,j+js,k+ks,nf,e1)
n = e1-c

! vertice più vicino scelto tra i patch selezionati
call strato_limite_4frame(igr,famn,nwv,e1,e2,c,n)

! se ha trovato un opportuno "nwv", cerca i donatori per le celle dello
! s.l. a partire da quel vertice, altrimenti ritorna
if (nwv(1).gt.0) then

   dw = blochi(ibl,igr)%cella(i,j,k)%dist
   pp = blochi(ibl,igr)%cella(i,j,k)%wall
   test_indici = (i-1)*(nii-i)*abs(is) &  ! se lo s.l. è spesso rispetto alle
               + (j-1)*(njj-j)*abs(js) &  ! dimensioni del blocco può succedere
               + (k-1)*(nkk-k)*abs(ks)    ! di sforare un indice

   do while ((dw.le.blt) .and. (pp.eq.p) .and. (test_indici.ge.0))

      ! se pp.ne.p vuol dire che la parete più vicina alla cella non è più
      ! quella per la quale sto calcolando lo strato limite (ad es. mi trovo in
      ! corrispondenza di un angolo)

      ! cerca iterativamente il vertice più vicino a "c" muovendosi prima lungo
      ! la coordinata uscente dalla parete e poi nel piano parallelo
      call centro_faccia_cornice(igr,ibl,i,j,k,nf,c)
      dq = 1
      qq = 1
      do while (dq.ne.0 .and. qq.lt.4)
         ! cerca il vertice con distanza da parete più prossima a partire da nwv
         ! lo memorizza in ip,jp,kp lasciando invariato nwv
         call ricerca_lungo_normale(igr,nwv,ip,jp,kp,dw)
         ! sul livello determinato da ip o jp o kp (dipende dalla faccia)
         ! determina la faccetta che contiene la proiezione di "c"
         ! memorizza il vertice della faccetta più vicino a "c" in ip,jp,kp
         ! aggiorna gli indici di nwv "paralleli" al piano (ad es. se la parete
         ! è sul piano J=0 aggiorna nwv(3) e nwv(5))
         iq = ip
         jq = jp
         kq = kp
         call ricerca_sul_piano(igr,nwv,ip,jp,kp,e1,e2,c,n)
         if (nwv(1).lt.1) return
         dq = abs(ip-iq)+abs(jp-jq)+abs(kp-kq)
         qq = qq+1
      end do

      ! individua donatori (attorno al vertice più vicino) e calcola pesi
      call trilinear_strato_limite(igr,nwv,ip,jp,kp,nd,md,dw,e1,e2,c)

      ! crea la cella chimera
      call newchi(igr,ibl,i,j,k,td,nd,md)
      call copychi(igr,ibl,i,j,k,i+ix,j+jx,k+kx)

      ! passa al punto successivo di s.l.
      i = i+is
      j = j+js
      k = k+ks
      dw = blochi(ibl,igr)%cella(i,j,k)%dist
      pp = blochi(ibl,igr)%cella(i,j,k)%wall
      test_indici = (i-1)*(nii-i)*abs(is) &
                  + (j-1)*(njj-j)*abs(js) &
                  + (k-1)*(nkk-k)*abs(ks)
   end do

end if

return
!5}}}
end subroutine sub2_framebl

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
! Equivalente di multidonor per le celle chimera di strato limite
! calcola i pesi dei donatori (le celle attorno al vertice ip,jp,kp) con una
! interpolazione trilineare nel riferimento locale
!===============================================================================
subroutine trilinear_strato_limite(igr,nwv,ip,jp,kp,nd,md,dw,e1,e2,c)
!{{{6
use prec
use moddef
use modpar
implicit none

integer (kind=I4P) :: igr,nwv(9),ip,jp,kp,nd
type (donor)     :: md(donmax)
real    (kind=R8P) :: dw
type (point)     :: e1,e2,c

integer (kind=I4P) :: i0,j0,k0,i1,j1,k1
integer (kind=I4P) :: i,j,k,l,m
integer (kind=I4P) :: bd
integer (kind=I4P) :: is,js,ks                 ,nplane
real (kind=R8P) :: u(donmax)                   ,ds(donmax),dn,sm,pm
real (kind=R8P) :: ph,pl
type (donor) :: a(donmax)
type (point) :: cp,aa,bb,cc
type (point) :: tri(donmax)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

call getvertex(igr,nwv(1),ip,jp,kp,cp)
bd = nwv(1)
is = nwv(6)
js = nwv(7)
ks = nwv(8)

! di default, nessun donatore
nd = 0

! gli indici delle 8 celle che insistono sul vertice trovato sopra
i0 = max(1,ip)
i1 = min(blochi(bd,igr)%ni,ip+1)
j0 = max(1,jp)
j1 = min(blochi(bd,igr)%nj,jp+1)
k0 = max(1,kp)
k1 = min(blochi(bd,igr)%nk,kp+1)

m = (i1-i0)*(j1-j0)*(k1-k0)
if (m.eq.0) go to 137

! ci sono 8 donatori: interpolazione trilineare

!  numero donatori
nd = 8

! coordinate cella chimera e donatori nel riferimento locale
aa = c-cp
bb%x = aa.dot.e1
bb%y = aa.dot.e2
bb%z = sqrt(dw)

m = 0
do k=k0,k1
do j=j0,j1
do i=i0,i1
   m = m+1
   call getcenter(igr,bd,i,j,k,cc)
   aa = cc-cp
   tri(m)%x = aa.dot.e1
   tri(m)%y = aa.dot.e2
   call getdis(igr,bd,i,j,k,dn)
   tri(m)%z = sqrt(dn)
   a(m)%b = bd
   a(m)%i = i
   a(m)%j = j
   a(m)%k = k
end do
end do
end do

! interpolazione trilineare
call trilinear(tri,bb,u)

! memorizza i donatori con i pesi
do m=1,nd
   md(m)%b = a(m)%b
   md(m)%i = a(m)%i
   md(m)%j = a(m)%j
   md(m)%k = a(m)%k
   md(m)%w = real(u(m),R8P)
end do

return

137 continue

! ci sono meno di 8 donatori: interpolazione mista

! ordina i donatori su due piani in base alla "distanza" dalla parete
m = 0
nplane = 0
if (is.ne.0) then
   l = (k1-k0+1)*(j1-j0+1)
   do i=i0,i1
      nplane = nplane+1
      do k=k0,k1
         do j=j0,j1
            m = m+1
            call getcenter(igr,bd,i,j,k,cc)
            u(m) = sqrt(dist(cc,c))
            a(m) = donor(bd,i,j,k,0.0)
            call getdis(igr,bd,i,j,k,ds(m))
         end do
      end do
   end do
else if (js.ne.0) then
   l = (k1-k0+1)*(i1-i0+1)
   do j=j0,j1
      nplane = nplane+1
      do k=k0,k1
         do i=i0,i1
            m = m+1
            call getcenter(igr,bd,i,j,k,cc)
            u(m) = sqrt(dist(cc,c))
            a(m) = donor(bd,i,j,k,0.0)
            call getdis(igr,bd,i,j,k,ds(m))
         end do
      end do
   end do
else if (ks.ne.0) then
   l = (j1-j0+1)*(i1-i0+1)
   do k=k0,k1
      nplane = nplane+1
      do j=j0,j1
         do i=i0,i1
            m = m+1
            call getcenter(igr,bd,i,j,k,cc)
            u(m) = sqrt(dist(cc,c))
            a(m) = donor(bd,i,j,k,0.0)
            call getdis(igr,bd,i,j,k,ds(m))
         end do
      end do
   end do
else
   stop 'Errore in trilinear_strato_limite: [ijk]skip nulli'
end if

!  numero donatori
nd = m

! calcola la somma dei pesi per la normalizzazione successiva
dn = 0.0d0
do m=1,nd
   u(m) = 1.0d0/max(u(m),1.0d-20)
   dn = dn + u(m)
end do
dn = 1.0d0/dn

if (nplane.eq.2) then
   ! il peso assegnato ai piani viene ridistribuito linearmente in base alla
   ! distanza della cella chimera dai due piani
   do m=1,l
      sm = u(m)+u(m+l)
      ph = abs(dw-ds(m))
      pl = abs(ds(m+l)-dw)
      pm = sm/(ph+pl)
      u(m)   = pm*pl
      u(m+l) = pm*ph
   end do
end if

! riempie la matrice md dei donatori
do m=1,nd
   md(m)%b = a(m)%b
   md(m)%i = a(m)%i
   md(m)%j = a(m)%j
   md(m)%k = a(m)%k
   md(m)%w = real(u(m)*dn,R8P)
end do

return
!6}}}
end subroutine trilinear_strato_limite

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
! partendo da nwv cerca il vertice a distanza dalla parete più prossima a dw
!===============================================================================
subroutine ricerca_lungo_normale(igr,nwv,ip,jp,kp,dw)
!{{{7
use prec
use modpar
use moddef, only: dist
implicit none

integer (kind=I4P) :: igr,nwv(9),ip,jp,kp
real    (kind=R8P) :: dw
integer (kind=I4P) :: bd,is,js,ks
real (kind=R8P)    :: ph,pl
type (point) :: pwall,pcur

! variabili utilizzate per comodita' (da non confondere is,js,ks con le
! omonime variabili nelle subroutine chiamanti)
bd = nwv(1)
ip = nwv(3)  ! inizializzazione
jp = nwv(4)  ! inizializzazione
kp = nwv(5)  ! inizializzazione
is = nwv(6)
js = nwv(7)
ks = nwv(8)

call getvertex(igr,bd,ip,jp,kp,pwall)

! inizia la ricerca (NOTA: lavorando con i vertici non posso utilizzare
! bufdis che e' definito nei centri cella)
ph = 0.0d0
do while (ph.lt.dw &
         .and. ((ip+is).ge.0 .and. (ip+is).le.blochi(bd,igr)%ni) &
         .and. ((jp+js).ge.0 .and. (jp+js).le.blochi(bd,igr)%nj) &
         .and. ((kp+ks).ge.0 .and. (kp+ks).le.blochi(bd,igr)%nk) &
         )
   ip = ip+is
   jp = jp+js
   kp = kp+ks
   call getvertex(igr,bd,ip,jp,kp,pcur)
   ph = dist(pwall,pcur)
end do

call getvertex(igr,bd,ip-is,jp-js,kp-ks,pcur)
pl = dist(pwall,pcur)
if ((ph-dw).gt.(dw-pl)) then
   ip = ip-is
   jp = jp-js
   kp = kp-ks
end if

return
!7}}}
end subroutine ricerca_lungo_normale

!===============================================================================
! cerca il vertice a distanza controvariante minima da "c"
!===============================================================================
subroutine ricerca_sul_piano(igr,nwv,ip,jp,kp,e1,e2,c,n)
!{{{8
use prec
use moddef
use modpar, only : blomet,blogrd,bp,blt
implicit none

integer (kind=I4P), intent(in) :: igr
integer (kind=I4P) :: nwv(9)
type (point), intent(in) :: c,n
!
logical (kind=I1P) :: coocontr,controllo_normali,controllo_distanza
!
integer (kind=I4P) :: nl(6),tgr
integer (kind=I4P) :: ip,jp,kp,pbl,nf
integer (kind=I4P) :: imin,jmin,kmin
integer (kind=I4P) :: is,js,ks
logical (kind=I1P) :: cs,ct
real (kind=R8P)    :: x1,x2,dis_norm,dis_plan
type (point)     :: q1,q2,e1,e2,cp,ry
type (point)     :: norm_face

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tgr  = 2**(igr-1)
pbl  = nwv(1)
nf   = nwv(2)

! nl   = bp(5:10,nwv(9))/tgr
nl(1) = max(bp( 5,nwv(9))/tgr , ip-2)
nl(2) = min(bp( 6,nwv(9))/tgr , ip+2)
nl(3) = max(bp( 7,nwv(9))/tgr , jp-2)
nl(4) = min(bp( 8,nwv(9))/tgr , jp+2)
nl(5) = max(bp( 9,nwv(9))/tgr , kp-2)
nl(6) = min(bp(10,nwv(9))/tgr , kp+2)

nwv(1) = 0

! NOTA: a seconda della faccia, ip o jp o kp indica il livello su cui fare la
! ricerca e NON va sovrascritto con nl(.)

! cerca la faccetta che contiene la proiezione
if (nf.lt.3) then
   is = 3-2*nf
   do kp=nl(5),nl(6)-1
   do jp=nl(3),nl(4)-1
      norm_face = real(is,kind(0d0)) * blomet(pbl,igr)%cella(ip,jp+1,kp+1)%sni
      ! se fallisce il controllo sulle normali passa direttamente alla prossima
      ! faccetta
      if (controllo_normali(n,norm_face)) then
         cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
         q1 = blogrd(pbl,igr)%nodo(ip,jp+1,kp)
         q2 = blogrd(pbl,igr)%nodo(ip,jp,kp+1)
         e1 = q1-cp
         e2 = q2-cp
         ry = c-cp
         !
         dis_norm = ry.dot.norm_face
         dis_plan = normq(ry-dis_norm*norm_face)
         cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
         ct = controllo_distanza(dis_norm**2,4d0*blt)
         !
         if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
            imin = ip
            jmin = jp+nint(x1)
            kmin = kp+nint(x2)
            nwv(4) = jmin
            nwv(5) = kmin
            go to 19
         end if
         cp = blogrd(pbl,igr)%nodo(ip,jp+1,kp+1)
         e1 = q2-cp  !!!
         e2 = q1-cp  !!!
         ry = c-cp
         !
         dis_norm = ry.dot.norm_face
         dis_plan = normq(ry-dis_norm*norm_face)
         cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
         ct = controllo_distanza(dis_norm**2,4d0*blt)
         !
         if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
            imin = ip
            jmin = jp+1-nint(x1)
            kmin = kp+1-nint(x2)
            nwv(4) = jmin
            nwv(5) = kmin
            go to 19
         end if
      end if
   end do
   end do

else if (nf.gt.4) then
   ks = 11-2*nf
   do jp=nl(3),nl(4)-1
   do ip=nl(1),nl(2)-1
      norm_face = real(ks,kind(0d0)) * blomet(pbl,igr)%cella(ip+1,jp+1,kp)%snk
      ! se fallisce il controllo sulle normali passa direttamente alla prossima
      ! faccetta
      if (controllo_normali(n,norm_face)) then
         cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
         q1 = blogrd(pbl,igr)%nodo(ip+1,jp,kp)
         q2 = blogrd(pbl,igr)%nodo(ip,jp+1,kp)
         e1 = q1-cp
         e2 = q2-cp
         ry = c-cp
         !
         dis_norm = ry.dot.norm_face
         dis_plan = normq(ry-dis_norm*norm_face)
         cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
         ct = controllo_distanza(dis_norm**2,4d0*blt)
         !
         if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
            imin = ip+nint(x1)
            jmin = jp+nint(x2)
            kmin = kp
            nwv(3) = imin
            nwv(4) = jmin
            go to 19
         end if
         cp = blogrd(pbl,igr)%nodo(ip+1,jp+1,kp)
         e1 = q2-cp  !!!
         e2 = q1-cp  !!!
         ry = c-cp
         !
         dis_norm = ry.dot.norm_face
         dis_plan = normq(ry-dis_norm*norm_face)
         cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
         ct = controllo_distanza(dis_norm**2,4d0*blt)
         !
         if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
            imin = ip+1-nint(x1)
            jmin = jp+1-nint(x2)
            kmin = kp
            nwv(3) = imin
            nwv(4) = jmin
            go to 19
         end if
      end if
   end do
   end do

else
   js = 7-2*nf
   do kp=nl(5),nl(6)-1
   do ip=nl(1),nl(2)-1
      norm_face = real(js,kind(0d0)) * blomet(pbl,igr)%cella(ip+1,jp,kp+1)%snj
      ! se fallisce il controllo sulle normali passa direttamente alla prossima
      ! faccetta
      if (controllo_normali(n,norm_face)) then
         cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
         q1 = blogrd(pbl,igr)%nodo(ip+1,jp,kp)
         q2 = blogrd(pbl,igr)%nodo(ip,jp,kp+1)
         e1 = q1-cp
         e2 = q2-cp
         ry = c-cp
         !
         dis_norm = ry.dot.norm_face
         dis_plan = normq(ry-dis_norm*norm_face)
         cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
         ct = controllo_distanza(dis_norm**2,4d0*blt)
         !
         if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
            imin = ip+nint(x1)
            jmin = jp
            kmin = kp+nint(x2)
            nwv(3) = imin
            nwv(5) = kmin
            go to 19
         end if
         cp = blogrd(pbl,igr)%nodo(ip+1,jp,kp+1)
         e1 = q2-cp  !!!
         e2 = q1-cp  !!!
         ry = c-cp
         !
         dis_norm = ry.dot.norm_face
         dis_plan = normq(ry-dis_norm*norm_face)
         cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
         ct = controllo_distanza(dis_norm**2,4d0*blt)
         !
         if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
            imin = ip+1-nint(x1)
            jmin = jp
            kmin = kp+1-nint(x2)
            nwv(3) = imin
            nwv(5) = kmin
            go to 19
         end if
      end if
   end do
   end do
end if

! se arriva a questa istruzione, non ha trovato un patch che contenga la
! proiezione di "c"
return

! ha trovato la faccetta in cui cade la proiezione di "c"
19 continue

nwv(1) = pbl
ip = imin
jp = jmin
kp = kmin
e1 = vers0(e1)
cp = e1.cross.e2
e2 = cp.cross.e1
e2 = vers0(e2)

return
!8}}}
end subroutine ricerca_sul_piano

!===============================================================================
! cerca il vertice a distanza controvariante minima da "c"
! praticamente identica a strato_limite_findonor ma adattata ai punti
! di strato limite nelle cornici dei blocchi
!===============================================================================
subroutine strato_limite_4frame(igr,famn,nwv,e1,e2,c,n)
!{{{9
use prec
use moddef
use modpar, only : blomet,blogrd,bp,blt,famw,ingopatch
implicit none

integer (kind=I4P), intent(in) :: igr,famn
integer (kind=I4P) :: nwv(9)
type (point), intent(in)     :: c,n
!
logical (kind=I1P) :: coocontr,controllo_normali,controllo_distanza
logical (kind=I1P) :: fuori_ingombro
!
integer (kind=I4P) :: p,q,nl(6),tgr
integer (kind=I4P) :: ip,jp,kp,pbl,nf
integer (kind=I4P) :: imin,jmin,kmin,pmin
integer (kind=I4P) :: is,js,ks
logical (kind=I4P) :: cs,ct
real (kind=R8P)    :: x1,x2,dis_norm,dis_plan
type (point)     :: q1,q2,e1,e2,cp,ry
type (point)     :: norm_face

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! cerca il vertice di parete più vicino
tgr = 2**(igr-1)
do q=1,famn
   p   = famw(q)
   if (fuori_ingombro(c,ingopatch(1,p))) go to 13
   pbl = bp(1,p)
   nf  = bp(2,p)
   nl  = bp(5:10,p)/tgr

   ! cerca sul patch più vicino la faccetta che contiene la proiezione
   if (nf.lt.3) then
      ip = nl(1)
      is = 3-2*nf
      js = 0
      ks = 0
      do kp=nl(5),nl(6)-1
      do jp=nl(3),nl(4)-1
         norm_face = real(is,kind(0d0)) * blomet(pbl,igr)%cella(ip,jp+1,kp+1)%sni
         ! se fallisce il controllo sulle normali passa direttamente alla
         ! prossima faccetta
         if (controllo_normali(n,norm_face)) then
            cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
            q1 = blogrd(pbl,igr)%nodo(ip,jp+1,kp)
            q2 = blogrd(pbl,igr)%nodo(ip,jp,kp+1)
            e1 = q1-cp
            e2 = q2-cp
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip
               jmin = jp+nint(x1)
               kmin = kp+nint(x2)
               go to 19
            end if
            cp = blogrd(pbl,igr)%nodo(ip,jp+1,kp+1)
            e1 = q2-cp  !!!
            e2 = q1-cp  !!!
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip
               jmin = jp+1-nint(x1)
               kmin = kp+1-nint(x2)
               go to 19
            end if
         end if
      end do
      end do

   else if (nf.gt.4) then
      kp = nl(5)
      is = 0
      js = 0
      ks = 11-2*nf
      do jp=nl(3),nl(4)-1
      do ip=nl(1),nl(2)-1
         norm_face = real(ks,kind(0d0)) * blomet(pbl,igr)%cella(ip+1,jp+1,kp)%snk
         ! se fallisce il controllo sulle normali passa direttamente alla
         ! prossima faccetta
         if (controllo_normali(n,norm_face)) then
            cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
            q1 = blogrd(pbl,igr)%nodo(ip+1,jp,kp)
            q2 = blogrd(pbl,igr)%nodo(ip,jp+1,kp)
            e1 = q1-cp
            e2 = q2-cp
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+nint(x1)
               jmin = jp+nint(x2)
               kmin = kp
               go to 19
            end if
            cp = blogrd(pbl,igr)%nodo(ip+1,jp+1,kp)
            e1 = q2-cp  !!!
            e2 = q1-cp  !!!
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+1-nint(x1)
               jmin = jp+1-nint(x2)
               kmin = kp
               go to 19
            end if
         end if
      end do
      end do

   else
      jp = nl(3)
      is = 0
      js = 7-2*nf
      ks = 0
      do kp=nl(5),nl(6)-1
      do ip=nl(1),nl(2)-1
         norm_face = real(js,kind(0d0)) * blomet(pbl,igr)%cella(ip+1,jp,kp+1)%snj
         ! se fallisce il controllo sulle normali passa direttamente alla
         ! prossima faccetta
         if (controllo_normali(n,norm_face)) then
            cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
            q1 = blogrd(pbl,igr)%nodo(ip+1,jp,kp)
            q2 = blogrd(pbl,igr)%nodo(ip,jp,kp+1)
            e1 = q1-cp
            e2 = q2-cp
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+nint(x1)
               jmin = jp
               kmin = kp+nint(x2)
               go to 19
            end if
            cp = blogrd(pbl,igr)%nodo(ip+1,jp,kp+1)
            e1 = q2-cp  !!!
            e2 = q1-cp  !!!
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+1-nint(x1)
               jmin = jp
               kmin = kp+1-nint(x2)
               go to 19
            end if
         end if
      end do
      end do
   end if

13 continue
end do

! se arriva qui non ha trovato donatori
nwv = 0
return

19 continue
nwv(1) = pbl
nwv(2) = nf
nwv(3) = imin
nwv(4) = jmin
nwv(5) = kmin
nwv(6) = is
nwv(7) = js
nwv(8) = ks
nwv(9) = pmin
e1 = vers0(e1)
cp = e1.cross.e2
e2 = cp.cross.e1
e2 = vers0(e2)

!9}}}
end subroutine strato_limite_4frame

!===============================================================================
! TRUE se le coordinate controvarianti di c nella base e1b,e2b sono 0<x1,x2<1
! Data la base covariante, quella controvariante è calcolata:
!   e1a = (e2b.cross.e1b).cross.e2b
!   e2a = (e1b.cross.e2b).cross.e1b
! tenendo conto che per 3 generici vettori u1,u2,u3 vale:
!   (u1.cross.u2).cross.u3 = -(u2.dot.u3)*u1 + (u1.dot.u3)*u2
!===============================================================================
function coocontr(ray,e1b,e2b,x1,x2)
!{{{10
use prec
use moddef
implicit none
logical (kind=I1P) :: coocontr
type (point), intent(in) :: ray,e1b,e2b
real (kind=R8P), parameter :: meps = -1.0d-3
real (kind=R8P) :: x1,x2,w
type (point)  :: e1a,e2a

! vettori controvarianti
w   = e1b.dot.e2b
e1a = (e2b.dot.e2b)*e1b - w*e2b
e2a = (e1b.dot.e1b)*e2b - w*e1b

! proiezioni del vettore ray sulla base controvariante
x1 = (ray.dot.e1a)/(e1b.dot.e1a)
x2 = (ray.dot.e2a)/(e2b.dot.e2a)

coocontr = (x1*(1d0-x1).gt.meps) .and. (x2*(1d0-x1-x2).gt.meps)
!10}}}
end function coocontr

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!{{{11
!===============================================================================
! cerca i patch nella stessa famiglia senza controllo sulla priorità
!===============================================================================
subroutine patch_relatives(p,famn)
use prec
use modpar, only : bp,family2patch,famw
implicit none
integer (kind=I4P) :: p,famn
integer (kind=I4P) :: famp,j

! famiglia del patch
famp = mod(bp(4,p),100)

! memorizza in variabili locali i parenti di "p"
j = 1
famn = 0
if (famp.gt.0) then
   do while (family2patch(j,famp).gt.0)
      if (family2patch(j,famp).ne.p) then
         famn = famn+1
         famw(famn) = family2patch(j,famp)
      end if
      j = j+1
   end do
end if

return
end subroutine patch_relatives

!===============================================================================

! inutilizzata, potrebbe servire
function edgepoint(nfm,nl,ipm,jpm,kpm)
use prec
implicit none
logical (kind=I1P) :: edgepoint,cor
integer (kind=I4P) :: nfm,nl(6),ipm,jpm,kpm

cor = .false.

if (nfm.lt.3) then
   cor = jpm.eq.nl(3).or.jpm.eq.nl(4) &
     .or.kpm.eq.nl(5).or.kpm.eq.nl(6)
else if (nfm.gt.4) then
   cor = jpm.eq.nl(3).or.jpm.eq.nl(4) &
     .or.ipm.eq.nl(1).or.ipm.eq.nl(2)
else
   cor = ipm.eq.nl(1).or.ipm.eq.nl(2) &
     .or.kpm.eq.nl(5).or.kpm.eq.nl(6)
end if

edgepoint = cor

return
end

!===============================================================================
! TRUE se il punto "c" cade fuori dall'angolo formato da (q2-q1) e (q3-q1)
! Calcolo aree tramite prodotti vettoriali
!===============================================================================
function fuori_angolo_vett(c,q1,q2,q3)
use prec
use moddef
implicit none
logical (kind=I1P) :: fuori_angolo_vett
type (point), intent(in) :: c,q1,q2,q3
type (point) :: r1,r2,r3
real (kind=R8P) :: d1,d2,rt

fuori_angolo_vett = .true.

r1 = c-q1
r2 = c-q2
r3 = c-q3

d1 = norma(r1.cross.r2) + norma(r2.cross.r3) + norma(r3.cross.r1)
d2 = norma((q1-q2).cross.(q1-q3))
rt = abs(d1-d2)/d2

if (rt.gt.1d-1) return

fuori_angolo_vett = .false.

return
end function fuori_angolo_vett

!===============================================================================
! TRUE se il punto "c" cade fuori dall'angolo formato dal triangolo q1-q2-q3
! (usa teorema Carnot)
!===============================================================================
function fuori_angolo_carnot(c,q1,q2,q3)
use prec
use moddef
implicit none
logical (kind=I1P) :: fuori_angolo_carnot
type (point), intent(in) :: c,q1,q2,q3
real (kind=R8P) :: r1,r2,r3, s1,s2,s3, pp
real (kind=R8P) :: d1,d2,rt

fuori_angolo_carnot = .true.

r1 = norma(c-q1)
r2 = norma(c-q2)
r3 = norma(c-q3)

s1 = norma(q2-q3)
s2 = norma(q1-q3)
s3 = norma(q1-q2)

pp = 0.5d0*(s1+s2+s3)
d1 = sqrt(pp*(pp-s1)*(pp-s2)*(pp-s3))

pp = 0.5d0*(r1+r2+s3)
d2 =      sqrt(pp*(pp-r1)*(pp-r2)*(pp-s3))
pp = 0.5d0*(r1+s2+r3)
d2 = d2 + sqrt(pp*(pp-r1)*(pp-s2)*(pp-r3))
pp = 0.5d0*(s1+r2+r3)
d2 = d2 + sqrt(pp*(pp-s1)*(pp-r2)*(pp-r3))

rt = abs(d1-d2)/(0.5d0*(d1+d2))

if (rt.gt.1d-1) return

fuori_angolo_carnot = .false.

return
end function fuori_angolo_carnot

!11}}}
!===============================================================================
! funzione per valutare la distanza massima del donatore di una cella di s.l.
!===============================================================================
function controllo_distanza(a,b)
use prec
implicit none
logical (kind=I1P) :: controllo_distanza
real (kind=R8P), intent(in) :: a,b
controllo_distanza = a .lt. b
end function

function controllo_normali(a,b)
use prec
use moddef, only : point,vers0,operator(.dot.)
implicit none
logical (kind=I1P) :: controllo_normali
type (point), intent(in) :: a,b
controllo_normali = (vers0(a).dot.vers0(b)) .gt. 0.7d0
end function
