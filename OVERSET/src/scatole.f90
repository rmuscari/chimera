!===============================================================================
subroutine alloca_scatole
use modpar
implicit none
allocate(lbox(nbox)) !.........................livello
allocate(pbox(nbox)) !.........................priorità
allocate(tbox(nbox)) !.........................tipo
allocate(gbox(nbox)) !.........................gruppo di appartenenza
allocate(ilex(nbox),jlex(nbox),klex(nbox)) !...dimensioni
allocate(obox(nbox)) !.........................centro della scatola
allocate(ibox(nbox),jbox(nbox),kbox(nbox)) !...versori (per orientamento)
allocate(vbox(8,nbox)) !.......................vertici esaedro
allocate(sbox(2,6,nbox)) !.....................normali facce esaedro
allocate(fbox(6,nbox)) !.......................centri facce esaedro
return
end

!===============================================================================

subroutine esaedro(q)
use prec
use moddef
use modpar, only : obox,fbox,vbox,sbox
implicit none
integer (kind=I4P) :: q

! centro scatola
obox(q) = mean(8,vbox(1:8,q))

! centri facce
fbox(1,q) = 0.25d0*(vbox(1,q)+vbox(7,q)+vbox(3,q)+vbox(5,q))
fbox(2,q) = 0.25d0*(vbox(2,q)+vbox(8,q)+vbox(4,q)+vbox(6,q))
fbox(3,q) = 0.25d0*(vbox(1,q)+vbox(6,q)+vbox(2,q)+vbox(5,q))
fbox(4,q) = 0.25d0*(vbox(3,q)+vbox(8,q)+vbox(4,q)+vbox(7,q))
fbox(5,q) = 0.25d0*(vbox(1,q)+vbox(4,q)+vbox(2,q)+vbox(3,q))
fbox(6,q) = 0.25d0*(vbox(5,q)+vbox(8,q)+vbox(6,q)+vbox(7,q))

! calcola le (due) normali per ciascuna faccia e ricalcola il "centro faccia"
! come punto di mezzo della diagonale
call calc_sbox(3,1,7,5,vbox(1,q),fbox(2,q),fbox(1,q),sbox(1,1,q))
call calc_sbox(2,4,6,8,vbox(1,q),fbox(1,q),fbox(2,q),sbox(1,2,q))
call calc_sbox(1,2,5,6,vbox(1,q),fbox(4,q),fbox(3,q),sbox(1,3,q))
call calc_sbox(4,3,8,7,vbox(1,q),fbox(3,q),fbox(4,q),sbox(1,4,q))
call calc_sbox(3,4,1,2,vbox(1,q),fbox(6,q),fbox(5,q),sbox(1,5,q))
call calc_sbox(5,6,7,8,vbox(1,q),fbox(5,q),fbox(6,q),sbox(1,6,q))

return
end

!===============================================================================
! Calcola le normali alle facce di una scatola
! Poiché le facce possono essere sghembe, divide ciascuna faccia lungo una
! diagonale in due triangoli e calcola le due normali
! Nella scelta della diagonale si assicura che i due triangoli risultanti
! formino un angolo < 180
!===============================================================================
subroutine calc_sbox(i1,i2,i3,i4,v,fop,fc,s)
use prec
use moddef
implicit none
real (kind=R8P), parameter :: neps=-1d-9
integer (kind=I4P), intent(in) :: i1,i2,i3,i4
type (point) :: v(8),fop,fc,s(2)
type (point) :: va,vb,vc
real (kind=R8P) :: chk

vc = 0.5d0*(v(i2)+v(i3))

va = (v(i3)-v(i2)).cross.(v(i1)-v(i2))
if ((va.dot.(vc-fop)).lt.0d0) va = -va
va = vers0(va)

chk = va.dot.vers0(vc-v(i4))
if (chk.gt.neps) then
   ! tutto bene, ho scelto la diagonale giusta
   ! definisco la seconda normale e torno
   vb = (v(i4)-v(i2)).cross.(v(i3)-v(i2))
   if ((vb.dot.(vc-fop)).lt.0d0) vb = -vb
   vb = vers0(vb)

   s(1) = va
   s(2) = vb
   fc   = vc
   return
end if

! ho scelto la diagonale sbagliata
!  write(*,'(a,e20.10)') 'chk1  = ',chk
!  write(*,'(a,3(e20.10))') 'va    = ',va
!  write(*,'(a,3(e20.10))') 'vc-v4 = ',vc-v(i4)

vc = 0.5d0*(v(i1)+v(i4))

va = (v(i4)-v(i1)).cross.(v(i3)-v(i1))
if ((va.dot.(vc-fop)).lt.0d0) va = -va
va = vers0(va)

chk = va.dot.vers0(vc-v(i2))
if (chk.gt.neps) then
   ! buona la seconda
   ! definisco la seconda normale e torno
   vb = (v(i2)-v(i1)).cross.(v(i4)-v(i1))
   if ((vb.dot.(vc-fop)).lt.0d0) vb = -vb
   vb = vers0(vb)

   s(1) = va
   s(2) = vb
   fc   = vc
   return
end if

! non va bene neanche la seconda diagonale ... c'è qualcosa che non va
!  write(*,'(a,e20.10)') 'chk2  = ',chk
!  write(*,'(a,3(e20.10))') 'va    = ',va
!  write(*,'(a,3(e20.10))') 'vc-v2 = ',vc-v(i2)
stop 'Errore in calc_sbox'

return
end

!===============================================================================
! controlla se il punto cade all'interno di una scatola appartenente ad un
! livello diverso dal proprio
!===============================================================================
subroutine check_parete_interna(ans,lvlblo,pryblo,lvbad,c)
use prec
use moddef
use modpar
implicit none
real (kind=R8P), parameter :: zero = 1d-10
integer (kind=I4P), intent(in) :: lvlblo,pryblo
integer (kind=I1P), intent(in) :: lvbad(0:lvmax)
integer (kind=I4P) :: ans,n,p
type (point) :: c,r
real (kind=R8P) :: s

ans = -1
do n=1,nbox

!1111111111111111111111111111111111111111111111111111111111111111111111111111111
!   VERSIONE 1
!  ! priorità scatola < priorità blocco
!  if (pbox(n).lt.pryblo) go to 11

!  ! priorità scatola = priorità blocco
!  if (pbox(n).eq.pryblo) then
!     if (lbox(n).eq.lvlblo .or. lvbad(lbox(n)).lt.0) go to 11
!  end if
!1111111111111111111111111111111111111111111111111111111111111111111111111111111

!2222222222222222222222222222222222222222222222222222222222222222222222222222222
!   VERSIONE 2
   ! scatola su un livello con donatori (accettabili o no)
   if (lvbad(lbox(n)).lt.0_I1P) go to 11
   ! livello scatola = livello blocco && priorità scatola <= priorità blocco
   if (lbox(n).eq.lvlblo .and. pbox(n).le.pryblo) go to 11
!2222222222222222222222222222222222222222222222222222222222222222222222222222222

   if (tbox(n).eq.1) then
      ! PARALLELEPIPEDO
      r = c-obox(n)
      s = r.dot.ibox(n)
      if (s.lt.(-ilex(n)) .or. s.gt.ilex(n)) go to 11
      s = r.dot.jbox(n)
      if (s.lt.(-jlex(n)) .or. s.gt.jlex(n)) go to 11
      s = r.dot.kbox(n)
      if (s.lt.(-klex(n)) .or. s.gt.klex(n)) go to 11
   else if (tbox(n).eq.2) then
      ! CILINDRO (OBSOLETO)
      r = c-obox(n)
      s = r.dot.ibox(n)
      if (s.lt.(-ilex(n)) .or. s.gt.ilex(n)) go to 11
      s = sqrt((r.dot.jbox(n))**2 + (r.dot.kbox(n))**2)
      if (s.gt.jlex(n)) go to 11
   else if (tbox(n).eq.3) then
      ! CILINDRI COASSIALI
      r = c-obox(n)
      s = r.dot.ibox(n)
      if (s.lt.(-ilex(n)) .or. s.gt.ilex(n)) go to 11
      s = sqrt((r.dot.jbox(n))**2 + (r.dot.kbox(n))**2)
      if ( (s.lt.jlex(n)) .or. (s.gt.klex(n)) ) go to 11
   else if (tbox(n).eq.9) then
      ! ESAEDRO
      do p=1,6
         r = c-fbox(p,n)
         s = r.dot.sbox(1,p,n)
         if (s.gt.zero) go to 11
         s = r.dot.sbox(2,p,n)
         if (s.gt.zero) go to 11
      end do
   end if

   ans = gbox(n)
   return

11 continue
end do

return
end
