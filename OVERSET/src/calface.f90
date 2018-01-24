! vim: set fdm=marker:
!===============================================================================
!  assegnazione c.c. sulle facce {{{1
!===============================================================================
subroutine calface(igr)
use prec
use modpar, only : ndeb, nccnat,offgen,offchi,estrap, &
                   npa,bp,patch_flags
integer (kind=I4P), intent(in) :: igr
integer (kind=I4P) :: tgr,ccf,p,q
integer (kind=I4P) :: bpp(10),bqq(10)
logical (kind=I1P) :: pflag(4)

!...............................................................................
!  stampa di controllo
!
write(*,'(a,i3,a)') 'inizio calface (igr =',igr,')'
tgr = 2**(igr-1)

!===============================================================================
!$omp parallel &
!$omp    default(none) &
!$omp    shared(igr,tgr,npa,bp,patch_flags) &
!$omp    private(p,q,ccf,bpp,bqq,pflag)
!===============================================================================

!$omp do

!  inizio loop su tutti i "Boundary Patches"
!
do p=1,npa

   bpp = bp(1:10,p)
   pflag = patch_flags(1:4,p)
   ccf = bpp(3)

   !............................................................................
   ! c.c. naturali
   !
   if (ccf.le.nccnat) then

      call naturali(igr,tgr,bpp)

   !............................................................................
   ! c.c. adiacenza uno a uno
   !
   else if (ccf.gt.100) then

      q = bpp(4)
      bqq = bp(1:10,q)
      call adiacen(igr,tgr,bpp,bqq)

   !............................................................................
   ! c.c. chimera con donatori per il centro cella
   !
   else if (ccf.ge.offgen .and. ccf.le.(offgen+nccnat)) then

      call chimcen(igr,tgr,bpp,pflag)

   !............................................................................
   ! c.c. chimera con donatori per centro faccia
   !
   else if (ccf.ge.offchi .and. ccf.le.(offchi+nccnat)) then

      call chimfac(igr,tgr,bpp,pflag)

   !............................................................................
   ! c.c. di tipo "estrap" (donatori per centro faccia oppure chimera con
   ! peso = -19.0)
   !
   else if (ccf.eq.estrap) then

      call chiestr(igr,tgr,bpp,pflag)

   !............................................................................
   ! c.c. sconosciuta (stop)
   !
   else

      write(*,*)    'Errore:',ccf,' c.c. sconosciuta in calface'
      write(ndeb,*) 'Errore:',ccf,' c.c. sconosciuta in calface'
      stop

   end if

!...............................................................................
!  fine loop sulle toppe
!
end do

!$omp end parallel

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
! subroutine nadia(l,m,n,nn,ia,ja,ka,nla) {{{1
subroutine nadia(l,m,n,nn,ia,ja,ka,nla)
use prec
implicit none
integer (kind=I4P) :: l,m,n,nn,ia,ja,ka,nla(6)

if      (l.eq.1) then
   ia = n - nn + m*(nn+nla(2))
else if (l.eq.2) then
   ia = 1 - n  + m*(nn+nla(2))
else if (l.eq.3) then
   ja = n - nn + m*(nn+nla(4))
else if (l.eq.4) then
   ja = 1 - n  + m*(nn+nla(4))
else if (l.eq.5) then
   ka = n - nn + m*(nn+nla(6))
else if (l.eq.6) then
   ka = 1 - n  + m*(nn+nla(6))
end if

end subroutine

!===============================================================================

! subroutine tadia(l,n,ia,ja,ka,nla) {{{1
subroutine tadia(l,n,ia,ja,ka,nla)
use prec
implicit none
integer (kind=I4P) :: l,n,ia,ja,ka,nla(6)

if      (l.eq.1) then
   ia = nla(1)+n
else if (l.eq.2) then
   ia = nla(2)+1-n
else if (l.eq.3) then
   ja = nla(3)+n
else if (l.eq.4) then
   ja = nla(4)+1-n
else if (l.eq.5) then
   ka = nla(5)+n
else if (l.eq.6) then
   ka = nla(6)+1-n
end if

end subroutine

!===============================================================================
! subroutine naturali(igr,tgr,bp) {{{1
subroutine naturali(igr,tgr,bp)
use prec
use modpar, only : blochi
implicit none
integer (kind=I4P), intent(in) :: igr,tgr,bp(10)
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: ibl,nf,ccf,nl(6),is,js,ks
integer (kind=I4P) :: i,j,k

ibl = bp(1)
nf  = bp(2)
ccf = bp(3)
nl  = bp(5:10)/tgr
is  = 0
js  = 0
ks  = 0

nii = blochi(ibl,igr)%ni
njj = blochi(ibl,igr)%nj
nkk = blochi(ibl,igr)%nk

! bp contiene gli indici dei vertici ma il loop per assegnare le con.cont. va
! eseguito sui centri faccia
nl( (/1,3,5/) ) = nl( (/1,3,5/) ) + 1

if (nf.eq.1) then
   nl(1) = 0
   nl(2) = 0
   is = -1
else if (nf.eq.2) then
   nl(1) = nii+1
   nl(2) = nii+1
   is = 1
else if (nf.eq.3) then
   nl(3) = 0
   nl(4) = 0
   js = -1
else if (nf.eq.4) then
   nl(3) = njj+1
   nl(4) = njj+1
   js = 1
else if (nf.eq.5) then
   nl(5) = 0
   nl(6) = 0
   ks = -1
else if (nf.eq.6) then
   nl(5) = nkk+1
   nl(6) = nkk+1
   ks = 1
end if

do k=nl(5),nl(6)
   do j=nl(3),nl(4)
      do i=nl(1),nl(2)
         blochi(ibl,igr)%cella(i,j,k)%t = ccf
         blochi(ibl,igr)%cella(i+is,j+js,k+ks)%t = ccf
      end do
   end do
end do

end subroutine

!===============================================================================
! subroutine adiacen(igr,tgr,bp,bq) {{{1
subroutine adiacen(igr,tgr,bp,bq)
use prec
use moddef, only : donor
use modpar, only : offbiu,donmax,blochi
implicit none
integer (kind=I4P), intent(in) :: igr,tgr,bp(10),bq(10)
integer (kind=I4P) :: ibl,nf,ccf
integer (kind=I4P) :: nii,njj,nkk,gcc(6)
integer (kind=I4P), dimension(6) :: nl,nlp,nlq
integer (kind=I4P) :: i,j,k,li,lj,lk,df,pf
integer (kind=I4P) :: td,nd,bd,id,jd,kd
logical (kind=I1P) :: il,jl,kl
type (donor) :: md(donmax)

ibl = bp(1)
nf  = bp(2)
ccf = bp(3)
nl  = bp(5:10)/tgr

nii = blochi(ibl,igr)%ni
njj = blochi(ibl,igr)%nj
nkk = blochi(ibl,igr)%nk
gcc = blochi(ibl,igr)%gc

! bp contiene gli indici dei vertici ma il loop per assegnare le con.cont. va
! eseguito sui centri faccia
nl( (/1,3,5/) ) = nl( (/1,3,5/) ) + 1

! estende i limiti del patch per includere le cornici esterne
if (nl(1).eq.1) nl(1) = -1
if (nl(2).eq.nii) nl(2) = nl(2)+2
if (nl(3).eq.1) nl(3) = -1
if (nl(4).eq.njj) nl(4) = nl(4)+2
if (nl(5).eq.1) nl(5) = -1
if (nl(6).eq.nkk) nl(6) = nl(6)+2

il  = .false.
jl  = .false.
kl  = .false.

df = mod(nf,2)
pf = 1-df
if (nf.lt.3) then
   nl(1) = pf*(nii+1) + df*(-gcc(1)+1)
   nl(2) = pf*(nii+gcc(2))
   il = .true.
else if (nf.gt.4) then
   nl(5) = pf*(nkk+1) + df*(-gcc(5)+1)
   nl(6) = pf*(nkk+gcc(6))
   kl = .true.
else
   nl(3) = pf*(njj+1) + df*(-gcc(3)+1)
   nl(4) = pf*(njj+gcc(4))
   jl = .true.
end if

nlp = bp(5:10)/tgr
nlq = bq(5:10)/tgr
li  =     ccf      / 100
lj  = mod(ccf,100) / 10
lk  = mod(ccf,10)
td  = offbiu+nf
nd  = 1
bd  = bq(1)

! direzione "I"
!
if (il) then
   pf = 1-mod(li+df,2)
   do k=nl(5),nl(6)
      call tadia(lk,k-nlp(5),id,jd,kd,nlq)
      do j=nl(3),nl(4)
         call tadia(lj,j-nlp(3),id,jd,kd,nlq)
         do i=nl(1),nl(2)
            call nadia(li,pf,i,nlp(2),id,jd,kd,nlq)
            md(1) = donor(bd,id,jd,kd,1.0)
            call newchi(igr,ibl,i,j,k,td,nd,md)
         end do
      end do
   end do

! direzione "J"
!
else if (jl) then
   pf = 1-mod(lj+df,2)
   do k=nl(5),nl(6)
      call tadia(lk,k-nlp(5),id,jd,kd,nlq)
      do j=nl(3),nl(4)
         call nadia(lj,pf,j,nlp(4),id,jd,kd,nlq)
         do i=nl(1),nl(2)
            call tadia(li,i-nlp(1),id,jd,kd,nlq)
            md(1) = donor(bd,id,jd,kd,1.0)
            call newchi(igr,ibl,i,j,k,td,nd,md)
         end do
      end do
   end do

! direzione "K"
!
else if (kl) then
   pf = 1-mod(lk+df,2)
   do k=nl(5),nl(6)
      call nadia(lk,pf,k,nlp(6),id,jd,kd,nlq)
      do j=nl(3),nl(4)
         call tadia(lj,j-nlp(3),id,jd,kd,nlq)
         do i=nl(1),nl(2)
            call tadia(li,i-nlp(1),id,jd,kd,nlq)
            md(1) = donor(bd,id,jd,kd,1.0)
            call newchi(igr,ibl,i,j,k,td,nd,md)
         end do
      end do
   end do

end if

end subroutine

!===============================================================================
! subroutine chimcen(igr,tgr,bp,pflag) {{{1
subroutine chimcen(igr,tgr,bp,pflag)
use prec
use moddef
use modpar, only : offgen,donmax,lvmax,blochi,blomet
implicit none
integer (kind=I4P), intent(in) :: igr,tgr,bp(10)
logical (kind=I1P), intent(in) :: pflag(4)
integer (kind=I4P) :: nii,njj,nkk,gcc(6)
integer (kind=I4P) :: ibl,nf,ccf
integer (kind=I4P), dimension(6) :: nl
integer (kind=I4P) :: i,j,k,df,pf
integer (kind=I4P) :: td,nd,bd
type (donor) :: md(donmax)

integer (kind=I4P) :: flag(10)
integer (kind=I1P) :: lvbad(0:lvmax)
integer (kind=I4P) :: q,xd
real (kind=R8P) :: d
type (point)     :: u

ibl = bp(1)
nf  = bp(2)
ccf = bp(3)
nl  = bp(5:10)/tgr

nii = blochi(ibl,igr)%ni
njj = blochi(ibl,igr)%nj
nkk = blochi(ibl,igr)%nk
gcc = blochi(ibl,igr)%gc

! bp contiene gli indici dei vertici ma il loop per assegnare le con.cont. va
! eseguito sui centri faccia
nl( (/1,3,5/) ) = nl( (/1,3,5/) ) + 1

df = mod(nf,2)
pf = 1-df
if (nf.lt.3) then
   nl(1) = pf*(nii+1) + df*(-gcc(1)+1)
   nl(2) = pf*(nii+gcc(2))
else if (nf.gt.4) then
   nl(5) = pf*(nkk+1) + df*(-gcc(5)+1)
   nl(6) = pf*(nkk+gcc(6))
else
   nl(3) = pf*(njj+1) + df*(-gcc(3)+1)
   nl(4) = pf*(njj+gcc(4))
end if

! tipo di c.c. così come verrà letto in Xnavis
bd = ibl
td = offgen+nf
q  = 1
if (ccf.gt.offgen .and. ccf.ne.(offgen+19)) q = 0
nd = 8
do k=nl(5),nl(6)
   do j=nl(3),nl(4)
      do i=nl(1),nl(2)

         ! ricerca donatori per il centro cella
         flag = (/ 1,1,1,q,q,q, 0,1,q,0 /)
         u = blomet(ibl,igr)%cella(i,j,k)%cen
         d = blochi(ibl,igr)%cella(i,j,k)%size
         call findonor(igr,ibl,d,flag,u,md,lvbad)

         if (flag(10).gt.1) then
            ! ha trovato un donatore
            call newchi(igr,ibl,i,j,k,td,nd,md)

         else if (ccf.eq.offgen) then
            ! non ha trovato un donatore e non ha un c.c. alternativa
            ! usa come donatore il punto più vicino di un altro blocco
            if (flag(10).ne.1) then
               write(*,*) 'Subroutine calcorn: pensavo di utilizzare ',&
               'come donatore il "punto più vicino"'
               write(*,*) 'ma non ricordo di averlo mai cercato !?!?'
               stop
            end if
            call newchi(igr,ibl,i,j,k,td,nd,md)

         else
            ! non ha trovato un donatore ma ha una c.c. alternativa
            blochi(ibl,igr)%cella(i,j,k)%t = ccf-offgen
         end if

      end do
   end do
end do

! ridefinisce, se può, le celle di cornice che si trovano nello s.l.
! SOLO QUELLE CHE NON HANNO UNA C.C. ALTERNATIVA (modifica del 22/05/2007
! necessaria per il timone della KVLCC suggerita da Andrea)
if (ccf.eq.offgen .and. pflag(4)) then
   xd = 0
   call strato_limite_cornice(igr,ibl,nf,nl,td,xd)
end if

end subroutine

!===============================================================================
! subroutine chimfac(igr,tgr,bp,pflag) {{{1
subroutine chimfac(igr,tgr,bp,pflag)
use prec
use moddef
use modpar, only : offchi,donmax,lvmax,blochi
implicit none
integer (kind=I4P), intent(in) :: igr,tgr,bp(10)
logical (kind=I1P), intent(in) :: pflag(4)
integer (kind=I4P) :: ibl,nf,ccf,pf
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P), dimension(6) :: nl
integer (kind=I4P) :: i,j,k,is,js,ks
integer (kind=I4P) :: td,nd,bd
type (donor) :: md(donmax)

integer (kind=I4P) :: flag(10)
integer (kind=I1P) :: lvbad(0:lvmax)
integer (kind=I4P) :: q,xd
real (kind=R8P) :: d
type (point)     :: u

ibl = bp(1)
nf  = bp(2)
ccf = bp(3)
nl  = bp(5:10)/tgr

nii = blochi(ibl,igr)%ni
njj = blochi(ibl,igr)%nj
nkk = blochi(ibl,igr)%nk

is  = 0
js  = 0
ks  = 0

! bp contiene gli indici dei vertici ma il loop per assegnare le con.cont. va
! eseguito sui centri faccia
nl( (/1,3,5/) ) = nl( (/1,3,5/) ) + 1

pf = 1-mod(nf,2)
if (nf.lt.3) then
   nl(1) = pf*(nii+1)
   nl(2) = nl(1)
   is = -1+2*pf
else if (nf.gt.4) then
   nl(5) = pf*(nkk+1)
   nl(6) = nl(5)
   ks = -1+2*pf
else
   nl(3) = pf*(njj+1)
   nl(4) = nl(3)
   js = -1+2*pf
end if

! tipo di c.c. così come verrà letto in Xnavis
bd = ibl
td = offchi+nf
q  = 1
if (ccf.gt.offchi .and. ccf.ne.(offchi+19)) q = 0
nd = 8
do k=nl(5),nl(6)
   do j=nl(3),nl(4)
      do i=nl(1),nl(2)

         ! ricerca donatori per il centro faccia
         flag = (/ 1,1,1,q,q,q, nf,1,q,0 /)
         d = blochi(ibl,igr)%cella(i,j,k)%size
         call centro_faccia_cornice(igr,ibl,i,j,k,nf,u)
         call findonor(igr,ibl,d,flag,u,md,lvbad)

         if (flag(10).gt.1) then
            ! ha trovato un donatore
            call newchi(igr,ibl,i,j,k,td,nd,md)

         else if (ccf.eq.offchi) then
            ! non ha trovato un donatore e non ha un c.c. alternativa
            ! usa come donatore il punto più vicino di un altro blocco
            if (flag(10).ne.1) then
               write(*,*) 'Subroutine chimfac: pensavo di utilizzare ',&
               'come donatore il "punto più vicino"'
               write(*,*) 'ma non ricordo di averlo mai cercato !?!?'
               stop
            end if
            call newchi(igr,ibl,i,j,k,td,nd,md)

         else
            ! non ha trovato un donatore ma ha una c.c. alternativa
            blochi(ibl,igr)%cella(i,j,k)%t = ccf-offchi

         end if

      end do
   end do
end do

! ridefinisce, se può, le celle di cornice che si trovano nello s.l.
! e non hanno una c.c. alternativa
if (ccf.eq.offchi .and. pflag(4)) then
   xd = nf
   call strato_limite_cornice(igr,ibl,nf,nl,td,xd)
end if

! duplica sulla corona più esterna (-1...n+2) le c.c. di quella più
! interna (0...n+1)
do k=nl(5),nl(6)
   do j=nl(3),nl(4)
      do i=nl(1),nl(2)
         call copychi(igr,ibl,i,j,k,i+is,j+js,k+ks)
      end do
   end do
end do

end subroutine

!===============================================================================

! subroutine chiestr(igr,tgr,bp,pflag) {{{1
subroutine chiestr(igr,tgr,bp,pflag)
use prec
use moddef
use modpar, only : offchi,donmax,lvmax,blochi
implicit none
integer (kind=I4P), intent(in) :: igr,tgr,bp(10)
logical (kind=I1P), intent(in) :: pflag(4)
integer (kind=I4P) :: ibl,nf,ccf,pf
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P), dimension(6) :: nl
integer (kind=I4P) :: i,j,k,is,js,ks
integer (kind=I4P) :: td,nd,bd
type (donor) :: md(donmax)

integer (kind=I4P) :: flag(10)
integer (kind=I1P) :: lvbad(0:lvmax)
integer (kind=I4P) :: xd
real (kind=R8P) :: d
type (point)     :: u

ibl = bp(1)
nf  = bp(2)
ccf = bp(3)
nl  = bp(5:10)/tgr

nii = blochi(ibl,igr)%ni
njj = blochi(ibl,igr)%nj
nkk = blochi(ibl,igr)%nk

is  = 0
js  = 0
ks  = 0

! bp contiene gli indici dei vertici ma il loop per assegnare le con.cont. va
! eseguito sui centri faccia
nl( (/1,3,5/) ) = nl( (/1,3,5/) ) + 1

pf = 1-mod(nf,2)
if (nf.lt.3) then
   nl(1) = pf*(nii+1)
   nl(2) = nl(1)
   is = -1+2*pf
else if (nf.gt.4) then
   nl(5) = pf*(nkk+1)
   nl(6) = nl(5)
   ks = -1+2*pf
else
   nl(3) = pf*(njj+1)
   nl(4) = nl(3)
   js = -1+2*pf
end if

! tipo di c.c. così come verrà letto in Xnavis
bd = ibl
td = offchi+nf
nd = 8
do k=nl(5),nl(6)
   do j=nl(3),nl(4)
      do i=nl(1),nl(2)

         ! ricerca donatori per il centro faccia
         flag = (/ 1,1,1,1,1,1, nf,1,0,0 /)
         d = blochi(ibl,igr)%cella(i,j,k)%size
         call centro_faccia_cornice(igr,ibl,i,j,k,nf,u)
         call findonor(igr,ibl,d,flag,u,md,lvbad)

         if (flag(10).gt.1) then
            ! ha trovato un donatore
            call newchi(igr,ibl,i,j,k,td,nd,md)

         else
            ! non ha trovato un donatore
            ! la definisce come cella chimera con donatrice se stessa e
            ! peso -19
            md(1) = donor(ibl,i,j,k,-19_R8P)
            call newchi(igr,ibl,i,j,k,td,1,md)

         end if

      end do
   end do
end do

! ridefinisce, se può, le celle di cornice che si trovano nello s.l.
if (pflag(4)) then
   xd = nf
   call strato_limite_cornice(igr,ibl,nf,nl,td,xd)
end if

! duplica sulla corona più esterna (-1...n+2) le c.c. di quella più
! interna (0...n+1)
do k=nl(5),nl(6)
   do j=nl(3),nl(4)
      do i=nl(1),nl(2)
         call copychi(igr,ibl,i,j,k,i+is,j+js,k+ks)
      end do
   end do
end do

end subroutine
