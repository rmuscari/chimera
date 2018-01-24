!===============================================================================
!  Calcolo dello spessore di strato limite per tutte le famiglie
!===============================================================================
subroutine calcolo_spessore_strato_limite
use prec
use moddef, only : point,dist
use modpar
implicit none
integer (kind=I4P) :: p,q,ibl,fam,nf,ip,jp,kp,mp
real (kind=R8P) :: diag_max,spes_min
type (point) :: c1,c2,c3,c4,c5

!...............................................................................

! alloca i vettori con gli spessori di strato limite
! boulap spessore caratteristico di ciascun patch di parete
! boulab (serve nella subroutine intpoint) spessore minimo tra i patch
!        di parete del blocco
allocate(boulap(npa))
allocate(boulaf(nfa))
allocate(boulab(nbl))

boulap = 0d0
boulaf = 0d0
boulab = 0d0
if (boulay.lt.1d-30) return

boulap = boulay
boulaf = boulay
boulab = boulay

! loop su tutti i patch
do p=1,npa

   if (.not.patch_flags(1,p) .and. .not.patch_flags(2,p)) cycle

   ibl = bp(1,p)
   nf  = bp(2,p)
   fam = mod(bp(4,p),100)
   diag_max = -1d0
   spes_min =  1d20

   ! NF = 1 o 2
   if (nf.lt.3) then
      ip = bp(5,p)
      mp = (2-nf)*blogrd(ibl,1)%ni
      do kp=bp(9,p),bp(10,p)-1
      do jp=bp(7,p),bp(8,p)-1
         call getvertex(1,ibl,ip,jp  ,kp  ,c1)
         call getvertex(1,ibl,ip,jp+1,kp  ,c2)
         call getvertex(1,ibl,ip,jp  ,kp+1,c3)
         call getvertex(1,ibl,ip,jp+1,kp+1,c4)
         call getvertex(1,ibl,mp,jp  ,kp  ,c5)
         diag_max = max(diag_max,dist(c1,c4),dist(c2,c3))
         spes_min = min(spes_min,dist(c1,c5))
      end do
      end do

   ! NF = 5 o 6
   else if (nf.gt.4) then
      kp = bp(9,p)
      mp = (6-nf)*blogrd(ibl,1)%nk
      do jp=bp(7,p),bp(8,p)-1
      do ip=bp(5,p),bp(6,p)-1
         call getvertex(1,ibl,ip  ,jp  ,kp,c1)
         call getvertex(1,ibl,ip  ,jp+1,kp,c2)
         call getvertex(1,ibl,ip+1,jp  ,kp,c3)
         call getvertex(1,ibl,ip+1,jp+1,kp,c4)
         call getvertex(1,ibl,ip  ,jp  ,mp,c5)
         diag_max = max(diag_max,dist(c1,c4),dist(c2,c3))
         spes_min = min(spes_min,dist(c1,c5))
      end do
      end do

   ! NF = 3 o 4
   else
      jp = bp(7,p)
      mp = (4-nf)*blogrd(ibl,1)%nj
      do kp=bp(9,p),bp(10,p)-1
      do ip=bp(5,p),bp(6,p)-1
         call getvertex(1,ibl,ip  ,jp,kp  ,c1)
         call getvertex(1,ibl,ip  ,jp,kp+1,c2)
         call getvertex(1,ibl,ip+1,jp,kp  ,c3)
         call getvertex(1,ibl,ip+1,jp,kp+1,c4)
         call getvertex(1,ibl,ip  ,mp,kp  ,c5)
         diag_max = max(diag_max,dist(c1,c4),dist(c2,c3))
         spes_min = min(spes_min,dist(c1,c5))
      end do
      end do

   end if

   boulap(p) = min(4d0*diag_max,0.9d0*spes_min,boulay)
   boulab(ibl) = min(boulab(ibl),boulap(p))
   ! se il patch appartiene ad una famiglia definisce
   ! boulaf come il minimo tra i patch che appartengono a quella famiglia e
   if (fam.gt.0) then
      boulaf(fam) = min(boulaf(fam),boulap(p))
   end if

end do

! ritocca fare un loop sui patch per definire
! boulab come il minimo tra le famiglie a cui appartengono i propri patch
do p=1,npa
   if (patch_flags(1,p) .or. patch_flags(2,p)) then
      ibl = bp(1,p)
      fam = mod(bp(4,p),100)
      if (fam.gt.0) then
         boulab(ibl) = min(boulab(ibl),boulaf(fam))
      end if
   end if
end do

!  stampa debug
write(ndeb,'(/,a)') 'Spessori di strato limite per famiglie {{{1'
do q=1,nfa
   if (family2patch(1,q).gt.0) then
      write(ndeb,'(a,i2,5x,a,e13.6)') 'Famiglia: ',q,&
                                      'spessore strato limite:',sqrt(boulaf(q))
   end if
end do

write(ndeb,'(/,a)') 'Spessori di strato limite per blocchi {{{1'
write(ndeb,'(a,e13.6)') 'Default:',sqrt(boulay)
do ibl=1,nbl
   if (boulab(ibl).lt.boulay) then
      write(ndeb,'(a,i4,5x,a,e13.6)') 'Blocco: ',ibl,&
                                     'spessore strato limite:',sqrt(boulab(ibl))
   end if
end do
write(ndeb,*)

end subroutine calcolo_spessore_strato_limite
