!===============================================================================
!  Calcolo della distanza minima dei centri cella da un patch di parete
!  (appartenente allo stesso blocco)
!===============================================================================
subroutine waldis(igr)
use prec
use moddef, only : point,dist
use modpar
implicit none
integer (kind=I4P) :: igr,ibl,nf,ip,jp,kp,ib,jb,kb,is,js,ks,im,jm,km
integer (kind=I4P) :: ipp,jpp,kpp,imin,jmin,kmin,nii,njj,nkk
integer (kind=I4P) :: i,j,k, ii,jj,kk
integer (kind=I4P) :: p
integer (kind=I4P) :: tgr,nl(6),ql(6,4)
real (kind=R8P)    :: dm
type (point)     :: cb,cp

!...............................................................................
!  stampa di controllo
write(*,'(a,i3,a)') 'inizio waldis  (igr =',igr,')'

! inizializzazione
do ibl=1,nbl
   blochi(ibl,igr)%cella(:,:,:)%dist = 1.0d20
   blochi(ibl,igr)%cella(:,:,:)%wall = 0
   blochi(ibl,igr)%cella(:,:,:)%boul = .false.
end do

! se c'è un solo blocco esce subito
if (nbl.eq.1 .or. boulay.lt.1d-30) return

! calcolo distanze nei vertici
tgr = 2**(igr-1)

! loop sui patches
do p=1,npa

   ! se non è un patch di parete salta al successivo
   if (.not.patch_flags(1,p) .and. .not.patch_flags(2,p)) cycle

   ibl = bp(1,p)
   nf  = bp(2,p)
   blt = boulap(p)
   nl  = bp(5:10,p)/tgr
   is  = 0
   js  = 0
   ks  = 0
   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   ! Calcolo distanza per i vertici "vicini" alla parete
   ! parte dai vertici sul patch e si muove lungo la linea coordinata che esce
   ! dal patch

   ! facce 1 e 2
   if (nf.lt.3) then
      ql(1:6,1) = (/ nl(nf),nl(nf), nl(3),nl(4), nl(5),nl(5) /)
      ql(1:6,2) = (/ nl(nf),nl(nf), nl(3),nl(4), nl(6),nl(6) /)
      ql(1:6,3) = (/ nl(nf),nl(nf), nl(3),nl(3), nl(5),nl(6) /)
      ql(1:6,4) = (/ nl(nf),nl(nf), nl(4),nl(4), nl(5),nl(6) /)
      is = 3-2*nf
      ip = nl(nf)
      do kp=nl(5),nl(6)
         do jp=nl(3),nl(4)
            blochi(ibl,igr)%cella(ip,jp,kp)%dist = 0.0d0
            blochi(ibl,igr)%cella(ip,jp,kp)%wall = p
            blochi(ibl,igr)%cella(ip,jp,kp)%boul = .true.
            jmin = jp
            kmin = kp
            ib = ip
            jm = jp
            km = kp
            dm = 0d0
            do while ( (dm.lt.5d0*blt) .and. (ib*(nii-ib).ge.0) )
               ib = ib+is
               call getvertex(igr,ibl,ib,jp,kp,cb)
               do kpp=max(nl(5),km-1),min(nl(6),km+1)
                  do jpp=max(nl(3),jm-1),min(nl(4),jm+1)
                     call getvertex(igr,ibl,ip,jpp,kpp,cp)
                     dm = dist(cp,cb)
                     if (dm.lt.blochi(ibl,igr)%cella(ib,jp,kp)%dist) then
                        blochi(ibl,igr)%cella(ib,jp,kp)%dist = dm
                        blochi(ibl,igr)%cella(ib,jp,kp)%wall = p
                        blochi(ibl,igr)%cella(ib,jp,kp)%boul = .true.
                        jmin = jpp
                        kmin = kpp
                     end if
                  end do
               end do
               dm = blochi(ibl,igr)%cella(ib,jp,kp)%dist
               jm = jmin
               km = kmin
            end do
         end do
      end do

   ! facce 5 e 6
   else if (nf.gt.4) then
      ql(1:6,1) = (/ nl(1),nl(2), nl(3),nl(3), nl(nf),nl(nf) /)
      ql(1:6,2) = (/ nl(1),nl(2), nl(4),nl(4), nl(nf),nl(nf) /)
      ql(1:6,3) = (/ nl(1),nl(1), nl(3),nl(4), nl(nf),nl(nf) /)
      ql(1:6,4) = (/ nl(2),nl(2), nl(3),nl(4), nl(nf),nl(nf) /)
      ks = 11-2*nf
      kp = nl(nf)
      do jp=nl(3),nl(4)
         do ip=nl(1),nl(2)
            blochi(ibl,igr)%cella(ip,jp,kp)%dist = 0.0d0
            blochi(ibl,igr)%cella(ip,jp,kp)%wall = p
            blochi(ibl,igr)%cella(ip,jp,kp)%boul = .true.
            imin = ip
            jmin = jp
            im = ip
            jm = jp
            kb = kp
            dm = 0d0
            do while ( (dm.lt.5d0*blt) .and. (kb*(nkk-kb).ge.0) )
               kb = kb+ks
               call getvertex(igr,ibl,ip,jp,kb,cb)
               do jpp=max(nl(3),jm-1),min(nl(4),jm+1)
                  do ipp=max(nl(1),im-1),min(nl(2),im+1)
                     call getvertex(igr,ibl,ipp,jpp,kp,cp)
                     dm = dist(cp,cb)
                     if (dm.lt.blochi(ibl,igr)%cella(ip,jp,kb)%dist) then
                        blochi(ibl,igr)%cella(ip,jp,kb)%dist = dm
                        blochi(ibl,igr)%cella(ip,jp,kb)%wall = p
                        blochi(ibl,igr)%cella(ip,jp,kb)%boul = .true.
                        imin = ipp
                        jmin = jpp
                     end if
                  end do
               end do
               dm = blochi(ibl,igr)%cella(ip,jp,kb)%dist
               im = imin
               jm = jmin
            end do
         end do
      end do

   ! facce 3 e 4
   else
      ql(1:6,1) = (/ nl(1),nl(2), nl(nf),nl(nf), nl(5),nl(5) /)
      ql(1:6,2) = (/ nl(1),nl(2), nl(nf),nl(nf), nl(6),nl(6) /)
      ql(1:6,3) = (/ nl(1),nl(1), nl(nf),nl(nf), nl(5),nl(6) /)
      ql(1:6,4) = (/ nl(2),nl(2), nl(nf),nl(nf), nl(5),nl(6) /)
      js = 7-2*nf
      jp = nl(nf)
      do kp=nl(5),nl(6)
         do ip=nl(1),nl(2)
            blochi(ibl,igr)%cella(ip,jp,kp)%wall = 0.0d0
            blochi(ibl,igr)%cella(ip,jp,kp)%wall = p
            blochi(ibl,igr)%cella(ip,jp,kp)%boul = .true.
            imin = ip
            kmin = kp
            im = ip
            jb = jp
            km = kp
            dm = 0d0
            do while ( (dm.lt.5d0*blt) .and. (jb*(njj-jb).ge.0) )
               jb = jb+js
               call getvertex(igr,ibl,ip,jb,kp,cb)
               do kpp=max(nl(5),km-1),min(nl(6),km+1)
                  do ipp=max(nl(1),im-1),min(nl(2),im+1)
                     call getvertex(igr,ibl,ipp,jp,kpp,cp)
                     dm = dist(cp,cb)
                     if (dm.lt.blochi(ibl,igr)%cella(ip,jb,kp)%dist) then
                        blochi(ibl,igr)%cella(ip,jb,kp)%dist = dm
                        blochi(ibl,igr)%cella(ip,jb,kp)%wall = p
                        blochi(ibl,igr)%cella(ip,jb,kp)%boul = .true.
                        imin = ipp
                        kmin = kpp
                     end if
                  end do
               end do
               dm = blochi(ibl,igr)%cella(ip,jb,kp)%dist
               im = imin
               km = kmin
            end do
         end do
      end do
   end if

end do

! loop sugli spigoli di parete
do p=1,nbe

   ibl = be(1,p)
   nl  = be(3:8,p)/tgr

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   do kb=0,nkk
      do jb=0,njj
         do ib=0,nii

            call getvertex(igr,ibl,ib,jb,kb,cb)

            dm = blochi(ibl,igr)%cella(ib,jb,kb)%dist
            do kp=nl(5),nl(6)
            do jp=nl(3),nl(4)
            do ip=nl(1),nl(2)
               call getvertex(igr,ibl,ip,jp,kp,cp)
               dm = min(dm,dist(cp,cb))
            end do
            end do
            end do
            blochi(ibl,igr)%cella(ib,jb,kb)%dist = dm

         end do
      end do
   end do

end do

! calcolo distanze nei centri cella
do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   do k=nkk,1,-1
      do j=njj,1,-1
         do i=nii,1,-1

            ! trasferisce l'informazione sul patch piu' vicino dai vertici
            ! ai centri cella
            dm = 1.0d20
            p = blochi(ibl,igr)%cella(i,j,k)%wall
            do kk=k-1,k
               do jj=j-1,j
                  do ii=i-1,i
                     if (blochi(ibl,igr)%cella(ii,jj,kk)%dist.lt.dm) then
                        dm = blochi(ibl,igr)%cella(ii,jj,kk)%dist
                        p  = blochi(ibl,igr)%cella(ii,jj,kk)%wall
                     end if
                  end do
               end do
            end do
            blochi(ibl,igr)%cella(i,j,k)%wall = p

            ! distanza nei centri cella come media di quella nei vertici
            ! (cioè, quadrato della media delle distanze lineari dei vertici)
            dm = 0.0d0
            do kk=k-1,k
               do jj=j-1,j
                  do ii=i-1,i
                     dm = dm + sqrt(blochi(ibl,igr)%cella(ii,jj,kk)%dist)
                  end do
               end do
            end do
            blochi(ibl,igr)%cella(i,j,k)%dist = (0.125d0*dm)**2

            ! una cella è di strato limite solo se lo sono i suoi 8 vertici
            if (.not.all(blochi(ibl,igr)%cella(i-1:i,j-1:j,k-1:k)%boul)) then
               blochi(ibl,igr)%cella(i,j,k)%boul = .false.
            end if

         end do
      end do
   end do

   ! distanza nei centri cella delle celle di cornice
   do k=1,nkk
      do j=1,njj
         blochi(ibl,igr)%cella(-1:0,j,k)%dist = &
         blochi(ibl,igr)%cella(1,j,k)%dist
         blochi(ibl,igr)%cella(nii+1:nii+2,j,k)%dist = &
         blochi(ibl,igr)%cella(nii,j,k)%dist
         !
         blochi(ibl,igr)%cella(-1:0,j,k)%wall = &
         blochi(ibl,igr)%cella(1,j,k)%wall
         blochi(ibl,igr)%cella(nii+1:nii+2,j,k)%wall = &
         blochi(ibl,igr)%cella(nii,j,k)%wall
         !
         blochi(ibl,igr)%cella(-1:0,j,k)%boul = &
         blochi(ibl,igr)%cella(1,j,k)%boul
         blochi(ibl,igr)%cella(nii+1:nii+2,j,k)%boul = &
         blochi(ibl,igr)%cella(nii,j,k)%boul
      end do
   end do

   do k=1,nkk
      do i=-1,nii+2
         blochi(ibl,igr)%cella(i,-1:0,k)%dist = &
         blochi(ibl,igr)%cella(i,1,k)%dist
         blochi(ibl,igr)%cella(i,njj+1:njj+2,k)%dist = &
         blochi(ibl,igr)%cella(i,njj,k)%dist
         !
         blochi(ibl,igr)%cella(i,-1:0,k)%wall = &
         blochi(ibl,igr)%cella(i,1,k)%wall
         blochi(ibl,igr)%cella(i,njj+1:njj+2,k)%wall = &
         blochi(ibl,igr)%cella(i,njj,k)%wall
         !
         blochi(ibl,igr)%cella(i,-1:0,k)%boul = &
         blochi(ibl,igr)%cella(i,1,k)%boul
         blochi(ibl,igr)%cella(i,njj+1:njj+2,k)%boul = &
         blochi(ibl,igr)%cella(i,njj,k)%boul
      end do
   end do

   do j=-1,njj+2
      do i=-1,nii+2
         blochi(ibl,igr)%cella(i,j,-1:0)%dist = &
         blochi(ibl,igr)%cella(i,j,1)%dist
         blochi(ibl,igr)%cella(i,j,nkk+1:nkk+2)%dist = &
         blochi(ibl,igr)%cella(i,j,nkk)%dist
         !
         blochi(ibl,igr)%cella(i,j,-1:0)%wall = &
         blochi(ibl,igr)%cella(i,j,1)%wall
         blochi(ibl,igr)%cella(i,j,nkk+1:nkk+2)%wall = &
         blochi(ibl,igr)%cella(i,j,nkk)%wall
         !
         blochi(ibl,igr)%cella(i,j,-1:0)%boul = &
         blochi(ibl,igr)%cella(i,j,1)%boul
         blochi(ibl,igr)%cella(i,j,nkk+1:nkk+2)%boul = &
         blochi(ibl,igr)%cella(i,j,nkk)%boul
      end do
   end do

end do

end subroutine
