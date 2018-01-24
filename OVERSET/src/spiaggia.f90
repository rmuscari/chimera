!===============================================================================
!  Cerca eventuali punti di spiaggia
!===============================================================================
subroutine spiaggia(igr)
use prec
use moddef, only : donor,dist,point
use modpar, only : stdflag,innerb,dispiaggia,donmax,&
                   bp,nbl,npa,patch_flags,blogrd,blochi
implicit none
integer (kind=I4P) :: igr
integer (kind=I4P) :: nii,njj,nkk
real (kind=R8P), parameter :: eps=1d-10
integer (kind=I4P) :: ibl,i,j,k,ip,jp,kp,p,nf
real (kind=R8P)    :: dd,disquare,dispatch
type (point)     :: cb,cp
type (donor)     :: md(donmax)

real (kind=R8P), allocatable, dimension(:,:,:) :: peso

!...............................................................................

!  stampa di controllo
write(*,'(a,i3,a)') 'inizio spiaggia  (igr =',igr,')'

! loop sui blocchi
do ibl=1,nbl

   nii = blogrd(ibl,igr)%ni
   njj = blogrd(ibl,igr)%nj
   nkk = blogrd(ibl,igr)%nk

   allocate(peso(0:nii,0:njj,0:nkk))
   peso = 0d0

   ! loop sui patches
   do p=1,npa

      ! se non è un patch di inflow/outflow salta al successivo
      if ((bp(1,p).ne.ibl).or.(bp(4,p).lt.1).or.(.not.patch_flags(3,p))) cycle

      dispatch = real(bp(4,p),kind(0d0))/dispiaggia
      disquare = dispatch**2

      nf = bp(2,p)

      do k=0,nkk
         do j=0,njj
            do i=0,nii

               ip = i
               jp = j
               kp = k
               select case (nf)
                  case(1) ; ip = 0
                  case(2) ; ip = nii
                  case(3) ; jp = 0
                  case(4) ; jp = njj
                  case(5) ; kp = 0
                  case(6) ; kp = nkk
               end select

               cb = blogrd(ibl,igr)%nodo(i,j,k)
               cp = blogrd(ibl,igr)%nodo(ip,jp,kp)
               dd = dist(cp,cb)

               if (dd.lt.disquare) then
                  peso(i,j,k) = max(peso(i,j,k),1d0-sqrt(dd)/dispatch)
               end if
            end do
         end do
      end do

   end do

   ! trasferisce i pesi nei centri cella
   do k=1,nkk
      do j=1,njj
         do i=1,nii
            p = blochi(ibl,igr)%cella(i,j,k)%t
            if (p.eq.stdflag) then
               dd = 0.125d0 * ( peso(i  ,j  ,k  ) + peso(i-1,j  ,k  ) &
                              + peso(i  ,j-1,k  ) + peso(i-1,j-1,k  ) &
                              + peso(i  ,j  ,k-1) + peso(i-1,j  ,k-1) &
                              + peso(i  ,j-1,k-1) + peso(i-1,j-1,k-1) )
               if (dd.gt.eps) then
                  md(1) = donor(ibl,i,j,k,dd)
                  call newchi(igr,ibl,i,j,k,innerb,1,md)
               end if
            end if
         end do
      end do
   end do

   deallocate(peso)
end do

end subroutine
