!===============================================================================
!  regolarizzazione dimensioni di cella
!===============================================================================
subroutine smovol(igr)
use prec
use modpar, only : nbl, blochi
implicit none
integer (kind=I4P) :: igr,ibl
integer (kind=I4P) :: nii,njj,nkk,nll
integer (kind=I4P) :: i,j,k,l
real    (kind=R8P) :: vtmp

do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   ! smoothing delle dimensioni di cella
   nll = max(min(nii,njj,nkk)/4,4)

   do l=1,nll

      do k=1,nkk
         do j=1,njj
            blochi(ibl,igr)%cella(-1:0,j,k)%size = &
            blochi(ibl,igr)%cella(1,j,k)%size
            blochi(ibl,igr)%cella(nii+1:nii+2,j,k)%size = &
            blochi(ibl,igr)%cella(nii,j,k)%size
         end do
      end do

      do k=1,nkk
         do i=-1,nii+2
            blochi(ibl,igr)%cella(i,-1:0,k)%size = &
            blochi(ibl,igr)%cella(i,1,k)%size
            blochi(ibl,igr)%cella(i,njj+1:njj+2,k)%size = &
            blochi(ibl,igr)%cella(i,njj,k)%size
         end do
      end do

      do j=-1,njj+2
         do i=-1,nii+2
            blochi(ibl,igr)%cella(i,j,-1:0)%size = &
            blochi(ibl,igr)%cella(i,j,1)%size
            blochi(ibl,igr)%cella(i,j,nkk+1:nkk+2)%size = &
            blochi(ibl,igr)%cella(i,j,nkk)%size
         end do
      end do

      do k=1,nkk
         do j=1,njj
            do i=1,nii
               vtmp = 6.0d0 * blochi(ibl,igr)%cella(i  ,j,k)%size &
                    +         blochi(ibl,igr)%cella(i-1,j,k)%size &
                    +         blochi(ibl,igr)%cella(i+1,j,k)%size &
                    +         blochi(ibl,igr)%cella(i,j-1,k)%size &
                    +         blochi(ibl,igr)%cella(i,j+1,k)%size &
                    +         blochi(ibl,igr)%cella(i,j,k-1)%size &
                    +         blochi(ibl,igr)%cella(i,j,k+1)%size
               blochi(ibl,igr)%cella(i,j,k)%size = vtmp/12.0d0
            end do
         end do
      end do

   end do

   ! ri-definisce la dimensione di cella come minimo tra la radice cubica del
   ! volume e la distanza dalla parete
   do k=1,nkk
      do j=1,njj
         do i=1,nii
            vtmp = blochi(ibl,igr)%cella(i,j,k)%size ** (1.0d0/3.0d0)
            blochi(ibl,igr)%cella(i,j,k)%size = &
            min(vtmp,sqrt(blochi(ibl,igr)%cella(i,j,k)%dist))
         end do
      end do
   end do

end do

end subroutine
