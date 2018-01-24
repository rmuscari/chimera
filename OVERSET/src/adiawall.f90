!===============================================================================
! Cambia tipo alle celle chimera di contorno adiacenti ad una parete interna
!===============================================================================
subroutine cambia_in_parete(igr)
use prec
use modpar, only : parete,innerw,offgen,offchi,blochi,nbl
implicit none
integer (kind=I4P) :: i,j,k,igr,ibl
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: tcc,m

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  stampa di controllo
write(*,'(a,i3,a)') 'inizio cambia_in_parete  (igr =',igr,')'

! Se una cella di cornice e' chimera e la cella interna adiacente e' una parete
! interna, cambia in parete la c.c. della cella di cornice
!
do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   i = 0
   m = 1
   do k=1,nkk
      do j=1,njj
         tcc = blochi(ibl,igr)%cella(m,j,k)%t
         if (tcc.eq.innerw) then
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if ((tcc.gt.offchi .and. tcc.lt.offchi+7) .or. &
                (tcc.gt.offgen .and. tcc.lt.offgen+7))     &
                blochi(ibl,igr)%cella(i,j,k)%t = parete
         end if
      end do
   end do
   i = nii+1
   m = nii
   do k=1,nkk
      do j=1,njj
         tcc = blochi(ibl,igr)%cella(m,j,k)%t
         if (tcc.eq.innerw) then
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if ((tcc.gt.offchi .and. tcc.lt.offchi+7) .or. &
                (tcc.gt.offgen .and. tcc.lt.offgen+7))     &
                blochi(ibl,igr)%cella(i,j,k)%t = parete
         end if
      end do
   end do

   j = 0
   m = 1
   do k=1,nkk
      do i=1,nii
         tcc = blochi(ibl,igr)%cella(i,m,k)%t
         if (tcc.eq.innerw) then
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if ((tcc.gt.offchi .and. tcc.lt.offchi+7) .or. &
                (tcc.gt.offgen .and. tcc.lt.offgen+7))     &
                blochi(ibl,igr)%cella(i,j,k)%t = parete
         end if
      end do
   end do
   j = njj+1
   m = njj
   do k=1,nkk
      do i=1,nii
         tcc = blochi(ibl,igr)%cella(i,m,k)%t
         if (tcc.eq.innerw) then
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if ((tcc.gt.offchi .and. tcc.lt.offchi+7) .or. &
                (tcc.gt.offgen .and. tcc.lt.offgen+7))     &
                blochi(ibl,igr)%cella(i,j,k)%t = parete
         end if
      end do
   end do

   k = 0
   m = 1
   do j=1,njj
      do i=1,nii
         tcc = blochi(ibl,igr)%cella(i,j,m)%t
         if (tcc.eq.innerw) then
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if ((tcc.gt.offchi .and. tcc.lt.offchi+7) .or. &
                (tcc.gt.offgen .and. tcc.lt.offgen+7))     &
                blochi(ibl,igr)%cella(i,j,k)%t = parete
         end if
      end do
   end do
   k = nkk+1
   m = nkk
   do j=1,njj
      do i=1,nii
         tcc = blochi(ibl,igr)%cella(i,j,m)%t
         if (tcc.eq.innerw) then
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if ((tcc.gt.offchi .and. tcc.lt.offchi+7) .or. &
                (tcc.gt.offgen .and. tcc.lt.offgen+7))     &
                blochi(ibl,igr)%cella(i,j,k)%t = parete
         end if
      end do
   end do

end do

return
end subroutine cambia_in_parete

!===============================================================================
! Se una cella di parete è adiacente ad una cella non standard, trasforma la
! cella di parete in una parete non attiva (tipo = 11)
! Per evitare di calcolare la forza nelle celle di parete inizialmente marcate
! chimera da patch a priorità più alta e ripristinate "standard" dall'opzione
! "more overlap"
!===============================================================================
subroutine parete_non_attiva(igr)
use prec
use modpar, only : parete,inawal,isowal,isoina,nbl,blochi
implicit none
integer (kind=I4P) :: igr,ibl
integer (kind=I4P) :: i,j,k
integer (kind=I4P) :: nii,njj,nkk

do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   ! controllo chimera-parete
   do k=1,nkk
      do j=1,njj
         ! pareti adiabatiche
         if ((blochi(ibl,igr)%cella(0,j,k)%t.eq.parete) .and. &
             (blochi(ibl,igr)%cella(1,j,k)%t.gt.0)) then
            blochi(ibl,igr)%cella( 0,j,k)%t = inawal
            blochi(ibl,igr)%cella(-1,j,k)%t = inawal
         end if
         if ((blochi(ibl,igr)%cella(nii+1,j,k)%t.eq.parete) .and. &
             (blochi(ibl,igr)%cella(nii,j,k)%t.gt.0)) then
            blochi(ibl,igr)%cella(nii+1,j,k)%t = inawal
            blochi(ibl,igr)%cella(nii+2,j,k)%t = inawal
         end if
         !
         ! pareti isoterme
         if ((blochi(ibl,igr)%cella(0,j,k)%t.eq.isowal) .and. &
             (blochi(ibl,igr)%cella(1,j,k)%t.gt.0)) then
            blochi(ibl,igr)%cella( 0,j,k)%t = isoina
            blochi(ibl,igr)%cella(-1,j,k)%t = isoina
         end if
         if ((blochi(ibl,igr)%cella(nii+1,j,k)%t.eq.isowal) .and. &
             (blochi(ibl,igr)%cella(nii,j,k)%t.gt.0)) then
            blochi(ibl,igr)%cella(nii+1,j,k)%t = isoina
            blochi(ibl,igr)%cella(nii+2,j,k)%t = isoina
         end if
      end do
   end do

   do k=1,nkk
      do i=1,nii
         ! pareti adiabatiche
         if ((blochi(ibl,igr)%cella(i,0,k)%t.eq.parete) .and. &
             (blochi(ibl,igr)%cella(i,1,k)%t.gt.0)) then
            blochi(ibl,igr)%cella(i, 0,k)%t = inawal
            blochi(ibl,igr)%cella(i,-1,k)%t = inawal
         end if
         if ((blochi(ibl,igr)%cella(i,njj+1,k)%t.eq.parete) .and.  &
             (blochi(ibl,igr)%cella(i,njj,k)%t.gt.0)) then
            blochi(ibl,igr)%cella(i,njj+1,k)%t = inawal
            blochi(ibl,igr)%cella(i,njj+2,k)%t = inawal
         end if
         ! pareti isoterme
         if ((blochi(ibl,igr)%cella(i,0,k)%t.eq.isowal) .and.  &
             (blochi(ibl,igr)%cella(i,1,k)%t.gt.0)) then
            blochi(ibl,igr)%cella(i, 0,k)%t = isoina
            blochi(ibl,igr)%cella(i,-1,k)%t = isoina
         end if
         if ((blochi(ibl,igr)%cella(i,njj+1,k)%t.eq.isowal) .and.  &
             (blochi(ibl,igr)%cella(i,njj,k)%t.gt.0)) then
            blochi(ibl,igr)%cella(i,njj+1,k)%t = isoina
            blochi(ibl,igr)%cella(i,njj+2,k)%t = isoina
         end if
      end do
   end do

   do j=1,njj
      do i=1,nii
         ! pareti adiabatiche
         if ((blochi(ibl,igr)%cella(i,j,0)%t.eq.parete) .and. &
             (blochi(ibl,igr)%cella(i,j,1)%t.gt.0)) then
            blochi(ibl,igr)%cella(i,j, 0)%t = inawal
            blochi(ibl,igr)%cella(i,j,-1)%t = inawal
         end if
         if ((blochi(ibl,igr)%cella(i,j,nkk+1)%t.eq.parete) .and. &
             (blochi(ibl,igr)%cella(i,j,nkk)%t.gt.0)) then
            blochi(ibl,igr)%cella(i,j,nkk+1)%t = inawal
            blochi(ibl,igr)%cella(i,j,nkk+2)%t = inawal
         end if
         ! pareti isoterme
         if ((blochi(ibl,igr)%cella(i,j,0)%t.eq.isowal) .and. &
             (blochi(ibl,igr)%cella(i,j,1)%t.gt.0)) then
            blochi(ibl,igr)%cella(i,j, 0)%t = isoina
            blochi(ibl,igr)%cella(i,j,-1)%t = isoina
         end if
         if ((blochi(ibl,igr)%cella(i,j,nkk+1)%t.eq.isowal) .and. &
             (blochi(ibl,igr)%cella(i,j,nkk)%t.gt.0)) then
            blochi(ibl,igr)%cella(i,j,nkk+1)%t = isoina
            blochi(ibl,igr)%cella(i,j,nkk+2)%t = isoina
         end if
      end do
   end do

end do

end subroutine
