!===============================================================================
!  postprocessore per XSHIP - scrittura file tecplot
!===============================================================================
subroutine printec
use prec
use modpar
implicit none

character (len=1), parameter :: nc = achar(0)

integer (kind=I4P) :: i,j,k,m,n

integer (kind=I4P) :: TECINI112,TECAUXSTR112,TECZNE112,TECDAT112,TECEND112
integer (kind=I4P), parameter :: fty=0, deb=0, dbl=1
integer (kind=I4P), parameter :: ZoneType=0
integer (kind=I4P), parameter :: zero=0, uno=1
integer (kind=I4P), dimension(50) :: ValueLocation
integer (kind=I4P), dimension(50), parameter :: PassiveVarList = 0
logical (kind=I1P) :: bplot

integer (kind=I4P), pointer :: ShareVarFromZone => null()

real (kind=R8P) :: SolutionTime
real (kind=R8P), allocatable :: var(:,:,:)
real (kind=R8P), allocatable :: tmp(:)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  Inizializzazione file Tecplot
i = TECINI112( trim(dset)//nc, &
               trim(vars)//nc, &
               trim(fout)//nc, &
               trim(work)//nc, &
               fty,deb,dbl)
print *,trim(vars)

write(work,'(f30.15)') uns
write(work,'(a)') adjustl(work)
write(temp,'(a)') trim(work) // nc
write(work,'(a)') "SolutionTime" // nc
i = TECAUXSTR112(work,temp);

!  Loop sui blocchi
do n=1,lb

   !  Allocazione memoria
   i = ni(n)+gc(1,n)+gc(2,n)+1
   j = nj(n)+gc(3,n)+gc(4,n)+1
   k = nk(n)+gc(5,n)+gc(6,n)+1
   allocate(var(i,j,k))
   allocate(tmp(i*j*k))

   bplot = n.ge.fb

   !  zona corrispondente al blocco "n"
   if (bplot) then
      write(temp,'(i3.3)') n
      temp = "blocco "//trim(temp)//nc
      SolutionTime = uns
      if (vertex) then
         i = ni(n)/is+1
         j = nj(n)/js+1
         k = nk(n)/ks+1
         ValueLocation = 1
         m = TECZNE112( trim(temp)//nc, &    ! zone title
                        ZoneType,       &    ! zone type = ordered
                        i,j,k,          &    ! I,J,K max
                        zero,zero,zero, &    ! non utilizzati
                        SolutionTime,   &    !
                        zero,zero,      &    ! non utilizzati
                        uno,            &    ! IsBlock deve essere = 1
                        zero,zero,zero, &    ! non utilizzati
                        zero,zero,      &    ! non utilizzati
                        PassiveVarList, &    ! 0 = attiva, 1 = passiva
                        ValueLocation,  &    ! node-centered = 1
                        ShareVarFromZone,zero)
!                       %var(0)         ,zero)
      else
         i = ni(n)/is+2
         j = nj(n)/js+2
         k = nk(n)/ks+2
!        i = ni(n)/is
!        j = nj(n)/js
!        k = nk(n)/ks
         ValueLocation = 1
         m = TECZNE112( trim(temp)//nc, &    ! zone title
                        ZoneType,       &    ! zone type = ordered
                        i,j,k,          &    ! I,J,K max
                        zero,zero,zero, &    ! non utilizzati
                        SolutionTime,   &    !
                        zero,zero,      &    ! non utilizzati
                        uno,            &    ! IsBlock deve essere = 1
                        zero,zero,zero, &    ! non utilizzati
                        zero,zero,      &    ! non utilizzati
                        PassiveVarList, &    ! 0 = attiva, 1 = passiva
                        ValueLocation,  &    ! cell-centered = 0
                        ShareVarFromZone,zero)
!                       %var(0)         ,zero)
      end if
   end if

   !  x, y, z
   do j=1,3
      call readxyz(var,tmp,ni(n),nj(n),nk(n),gc(1,n),m,bplot)
      if (bplot)  i = TECDAT112(m,tmp,dbl)
   end do

   !  variabili dipendenti
   do j=1,nv
      call readvar(var,tmp,ni(n),nj(n),nk(n),gc(1,n),m,bplot)
      if (bplot)  i = TECDAT112(m,tmp,dbl)
   end do

   !  icc
   if (cc) then
      call readicc(n,tmp,m,bplot)
      if (bplot)  i = TECDAT112(m,tmp,dbl)
   end if

   !  Libera memoria allocata
   deallocate(var)
   deallocate(tmp)

end do

13 continue

!  Chiusura file Tecplot
i = TECEND112()

end subroutine

!===============================================================================
!  legge le coordinate dal file di reticolo e le memorizza in "tmp"
!===============================================================================
subroutine readxyz(var,tmp,ni,nj,nk,gc,m,p)
use prec
use modpar, only : nret, &
                   is,js,ks, &
                   vertex
implicit none
integer (kind=I4P), intent(IN) :: ni,nj,nk,gc(6)
integer (kind=I4P) :: m
logical (kind=I1P) :: p
integer (kind=I4P) :: i,j,k,cl,cr
real (kind=R8P) :: var(-gc(1):ni+gc(2),-gc(3):nj+gc(4),-gc(5):nk+gc(6))
real (kind=R8P) :: tmp(*)

!  lettura coordinate
read(nret) (((var(i,j,k),i=-gc(1),ni+gc(2)),&
                         j=-gc(3),nj+gc(4)),&
                         k=-gc(5),nk+gc(6))
if (.not.p) return

!  numero di dati da scrivere (usato da TECDAT112)
m = 0

!  valori nei vertici delle celle
if (vertex) then
   do k=0,nk,ks
      do j=0,nj,js
         do i=0,ni,is
            m = m+1
            tmp(m) = var(i,j,k)
         end do
      end do
   end do

!  valori nei centri cella
else
   cr = 1 !!!
   cl = 1-cr
   do k=cl,nk+cr,ks
      do j=cl,nj+cr,js
         do i=cl,ni+cr,is
            m = m+1
            tmp(m) = var(i-1,j-1,k-1) + var(i,j-1,k-1) &
                   + var(i-1,j,k-1)   + var(i,j,k-1)   &
                   + var(i-1,j-1,k)   + var(i,j-1,k)   &
                   + var(i-1,j,k)     + var(i,j,k)
            tmp(m) = 0.125_R8P*tmp(m)
         end do
      end do
   end do
end if

end subroutine

!...............................................................................

subroutine readvar(var,tmp,ni,nj,nk,gc,m,p)
use prec
use modpar, only : nsol, &
                   is,js,ks, &
                   vertex
implicit none
integer (kind=I4P), intent(IN) :: ni,nj,nk,gc(6)
integer (kind=I4P) :: m
logical (kind=I1P) :: p
integer (kind=I4P) :: i,j,k,cl,cr
real (kind=R8P) :: var(-gc(1)+1:ni+gc(2),-gc(3)+1:nj+gc(4),-gc(5)+1:nk+gc(6))
real (kind=R8P) :: tmp(*)

read(nsol) (((var(i,j,k),i=-gc(1)+1,ni+gc(2)),&
                         j=-gc(3)+1,nj+gc(4)),&
                         k=-gc(5)+1,nk+gc(6))
if (.not.p) return

m = 0

if (vertex) then
   do k=0,nk,ks
      do j=0,nj,js
         do i=0,ni,is
            m = m+1
            tmp(m) = var(i+1,j+1,k+1) + var(i,j+1,k+1) &
                   + var(i+1,j,k+1)   + var(i,j,k+1)   &
                   + var(i+1,j+1,k)   + var(i,j+1,k)   &
                   + var(i+1,j,k)     + var(i,j,k)
            tmp(m) = 0.125_R8P*tmp(m)
         end do
      end do
   end do
else
   cr = 1 !!!
   cl = 1-cr
   do k=cl,nk+cr,ks
      do j=cl,nj+cr,js
         do i=cl,ni+cr,is
            m = m+1
            tmp(m) = var(i,j,k)
         end do
      end do
   end do
end if

end subroutine

!...............................................................................

subroutine readicc(n,tmp,m,p)
use prec
use modpar, only : nicc, is,js,ks, vertex,anycc, blochi
implicit none
integer (kind=I4P), intent(IN) :: n
integer (kind=I4P) :: m
logical (kind=I1P) :: p
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: i,j,k,q,cl,cr
real (kind=R8P) :: tmp(*)

read(nicc) (((blochi(n)%cella(i,j,k)%t,&
              i = -blochi(n)%gc(1)+1, blochi(n)%ni+blochi(n)%gc(2)),&
              j = -blochi(n)%gc(3)+1, blochi(n)%nj+blochi(n)%gc(4)),&
              k = -blochi(n)%gc(5)+1, blochi(n)%nk+blochi(n)%gc(6))
if (.not.p) return

nii = blochi(n)%ni
njj = blochi(n)%nj
nkk = blochi(n)%nk

m = 0

if (vertex) then

   if (anycc) then
      ! un nodo è chimera se una sola delle celle limitrofe lo è
      blochi(n)%cella(0    ,:,:)%t = blochi(n)%cella(1  ,:,:)%t
      blochi(n)%cella(nii+1,:,:)%t = blochi(n)%cella(nii,:,:)%t
      blochi(n)%cella(:,0    ,:)%t = blochi(n)%cella(:,1  ,:)%t
      blochi(n)%cella(:,njj+1,:)%t = blochi(n)%cella(:,njj,:)%t
      blochi(n)%cella(:,:,0    )%t = blochi(n)%cella(:,:,1  )%t
      blochi(n)%cella(:,:,nkk+1)%t = blochi(n)%cella(:,:,nkk)%t
      do k=0,nkk,ks
         do j=0,njj,js
            do i=0,nii,is
               m = m+1
               tmp(m) = 0.0_R8P
               if (any(blochi(n)%cella(i:i+1,j:j+1,k:k+1)%t.gt.19)) &
                   tmp(m) = 1.0_R8P
            end do
         end do
      end do
   else
      ! un nodo è chimera se lo sono tutte le celle limitrofe
      ! opzione di default
      do k=0,nkk,ks
         do j=0,njj,js
            do i=0,nii,is
               m = m+1
               tmp(m) = 0.0_R8P
               if (all(blochi(n)%cella(i:i+1,j:j+1,k:k+1)%t.gt.0)) &
                   tmp(m) = 1.0_R8P
               q = min(i*(nii-i),1) +  min(j*(njj-j),1) + min(k*(nkk-k),1)
               ! q = 0  il nodo è su un vertice
               ! q = 1  il nodo è su uno spigolo
               ! q = 2  il nodo è su una faccia
               ! q = 3  il nodo è all'interno del blocco
               if (q.eq.2 .and. tmp(m).lt.0.6_R8P) then
                  if (i.eq.0) then
                     if (all(blochi(n)%cella(i,j:j+1,k:k+1)%t.gt.19)) &
                     tmp(m) = 0.5_R8P
                  else if (i.eq.nii) then
                     if (all(blochi(n)%cella(i+1,j:j+1,k:k+1)%t.gt.19)) &
                     tmp(m) = 0.5_R8P
                  end if
                  if (j.eq.0) then
                     if (all(blochi(n)%cella(i:i+1,j  ,k:k+1)%t.gt.19)) &
                     tmp(m) = 0.5_R8P
                  else if (j.eq.njj) then
                     if (all(blochi(n)%cella(i:i+1,j+1,k:k+1)%t.gt.19)) &
                     tmp(m) = 0.5_R8P
                  end if
                  if (k.eq.0) then
                     if (all(blochi(n)%cella(i:i+1,j:j+1,k  )%t.gt.19)) &
                     tmp(m) = 0.5_R8P
                  else if (k.eq.nkk) then
                     if (all(blochi(n)%cella(i:i+1,j:j+1,k+1)%t.gt.19)) &
                     tmp(m) = 0.5_R8P
                  end if
               end if
            end do
         end do
      end do
   !
   end if

else
   cr = 1 !!!
   cl = 1-cr
   do k=cl,nkk+cr,ks
      do j=cl,njj+cr,js
         do i=cl,nii+cr,is
            m = m+1
            tmp(m) = 0.0_R8P
            if (blochi(n)%cella(i,j,k)%t.gt.0) tmp(m) = 1.0_R8P
         end do
      end do
   end do
end if

end subroutine
