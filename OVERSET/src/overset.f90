!===============================================================================
!  preprocessore CHIMERA per XSHIP - main
!===============================================================================
program overset
use prec
use modpar
implicit none

character (len=50) :: ftmp
integer (kind=I4P) :: hms(8),igr
integer (kind=I4P) :: omp_get_num_threads,omp_get_thread_num

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  apre file di debug
open(ndeb,file='overset.deb',form='formatted',status='unknown')
rewind(ndeb)

! check sui thread omp
igr = omp_get_num_threads()
!$omp parallel   &
!$omp default(none)
write(*,'(a,i2,a,i2,a)')  ' Thread: ', omp_get_thread_num(), &
        ' of ', omp_get_num_threads(), ' is alive'
!$omp end parallel

!  lettura parametri
call readpar

!  lettura reticoli [ generazione cornici, scrittura reticoli con cornici ]
call iogrid

!  ...
call calcolo_spessore_strato_limite

!  inizializzazione matrici icc e rcc e puntatore prc
call initcc

!  inizio loop sui reticoli
do igr=ngr,fgr,-1

   !  chiamata sistema per stampa ora
   call date_and_time(ftmp,ftmp,ftmp,hms)
   write(ftmp,'(i2.2,a,i2.2,a,i2.2)') hms(5),':',hms(6),':',hms(7)
   write(*   ,'(/,3x,a,4x,a,i3,a)') trim(ftmp),'(igr =',igr,')'
   write(ndeb,'(/,3x,a,4x,a,i3,a)') trim(ftmp),'(igr =',igr,')'

   !  calcolo metrica
   call metrix(igr)

   ! calcolo distanza dalla parete
   call waldis(igr)

   ! regolarizzazione dimensioni di cella
   call smovol(igr)

   !  prima le celle chimera all'interno dei blocchi ...
   call intpoint(igr)

   !  c.c. sulle facce di contorno dei blocchi
   call calface(igr)

   !  c.c. sugli spigoli dei blocchi
   call calspig(igr)

   !  c.c. nei vertici dei blocchi
   call calvert(igr)

   ! trasforma in pareti non attive le pareti adiacenti a celle non standard
   ! (da eseguire prima di un eventuale "more_overlap")
   call parete_non_attiva(igr)

   !  aumenta la sovrapposizione delle griglie
   if (mover) call more_overlap(igr)

   !  per le celle standard superstiti, vede se finiscono dentro una spiaggia
   if (dispiaggia.gt.0d0) call spiaggia(igr)

   !  cambia tipo alle celle chimera di contorno adiacenti ad una parete interna
   if (adiaw) call cambia_in_parete(igr)

   !  scrittura file tecplot per debugging
   call printdis(igr)

   !  verifica che tutte le celle abbiano una valida dipendenza
   !  calcola la funzione distanza
   call xcheck(igr)

   !  scrittura matrice ICC e vettore RCC
   call printcc(igr)

   !  chiamata sistema per stampa ora
   call date_and_time(ftmp,ftmp,ftmp,hms)
   write(ftmp,'(i2.2,a,i2.2,a,i2.2)') hms(5),':',hms(6),':',hms(7)
   write(*   ,'(3x,a,4x,a,i3,a)') trim(ftmp),'(igr =',igr,')'
   write(ndeb,'(3x,a,4x,a,i3,a)') trim(ftmp),'(igr =',igr,')'

!  fine loop sui reticoli
end do

!  chiude il file di debug
close(ndeb)

end program

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  inizializzazione matrice ICC e RCC
!===============================================================================
subroutine initcc
use prec
use modpar, only : stdflag,blochi,nbl,ngr
implicit none
integer (kind=I4P) :: igr,ibl

! inizializzazione "tipo" per tutte le celle
do igr=1,ngr
   do ibl=1,nbl
      blochi(ibl,igr)%cella(:,:,:)%t = stdflag
   end do
end do

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  aumenta la sovrapposizione delle griglie (riazzera icc per le celle chimera
!  confinanti con almeno una cella standard)
!===============================================================================
subroutine more_overlap(igr)
use prec
use modpar, only : innerx,innerw,stdflag, nbl, blochi
implicit none
integer (kind=I4P), allocatable, dimension(:,:,:) :: mask
integer (kind=I4P) :: igr,ibl,bd
integer (kind=I4P) :: nii,njj,nkk, i,j,k, id,jd,kd
real    (kind=R8P) :: dimchi,dimdon

do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   allocate(mask(0:nii+1,0:njj+1,0:nkk+1))

   ! inizializzazione di mask
   do k=1,nkk
      do j=1,njj
         do i=1,nii
            ! mask può essere innerx=27, innerw=28 o innerb=29
            mask(i,j,k) = blochi(ibl,igr)%cella(i,j,k)%t
         end do
      end do
   end do

   ! trasferisce la maschera sulle cornici
   do k=1,nkk
      do j=1,njj
         mask(0    ,j,k) = mask(1  ,j,k)
         mask(nii+1,j,k) = mask(nii,j,k)
      end do
   end do

   do k=1,nkk
      do i=0,nii+1
         mask(i,0    ,k) = mask(i,1  ,k)
         mask(i,njj+1,k) = mask(i,njj,k)
      end do
   end do

   do j=0,njj+1
      do i=0,nii+1
         mask(i,j,0    ) = mask(i,j,1  )
         mask(i,j,nkk+1) = mask(i,j,nkk)
      end do
   end do

   ! loop di riciclo
   do k=1,nkk
      do j=1,njj
         do i=1,nii

            ! se non è chimera, prosegue
            if (mask(i,j,k).ne.innerx) cycle

            ! se è adiacente lungo qualsiasi direzione ad una cella di parete
            ! interna (innerw), prosegue (per evitare adiacenze
            ! (cella standard)<-->(parete interna))
            if (any(mask(i-1:i+1,j-1:j+1,k-1:k+1).eq.innerw)) cycle

            ! se la sua dimensione è "molto" maggiore della prima donatrice,
            ! prosegue (per evitare sovrapposizioni eccessive)
            dimchi = blochi(ibl,igr)%cella(i,j,k)%size
            !
            bd = blochi(ibl,igr)%cella(i,j,k)%q(1)%b
            id = blochi(ibl,igr)%cella(i,j,k)%q(1)%i
            jd = blochi(ibl,igr)%cella(i,j,k)%q(1)%j
            kd = blochi(ibl,igr)%cella(i,j,k)%q(1)%k
            dimdon = blochi(bd,igr)%cella(id,jd,kd)%size
            !
            if (dimchi.gt.(1.5d0*dimdon)) cycle

            ! se tutte le celle adiacenti sono standard, resetta il tipo
            if ( (mask(i-1,j,k).eq.stdflag) .or. &
                 (mask(i+1,j,k).eq.stdflag) .or. &
                 (mask(i,j-1,k).eq.stdflag) .or. &
                 (mask(i,j+1,k).eq.stdflag) .or. &
                 (mask(i,j,k-1).eq.stdflag) .or. &
                 (mask(i,j,k+1).eq.stdflag)    ) &
                  blochi(ibl,igr)%cella(i,j,k)%t = stdflag

         end do
      end do
   end do

   deallocate(mask)

end do

end subroutine
