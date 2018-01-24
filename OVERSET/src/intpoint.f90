!===============================================================================
!  determinazione celle donatrici per celle chimera
!===============================================================================
subroutine intpoint(igr)
use prec
use moddef
use modpar, only : nbl,pri,lv,lvmax, &
                   innerx,innerw,donmax, &
                   npa,patch_flags,boulay,boulab,&
                   blochi,blomet
implicit none
integer (kind=I4P) :: igr,ibl,jbl
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: i,j,k,p,flag(10)
logical (kind=I1P) :: test,csl
real (kind=R8P) :: d,s
type (donor)     :: md(donmax)
type (point)     :: c
integer (kind=I1P) :: lvbad(0:lvmax)

interface
   subroutine nelblocco(igr,ibl,t,c)
   use prec
   use moddef, only : point
   integer (kind=I4P) :: igr,ibl
   logical (kind=I1P) :: t
   type (point)     :: c
   end subroutine
end interface

!...............................................................................
!  stampa di controllo
write(*,'(a,i3,a)') 'inizio punti fuori s.l. (igr =',igr,')'

! punti fuori dallo strato limite
do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

!===============================================================================
!$omp parallel &
!$omp    default(none) &
!$omp    shared(igr,ibl,nbl,nii,njj,nkk,lv,pri,boulab,blochi,blomet) &
!$omp    private(i,j,k,p,d,s,c,flag,lvbad,md,test,csl)
!===============================================================================

!$omp do
   do k=1,nkk
   do j=1,njj
   do i=1,nii

      c = blomet(ibl,igr)%cella(i,j,k)%cen
      s = blochi(ibl,igr)%cella(i,j,k)%size
      d = blochi(ibl,igr)%cella(i,j,k)%dist
      csl = blochi(ibl,igr)%cella(i,j,k)%boul
      if (d.gt.boulab(ibl) .or. .not.csl) then
         flag = (/ 1,1,0,1,1,0, 0,0,0,0 /)

         ! inizializza lvbad per la marcatura delle pareti interne
         lvbad(0:lvmax) = 1_I1P

         ! cerca possibili donatori (memorizza tutti i livelli per i quali trova
         ! una cella sovrapposta ponendo lvbad = -1)
         call findonor(igr,ibl,s,flag,c,md,lvbad)

         ! cerca possibili donatori da blocchi a priorità più bassa
         ! (serve solo per non far marcare come pareti interne i punti che
         ! cascano nelle scatole a priorità più bassa)
         do jbl=1,nbl
            if (lvbad(lv(jbl)).gt.0_I1P) then
               if (pri(jbl).lt.pri(ibl)) then
                  call nelblocco(igr,jbl,test,c)
                  if (test) lvbad(lv(jbl)) = -1_I1P
               end if
            end if
         end do
         ! controlla se casca in una scatola di un livello con il quale non ha
         ! sovrapposizione
         call check_parete_interna(p,lv(ibl),pri(ibl),lvbad,c)

         if (p.ge.0) then
            ! casca nella scatola: viene marcata come parete interna
            md(1) = donor(ibl,i,j,k,real(p,R8P))
            call newchi(igr,ibl,i,j,k,innerw,1,md)
         else if (flag(10).gt.2) then
            ! non casca in nessuna scatola
            ! se ha trovato una donatrice a priorità più alta o di dimensioni
            ! più piccole diventa una cella chimera
            call newchi(igr,ibl,i,j,k,innerx,8,md)
         end if

      end if

   end do
   end do
   end do

!$omp end parallel

end do

if (boulay.lt.1d-30) return

! punti nello strato limite
!

!  stampa di controllo
write(*,'(a,i3,a)') 'inizio punti nello s.l. (igr =',igr,')'

! loop sui patches
do p=1,npa
   ! considera solo patch di parete
   if (patch_flags(1,p)) then
      call strato_limite_main(igr,p)
   end if
end do

return
end subroutine intpoint

!===============================================================================
!
!===============================================================================
subroutine nelblocco(igr,ibl,t,c)
use prec
use modpar
implicit none

integer (kind=I4P) :: igr,ibl
logical (kind=I1P) :: t
logical (kind=I1P) :: fuori_ingombro
type (point)     :: c
integer (kind=I4P) :: n,m

!...............................................................................

! inizializza t
t = .false.

!  se "c" esce dai limiti del blocco torna subito alla subroutine chiamante
if (fuori_ingombro(c,ingombro(1,0,ibl))) return

! ripeto il controllo considerando i sotto-blocchi in cui ho diviso ibl
do n=1,nr_sottoblocchi
   if (.not.fuori_ingombro(c,ingombro(1,n,ibl))) go to 29
end do
! "c" esce dai limiti di TUTTI i sotto-blocchi di ibl
return

! ha trovato un sotto-blocco nei cui limiti casca "c"
29 continue

!  conta numero di intersezioni raggio x=xc,y=yc,z>zc  con i confini del blocco
call intersec(m,igr,ibl,ingofax(1,ibl),ingofay(1,ibl),ingofaz(1,ibl),&
              ingofacc(1,1,ibl),c)
if (m.eq.0) return

t = .true.

return
end subroutine nelblocco
