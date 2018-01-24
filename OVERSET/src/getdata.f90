!===============================================================================
! routine per l'estrazione di informazioni dalle aree di buffer
!===============================================================================

!===============================================================================
!  estrae il valore di icc nella cella i,j,k del blocco ibl del reticolo igr
!===============================================================================
subroutine geticc(igr,ibl,i,j,k,ic)
use prec
use modpar, only : blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
integer (kind=I4P) :: ic
ic = blochi(ibl,igr)%cella(i,j,k)%t
end subroutine

!===============================================================================
!  assegna il valore di icc alla cella i,j,k del blocco ibl del reticolo igr
!===============================================================================
subroutine seticc(igr,ibl,i,j,k,ic)
use prec
use modpar, only : blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k,ic
blochi(ibl,igr)%cella(i,j,k)%t = ic
end subroutine

!===============================================================================
!  torna il numero di patch di parete più vicino
!===============================================================================
subroutine getwal(igr,ibl,i,j,k,p)
use prec
use modpar, only : blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
integer (kind=I4P) :: p
p = blochi(ibl,igr)%cella(i,j,k)%wall
end subroutine

!===============================================================================
!  valore della distanza da parete della cella i,j,k del blocco ibl
!  del reticolo igr
!===============================================================================
subroutine getdis(igr,ibl,i,j,k,d)
use prec
use modpar, only : blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
real (kind=R8P) :: d
d = blochi(ibl,igr)%cella(i,j,k)%dist
end subroutine

!===============================================================================
!  calcola una dimensione caratteristica della cella
!===============================================================================
subroutine getsize(igr,ibl,i,j,k,d)
use prec
use modpar, only : blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
real (kind=R8P) :: d
d = blochi(ibl,igr)%cella(i,j,k)%size
end subroutine

!===============================================================================
!  ritorna il vertice i,j,k del blocco ibl del reticolo igr
!===============================================================================
subroutine getvertex(igr,ibl,i,j,k,c)
use prec
use moddef, only : point
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
type (point) :: c
c = blogrd(ibl,igr)%nodo(i,j,k)
end subroutine

!===============================================================================
!  calcola il centro della cella i,j,k del blocco ibl del reticolo igr
!===============================================================================
subroutine getcenter(igr,ibl,i,j,k,c)
use prec
use moddef, only : point
use modpar, only : blomet
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
type (point) :: c
c = blomet(ibl,igr)%cella(i,j,k)%cen
end subroutine

!===============================================================================
!  estrae i 4 vertici dell'interfaccia tra due celle
!  nota: per nf=1 (analogamente per nf=3,5) la "i" in ingresso è "0"
!        per nf=2 (analogamente per nf=4,6) la "i" in ingresso è "n+1"
!===============================================================================
subroutine getintface(igr,ibl,i,j,k,nf,w)
use prec
use moddef, only : point
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k,nf
type (point) :: w(2,2)

select case (nf)
   case (1)
      w(1,1) = blogrd(ibl,igr)%nodo(i,j-1,k-1)
      w(2,1) = blogrd(ibl,igr)%nodo(i,j  ,k-1)
      w(1,2) = blogrd(ibl,igr)%nodo(i,j-1,k  )
      w(2,2) = blogrd(ibl,igr)%nodo(i,j  ,k  )
   case (2)
      ! vedi nota in testa alla routine
      w(1,1) = blogrd(ibl,igr)%nodo(i-1,j-1,k-1)
      w(2,1) = blogrd(ibl,igr)%nodo(i-1,j  ,k-1)
      w(1,2) = blogrd(ibl,igr)%nodo(i-1,j-1,k  )
      w(2,2) = blogrd(ibl,igr)%nodo(i-1,j  ,k  )
   case (3)
      w(1,1) = blogrd(ibl,igr)%nodo(i-1,j,k-1)
      w(2,1) = blogrd(ibl,igr)%nodo(i  ,j,k-1)
      w(1,2) = blogrd(ibl,igr)%nodo(i-1,j,k  )
      w(2,2) = blogrd(ibl,igr)%nodo(i  ,j,k  )
   case (4)
      ! vedi nota in testa alla routine
      w(1,1) = blogrd(ibl,igr)%nodo(i-1,j-1,k-1)
      w(2,1) = blogrd(ibl,igr)%nodo(i  ,j-1,k-1)
      w(1,2) = blogrd(ibl,igr)%nodo(i-1,j-1,k  )
      w(2,2) = blogrd(ibl,igr)%nodo(i  ,j-1,k  )
   case (5)
      w(1,1) = blogrd(ibl,igr)%nodo(i-1,j-1,k)
      w(2,1) = blogrd(ibl,igr)%nodo(i  ,j-1,k)
      w(1,2) = blogrd(ibl,igr)%nodo(i-1,j  ,k)
      w(2,2) = blogrd(ibl,igr)%nodo(i  ,j  ,k)
   case (6)
      ! vedi nota in testa alla routine
      w(1,1) = blogrd(ibl,igr)%nodo(i-1,j-1,k-1)
      w(2,1) = blogrd(ibl,igr)%nodo(i  ,j-1,k-1)
      w(1,2) = blogrd(ibl,igr)%nodo(i-1,j  ,k-1)
      w(2,2) = blogrd(ibl,igr)%nodo(i  ,j  ,k-1)
   case default
      stop 'Errore in getintface'
end select

end subroutine

!===============================================================================
! Calcola il centro della faccia della cella i,j,k CHE SI TROVA SULLA
! CORNICE nf DEL BLOCCO ibl (assume che la cella sia una cella di cornice)
! Esempio: se nf=1 ignora il valore di "i" e calcola la media dei vertici
!          del reticolo con indici (0,j-1:j,k-1:k)
!          se nf=4 ignora il valore di "j" e calcola la media dei vertici
!          del reticolo con indici (i-1:i,nj,k-1:k)
!===============================================================================
subroutine centro_faccia_cornice(igr,ibl,i,j,k,nf,c)
use prec
use moddef, only : point,mean
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k,nf
integer (kind=I4P) :: nii,njj,nkk
type (point) :: c
nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk
select case (nf)
   case (1) ; c = mean(blogrd(ibl,igr)%nodo(0    ,j-1:j,k-1:k))
   case (2) ; c = mean(blogrd(ibl,igr)%nodo(nii  ,j-1:j,k-1:k))
   case (3) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,0    ,k-1:k))
   case (4) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,njj  ,k-1:k))
   case (5) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1:j,0    ))
   case (6) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1:j,nkk  ))
end select
end subroutine

!===============================================================================
! Calcola il centro della faccia nf della cella i,j,k
!===============================================================================
subroutine centro_faccia_cella(igr,ibl,i,j,k,nf,c)
use prec
use moddef, only : point,mean
use modpar, only : blogrd
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k,nf
type (point) :: c
select case (nf)
   case (1) ; c = mean(blogrd(ibl,igr)%nodo(i-1  ,j-1:j,k-1:k))
   case (2) ; c = mean(blogrd(ibl,igr)%nodo(i    ,j-1:j,k-1:k))
   case (3) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1  ,k-1:k))
   case (4) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,j    ,k-1:k))
   case (5) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k-1  ))
   case (6) ; c = mean(blogrd(ibl,igr)%nodo(i-1:i,j-1:j,k    ))
end select
end subroutine
