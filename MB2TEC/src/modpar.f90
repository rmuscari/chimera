module modpar

   use prec, only : R8P,R4P,I4P,I1P
   use moddef

   implicit none

   private :: I1P,I4P,R4P,R8P

   integer (kind=I4P), parameter :: nlog = 51
   integer (kind=I4P), parameter :: nret = 52
   integer (kind=I4P), parameter :: nsol = 53
   integer (kind=I4P), parameter :: nicc = 54

   real (kind=R8P), parameter  :: eps = 1d-10
   real (kind=R8P) :: uns

   character (len=80) :: fsol,fret,ficc,fout,work,dset
   character (len=200) :: vars,temp

   logical (kind=I1P) :: vertex,cc,compressible,anycc

   ! nb = Numero Blocchi
   ! fb = First Block
   ! lb = Last Block
   ! nv = Numero Variabili
   ! is = I Skip
   ! js = J Skip
   ! ks = K Skip
   integer (kind=I4P) :: nb,fb,lb,nv
   integer (kind=I4P) :: is,js,ks

   integer (kind=I4P), allocatable, dimension(:)   :: ni,nj,nk
   integer (kind=I4P), allocatable, dimension(:,:) :: gc

   type (blocco_chimera),  allocatable, dimension(:)  :: blochi

end module modpar
