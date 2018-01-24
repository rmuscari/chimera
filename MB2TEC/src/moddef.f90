module moddef

   use prec, only : I1P,I2P,I4P,R4P,R8P
   implicit none
   save

   private :: I1P,I2P,I4P,R4P,R8P

! TIPI DERIVATI PER IL RETICOLO

   type point
     sequence
     real (kind=R8P) :: x
     real (kind=R8P) :: y
     real (kind=R8P) :: z
   end type

   type metrica
     sequence
     type (point) :: cen
     type (point) :: sni
     type (point) :: snj
     type (point) :: snk
   end type

! TIPI DERIVATI PER IL CHIMERA

   type donor
     sequence
     integer (kind=I4P) :: b
     integer (kind=I4P) :: i
     integer (kind=I4P) :: j
     integer (kind=I4P) :: k
     real    (kind=R8P) :: w
   end type

   type cella_chimera
     sequence
     integer (kind=I2P) :: t     ! tipo cella (se =0 è una cella standard)
     integer (kind=I1P) :: n     ! numero donatori
     type (donor), allocatable :: q(:)
   end type

! DEFINIZIONE BLOCCHI

   type blocco_reali
     sequence
     integer (kind=I4P) :: ni
     integer (kind=I4P) :: nj
     integer (kind=I4P) :: nk
     integer (kind=I4P) :: gc(6)
     real (kind=R8P), allocatable :: cella(:,:,:)
   end type

   type blocco_reticolo
     sequence
     integer (kind=I4P) :: ni
     integer (kind=I4P) :: nj
     integer (kind=I4P) :: nk
     integer (kind=I4P) :: gc(6)
     type (point), allocatable :: nodo(:,:,:)
   end type

   type blocco_metrica
     sequence
     integer (kind=I4P) :: ni
     integer (kind=I4P) :: nj
     integer (kind=I4P) :: nk
     integer (kind=I4P) :: gc(6)
     type (metrica), allocatable :: cella(:,:,:)
   end type

   type blocco_chimera
     sequence
     integer (kind=I4P) :: ni
     integer (kind=I4P) :: nj
     integer (kind=I4P) :: nk
     integer (kind=I4P) :: gc(6)
     type (cella_chimera), allocatable :: cella(:,:,:)
   end type

end module moddef
