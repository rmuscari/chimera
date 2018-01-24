module moddef

   use prec, only : I1P,I2P,I4P,R4P,R8P
   implicit none
   save

   private :: I1P,I2P,I4P,R4P,R8P

   private :: point_from_real
   private :: add_point_point
   private :: sub_point_point, neg_point
   private :: mul_real_point, mul_point_real, mul_point_point
   private :: div_point_real
   private :: dotproduct
   private :: crossproduct
   private :: pler
   private :: abs_point
   private :: dist_point, dist_proj
   private :: mean_pnt, mean_pnt3d, mean_pnt2d

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
     real    (kind=R8P) :: dist  ! distanza dalla parete più vicina
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

   !  blocchi
   type (blocco_chimera),  allocatable :: blochi(:)
   type (blocco_reticolo), allocatable :: blogrd(:)

! INTERFACCE

   interface assignment (=)
     module procedure point_from_real
   end interface

   interface operator (+)
     module procedure add_point_point
   end interface

   interface operator (-)
     module procedure sub_point_point
     module procedure neg_point
   end interface

   interface operator (*)
     module procedure mul_real_point
     module procedure mul_point_real
     module procedure mul_point_point
   end interface

   interface operator (/)
     module procedure div_point_real
   end interface

   interface operator (.dot.)
     module procedure dotproduct
   end interface

   interface operator (.cross.)
     module procedure crossproduct
   end interface

   interface operator (.le.)
     module procedure pler
   end interface

   interface abs
     module procedure abs_point
   end interface

   interface dist
     module procedure dist_point
     module procedure dist_proj
   end interface

   interface mean
     module procedure mean_pnt
     module procedure mean_pnt3d
     module procedure mean_pnt2d
   end interface

contains

   !............................................................................
   ! ASSEGNAZIONE TRA TYPE DIFFERENTI
   subroutine point_from_real(a,b)
   type (point),  intent(out) :: a
   real (kind=R8P), intent(in)  :: b
   a%x = b
   a%y = b
   a%z = b
   return
   end subroutine

   !............................................................................
   ! SOMMA DI DUE VETTORI
   function add_point_point(a,b)
   type (point)             :: add_point_point
   type (point), intent(in) :: a,b
   add_point_point%x = a%x + b%x
   add_point_point%y = a%y + b%y
   add_point_point%z = a%z + b%z
   end function

   !............................................................................
   ! DIFFERENZA TRA DUE VETTORI
   function sub_point_point(a,b)
   type (point)             :: sub_point_point
   type (point), intent(in) :: a,b
   sub_point_point%x = a%x - b%x
   sub_point_point%y = a%y - b%y
   sub_point_point%z = a%z - b%z
   end function

   !............................................................................
   ! CAMBIO SEGNO DI UN VETTORE
   function neg_point(a)
   type (point)             :: neg_point
   type (point), intent(in) :: a
   neg_point%x = -a%x
   neg_point%y = -a%y
   neg_point%z = -a%z
   end function

   !............................................................................
   ! MOLTIPLICAZIONE DI UNO SCALARE PER UN VETTORE
   function mul_real_point(a,b)
   type (point)              :: mul_real_point
   real (kind=R8P), intent(in) :: a
   type (point),  intent(in) :: b
   mul_real_point%x = a * b%x
   mul_real_point%y = a * b%y
   mul_real_point%z = a * b%z
   end function

   function mul_point_real(a,b)
   type (point)              :: mul_point_real
   type (point),  intent(in) :: a
   real (kind=R8P), intent(in) :: b
   mul_point_real%x = a%x * b
   mul_point_real%y = a%y * b
   mul_point_real%z = a%z * b
   end function

   !............................................................................
   ! MOLTIPLICAZIONE DI DUE VETTORI (COMPONENTE PER COMPONENTE)
   function mul_point_point(a,b)
   type (point)             :: mul_point_point
   type (point), intent(in) :: a,b
   mul_point_point%x = a%x * b%x
   mul_point_point%y = a%y * b%y
   mul_point_point%z = a%z * b%z
   end function

   !............................................................................
   ! DIVISIONE DI UN VETTORE PER UNO SCALARE
   function div_point_real(a,b)
   type (point)              :: div_point_real
   type (point),  intent(in) :: a
   real (kind=R8P), intent(in) :: b
   div_point_real%x = a%x / b
   div_point_real%y = a%y / b
   div_point_real%z = a%z / b
   end function

   !............................................................................
   ! PRODOTTO SCALARE DI DUE VETTORI
   function dotproduct(a,b)
   real (kind=R8P)            :: dotproduct
   type (point), intent(in) :: a,b
   dotproduct = a%x * b%x + a%y * b%y + a%z * b%z
   end function

   !............................................................................
   ! PRODOTTO VETTORIALE DI DUE VETTORI
   function crossproduct(a,b)
   type (point)             :: crossproduct
   type (point), intent(in) :: a,b
   crossproduct%x = a%y * b%z - a%z * b%y
   crossproduct%y = a%z * b%x - a%x * b%z
   crossproduct%z = a%x * b%y - a%y * b%x
   end function

   !............................................................................
   ! MODULO DI UN VETTORE
   function norma(a)
   real (kind=R8P)            :: norma
   type (point), intent(in) :: a
   norma = sqrt ( a%x * a%x + a%y * a%y + a%z * a%z )
   end function

   function normq(a)
   real (kind=R8P)            :: normq
   type (point), intent(in) :: a
   normq = a%x * a%x + a%y * a%y + a%z * a%z
   end function

   ! VERSORE DI UN VETTORE: se il vettore è nullo ==> vers=NULL
   function vers0(a)
   type (point) :: vers0
   type (point), intent(in) :: a
   real (kind=R8P) :: d
   d = a%x * a%x + a%y * a%y + a%z * a%z
   if (d.gt.1.0e-20_R8P) then
      vers0 = a / sqrt(d)
   else
      vers0 = 0.0_R8P
   end if
   end function

   ! VERSORE DI UN VETTORE: se il vettore è nullo ==> vers=HUGE
   function vers1(a)
   type (point) :: vers1
   type (point), intent(in) :: a
   real (kind=R8P) :: d
   d = a%x * a%x + a%y * a%y + a%z * a%z
   if (d.gt.1.0e-20_R8P) then
      vers1 = a / sqrt(d)
   else
      vers1 = 1.0e20_R8P
   end if
   end function

   ! VERSORE DI UN VETTORE: versore orientato di un prodotto vettoriale
   function versore_orientato(a,b,c)
   type (point) :: versore_orientato
   type (point), intent(in) :: a,b,c
   type (point)  :: v
   real (kind=R8P) :: d
   v = crossproduct(a,b)
   d = dotproduct(v,c)
   if (d.gt.0.0_R8P) then
      versore_orientato = vers0(v)
   else
      versore_orientato = -vers0(v)
   end if
   end function

   !............................................................................
   ! DISTANZA TRA DUE PUNTI AL QUADRATO
   function dist_point(a,b)
   real (kind=R8P)            :: dist_point
   type (point), intent(in) :: a,b
   dist_point = (a%x-b%x)**2 + (a%y-b%y)**2 + (a%z-b%z)**2
   end function

   function dist_proj(a,b,n)
   real (kind=R8P)            :: dist_proj
   type (point), intent(in) :: a,b,n
   dist_proj = ((b-a).dot.n)**2
   end function

   !............................................................................
   ! ALTRO ......

   function abs_point(a)
   type (point) :: abs_point
   type (point), intent(in) :: a
   abs_point%x = abs(a%x)
   abs_point%y = abs(a%y)
   abs_point%z = abs(a%z)
   end function

   function pler(a,b)
   logical (kind=I1P)          :: pler
   real (kind=R8P), intent(in) :: b
   type (point),  intent(in) :: a
   pler = a%x.le.b .or. a%y.le.b .or. a%z.le.b
   end function

   function mean_pnt(n,a)
   type (point) :: mean_pnt
   integer (kind=I4P), intent(in) :: n
   type (point), intent(in) :: a(n)
   integer (kind=I4P) :: i
   type (point) :: s
   s = 0.0_R8P
   do i=1,n
      s%x = s%x+a(i)%x
      s%y = s%y+a(i)%y
      s%z = s%z+a(i)%z
   end do
   mean_pnt = s / real(n,8)
   end function

   function mean_pnt3d(a)
   type (point)             :: mean_pnt3d
   type (point), intent(in) :: a(:,:,:)
   real (kind=R8P)            :: isize
   isize = 1.0_R8P/size(a)
   mean_pnt3d%x = isize * sum(a(:,:,:)%x)
   mean_pnt3d%y = isize * sum(a(:,:,:)%y)
   mean_pnt3d%z = isize * sum(a(:,:,:)%z)
   end function

   function mean_pnt2d(a)
   type (point)             :: mean_pnt2d
   type (point), intent(in) :: a(:,:)
   real (kind=R8P)            :: isize
   isize = 1.0_R8P/size(a)
   mean_pnt2d%x = isize * sum(a(:,:)%x)
   mean_pnt2d%y = isize * sum(a(:,:)%y)
   mean_pnt2d%z = isize * sum(a(:,:)%z)
   end function

!===============================================================================

end module moddef
