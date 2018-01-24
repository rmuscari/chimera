!< Xnavis *point* class definition
module xnavis_point_t

   use xnavis_prec_t

   implicit none
   private
   public :: point

   real(R8P), parameter :: smallR8P = tiny(1._R8P)

   type point
      real(R8P) :: x = 0._R8P
      real(R8P) :: y = 0._R8P
      real(R8P) :: z = 0._R8P
      contains
         generic :: assignment(=)     => point_from_real
         generic :: operator(+)       => add_point_point
         generic :: operator(-)       => sub_point_point,neg_point
         generic :: operator(*)       => mul_real_point,mul_point_real,&
                                         mul_point_point
         generic :: operator(/)       => div_point_real
         generic :: operator(.dot.)   => dotproduct
         generic :: operator(.cross.) => crossproduct

         procedure :: point_from_real
         procedure :: add_point_point
         procedure :: sub_point_point
         procedure :: neg_point
         procedure, pass(b) :: mul_real_point
         procedure :: mul_point_real
         procedure :: mul_point_point
         procedure :: div_point_real
         procedure :: dotproduct
         procedure :: crossproduct
         procedure :: norm
         procedure :: sq_norm
         procedure :: vers0
         procedure :: vers1
   endtype point

   contains

   !..............................
   !..............................
   subroutine point_from_real(a,b)
   class(point), intent(inout) :: a
   real(R8P), intent(in)  :: b
   a%x = b
   a%y = b
   a%z = b
   end subroutine

   !............................
   !............................
   function add_point_point(a,b)
   class (point), intent(in) :: a,b
   type(point) :: add_point_point
   add_point_point%x = a%x + b%x
   add_point_point%y = a%y + b%y
   add_point_point%z = a%z + b%z
   end function

   !............................
   !............................
   function sub_point_point(a,b)
   class(point), intent(in) :: a,b
   type(point) :: sub_point_point
   sub_point_point%x = a%x - b%x
   sub_point_point%y = a%y - b%y
   sub_point_point%z = a%z - b%z
   end function

   function neg_point(a)
   class(point), intent(in) :: a
   type(point) :: neg_point
   neg_point%x = -a%x
   neg_point%y = -a%y
   neg_point%z = -a%z
   end function

   !...........................
   !...........................
   function mul_real_point(a,b)
   real(R8P), intent(in) :: a
   class(point), intent(in) :: b
   type(point) :: mul_real_point
   mul_real_point%x = a * b%x
   mul_real_point%y = a * b%y
   mul_real_point%z = a * b%z
   end function

   function mul_point_real(a,b)
   class(point), intent(in) :: a
   real(R8P), intent(in) :: b
   type(point) :: mul_point_real
   mul_point_real%x = a%x * b
   mul_point_real%y = a%y * b
   mul_point_real%z = a%z * b
   end function

   function mul_point_point(a,b)
   class(point), intent(in) :: a,b
   type(point) :: mul_point_point
   mul_point_point%x = a%x * b%x
   mul_point_point%y = a%y * b%y
   mul_point_point%z = a%z * b%z
   end function

   !...........................
   !...........................
   function div_point_real(a,b)
   class(point), intent(in) :: a
   real(R8P), intent(in) :: b
   type(point) :: div_point_real
   div_point_real%x = a%x / b
   div_point_real%y = a%y / b
   div_point_real%z = a%z / b
   end function

   !.......................
   !.......................
   function dotproduct(a,b)
   class(point), intent(in) :: a,b
   real(R8P) :: dotproduct
   dotproduct = a%x * b%x + a%y * b%y + a%z * b%z
   end function

   function crossproduct(a,b)
   class(point), intent(in) :: a,b
   type(point) :: crossproduct
   crossproduct%x = a%y * b%z - a%z * b%y
   crossproduct%y = a%z * b%x - a%x * b%z
   crossproduct%z = a%x * b%y - a%y * b%x
   end function

   !................
   !................
   function norm(a)
   class(point), intent(in) :: a
   real(R8P) :: norm
   norm = sqrt ( a%x * a%x + a%y * a%y + a%z * a%z )
   end function

   function sq_norm(a)
   class(point), intent(in) :: a
   real(R8P) :: sq_norm
   sq_norm = a%x * a%x + a%y * a%y + a%z * a%z
   end function

   ! VERSORE DI UN VETTORE: se la norma è nulla ==> vers=NULL
   function vers0(a)
   class(point), intent(in) :: a
   type(point) :: vers0
   real(R8P) :: d
   d = sqrt(a%x * a%x + a%y * a%y + a%z * a%z)
   if (d.lt.smallR8P) then
      vers0 = 0._R8P
      return
   end if
   vers0 = a / d
   end function

   ! VERSORE DI UN VETTORE: se la norma è nulla usa smallR8P
   function vers1(a)
   class(point), intent(in) :: a
   type(point) :: vers1
   real(R8P) :: d
   d = sqrt(a%x * a%x + a%y * a%y + a%z * a%z)
   if (d.lt.smallR8P) d = smallR8P
   vers1 = a / d
   end function

   ! VERSORE DI UN VETTORE: versore orientato di un prodotto vettoriale
   function versore_orientato(a,b,c)
   class(point), intent(in) :: a,b,c
   type(point) :: versore_orientato
   type(point)  :: v
   real(R8P) :: d
   v = crossproduct(a,b)
   d = dotproduct(v,c)
   if (d.gt.0._R8P) then
      versore_orientato = vers0(v)
   else
      versore_orientato = -vers0(v)
   end if
   end function

   ! DISTANZA TRA DUE PUNTI AL QUADRATO
   function dist_point(a,b)
   class(point), intent(in) :: a,b
   real(R8P) :: dist_point
   dist_point = (a%x-b%x)**2 + (a%y-b%y)**2 + (a%z-b%z)**2
   end function

   function dist_proj(a,b,n)
   class(point), intent(in) :: a,b,n
   real(R8P) :: dist_proj
   dist_proj = ((b-a).dot.n)**2
   end function

   function abs_point(a)
   class(point), intent(in) :: a
   type(point) :: abs_point
   abs_point%x = abs(a%x)
   abs_point%y = abs(a%y)
   abs_point%z = abs(a%z)
   end function

   function pler(a,b)
   class(point), intent(in) :: a
   real(R8P), intent(in) :: b
   logical (I1P) :: pler
   pler = a%x.le.b .or. a%y.le.b .or. a%z.le.b
   end function

   function mean_pnt(n,a)
   integer(I4P), intent(in) :: n
   class(point), intent(in) :: a(n)
   type(point) :: mean_pnt
   integer(I4P) :: i
   type(point) :: s
   s = 0.0_R8P
   do i=1,n
      s%x = s%x+a(i)%x
      s%y = s%y+a(i)%y
      s%z = s%z+a(i)%z
   end do
   mean_pnt = s / real(n,8)
   end function

   function mean_pnt3d(a)
   class(point), intent(in) :: a(:,:,:)
   real(R8P) :: isize
   type(point) :: mean_pnt3d
   isize = 1.0_R8P/size(a)
   mean_pnt3d%x = isize * sum(a(:,:,:)%x)
   mean_pnt3d%y = isize * sum(a(:,:,:)%y)
   mean_pnt3d%z = isize * sum(a(:,:,:)%z)
   end function

   function mean_pnt2d(a)
   class(point), intent(in) :: a(:,:)
   real(R8P) :: isize
   type(point) :: mean_pnt2d
   isize = 1.0_R8P/size(a)
   mean_pnt2d%x = isize * sum(a(:,:)%x)
   mean_pnt2d%y = isize * sum(a(:,:)%y)
   mean_pnt2d%z = isize * sum(a(:,:)%z)
   end function

!........................
!........................
end module xnavis_point_t
