module xnavis_prec_t

   implicit none
   save

   public

   private :: castr_R8P
   private :: castr_R4P
   private :: castr_I4P
   private :: castr_I1P
   private :: nfig_I4P

   integer, parameter:: R8P  = selected_real_kind(15,307)
   integer, parameter:: R4P  = selected_real_kind( 6,37)
   integer, parameter:: I8P  = selected_int_kind(18)
   integer, parameter:: I4P  = selected_int_kind(9)
   integer, parameter:: I2P  = selected_int_kind(4)
   integer, parameter:: I1P  = selected_int_kind(1)

   ! funzioni per la stampa dei vari tipi
   interface castr
      module procedure castr_R8P
      module procedure castr_R4P
      module procedure castr_I4P
      module procedure castr_I2P
      module procedure castr_I1P
   end interface

   interface nfig
      module procedure nfig_I4P
   end interface

contains

! ESEMPIO PER L'UTILIZZO DI ARGOMENTI OPZIONALI
!
! elemental function castr_R8P(no_sign,n) result(str)
! implicit none
! logical (kind=I1P), intent(IN), optional:: no_sign
! real (kind=R8P), intent(IN) :: n
! character (len=23) :: str
! write(str,'(E23.15E3)') n
! if (n>0.0_R8P) str(1:1) = '+'
! if (present(no_sign)) str=str(2:)
! end function

   elemental function castr_R8P(n) result(str)
   implicit none
   real (kind=R8P), intent(IN) :: n
   character (len=23) :: str
   write(str,'(E23.15E3)') n
   end function

   elemental function castr_R4P(n) result(str)
   implicit none
   real (kind=R4P), intent(IN) :: n
   character (len=13) :: str
   write(str,'(E13.6E2)') n
   end function

   elemental function castr_I4P(n) result(str)
   implicit none
   integer (kind=I4P), intent(IN) :: n
   character (len=12) :: str
   write(str,'(I12)') n
   str = adjustl(str)
   end function

   elemental function castr_I2P(n) result(str)
   implicit none
   integer (kind=I2P), intent(IN) :: n
   character (len=6) :: str
   write(str,'(I6)') n
   str = adjustl(str)
   end function

   elemental function castr_I1P(n) result(str)
   implicit none
   integer (kind=I1P), intent(IN) :: n
   character (len=4) :: str
   write(str,'(I4)') n
   str = adjustl(str)
   end function

   ! ritorna il numero di cifre necessario per rappresentare un intero
   elemental function nfig_I4P(n) result(fig)
   implicit none
   integer (kind=I4P), intent(IN) :: n
   integer (kind=I4P) :: fig,m
   fig = 1
   m = n/10
   do while (m.ne.0_I4P)
      fig = fig+1
      m = m/10
   end do
   if (n.lt.0_I4P) fig=fig+1
   end function

end module
