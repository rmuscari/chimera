!< Xnavis *chimera* class definition.
module xnavis_chimera
!< Xnavis *chimera* class definition.
   use xnavis_prec

   implicit none
   private
   public :: donor
   public :: chimera

   type :: donor
      !< Donor class for handling chimera donors informations.
      integer(I4P) :: b=0_I4P  !< Block index.
      integer(I4P) :: i=0_I4P  !< Index along i-rank.
      integer(I4P) :: j=0_I4P  !< Index along j-rank.
      integer(I4P) :: k=0_I4P  !< Index along k-rank.
      real(R8P)    :: w=0._R8P !< Weight of interpolation.
   endtype donor

   type :: chimera
      !< Chimera class for handling chimera informations.
      integer(I4P)             :: donors_number=0_I4P !< Numbers of donors.
      type(donor), allocatable :: donors(:)           !< Donors list.
      contains
         procedure, pass(self) :: alloc          !< Allocate chimera.
         procedure, pass(self) :: destroy        !< Destroy chimera.
         procedure, pass(self) :: realloc        !< Reallocate chimera.
         procedure, pass(self) :: load_from_file !< Load from file.
         procedure, pass(self) :: save_into_file !< Save into file.
   endtype chimera

   contains

      ! destroy chimera
      elemental subroutine destroy(self)
         class(chimera), intent(inout) :: self
         !
         if (allocated(self%donors)) deallocate(self%donors)
         self%donors_number=0_I4P
      end subroutine destroy

      ! allocate chimera
      subroutine alloc(self, donors_number)
         class(chimera), intent(inout) :: self
         integer(I4P),   intent(in)    :: donors_number
         !
         if (allocated(self%donors)) call errore('Cella chimera già allocata',1)
         allocate(self%donors(1:donors_number))
         self%donors_number = donors_number
      end subroutine alloc

      ! REALLOCATE chimera (does not stop if already allocated!)
      elemental subroutine realloc(self, donors_number)
         class(chimera), intent(inout) :: self
         integer(I4P),   intent(in)    :: donors_number
         !
         if (allocated(self%donors)) call self%destroy
         allocate(self%donors(1:donors_number))
         self%donors_number = donors_number
      end subroutine realloc

      subroutine load_from_file(self, file_unit)
         class(chimera), intent(inout) :: self
         integer(I4P),   intent(in)    :: file_unit
         integer(I4P)                  :: d
         !
         read(file_unit) self%donors_number,&
            (self%donors(d)%b,self%donors(d)%i,self%donors(d)%j,&
             self%donors(d)%k,self%donors(d)%w,d=1,self%donors_number)
      endsubroutine load_from_file

      subroutine save_into_file(self, file_unit)
         class(chimera), intent(in) :: self
         integer(I4P),   intent(in) :: file_unit
         integer(I4P)               :: d
         !
         write(file_unit) self%donors_number,&
            (self%donors(d)%b,self%donors(d)%i,self%donors(d)%j,&
             self%donors(d)%k,self%donors(d)%w,d=1,self%donors_number)
      endsubroutine save_into_file

      subroutine errore(str,ierr)
         character(*),intent(in) :: str
         integer(I4P),intent(in) :: ierr
         !
         if (ierr==0) return
         write(*,'(a)') trim(str)
         stop
      end subroutine errore

endmodule xnavis_chimera
