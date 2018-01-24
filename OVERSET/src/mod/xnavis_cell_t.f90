module xnavis_cell_t

   use xnavis_prec_t
   use xnavis_point_t
   use xnavis_chimera_t

   implicit none
   private
   public :: cell

   ! cell definition 
   type :: cell

      type(point) :: vertex
      type(point) :: center

      type(point) :: sni   ! normal i-face
      type(point) :: snj   ! normal j-face
      type(point) :: snk   ! normal k-face

      type(point) :: jac(3)

      real(R8P) :: ai = 0._R8P    ! area i-face
      real(R8P) :: aj = 0._R8P    ! area j-face
      real(R8P) :: ak = 0._R8P    ! area k-face
      real(R8P) :: vol = 0._R8P   ! volume
      real(R8P) :: dist = 0._R8P  ! wall distance

      ! start chimera data

      integer(I4P) :: tycc=regular_cell
      type(chimera), allocatable :: xdata  ! Allocated only if tycc != standard

      contains
         procedure, pass(self) :: alloc    !< Allocate chimera data.
         procedure, pass(self) :: destroy  !< Destroy chimera data.
         ! procedure, pass(self) :: load_from_file !< Load from file.
         ! procedure, pass(self) :: save_into_file !< Save into file.
   endtype cell

   contains

      elemental subroutine alloc(self)
      class(cell), intent(inout) :: self !< The cell.
      !
      if (self%tycc>0) then
         call self%destroy
         allocate(self%xdata)
      endif
      end subroutine alloc

      elemental subroutine destroy(self)
      class(cell), intent(inout) :: self !< The cell.
      !
      if (allocated(self%xdata)) then
         deallocate(self%xdata)
      endif
      end subroutine destroy

end module xnavis_cell_t
