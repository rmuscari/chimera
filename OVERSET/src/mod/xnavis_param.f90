!===============================================================================
!  preprocessore CHIMERA per XSHIP -  modulo con tutti i parametri
!===============================================================================
module xnavis_param

   use xnavis_prec

   implicit none
   public

   ! parametri 
   integer (kind=I4P), parameter :: donmax = 8       ! nr max donatori
   integer (kind=I4P), parameter :: nccnat = 19      ! nr max c.c. naturali
   integer (kind=I4P), parameter :: lvmax  = 100     ! nr mx livelli

   integer (kind=I4P), parameter :: npar = 51
   integer (kind=I4P), parameter :: ngrd = 52
   integer (kind=I4P), parameter :: nout = 53
   integer (kind=I4P), parameter :: ndeb = 54
   integer (kind=I4P), parameter :: nbod = 55

   ! etichette per celle nulle o chimera

   ! cella standard
   integer (kind=I4P), parameter :: regular_cell = 0

   ! parete adiabatica
   integer (kind=I4P), parameter :: adiabatic_wall_cell = 1

   ! parete adiabatica senza calcolo forze
   integer (kind=I4P), parameter :: fake_adiabatic_wall_cell = 11

   ! parete isoterma
   integer (kind=I4P), parameter :: isothermal_wall_cell = 12

   ! parete isoterma senza calcolo forze
   integer (kind=I4P), parameter :: fake_isothermal_cell = 13

   ! offset per c.c. chimera (21,..,26, celle chimera di contorno
   ! con donatori per centro faccia)
   integer (kind=I4P), parameter :: offchi = 20

   ! cella chimera standard (interna ad un blocco di calcolo)
   integer (kind=I4P), parameter :: inner_chimera_cell = 27 

   ! cella di parete interna
   integer (kind=I4P), parameter :: inner_wall_cell = 28

   ! cella di spiaggia
   integer (kind=I4P), parameter :: beach_cell = 29

   ! offset per c.c. chimera (41,..,46, celle chimera di contorno
   ! con donatori per centro cella)
   integer (kind=I4P), parameter :: offgen = 40

   ! offset per c.c. di adiacenza
   integer (kind=I4P), parameter :: offbiu = 60

   ! estrapolazione per corpi in movimento
   integer (kind=I4P), parameter :: estrap = 70

   ! cella di spigolo
   integer (kind=I4P), parameter :: xedge = 80

end module

