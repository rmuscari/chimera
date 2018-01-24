!===============================================================================
!  preprocessore CHIMERA per XSHIP -  modulo con tutti i parametri
!===============================================================================
module xnavis_param

   use xnavis_prec_t

   implicit none
   private

   ! parametri 
   integer (kind=I4P), parameter :: donmax = 8       ! nr max donatori
   integer (kind=I4P), parameter :: nccnat = 19      ! nr max c.c. naturali
   integer (kind=I4P), parameter :: lvmax  = 100     ! nr mx livelli

   integer (kind=I4P), parameter :: npar = 51
   integer (kind=I4P), parameter :: ngrd = 52
   integer (kind=I4P), parameter :: nout = 53
   integer (kind=I4P), parameter :: ndeb = 54
   integer (kind=I4P), parameter :: nbod = 55

   !  etichette per celle nulle o chimera
   integer (kind=I4P), parameter :: regular_cell  = 0  ! cella standard
   integer (kind=I4P), parameter :: adiawall_cell = 1  ! parete adiabatica
   integer (kind=I4P), parameter :: fakewall_cell = 11 ! parete adiab inerte
   integer (kind=I4P), parameter :: isotwall_cell = 12 ! parete isoterma
   integer (kind=I4P), parameter :: fakeisot_cell = 13 ! parete isoterma inerte
   integer (kind=I4P), parameter :: offchi  = 20  ! offset per c.c. chimera
                                       ! (21,..,26, celle chimera di contorno
                                       ! con donatori per centro faccia)
   integer (kind=I4P), parameter :: inner_chimera_cell = 27 
   integer (kind=I4P), parameter :: inner_wall_cell    = 28  ! parete interna
   integer (kind=I4P), parameter :: beach_cell         = 29  ! spiaggia
   integer (kind=I4P), parameter :: offgen  = 40  ! offset per c.c. chimera
                                       ! (41,..,46, celle chimera di contorno
                                       ! con donatori per centro cella)
   integer (kind=I4P), parameter :: offbiu  = 60  ! offset per c.c. di adiacenza
   integer (kind=I4P), parameter :: estrap  = 70  ! estrapolazione per corpi in
                                                  ! movimento
   integer (kind=I4P), parameter :: xedge   = 80  ! cella di spigolo

   !  ngr ............... numero livelli multigrid
   !  fgr ............... livello più fino calcolato
   !  nbl ............... numero blocchi
   !  ingombro .......... ingombro massimo di blocchi e sottoblocchi
   !  i_centro
   !  j_centro
   !  k_centro .......... indici del centro di blocchi e sottoblocchi
   !  ingofacc .......... ingombro massimo delle facce dei blocchi
   !  ingofa{x,y,z} ..... valore massimo della componente {x,y,z} delle normali
   !                      alle faccette per ciascun confine di blocco
   !  ghostadj .......... numero di celle di cornice sulle facce
   !                      con c.c. di adiacenza
   integer (kind=I4P) :: ngr,fgr,nbl
   integer (kind=I4P) :: nr_sottoblocchi
   integer (kind=I4P), allocatable, dimension(:,:) :: i_centro,j_centro,k_centro
   type (point), allocatable, dimension(:,:,:)  :: ingombro,ingofacc
   real (kind=R8P), allocatable, dimension(:,:)  :: ingofax,ingofay,ingofaz
   type (point), allocatable, dimension(:,:)  :: ingopatch
   integer (kind=I4P) :: ghostadj
   type (blocco_chimera),  allocatable, dimension(:,:)  :: blochi
   type (blocco_metrica),  allocatable, dimension(:,:)  :: blomet
   type (blocco_reticolo), allocatable, dimension(:,:)  :: blogrd

   !  lv  = livello del blocco
   !  pri = priorità
   !  grp = gruppo
   integer (kind=I4P), allocatable, dimension(:) :: lv,pri,grp
   integer (kind=I4P), allocatable, dimension(:,:) :: nlv,blv

   !  npa = numero patch
   !  bp = Boundary Patches (blocco, faccia, tipo, jolly, indici)
   !  famw,famd = vettori di appoggio
   !  family2patch = prende un contatore e una famiglia e tira fuori un patch
   !  patch_flags = parete?, perforante?, ...
   !  blt = boundary layer thickness, variabile di lavoro
   integer (kind=I4P) :: npa,nfa
   integer (kind=I4P), allocatable, dimension(:,:) :: bp
   integer (kind=I4P), allocatable, dimension(:,:) :: family2patch
   logical (kind=I1P), allocatable, dimension(:,:) :: patch_flags
   integer (kind=I4P), allocatable, dimension(:) :: famw
   real (kind=R8P), allocatable, dimension(:) :: famd
   real (kind=R8P), allocatable, dimension(:) :: boulap,boulaf,boulab
   real (kind=R8P) :: blt

   !  be = Boundary Edges (blocco, direzione (1=i, 2=j, 3=k), indici)
   !  nbe = numero spigoli di parete
   integer (kind=I4P) :: nbe
   integer (kind=I4P), allocatable, dimension(:,:) :: be

   !  scatole per individuazione pareti interne
   ! nbox ............. numero scatole
   ! tbox ............. tipo
   ! lbox ............. livello
   ! pbox ............. priorità
   ! gbox ............. gruppo di appartenenza
   ! ilex,jlex,klex ... dimensioni
   ! obox ............. centro della scatola
   ! ibox,jbox,kbox ... versori
   ! vbox ............. vertici esaedro
   ! fbox ............. centri facce esaedro
   integer (kind=I4P) :: nbox
   integer (kind=I4P), allocatable, dimension(:) :: tbox,lbox,pbox,gbox
   real (kind=R8P), allocatable, dimension(:) :: ilex,jlex,klex
   type (point), allocatable, dimension(:) :: obox,ibox,jbox,kbox
   type (point), allocatable, dimension(:,:) :: vbox,fbox
   type (point), allocatable, dimension(:,:,:) :: sbox

   ! spessore strato limite
   real (kind=R8P) :: boulay

   ! ampiezza della spiaggia
   real (kind=R8P) :: dispiaggia

end module xnavis_param

