
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

