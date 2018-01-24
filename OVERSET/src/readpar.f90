!===============================================================================
!  preprocessore CHIMERA per XSHIP - lettura parametri
!===============================================================================
subroutine readpar
use prec
use moddef
use modpar
implicit none

integer (kind=I4P) :: i,j,k,p,q
integer (kind=I4P) :: mgr,ibl,jbl,igr,nf
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: ii0,iin,jj0,jjn,kk0,kkn
integer (kind=I4P) :: numero_pareti
character (len=80) :: str

type (point) :: va,vb

! per la lettura di body.dat
integer (kind=I4P) :: itmp(8)
type (point), allocatable, dimension(:) :: ptmp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  apertura file parametri
open(npar,file='cc.par',form='formatted',status='old')
rewind(npar)

!  nome ed apertura file reticolo
read(npar,*) fgrd
write(ndeb,'(a)') '# vim: set fdm=marker:'
write(ndeb,'(a)') '#'
write(ndeb,'(a,a)') 'file reticolo  = ',trim(fgrd)

!  nome file output
read(npar,*) fout
write(ndeb,'(a,a,a)') 'file output    = ',trim(fout),'...'

!  opzione per la scrittura dei file di reticolo
read(npar,*) lgrid

! aumenta la sovrapposizione tra reticoli
read(npar,*) mover

! cambia tipo alle celle chimera di contorno adiacenti a pareti interne
read(npar,*) adiaw

!  numero di livelli multigrid
read(npar,*) ngr,fgr
write(ndeb,'(a,i3)') 'numero livelli = ',ngr
write(ndeb,'(a,i3)') 'LIVELLO FINO   = ',fgr
mgr = 2**(ngr-1)
nr_sottoblocchi = mgr**3

!  numero di celle di cornice sulle facce dei blocchi con c.c. di adiacenza
read(npar,*) ghostadj
ghostadj = max(2,ghostadj)
write(ndeb,'(a,i3)') 'celle di cornice sulle adiacenze = ',ghostadj

!  spessore di strato limite stimato
read(npar,*) boulay
write(ndeb,'(a,g12.3)') 'spessore di strato limite (stimato) = ',boulay
if (boulay.gt.0d0) then
   boulay = boulay**2
else
   boulay = 0d0
end if

!  ampiezza della spiaggia numerica
read(npar,*) dispiaggia
write(ndeb,'(a,g12.3)') 'ampiezza della spiaggia numerica = ',dispiaggia

!  numero e priorità dei blocchi
read(npar,*) nbl

! questi andranno integrati in blocco_chimera
allocate(blogrd(nbl,ngr))
allocate(blomet(nbl,ngr))
allocate(blochi(nbl,ngr))

allocate(lv(nbl))
allocate(nlv(0:6,nbl))
allocate(blv(nbl-1,nbl))
allocate(grp(nbl))
allocate(pri(nbl))
allocate(ingombro(2,0:nr_sottoblocchi,nbl))
allocate(ingofacc(2,6,nbl))
allocate(ingofax(6,nbl))
allocate(ingofay(6,nbl))
allocate(ingofaz(6,nbl))
allocate(i_centro(nr_sottoblocchi,nbl))
allocate(j_centro(nr_sottoblocchi,nbl))
allocate(k_centro(nr_sottoblocchi,nbl))

! livelli di appartenenza dei blocchi
q = 0
do ibl=1,nbl
   read(npar,*) lv(ibl),grp(ibl),pri(ibl)
   q = max(q,lv(ibl))
end do
if (q.gt.lvmax) then
   write(*,'(a)') 'Il parametro lvmax in mod/modpar.f90 è troppo piccolo'
   write(*,'(a)') 'Va impostato almeno pari a ' // castr(q)
   write(*,'(a)') 'Dopodiché ricompila e rilancia'
   stop
end if

   ! * ri-lettura numero blocchi dal file reticolo per verifica consistenza
   ! * lettura dimensioni dei blocchi
   !
   ! se il numero di blocchi nel file reticolo (ibl) è maggiore di quello
   ! specificato nel file parametri (nbl), invece di fermarsi legge solo i
   ! primi nbl blocchi
   open(ngrd,file=fgrd,form='unformatted',status='old',position='rewind')
   read(ngrd) ibl
   if (ibl.ne.nbl) then
      write(*,'(/,9x,a,i4)') 'File reticolo inconsistente con file parametri'
      write(*,'(9x,a,i4  )') 'numero blocchi file parametri = ',nbl
      write(*,'(9x,a,i4,/)') 'numero blocchi file reticolo  = ',ibl
      if (ibl.lt.nbl) stop
   end if
   nimax = 0
   njmax = 0
   nkmax = 0
   do jbl=1,nbl
      read(ngrd) i,j,k
      blogrd(jbl,1)%ni = i
      blogrd(jbl,1)%nj = j
      blogrd(jbl,1)%nk = k
      nimax = max(nimax,i)
      njmax = max(njmax,j)
      nkmax = max(nkmax,k)
   end do
   namax = max(nimax,njmax,nkmax)
   do jbl=nbl+1,ibl
      read(ngrd) i,j,k
   end do
   close(ngrd)

!  condizioni al contorno
!  si assicura che i limiti sugli indici abbiano valori ragionevoli
!  verifica compatibilità dimensioni patches con numero reticoli
read(npar,*) npa

! alloca le matrici che dipendono dal numero di patches
allocate( bp(10,npa) )
allocate( patch_flags(4,npa) )
allocate( famw(npa),famd(npa) )
allocate( ingopatch(2,npa))

ibl = 0
numero_pareti = 0
patch_flags = .false.
nfa = 0
do p=1,npa
   read(npar,*) bp(1:10,p)

   if (bp(3,p).gt.99 .and. bp(4,p).lt.1) then
      write(*,'(/,a)') 'patch negativo adiacente ad un patch non-chimera'
      write(*,'(11(i6),/)') p,bp(:,p)
      stop
   end if

   ! contrassegna i patch di parete e ne assegna la priorità
   if (    (bp(3,p).eq.parete) & ! parete standard adiabatica
      .or. (bp(3,p).eq.inawal) & ! parete non attiva (niente forze) adiabatica
      .or. (bp(3,p).eq.isowal) & ! parete standard isoterma
      .or. (bp(3,p).eq.isoina) & ! parete non attiva isoterma
      ) then
      patch_flags(1,p) = .true.
      ! la parete appartiene ad una famiglia
      if (bp(4,p).gt.0) then
         numero_pareti = numero_pareti + 1
         ! se viene indicata solo la famiglia, aggiunge la priorità del blocco
         if (bp(4,p).lt.100) then
            jbl = bp(1,p)
            bp(4,p) = 100*pri(jbl) + bp(4,p)
         end if
         nfa = max(nfa,mod(bp(4,p),100))
      end if
   end if

   ! contrassegna i patch con c.c. alternativa di parete
   if ( (bp(3,p).eq.(offchi+parete)) .or. &
        (bp(3,p).eq.(offchi+inawal)) .or. &
        (bp(3,p).eq.(offchi+isowal)) .or. &
        (bp(3,p).eq.(offchi+isoina)) .or. &
        (bp(3,p).eq.(offgen+parete)) .or. &
        (bp(3,p).eq.(offgen+inawal)) .or. &
        (bp(3,p).eq.(offgen+isowal)) .or. &
        (bp(3,p).eq.(offgen+isoina))   )  &
        patch_flags(2,p) = .true.

   ! contrassegna i patch di inflow,outflow
   if ((bp(3,p).eq.3) .or. (bp(3,p).eq.4)) patch_flags(3,p) = .true.

   i = blogrd(bp(1,p),1)%ni
   j = blogrd(bp(1,p),1)%nj
   k = blogrd(bp(1,p),1)%nk

   bp(5,p)  = max(bp(5,p) ,0)
   bp(6,p)  = min(bp(6,p) ,i)
   bp(7,p)  = max(bp(7,p) ,0)
   bp(8,p)  = min(bp(8,p) ,j)
   bp(9,p)  = max(bp(9,p) ,0)
   bp(10,p) = min(bp(10,p),k)

   if      (bp(2,p).eq.1) then
      bp(5,p) = 0 ; bp(6,p) = 0
   else if (bp(2,p).eq.2) then
      bp(5,p) = i ; bp(6,p) = i
   else if (bp(2,p).eq.3) then
      bp(7,p) = 0 ; bp(8,p) = 0
   else if (bp(2,p).eq.4) then
      bp(7,p) = j ; bp(8,p) = j
   else if (bp(2,p).eq.5) then
      bp(9,p) = 0 ; bp(10,p) = 0
   else if (bp(2,p).eq.6) then
      bp(9,p) = k ; bp(10,p) = k
   end if

   if (any(mod(bp(5:10,p),mgr).ne.0)) then
       write(*,'(/,a)') 'dimensioni patch incompatibili con numero livelli'
       write(*,'(a,i4,a,6(i4),/)') 'p =',p,'      dim = ',bp(5:10,p)
       ibl = ibl+1
   end if
end do
if (ibl.gt.0)  stop

!  controllo consistenza: se p è adiacente a q, q deve essere adiacente a p :-)
ibl = 0
do p=1,npa
   if (bp(3,p).gt.99) then
      q = bp(4,p)
      if (bp(4,q).ne.p) then
         write(*,'(i3,a,i3,a)') p,' è adiacente a ',q,' ma non viceversa'
         ibl = ibl+1
      end if
   end if
end do
if (ibl.gt.0)  stop

! marca i patch adiacenti a patch di parete per il calcolo dello strato limite
! sulle facce chimera
do p=1,npa
   ibl = bp(1,p)
   nf  = bp(2,p)
   do q=1,npa
      if (bp(1,q).eq.ibl .and. patch_flags(1,q)) then
         if (nf.lt.3) then
            if (bp(2,q).gt.2) patch_flags(4,p) = .true.
         else if (nf.gt.4) then
            if (bp(2,q).lt.5) patch_flags(4,p) = .true.
         else
            if (bp(2,q).lt.3 .or. bp(2,q).gt.4) patch_flags(4,p) = .true.
         end if
      end if
   end do
end do

! ordina i patch di parete in base alla priorità
allocate(family2patch(numero_pareti+1,nfa))
family2patch = 0
! loop su tutte le possibili famiglie di patch di parete
do q=1,nfa
   i = 0
   do p=1,npa
      if (patch_flags(1,p) .and. (mod(bp(4,p),100).eq.q)) then
         i = i+1
         family2patch(i,q) = p
      end if
   end do
   ! ordina le righe (o colonne?) di family2patch in base alla priorità
   do j=1,i-1
   do k=1,i-1
      if (bp(4,family2patch(k,q)).lt.bp(4,family2patch(k+1,q))) then
         p = family2patch(k,q)
         family2patch(k,q) = family2patch(k+1,q)
         family2patch(k+1,q) = p
      end if
   end do
   end do
end do

!  spigoli di parete
read(npar,*) nbe
allocate (be(8,nbe))
do p=1,nbe
   read(npar,*) be(1:8,p)
end do

! lettura scatole per definizione punti di parete interni
! nbox .... numero scatole
! pbox .... priorità
! tbox .... tipo
! gbox .... gruppo di appartenenza
! ibox .... versori per definire l'orientamento
! ilex .... lunghezze lungo le direzioni dei versori
! vbox .... vertici esaedro
! sbox .... normali facce esaedro
read(npar,*) nbox

! LETTURA SCATOLE DA FILE CC.PAR
if (nbox.gt.0) then
   call alloca_scatole
   do q=1,nbox
      read(npar,*) tbox(q),ibl    ! tipo e blocco associato
      lbox(q) = lv(ibl)
      gbox(q) = grp(ibl)
      pbox(q) = pri(ibl)

      if (tbox(q).eq.1) then
         ! PARALLELEPIPEDO
         read(npar,*) obox(q)   ! estremo diagonale
         read(npar,*) ibox(q)   ! estremo diagonale
         read(npar,*) jbox(q)   ! jbox è normale (quindi, definisce) una faccia
         read(npar,*) kbox(q)   ! jbox e kbox definiscono una seconda faccia

         va = 0.5d0*(ibox(q)+obox(q))
         vb = 0.5d0*(ibox(q)-obox(q))

         obox(q) = va
         ibox(q) = jbox(q).cross.kbox(q)
         kbox(q) = ibox(q).cross.jbox(q)
         ibox(q) = ibox(q)/norma(ibox(q))
         jbox(q) = jbox(q)/norma(jbox(q))
         kbox(q) = kbox(q)/norma(kbox(q))

         ilex(q) = abs(ibox(q).dot.vb)
         jlex(q) = abs(jbox(q).dot.vb)
         klex(q) = abs(kbox(q).dot.vb)
      else if (tbox(q).eq.2) then
         ! CILINDRO (OBSOLETO)
         read(npar,*) obox(q)   ! centro
         read(npar,*) ibox(q)   ! vettore parallelo all'asse
         read(npar,*) ilex(q)   ! semi-altezza
         read(npar,*) jlex(q)   ! raggio

         ! vettore a cazzo (se è parallelo a ibox vado a Lourdes)
         va = point(0.2314d0,0.4117d0,0.7567d0)

         jbox(q) = ibox(q)+va
         kbox(q) = ibox(q).cross.jbox(q)
         jbox(q) = kbox(q).cross.ibox(q)

         ibox(q) = ibox(q)/norma(ibox(q))
         jbox(q) = jbox(q)/norma(jbox(q))
         kbox(q) = kbox(q)/norma(kbox(q))
      else if (tbox(q).eq.3) then
         ! Coppia di cilindri coassiali
         read(npar,*) va        ! un estremo dell'asse
         read(npar,*) vb        ! l'altro estremo dell'asse
         read(npar,*) jlex(q)   ! raggio cilindro interno
         read(npar,*) klex(q)   ! raggio cilindro esterno

         obox(q) = 0.5d0*(va+vb)        ! centro dell'asse
         ilex(q) = 0.5d0*norma(vb-va)   ! semi-altezza
         ibox(q) = vers0(vb-va)         ! versore parallelo all'asse
         ! tento di completare la terna
         va = ibox(q)
         if (va%x.le.va%y .and. va%x.le.va%z) then
            vb%x = 0d0
            vb%y =  va%z
            vb%z = -va%y
         else if (va%y.le.va%x .and. va%y.le.va%z) then
            vb%x =  va%z
            vb%y = 0d0
            vb%z = -va%x
         else
            vb%x =  va%y
            vb%y = -va%x
            vb%z = 0d0
         end if
         va = ibox(q).cross.vb
         vb = ibox(q).cross.va

         jbox(q) = vers0(va)
         kbox(q) = vers0(vb)
      else if (tbox(q).eq.9) then
         ! ESAEDRO
         do p=1,8
            read(npar,*) vbox(p,q)
         end do
         call esaedro(q)
      else
         write(*,'(i2,a)') tbox(q),': tipo di scatola non riconosciuto'
         stop
      end if
   end do

   ! a scanso di equivoci, scrive il file body.dat, tanto per essere
   ! sicuri che quello generato da infocc (o rotbox o qualche altro
   ! tool) sia corretto
   open(nbod,file='body.dat',form='formatted',status='unknown')
   rewind(nbod)
   write(nbod,*) nbox*8
   do q=1,nbox
      write(nbod,'(7(3(e22.12),/),3(e22.12))') (vbox(p,q),p=1,8)
   end do
   write(nbod,*) nbox
   do q=1,nbox
      write(nbod,'(9i6)') gbox(q),(q-1)*8+(/1,2,4,3,5,6,8,7/)
   end do
   close(nbod)

! LETTURA FILE BODY.DAT
else if (nbox.lt.0) then

   open(nbod,file='body.dat',form='formatted',status='old')
   rewind(nbod)
   read(nbod,*) q
   allocate(ptmp(q))
   do i=1,q
      read(nbod,*) ptmp(i)
   end do
   read(nbod,*) nbox
   if ((q/8).ne.nbox) then
      write(*,*) 'Non ho capito quante scatole vengono definite in body.dat'
      stop
   end if
   call alloca_scatole
   do q=1,nbox
      tbox(q) = 9    ! tipo di scatola: esaedro
      ! livello, gruppo, priorità e vertici della scatola
      ! NOTA: nel file body.dat letto da Xnavis la prima e terza colonna vanno
      !       cancellate
      read(nbod,'(a80)') str
      read(str,*,err=101,end=101) lbox(q),gbox(q),pbox(q),itmp(1:8)
      vbox(1,q) = ptmp(itmp(1))
      vbox(2,q) = ptmp(itmp(2))
      vbox(3,q) = ptmp(itmp(4))
      vbox(4,q) = ptmp(itmp(3))
      vbox(5,q) = ptmp(itmp(5))
      vbox(6,q) = ptmp(itmp(6))
      vbox(7,q) = ptmp(itmp(8))
      vbox(8,q) = ptmp(itmp(7))
      call esaedro(q)   ! calcola proprietà della scatola
   end do
   close(nbod)
   deallocate(ptmp)
end if

!  chiusura file parametri
close(npar)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  verifica compatibilità dimensioni blocchi con numero reticoli
i = 0
do ibl=1,nbl
   if (mod(blogrd(ibl,1)%ni,mgr).ne.0 .or. &
       mod(blogrd(ibl,1)%nj,mgr).ne.0 .or. &
       mod(blogrd(ibl,1)%nk,mgr).ne.0 ) then

       write(*,'(/,a)') 'dimensioni blocco incompatibili con numero livelli'
       write(*,'(a,i3,a,3(i5))') 'ibl =',ibl,'      ni,nj,nk = ',&
                         blogrd(ibl,1)%ni,blogrd(ibl,1)%nj,blogrd(ibl,1)%nk
       i = i+1
   end if
end do
if (i.gt.0)  stop

!  dimensioni dei blocchi dei reticoli radi
do igr=2,ngr
   do ibl=1,nbl
      blogrd(ibl,igr)%ni = blogrd(ibl,igr-1)%ni / 2
      blogrd(ibl,igr)%nj = blogrd(ibl,igr-1)%nj / 2
      blogrd(ibl,igr)%nk = blogrd(ibl,igr-1)%nk / 2
   end do
end do

! definizione di nlv e blv
! dentro blv sono ordinati i possibili blocchi donatori secondo l'ordine:
! prima quelli dello stesso livello in funzione della priorità
! poi quelli degli altri livelli in funzione della priorità
! BLOCCHI A PRIORITÀ PIÙ BASSA (SUL PROPRIO O SU ALTRI LIVELLI) NON VENGONO
! CONSIDERATI
do ibl=1,nbl
   p = 0
   nlv(0,ibl) = 0

   ! stesso livello, priorità più alta
   do jbl=1,nbl
      if (lv(ibl).eq.lv(jbl) .and. pri(jbl).gt.pri(ibl)) then
         p = p+1
         blv(p,ibl) = jbl
      end if
   end do
   nlv(1,ibl) = p

   ! stesso livello, stessa priorità
   do jbl=1,nbl
      if (lv(ibl).eq.lv(jbl) .and. pri(jbl).eq.pri(ibl) .and. ibl.ne.jbl) then
         p = p+1
         blv(p,ibl) = jbl
      end if
   end do
   nlv(2,ibl) = p

   ! stesso livello, priorità più bassa
   do jbl=1,nbl
      if (lv(ibl).eq.lv(jbl) .and. pri(jbl).lt.pri(ibl)) then
         p = p+1
         blv(p,ibl) = jbl
      end if
   end do
   nlv(3,ibl) = p

   ! livello diverso, priorità più alta
   do jbl=1,nbl
      if (lv(ibl).ne.lv(jbl) .and. pri(jbl).gt.pri(ibl)) then
         p = p+1
         blv(p,ibl) = jbl
      end if
   end do
   nlv(4,ibl) = p

   ! livello diverso, stessa priorità
   do jbl=1,nbl
      if (lv(ibl).ne.lv(jbl) .and. pri(jbl).eq.pri(ibl)) then
         p = p+1
         blv(p,ibl) = jbl
      end if
   end do
   nlv(5,ibl) = p

   ! livello diverso, priorità più bassa
   do jbl=1,nbl
      if (lv(ibl).ne.lv(jbl) .and. pri(jbl).lt.pri(ibl)) then
         p = p+1
         blv(p,ibl) = jbl
      end if
   end do
   nlv(6,ibl) = p
end do

!  stampa debug
write(ndeb,'(/,a)') 'Livello, gruppo e priorità dei blocchi {{{1'
write(ndeb,'(/,a,i3)') 'numero blocchi = ',nbl
do ibl=1,nbl
   write(ndeb,'("blocco",i4,":")',advance='no') ibl
   write(ndeb,'(5x,"lvl =",i3)',advance='no') lv(ibl)
   write(ndeb,'(5x,"grp =",i3)',advance='no') grp(ibl)
   write(ndeb,'(5x,"pri =",i3)') pri(ibl)
end do

! write(ndeb,*)
! write(ndeb,'(/,a)') "ordine di ricerca di blocchi donatori"
! do ibl=1,nbl
!    write(ndeb,'(/,a,i3.3)') "blocco #",ibl
!    do p=1,6
!       if (nlv(p,ibl).gt.nlv(p-1,ibl)) then
!          write(ndeb,'(3x,15(i4))') blv(nlv(p-1,ibl)+1:nlv(p,ibl),ibl)
!       else
!          write(ndeb,'(3x,a)') 'nessuno'
!       end if
!    end do
! end do

write(ndeb,*)
write(ndeb,'(a)') 'parametri utilizzati per i patches (tutti) {{{1'
write(ndeb,*) '  ibl   nf   icc  adiac'//&
              '  Imin  Imax  Jmin  Jmax  Kmin  Kmax  patch'
write(str,'(a,i1,a)') '(10i5,a,i',nfig(npa)+1,')'
do p=1,npa
   write(ndeb,str) bp(:,p),'     !',p
end do

write(ndeb,'(/,a)') 'ordinamento patch di parete per famiglia e priorità {{{1'
write(str,'(a,i1,a)') '(i',nfig(npa)+1,')'
do q=1,nfa
   if (family2patch(1,q).gt.0) then
      i=1
      write(ndeb,'(a,i2,a)',advance='no') 'Famiglia #',q,' :'
      do while (family2patch(i,q).gt.0)
         write(ndeb,str,advance='no') family2patch(i,q)
         i = i+1
      end do
      write(ndeb,*)
   end if
end do

if (nbe.gt.0) then
   write(ndeb,*)
   write(ndeb,'(a)') 'parametri utilizzati per gli spigoli di parete {{{1'
   write(ndeb,*) '  ibl   dir    Imin  Imax  Jmin  Jmax  Kmin  Kmax  spigolo'
   do p=1,nbe
      write(ndeb,'(8i5,a,i4)') be(:,p),'     !',p
   end do
end if

if (nbox.gt.0) then
   write(ndeb,*)
   write(ndeb,'(a,/)') 'scatole {{{1'
   do p=1,nbox
      write(ndeb,'(a,i3.3)',advance='no') '   #',p
      write(ndeb,'(3x,a,i2)',advance='no') '   lvl = ',lbox(p)
      write(ndeb,'(3x,a,i2)',advance='no') '   grp = ',gbox(p)
      write(ndeb,'(3x,a,i2)',advance='no') '   pri = ',pbox(p)
      write(ndeb,'(3x,a,i1)') '   tipo = ',tbox(p)
!     write(ndeb,'(a,3(e15.6))') '   centro :',obox(p)
!     if (tbox(p).lt.3) then
!        write(ndeb,'(a,3(e15.6))') '   vers i :',ibox(p)
!        write(ndeb,'(a,3(e15.6))') '   vers j :',jbox(p)
!        write(ndeb,'(a,3(e15.6))') '   vers k :',kbox(p)
!        write(ndeb,'(a,e15.6)') '   leng i :',ilex(p)
!        write(ndeb,'(a,e15.6)') '   leng j :',jlex(p)
!        write(ndeb,'(a,e15.6)') '   leng k :',klex(p)
!     else if (tbox(p).eq.9) then
!        do q=1,6
!           write(ndeb,'(a,i1.1,a,3(e15.6))') '   nor1 ',q,' :',sbox(1,q,p)
!           write(ndeb,'(a,i1.1,a,3(e15.6))') '   nor2 ',q,' :',sbox(2,q,p)
!           write(ndeb,'(a,i1.1,a,3(e15.6))') '   face ',q,' :',fbox(q,p)
!        end do
!     end if
   end do
end if
call printbox

write(ndeb,*)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ridefinizione delle celle di cornice a seconda del tipo di c.c. sulle facce di
! ciascun blocco
do igr=1,ngr
   do ibl=1,nbl
      blogrd(ibl,igr)%gc(:) = 2
   end do
end do
do igr=1,ngr
   do p=1,npa
      if ( (bp(3,p).gt.99) ) then
         ibl = bp(1,p)
         nf  = bp(2,p)
         q   = bp(4,p) ! patch adiacente a quello considerato
         jbl = bp(1,q) ! blocco adiacente al patch considerato
         if (bp(2,q).lt.3) then
            blogrd(ibl,igr)%gc(nf) = min(ghostadj,blogrd(jbl,igr)%ni)
         else if (bp(2,q).gt.4) then
            blogrd(ibl,igr)%gc(nf) = min(ghostadj,blogrd(jbl,igr)%nk)
         else
            blogrd(ibl,igr)%gc(nf) = min(ghostadj,blogrd(jbl,igr)%nj)
         end if
      end if
   end do
end do

do igr=1,ngr
   do ibl=1,nbl
      blomet(ibl,igr)%ni = blogrd(ibl,igr)%ni
      blomet(ibl,igr)%nj = blogrd(ibl,igr)%nj
      blomet(ibl,igr)%nk = blogrd(ibl,igr)%nk
      blomet(ibl,igr)%gc = blogrd(ibl,igr)%gc
      blochi(ibl,igr)%ni = blogrd(ibl,igr)%ni
      blochi(ibl,igr)%nj = blogrd(ibl,igr)%nj
      blochi(ibl,igr)%nk = blogrd(ibl,igr)%nk
      blochi(ibl,igr)%gc = blogrd(ibl,igr)%gc
      ! A questo punto può allocare le matrici che dipendono dalle dimensioni
      ! dei blocchi (i nodi del reticolo sono definiti, con due celle di
      ! cornice, tra -2:n+2 mentre tutto il resto, compreso le normali alle
      ! interfacce, è definito tra -1:n+2
      ii0 = blogrd(ibl,igr)%gc(1)
      iin = blogrd(ibl,igr)%ni + blogrd(ibl,igr)%gc(2)
      jj0 = blogrd(ibl,igr)%gc(3)
      jjn = blogrd(ibl,igr)%nj + blogrd(ibl,igr)%gc(4)
      kk0 = blogrd(ibl,igr)%gc(5)
      kkn = blogrd(ibl,igr)%nk + blogrd(ibl,igr)%gc(6)
      allocate(blogrd(ibl,igr)%nodo (-ii0  :iin,-jj0  :jjn,-kk0  :kkn))
      allocate(blomet(ibl,igr)%cella(-ii0+1:iin,-jj0+1:jjn,-kk0+1:kkn))
      allocate(blochi(ibl,igr)%cella(-ii0+1:iin,-jj0+1:jjn,-kk0+1:kkn))
   end do
end do

write(ndeb,'(/,a)') 'Dimensioni dei blocchi e ghost cells {{{1'
write(ndeb,*)
do ibl=1,nbl
   write(ndeb,'("blocco",i4,":",3(i5),6(i4))') ibl,&
      blogrd(ibl,1)%ni,blogrd(ibl,1)%nj,blogrd(ibl,1)%nk,blogrd(ibl,1)%gc(1:6)
end do

!  conteggio punti effettivi di calcolo (senza cornici)
p = 0
j = 0
do ibl=1,nbl
   nii = blogrd(ibl,1)%ni
   njj = blogrd(ibl,1)%nj
   nkk = blogrd(ibl,1)%nk
   ii0 = blogrd(ibl,1)%gc(1)
   iin = blogrd(ibl,1)%gc(2)
   jj0 = blogrd(ibl,1)%gc(3)
   jjn = blogrd(ibl,1)%gc(4)
   kk0 = blogrd(ibl,1)%gc(5)
   kkn = blogrd(ibl,1)%gc(6)
   p = p + nii*njj*nkk
   j = j + (ii0+nii+iin)*(jj0+njj+jjn)*(kk0+nkk+kkn)
end do

write(*,'(a,i8)')   'numero punti reticolo fitto (no cornici)   : ', p
write(*,'(a,i8,/)') 'numero punti complessivo (incluse cornici) : ', j

return

101 continue
write(*,'(/,80(a1))') ('*',i=1,80)
write(*,'(a)') '                     Errore nella lettura di body.dat'
write(*,'(a)',advance='no') 'Ricorda che il file body.dat letto da overset NON '
write(*,'(a)') 'è uguale a quello letto da'
write(*,'(a)',advance='no') 'Xnavis. In particolare, nella definizione delle '
write(*,'(a)') 'scatole, overset si aspetta:'
write(*,'(a)') 'livello, GRUPPO, priorità, GLI 8 INDICI DEI VERTICI'
write(*,'(80(a1),/)') ('*',i=1,80)
stop

end subroutine

!===============================================================================

! inutilizzata, potrebbe servire
subroutine sort(n,blv,pri)
use prec
implicit none
integer (kind=I4P), intent(in) :: n
integer (kind=I4P) :: blv(n)
integer (kind=I4P), intent(in) :: pri(*)
integer (kind=I4P) :: p,q,tmp

if (n.lt.2) return

do p=1,n-1
do q=1,n-1
   if (pri(blv(q)).lt.pri(blv(q+1))) then
      tmp = blv(q)
      blv(q) = blv(q+1)
      blv(q+1) = tmp
   end if
end do
end do

end subroutine
