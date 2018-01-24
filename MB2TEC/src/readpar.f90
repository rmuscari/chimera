!===============================================================================
!  postprocessore per XNAVIS & codice compressibile - lettura parametri
!===============================================================================
subroutine readpar
use prec, only : R8P,R4P,I4P,I1P
use modpar
implicit none

interface assopz
   subroutine setstr (n,args,opz)
      use prec, only : R8P,R4P,I4P,I1P
      implicit none
      integer (kind=I4P), intent(inout) :: n
      integer (kind=I4P), intent(in)    :: args
      character (len=*),  intent(out)   :: opz
   end subroutine
   subroutine setint (n,args,opz)
      use prec, only : R8P,R4P,I4P,I1P
      implicit none
      integer (kind=I4P), intent(inout) :: n
      integer (kind=I4P), intent(in)    :: args
      integer (kind=I4P), intent(out)   :: opz
   end subroutine
end interface

character (len=80) :: opz,flog
integer (kind=I4P) :: iargc
integer (kind=I4P) :: i,ivar
logical (kind=I1P) :: test
character (len=6)  :: tvar
character (len=3)  :: tgrd

!...............................................................................
!  Inizializza variabili
fsol = ''
fret = ''
ficc = ''
fout = 'plot.plt'
work = '/tmp'

vertex       = .true.
cc           = .false.
compressible = .false.
anycc        = .false.

fb = 1   ! First Block
lb = 0   ! Last Block
nv = 0   ! Numero Variabili
is = 1   ! I Skip
js = 1   ! J Skip
ks = 1   ! K Skip

!...............................................................................
!  Cerca un file log con le ultime opzioni usate
call getenv('USER',temp)
flog = trim(work)//'/gh2'//trim(temp)// '.log'
open(nlog,file=flog,form='formatted',status='unknown',err=10)
rewind(nlog)
read(nlog,'(a80)',err=10,end=10) fsol
read(nlog,*,err=10,end=10) nv
close(nlog)

10 continue

!...............................................................................
!  lettura delle opzioni dalla riga di comando
i=0
do while (i.lt.iargc())
   i = i+1
   call getarg(i,opz)
   if      (opz.eq."-i") then
      call assopz(i,iargc(),fsol)
   else if (opz.eq."-g") then
      call assopz(i,iargc(),fret)
   else if (opz.eq."-x") then
      call assopz(i,iargc(),ficc)
      cc = .true.
   else if (opz.eq."-o") then
      call assopz(i,iargc(),fout)
   else if (opz.eq."-fb") then
      call assopz(i,iargc(),fb)
   else if (opz.eq."-lb") then
      call assopz(i,iargc(),lb)
   else if (opz.eq."-v") then
      call assopz(i,iargc(),nv)
   else if (opz.eq."-any") then
      anycc = .true.
   else if (opz.eq."-C") then
      vertex = .false.
   else if (opz.eq."-c") then
      compressible = .true.
   else if (opz.eq."-t") then
      call assopz(i,iargc(),work)
   else if (opz.eq."-s") then
      call assopz(i,iargc(),is)
      js = is
      ks = is
   else if (opz.eq."-is") then
      call assopz(i,iargc(),is)
   else if (opz.eq."-js") then
      call assopz(i,iargc(),js)
   else if (opz.eq."-ks") then
      call assopz(i,iargc(),ks)
   else
      write(*,*)
      write(*,30) '    Opzioni valide.........  default'
      write(*,30) '-i  file.di.input..........  '//trim(fsol)
      write(*,30) '-g  file.di.reticolo.......  '//trim(fret)
      write(*,30) '-x  file.icc...............  '//trim(ficc)
      write(*,30) '-o  file.di.output.........  '//trim(fout)
      write(*,40) '-fb primo.blocco...........  ',fb
      write(*,40) '-lb ultimo.blocco..........  ',lb
      write(*,40) '-v  numero.variabili.......  ',nv
      write(*,30) '-t  nome.directory.lavoro..  '//trim(work)
      write(*,40) '-s  skip.globale...........  ',is
      write(*,50) '-[ijk]s skip.[ijk].........  ',is,js,ks
      write(*,*)
      write(*,30) '    Flag utili.............   effetto'
      write(*,30) '-C  .......................   variabili nei centri cella'
      write(*,30) '-c  .......................   soluzione compressibile'
      write(*,30) '-any ......................   cc calcolato con ANY'
      stop
   endif
enddo

!...............................................................................
!  Cerca di indovinare il nome del file di c.c. (se non assegnato)
i = index(fsol,".0",back=.true.)
tgrd = fsol(i:i+2)
if (ficc.eq.'') then
   opz='CC/cc'//tgrd
   inquire(file=opz,exist=test)
   if (test) then
      ficc=opz
      go to 11
   end if
   opz='../CC/cc'//tgrd
   inquire(file=opz,exist=test)
   if (test) then
      ficc=opz
      go to 11
   end if
end if
11 continue

!...............................................................................
!  Cerca di indovinare il nome del file di griglia (se non assegnato)
if (fret.eq.'') then
   opz=trim(ficc)//'.grd'
   inquire(file=opz,exist=test)
   if (test) then
      fret=opz
      go to 13
   end if
   opz=trim(fsol)//'.grd'
   inquire(file=opz,exist=test)
   if (test) then
      fret=opz
      go to 13
   end if
end if
13 continue

!...............................................................................
!  Scrive in ``vars'' il nome delle variabili
temp = ' '
ivar = 1
if (compressible) then
   ! nome delle variabili per casi compressibili
   if (nv.eq.5) then
      temp = ' r ru rv rw e'
      ivar = 6
   else if (nv.eq.6) then
      temp = ' r1 r2 ru rv rw e'
      ivar = 7
   end if
else
   ! nome delle variabili per casi incompressibili
   if (nv.eq.1) then
      temp = ' d'
      ivar = nv+1
   else if (nv.eq.2) then
      temp = ' TKEmod TKEris'
      ivar = nv+1
   else if (nv.eq.4) then
      temp = ' u v w p'
      ivar = nv+1
   else if (nv.eq.5) then
      temp = ' TKEmod TKEris uuav vvav wwav'
      ivar = nv+1
   else if (nv.eq.7) then
      temp = ' u v w p ni_tot ni_tilde lam_2'
      ivar = nv+1
   else if (nv.eq.8) then
      temp = ' u v w p f f0 ni_tot ni_tilde'
      ivar = nv+1
   else if (nv.eq.9) then
      temp = ' u v w p f f0 ni_tot ni_tilde lam_2'
      ivar = nv+1
   else if (nv.eq.11) then
      temp = ' uav vav wav pav turbav uuav vvav wwav uvav uwav vwav'
      ivar = nv+1
   else if (nv.gt.4) then
      temp = ' u v w p'
      ivar = 5
   end if
end if

do i=ivar,nv
   write(tvar,'(a,i2.2)') ' var',i
   temp = trim(temp)//tvar
end do

vars = 'x y z' // trim(temp)
if (cc) vars = 'x y z' // trim(temp) // ' icc'

!...............................................................................
!  Stampa su file log
open(nlog,file=flog,form='formatted',status='unknown',err=20)
rewind(nlog)
write(nlog,'(a)') trim(fsol)
write(nlog,'(i4)') nv
close(nlog)

20 continue

!...............................................................................
!  Definisce il nome completo del file output
call getenv('tecp',temp)
if (temp .eq. '') then
   temp = '.'
end if

! En passant, definisce il nome del DataSet in TECPLOT
i = index(fout,'.plt')
dset = fout(1:i-1)

! if (fout(1:1).ne.'/') fout = trim(temp) // '/' // trim(fout)
if (index(fout,'/').eq.0) fout = trim(temp) // '/' // trim(fout)

!...............................................................................
!  Per la stampa di controllo
if (nv.eq.0)  fsol = ""

30 format(a)
40 format(a,i3)
50 format(a,"{",i2,",",i2,",",i2,"}")
return
end

!===============================================================================
!===============================================================================

subroutine setstr(n,args,opz)
use prec, only : R8P,R4P,I4P,I1P
implicit none
integer (kind=I4P), intent(inout) :: n
integer (kind=I4P), intent(in)    :: args
character (len=*), intent(out)  :: opz
if (n.ge.args)  return
n = n+1
call getarg(n,opz)
end subroutine

!...............................................................................
subroutine setint(n,args,opz)
use prec, only : R8P,R4P,I4P,I1P
implicit none
integer (kind=I4P), intent(inout) :: n
integer (kind=I4P), intent(in)    :: args
integer (kind=I4P), intent(out)   :: opz
character (len=50) :: tmp
if (n.ge.args)  return
n = n+1
call getarg(n,tmp)
read(tmp,*) opz
end subroutine
