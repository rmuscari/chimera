!===============================================================================
!  postprocessore per XNAVIS
!===============================================================================
program mbvtec
use prec
use modpar
implicit none

character (len=30) :: form1,form2,form3
integer (kind=I4P) :: i,j,k,n
integer (kind=I4P) :: ii0,iin,jj0,jjn,kk0,kkn
logical (kind=I1P) :: test

!  Lettura parametri
call readpar

!  Apre file reticolo
call myfopen(nret,fret,'unformatted','old')

!  Apre file soluzione
if (nv.gt.0) call myfopen(nsol,fsol,'unformatted','old')

!  Apre file condizioni al contorno
if (cc) call myfopen(nicc,ficc,'unformatted','old')

!  Numero e dimensioni blocchi da file reticolo
read(nret) nb
allocate(ni(nb))
allocate(nj(nb))
allocate(nk(nb))
allocate(gc(6,nb))
do n=1,nb
   read(nret) ni(n),nj(n),nk(n),gc(1:6,n)
end do
!
allocate(blochi(nb))
do n=1,nb
   blochi(n)%ni = ni(n)
   blochi(n)%nj = nj(n)
   blochi(n)%nk = nk(n)
   blochi(n)%gc = gc(1:6,n)
   ii0 = gc(1,n)
   iin = ni(n) + gc(2,n)
   jj0 = gc(3,n)
   jjn = nj(n) + gc(4,n)
   kk0 = gc(5,n)
   kkn = nk(n) + gc(6,n)
   allocate(blochi(n)%cella(-ii0+1:iin,-jj0+1:jjn,-kk0+1:kkn))
end do

!  Controllo sull'intervallo di blocchi da visualizzare
!  deve essere:   1 <= fb <= lb <= nb
fb = max(1,fb)
if (lb.lt.fb) lb = nb
lb = min(nb,lb)

!  Stampa di controllo
write(form1,'(4(a,i1))') '(a,i',nfig(fb),',a,i',nfig(lb),',a,i',nfig(nb),')'
write(form2,'(a)') '(a,"{",i2,",",i2,",",i2,"}")'
write(form3,'(a,i1,a)') '("blocco ",i',nfig(nb),'," :",i5,i5,i5)'

write(*,'(a)') 'file.di.input................ '//trim(fsol)
write(*,'(a)') 'file.di.reticolo............. '//trim(fret)
write(*,'(a)') 'file.icc..................... '//trim(ficc)
write(*,'(a)') 'file.di.output............... '//trim(fout)
write(*,'(a,i3)') 'numero.variabili............. ',nv
write(*,'(a)') 'nome.directory.lavoro........ '//trim(work)
write(*,form2) '{I,J,K}Skip.................. ',is,js,ks
write(*,form1) 'blocchi ..................... ',fb,':',lb,' di ',nb
write(*,*)
if (vertex) then
   write(*,'(a)') 'variabili.prese.............. nei vertici'
else
   write(*,'(a)') 'variabili.prese.............. nei centri cella'
end if
write(*,'(a,l1)') 'compressibile................ ',compressible
if (cc) then
   if (anycc) then
      write(*,'(a)') 'opzione per icc ............. ANY'
   else
      write(*,'(a)') 'opzione per icc ............. ALL'
   end if
end if

! Intestazione per non stazionari (tsol = 0 per i conti stazionari),
! dimensioni blocchi da file dati e verifica consistenza con file reticolo
if (nv.gt.0) then
   read(nsol) uns
   if (uns.gt.eps) then
      write(temp,'(e13.6)') uns
      write(*,'(a,a)') 'Soluzione.al.tempo........... ',adjustl(temp)
   end if
   do n=1,nb
      read(nsol) i,j,k
      if (i.ne.ni(n) .or. j.ne.nj(n) .or. k.ne.nk(n)) then
         write(*,*) 'blocco ',n
         write(*,*) 'dimensioni reticolo: ',ni(n),nj(n),nk(n)
         write(*,*) 'file dati          : ',i,j,k
         stop
      end if
   end do
end if

if (cc) then
   read(nicc) n
   if (n.ne.nb) then
      stop 'file condizioni contorno inconsistente con file reticolo'
   end if
   do n=1,nb
      read(nicc) i,j,k
      if (i.ne.ni(n) .or. j.ne.nj(n) .or. k.ne.nk(n)) then
         stop 'file condizioni contorno inconsistente con file reticolo'
      end if
   end do
end if

!  Scrittura file tecplot
call printec

!  Uscita
close(nret)

inquire(nsol,opened=test)
if (test) close(nsol)

inquire(nicc,opened=test)
if (test) close(nicc)

if (cc) then
   write(*,'(/,a)') '             ATTENZIONE!'
   if (anycc) then
      write(*,'(a,/)') 'MB2TEC = ANY   =====>   TECPLOT = ALL'
   else
      write(*,'(a,/)') 'MB2TEC = ALL   =====>   TECPLOT = ANY'
   end if
end if

end program
