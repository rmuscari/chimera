!===============================================================================
!  legge un file di output di "overset", chiede bl,i,j,k e ritorna tutte
!  le informazioni disponibili sulla cella
!===============================================================================
program getchi
use prec
use moddef
implicit none

integer (kind=I4P), parameter :: idcc = 52
integer (kind=I4P), parameter :: idgr = 53
integer (kind=I4P), parameter :: idis = 54

integer (kind=I4P), parameter :: nccnat = 19

!  nomi file
character (len=50) :: fcc,fgr,fdi
character (len=3)  :: ngr

integer (kind=I4P) :: iargc
integer (kind=I4P) :: nb,bl
integer (kind=I4P) :: m,n,i,j,k,p,t
integer (kind=I4P) :: ii0,iin,jj0,jjn,kk0,kkn
logical (kind=I1P) :: ldis
type (point) :: c
real (kind=R8P) :: d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! legge nome file dalla riga di comando
!
if (iargc().ne.1) then
   write(*,'(/,a,/)') 'Usage: getchi nome_file'
   stop
end if
call getarg(1,fcc)

! apertura file
!
open(idcc,file=fcc,status='old',form='unformatted',position='rewind')

fgr = trim(fcc)//'.grd'  ! file reticolo standard
j = index(fcc,'CC')
i = index(fcc,'.p0')
if (j.gt.0 .and. i.gt.0 .and. (j+1).lt.(i-1)) then
   fgr = fcc(1:j-1)//fcc(j+2:i-1)//'.grd'//trim(fcc(i:))
else if (i.gt.0 .and. j.lt.1) then
   fgr = fcc(1:i-1)//'.grd'//trim(fcc(i:))
else if (j.gt.0 .and. i.lt.1) then
   fgr = fcc(1:j-1)//trim(fcc(j+2:))//'.grd'
end if
open(idgr,file=trim(fgr),status='old',form='unformatted')
rewind(idgr)

i = index(fcc,'.0')
ngr = fcc(i:i+2)
fdi = 'dist'//ngr
inquire(file=fdi,exist=ldis)
if (ldis) then
   open(idis,file=fdi,status='old',form='unformatted',position='rewind')
   read(idis) (d,i=1,7)
end if

! lettura numero e dimensioni blocchi
!
read(idcc) nb
read(idgr) m
call readcheck(m,nb,1)

allocate(blochi(nb))
allocate(blogrd(nb))

do bl=1,nb
   read(idcc) blochi(bl)%ni,blochi(bl)%nj,blochi(bl)%nk,blochi(bl)%gc(1:6)
   read(idgr) blogrd(bl)%ni,blogrd(bl)%nj,blogrd(bl)%nk,blogrd(bl)%gc(1:6)
   call readcheck(blochi(bl)%ni,blogrd(bl)%ni,2)
   call readcheck(blochi(bl)%nj,blogrd(bl)%nj,2)
   call readcheck(blochi(bl)%nk,blogrd(bl)%nk,2)
   do i=1,6
      call readcheck(blochi(bl)%gc(i),blogrd(bl)%gc(i),3)
   end do
   if (ldis) then
      read(idis) i,j,k
      call readcheck(i,blochi(bl)%ni,5)
      call readcheck(j,blochi(bl)%nj,6)
      call readcheck(k,blochi(bl)%nk,7)
   end if
   ! allocazione blocchi
   ii0 = - blochi(bl)%gc(1)
   iin = blochi(bl)%ni + blochi(bl)%gc(2)
   jj0 = - blochi(bl)%gc(3)
   jjn = blochi(bl)%nj + blochi(bl)%gc(4)
   kk0 = - blochi(bl)%gc(5)
   kkn = blochi(bl)%nk + blochi(bl)%gc(6)
   allocate(blochi(bl)%cella(-ii0+1:iin,-jj0+1:jjn,-kk0+1:kkn))
   allocate(blogrd(bl)%nodo (-ii0  :iin,-jj0  :jjn,-kk0  :kkn))
end do

write(*,'(/,a)') 'file cc  : ' // trim(fcc)
write(*,'(  a)') 'file grd : ' // trim(fgr)
if (ldis) write(*,'(a)') 'file dis : ' // trim(fdi)

do bl=1,nb
   ii0 = - blochi(bl)%gc(1) + 1
   iin = blochi(bl)%ni + blochi(bl)%gc(2)
   jj0 = - blochi(bl)%gc(3) + 1
   jjn = blochi(bl)%nj + blochi(bl)%gc(4)
   kk0 = - blochi(bl)%gc(5) + 1
   kkn = blochi(bl)%nk + blochi(bl)%gc(6)
   write(*,'(7(i5))') bl,ii0,iin,jj0,jjn,kk0,kkn
end do

! lettura tipo
do bl=1,nb
   ii0 = - blochi(bl)%gc(1) + 1
   iin = blochi(bl)%ni + blochi(bl)%gc(2)
   jj0 = - blochi(bl)%gc(3) + 1
   jjn = blochi(bl)%nj + blochi(bl)%gc(4)
   kk0 = - blochi(bl)%gc(5) + 1
   kkn = blochi(bl)%nk + blochi(bl)%gc(6)
   write(*,'(a,10(i5))') 'pippo 1',bl,blochi(bl)%ni,blochi(bl)%nj,&
                          blochi(bl)%nk,blochi(bl)%gc
   read(idcc) (((blochi(bl)%cella(i,j,k)%t,&
                  i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
end do
print*,'pippo 2'

!  per le celle chimera (tipo>nccnat) lettura dei dati
do bl=1,nb
   ii0 = - blochi(bl)%gc(1) + 1
   iin = blochi(bl)%ni + blochi(bl)%gc(2)
   jj0 = - blochi(bl)%gc(3) + 1
   jjn = blochi(bl)%nj + blochi(bl)%gc(4)
   kk0 = - blochi(bl)%gc(5) + 1
   kkn = blochi(bl)%nk + blochi(bl)%gc(6)
   write(*,'(a,10(i5))') 'pippo 3',bl,blochi(bl)%ni,blochi(bl)%nj,&
                          blochi(bl)%nk,blochi(bl)%gc
   do k=kk0,kkn
      do j=jj0,jjn
         do i=ii0,iin
            if (blochi(bl)%cella(i,j,k)%t.gt.nccnat) then
               m = blochi(bl)%cella(i,j,k)%n
               allocate(blochi(bl)%cella(i,j,k)%q(m))
               read(idcc)  blochi(bl)%cella(i,j,k)%n
               read(idcc) (blochi(bl)%cella(i,j,k)%q(n)%b,&
                            blochi(bl)%cella(i,j,k)%q(n)%i,&
                            blochi(bl)%cella(i,j,k)%q(n)%j,&
                            blochi(bl)%cella(i,j,k)%q(n)%k,&
                            blochi(bl)%cella(i,j,k)%q(n)%w,n=1,m)
            end if
         end do
      end do
   end do
end do
print*,'pippo 4'

! lettura grd
!
do bl=1,nb
   ii0 = - blogrd(bl)%gc(1)
   iin = blogrd(bl)%ni + blogrd(bl)%gc(2)
   jj0 = - blogrd(bl)%gc(3)
   jjn = blogrd(bl)%nj + blogrd(bl)%gc(4)
   kk0 = - blogrd(bl)%gc(5)
   kkn = blogrd(bl)%nk + blogrd(bl)%gc(6)
   read(idgr) (((blogrd(bl)%nodo(i,j,k)%x,&
                 i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
   read(idgr) (((blogrd(bl)%nodo(i,j,k)%y,&
                 i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
   read(idgr) (((blogrd(bl)%nodo(i,j,k)%z,&
                  i=ii0,iin),j=jj0,jjn),k=kk0,kkn)
end do

! lettura dis
!
if (ldis) then
   do bl=1,nb
      iin = blochi(bl)%ni
      jjn = blochi(bl)%nj
      kkn = blochi(bl)%nk
      read(idis) (((blochi(bl)%cella(i,j,k)%dist,&
                    i=0,iin+1),j=0,jjn+1),k=0,kkn+1)
   end do
else
   do bl=1,nb
      blochi(bl)%cella(:,:,:)%dist = -1d0
   end do
end if

! chiusura file
!
close(idcc)
close(idgr)
if (ldis) close(idis)

!...............................................................................

!  loop infinito
do while (.true.)

   write(*,'(/,a,$)') 'bl,i,j,k : '
   read(*,*) bl,i,j,k

   ii0 = - blogrd(bl)%gc(1)+1
   iin = blogrd(bl)%ni + blogrd(bl)%gc(2)
   jj0 = - blogrd(bl)%gc(3)+1
   jjn = blogrd(bl)%nj + blogrd(bl)%gc(4)
   kk0 = - blogrd(bl)%gc(5)+1
   kkn = blogrd(bl)%nk + blogrd(bl)%gc(6)

   m = 0
   if (bl.lt.1 .or. bl.gt.nb) then
      write(*,*) 'Valore di bl non valido'
      write(*,'(a,i4)') '1 < bl <',nb
      bl = max(1,min(bl,nb))
      m = 1
   end if
   if (i.lt.ii0 .or. i.gt.iin) then
      write(*,'(a)') 'I supera i limiti : ' // castr(ii0) &
                     // ' <= I <='          // castr(iin)
      i  = max(0,min(i,blochi(bl)%ni+1))
      m = 1
   end if
   if (j.lt.jj0 .or. j.gt.jjn) then
      write(*,'(a)') 'J supera i limiti : ' // castr(jj0) &
                     // ' <= J <='          // castr(jjn)
      j  = max(0,min(j,blochi(bl)%nj+1))
      m = 1
   end if
   if (k.lt.kk0 .or. k.gt.kkn) then
      write(*,'(a)') 'K supera i limiti : ' // castr(kk0) &
                     // ' <= K <='          // castr(kkn)
      k  = max(0,min(k,blochi(bl)%nk+1))
      m = 1
   end if
   if (m.eq.1) then
      write(*,'(a,4(i5))') ' bl,i,j,k : ',bl,i,j,k
   end if

   t = blochi(bl)%cella(i,j,k)%t
   write(*,41) '   tipo.cella.......... :',t

   if (t.gt.nccnat) then
      m = blochi(bl)%cella(i,j,k)%n
      write(*,51) '   numero.donatori..... :',m
      write(*,'(a)') '   nr.   bd    id    jd    kd    peso'
      do p=1,m
         write(*,61) p,blochi(bl)%cella(i,j,k)%q(p)%b, &
                       blochi(bl)%cella(i,j,k)%q(p)%i, &
                       blochi(bl)%cella(i,j,k)%q(p)%j, &
                       blochi(bl)%cella(i,j,k)%q(p)%k, &
                       blochi(bl)%cella(i,j,k)%q(p)%w
      end do
   end if

   if (t.ge.21 .and. t.le.26) then
      call centrofaccia(bl,i,j,k,t-20,c)
   else
      call centrocella(bl,i,j,k,c)
   end if

   if (ldis) then
      d = blochi(bl)%cella(i,j,k)%dist
      write(*,'(a,e15.8)') '   distanza chimera-baricentro: ',d
   end if

end do

41 format(/,a,i4)
51 format(a,i4,/)
61 format(2x,i3,4(i6),f10.4)

end program

!===============================================================================
subroutine readcheck(a,b,err)
use prec
implicit none
integer (kind=I4P), intent(in) :: a,b,err
if (a.ne.b) then
   write(*,*) 'LETTURE INCONSISTENTI'
   write(*,*) 'VALORI LETTI :',a,'VS.',b
   write(*,*) 'ERRORE = ',err
   stop
end if
return
end subroutine

!===============================================================================
subroutine centrocella(b,i,j,k,c)
use prec
use moddef, only : point,blogrd, operator(+), operator(*)
implicit none
integer (kind=I4P), intent(in) :: b,i,j,k
type (point) :: c

c = blogrd(b)%nodo(i-1,j-1,k-1) &
  + blogrd(b)%nodo(i-1,j-1,k  ) &
  + blogrd(b)%nodo(i-1,j  ,k-1) &
  + blogrd(b)%nodo(i  ,j-1,k-1) &
  + blogrd(b)%nodo(i  ,j  ,k-1) &
  + blogrd(b)%nodo(i  ,j-1,k  ) &
  + blogrd(b)%nodo(i-1,j  ,k  ) &
  + blogrd(b)%nodo(i  ,j  ,k  )

c = 0.125d0*c

!  write(*,'(/,3(a,f14.8,/))') '   centro cella: ',c%x,&
!                              '                 ',c%y,&
!                              '                 ',c%z
write(*,'(/,a,e17.9,e17.9,e17.9,/)') '   centro cella:',c%x,c%y,c%z

end subroutine

!===============================================================================
subroutine centrofaccia(b,i,j,k,f,c)
use prec
use moddef, only : point,blogrd, operator(+), operator(*)
implicit none
integer (kind=I4P), intent(in) :: b,i,j,k,f
type (point) :: c
integer (kind=I4P) :: iin,jjn,kkn

iin = blogrd(b)%ni
jjn = blogrd(b)%nj
kkn = blogrd(b)%nk

select case (f)
   case (1)
      c = blogrd(b)%nodo(0 ,j-1,k-1) &
        + blogrd(b)%nodo(0 ,j-1,k  ) &
        + blogrd(b)%nodo(0 ,j  ,k-1) &
        + blogrd(b)%nodo(0 ,j  ,k  )
   case (2)
      c = blogrd(b)%nodo(iin,j-1,k-1) &
        + blogrd(b)%nodo(iin,j-1,k  ) &
        + blogrd(b)%nodo(iin,j  ,k-1) &
        + blogrd(b)%nodo(iin,j  ,k  )
   case (3)
      c = blogrd(b)%nodo(i-1,0 ,k-1) &
        + blogrd(b)%nodo(i-1,0 ,k  ) &
        + blogrd(b)%nodo(i  ,0 ,k-1) &
        + blogrd(b)%nodo(i  ,0 ,k  )
   case (4)
      c = blogrd(b)%nodo(i-1,jjn,k-1) &
        + blogrd(b)%nodo(i-1,jjn,k  ) &
        + blogrd(b)%nodo(i  ,jjn,k-1) &
        + blogrd(b)%nodo(i  ,jjn,k  )
   case (5)
      c = blogrd(b)%nodo(i-1,j-1,0 ) &
        + blogrd(b)%nodo(i-1,j  ,0 ) &
        + blogrd(b)%nodo(i  ,j-1,0 ) &
        + blogrd(b)%nodo(i  ,j  ,0 )
   case (6)
      c = blogrd(b)%nodo(i-1,j-1,kkn) &
        + blogrd(b)%nodo(i-1,j  ,kkn) &
        + blogrd(b)%nodo(i  ,j-1,kkn) &
        + blogrd(b)%nodo(i  ,j  ,kkn)
end select

c = 0.25d0*c

write(*,'(/,3(a,f14.8,/))') '   centro faccia: ',c%x,&
                            '                  ',c%y,&
                            '                  ',c%z

end subroutine
