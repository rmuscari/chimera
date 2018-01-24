! vim: set fdm=marker:
!===============================================================================
!  assegnazione c.c. sugli spigoli {{{1
!===============================================================================
subroutine calspig(igr)
use prec
use moddef
use modpar, only : nccnat,offbiu,blochi,xedge,donmax,nbl
implicit none
integer (kind=I4P), intent(in) :: igr
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: s(14,12)
integer (kind=I4P) :: ibl,i,j,k,n
integer (kind=I4P) :: ii,jj,kk,id,jd,kd,of,oi,oj,ok,ois,ojs,oks
integer (kind=I4P) :: pic,sic,tcc,nd
real (kind=R8P)    :: w1,w2
type (donor) :: md(donmax)

!  stampa di controllo
write(*,'(a,i3,a)') 'inizio calspig (igr =',igr,')'

!  i punti di spigolo sono SEMPRE chimera
tcc = xedge

do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   !............................................................................
   !  indici pertinenti ai 12 spigoli; un po' criptico :-)
   !  estremi loop .. si commenta da sé
   !  ia,ja,ka ...... offset rispetto (i,j,k) della cella adiacente (primaria)
   !  ib,jb,kb ...... offset rispetto (i,j,k) della cella adiacente (secondaria)
   !  fa,fb ......... numero di faccia attraverso cui applicare le b.c.
   !
   !                   estremi loop               ia,ja,ka  ib,jb,kb  fa,fb
   !            <----------------------------->   <------>  <------>  <--->
   s(:,1)  = (/ 1,nii, 0,0, 0,0,                   0, 1, 0,  0, 0, 1,  5, 3 /)
   s(:,2)  = (/ 1,nii, njj+1,njj+1, 0,0,           0,-1, 0,  0, 0, 1,  5, 4 /)
   s(:,3)  = (/ 1,nii, 0,0, nkk+1,nkk+1,           0, 1, 0,  0, 0,-1,  6, 3 /)
   s(:,4)  = (/ 1,nii, njj+1,njj+1, nkk+1,nkk+1,   0,-1, 0,  0, 0,-1,  6, 4 /)
   s(:,5)  = (/ 0,0, 1,njj, 0,0,                   1, 0, 0,  0, 0, 1,  5, 1 /)
   s(:,6)  = (/ nii+1,nii+1, 1,njj, 0,0,          -1, 0, 0,  0, 0, 1,  5, 2 /)
   s(:,7)  = (/ 0,0, 1,njj, nkk+1,nkk+1,           1, 0, 0,  0, 0,-1,  6, 1 /)
   s(:,8)  = (/ nii+1,nii+1, 1,njj, nkk+1,nkk+1,  -1, 0, 0,  0, 0,-1,  6, 2 /)
   s(:,9)  = (/ 0,0, 0,0, 1,nkk,                   1, 0, 0,  0, 1, 0,  3, 1 /)
   s(:,10) = (/ nii+1,nii+1, 0,0, 1,nkk,          -1, 0, 0,  0, 1, 0,  3, 2 /)
   s(:,11) = (/ 0,0, njj+1,njj+1, 1,nkk,           1, 0, 0,  0,-1, 0,  4, 1 /)
   s(:,12) = (/ nii+1,nii+1, njj+1,njj+1, 1,nkk,  -1, 0, 0,  0,-1, 0,  4, 2 /)
   !............................................................................

   !  loop sui 12 spigoli del blocco
   !
   do n=1,12

      do k=s(5,n),s(6,n)
         do j=s(3,n),s(4,n)
            do i=s(1,n),s(2,n)

               oi = -s(7,n)-s(10,n)
               oj = -s(8,n)-s(11,n)
               ok = -s(9,n)-s(12,n)

               ois = oi
               if (ois.eq.0) ois = 1
               ojs = oj
               if (ojs.eq.0) ojs = 1
               oks = ok
               if (oks.eq.0) oks = 1

               !  tipo c.c. sulle due facce che insistono sullo spigolo
               !
               pic = blochi(ibl,igr)%cella(i+s( 7,n),j+s( 8,n),k+s( 9,n))%t
               sic = blochi(ibl,igr)%cella(i+s(10,n),j+s(11,n),k+s(12,n))%t
               if (pic.lt.1 .or. sic.lt.1) then
                  stop 'Errore in calspig: icc sulle facce deve essere > 0'
               end if

               ! scommentare per dare la priorità alle adiacenze
               ! if (pic.gt.offbiu .or. sic.gt.offbiu) cycle

               ! c.c. naturali
               ! il punto viene definito come chimera con un solo donatore MA ...
               ! ... al posto del blocco (ibl) scrivo il numero di faccia ...
               ! ... al posto del peso (1.0) scrivo il tipo di c.c. naturale
               ! NOTA: gli indici del donatore non sono quelli del punto da dove
               !       viene presa la condizione al contorno!!
               !
               if (pic.le.nccnat .or. sic.le.nccnat) then

                  nd  = 1
                  if (sic.gt.pic .or. (sic.eq.5 .and. pic.le.nccnat)) then
                     of = s(13,n)
                     w1 = pic
                  else
                     of = s(14,n)
                     w1 = sic
                  end if

                  do kk=k,k+ok,oks
                     do jj=j,j+oj,ojs
                        do ii=i,i+oi,ois
                           call mirror(of,nii,njj,nkk,ii,jj,kk,id,jd,kd)
                           md(1) = donor(of,id,jd,kd,w1)
                           call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)
                        end do
                     end do
                  end do

               ! assegna le c.c. chimera solo se entrambe le facce hanno una
               ! c.c. chimera (NON adiacenza)
               !
               else if (pic.lt.offbiu .and.  sic.lt.offbiu) then

                  nd = 4
                  w1 = 1.0_R8P
                  w2 = -0.5_R8P

                  md(1) = donor(ibl,i+  s( 7,n),j+  s( 8,n),k+  s( 9,n),w1)
                  md(2) = donor(ibl,i+  s(10,n),j+  s(11,n),k+  s(12,n),w1)
                  md(3) = donor(ibl,i+2*s( 7,n),j+2*s( 8,n),k+2*s( 9,n),w2)
                  md(4) = donor(ibl,i+2*s(10,n),j+2*s(11,n),k+2*s(12,n),w2)
                  call newchi(igr,ibl,i,j,k,tcc,nd,md)

                  ii = i+oi
                  jj = j+oj
                  kk = k+ok
                  md(1) = donor(ibl,ii+2*s( 7,n),jj+2*s( 8,n),kk+2*s( 9,n),w1)
                  md(2) = donor(ibl,ii+2*s(10,n),jj+2*s(11,n),kk+2*s(12,n),w1)
                  md(3) = donor(ibl,ii+3*s( 7,n),jj+3*s( 8,n),kk+3*s( 9,n),w2)
                  md(4) = donor(ibl,ii+3*s(10,n),jj+3*s(11,n),kk+3*s(12,n),w2)
                  call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

                  nd = 2
                  w1 = 2.0_R8P
                  w2 = -1.0_R8P

                  ii = i-s(7,n)
                  jj = j-s(8,n)
                  kk = k-s(9,n)
                  md(1) = donor(ibl,ii+  s(10,n),jj+  s(11,n),kk+  s(12,n),w1)
                  md(2) = donor(ibl,ii+2*s(10,n),jj+2*s(11,n),kk+2*s(12,n),w2)
                  call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

                  ii = i-s(10,n)
                  jj = j-s(11,n)
                  kk = k-s(12,n)
                  md(1) = donor(ibl,ii+  s(7,n),jj+  s(8,n),kk+  s(9,n),w1)
                  md(2) = donor(ibl,ii+2*s(7,n),jj+2*s(8,n),kk+2*s(9,n),w2)
                  call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

               end if

            end do
         end do
      end do

   !............................................................................
   !  fine loop sugli spigoli
   end do

end do

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!  assegnazione c.c. sui vertici {{{1
!===============================================================================
subroutine calvert(igr)
use prec
use moddef
use modpar, only : donmax,offbiu,blochi,nccnat,xedge,nbl
implicit none
integer (kind=I4P), intent(in) :: igr
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: li,lj,lk, i,j,k, iin,jin,kin, ii,jj,kk, id,jd,kd
integer (kind=I4P) :: of,oi,oj,ok,nd
integer (kind=I4P) :: ibl,tcc,pic,sic,tic,mic
real (kind=R8P)    :: w1,w2
type (donor) :: md(donmax)

!...............................................................................
!  stampa di controllo
write(*,'(a,i3,a)') 'inizio calvert (igr =',igr,')'

!  i punti di vertice sono SEMPRE chimera
tcc = xedge

do ibl=1,nbl

   nii = blochi(ibl,igr)%ni
   njj = blochi(ibl,igr)%nj
   nkk = blochi(ibl,igr)%nk

   !  loop sugli 8 vertici del blocco
   do lk=0,1
      do lj=0,1
         do li=0,1

            i = li*(nii+1)
            j = lj*(njj+1)
            k = lk*(nkk+1)

            iin = li*(nii-1) + 1
            jin = lj*(njj-1) + 1
            kin = lk*(nkk-1) + 1

            oi = i-iin
            oj = j-jin
            ok = k-kin

            !  tipo c.c. sulle due facce che insistono sullo vertice
            !
            pic = blochi(ibl,igr)%cella(iin,jin,k)%t
            sic = blochi(ibl,igr)%cella(iin,j,kin)%t
            tic = blochi(ibl,igr)%cella(i,jin,kin)%t
            if (pic.lt.1 .or. sic.lt.1 .or. tic.lt.1) then
               stop 'Errore in calvert: tipo sugli spigoli deve essere > 0'
            end if

            ! scommentare per dare la priorità alle adiacenze
            ! if (pic.gt.offbiu .or. sic.gt.offbiu .or. tic.gt.offbiu) cycle

            !  c.c. naturali
            !
            if (pic.le.nccnat .or. sic.le.nccnat .or. tic.le.nccnat) then

               nd  = 1
               mic = min(pic,sic,tic)
               w1  = mic
               if (mic.eq.pic) then
                  of = 5+lk
               else if (mic.eq.sic) then
                  of = 3+lj
               else
                  of = 1+li
               end if

               do kk=k,k+ok,ok
                  do jj=j,j+oj,oj
                     do ii=i,i+oi,oi
                        call mirror(of,nii,njj,nkk,ii,jj,kk,id,jd,kd)
                        md(1) = donor(of,id,jd,kd,w1)
                        call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)
                     end do
                  end do
               end do

            ! assegna le c.c. sui vertici solo se le tre facce hanno una
            ! c.c. chimera (NON adiacenza)
            !
            else if (pic.lt.offbiu .and. sic.lt.offbiu .and. tic.lt.offbiu) then

               nd  = 6
               w1 =  0.6666666666_R8P
               w2 = -0.3333333333_R8P

               md(1) = donor(ibl,i-oi,j,k,w1)
               md(2) = donor(ibl,i,j-oj,k,w1)
               md(3) = donor(ibl,i,j,k-ok,w1)
               md(4) = donor(ibl,i-2*oi,j,k,w2)
               md(5) = donor(ibl,i,j-2*oj,k,w2)
               md(6) = donor(ibl,i,j,k-2*ok,w2)
               call newchi(igr,ibl,i,j,k,tcc,nd,md)

               ii = i+oi
               jj = j+oj
               kk = k+ok
               md(1) = donor(ibl,i-oi  ,jj    ,kk    ,w1)
               md(2) = donor(ibl,ii    ,j-oj  ,kk    ,w1)
               md(3) = donor(ibl,ii    ,jj    ,k-ok  ,w1)
               md(4) = donor(ibl,i-2*oi,jj    ,kk    ,w2)
               md(5) = donor(ibl,ii    ,j-2*oj,kk    ,w2)
               md(6) = donor(ibl,ii    ,jj    ,k-2*ok,w2)
               call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

               nd = 4
               w1 = 1.0_R8P
               w2 = -0.5_R8P

               ii = i+oi
               jj = j
               kk = k
               md(1) = donor(ibl,ii,jj-oj  ,kk     ,w1)
               md(2) = donor(ibl,ii,jj     ,kk-ok  ,w1)
               md(3) = donor(ibl,ii,jj-2*oj,kk     ,w2)
               md(4) = donor(ibl,ii,jj     ,kk-2*ok,w2)
               call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

               ii = i
               jj = j+oj
               kk = k
               md(1) = donor(ibl,ii-oi  ,jj,kk     ,w1)
               md(2) = donor(ibl,ii     ,jj,kk-ok  ,w1)
               md(3) = donor(ibl,ii-2*oi,jj,kk     ,w2)
               md(4) = donor(ibl,ii     ,jj,kk-2*ok,w2)
               call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

               ii = i
               jj = j
               kk = k+ok
               md(1) = donor(ibl,ii-oi  ,jj     ,kk,w1)
               md(2) = donor(ibl,ii     ,jj-oj  ,kk,w1)
               md(3) = donor(ibl,ii-2*oi,jj     ,kk,w2)
               md(4) = donor(ibl,ii     ,jj-2*oj,kk,w2)
               call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

               nd = 2
               w1 = 2.0_R8P
               w2 = -1.0_R8P

               ii = i
               jj = j+oj
               kk = k+ok
               md(1) = donor(ibl,ii-oi  ,jj,kk,w1)
               md(2) = donor(ibl,ii-2*oi,jj,kk,w2)
               call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

               ii = i+oi
               jj = j
               kk = k+ok
               md(1) = donor(ibl,ii,jj-oj  ,kk,w1)
               md(2) = donor(ibl,ii,jj-2*oj,kk,w2)
               call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

               ii = i+oi
               jj = j+oj
               kk = k
               md(1) = donor(ibl,ii,jj,kk-ok  ,w1)
               md(2) = donor(ibl,ii,jj,kk-2*ok,w2)
               call newchi(igr,ibl,ii,jj,kk,tcc,nd,md)

            end if

         end do
      end do
   end do

end do

end subroutine

!=/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\==/\
!/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==\/==

!===============================================================================
!===============================================================================
subroutine mirror(f,ni,nj,nk,i,j,k,im,jm,km)
use prec
implicit none
integer (kind=I4P) :: f,ni,nj,nk,i,j,k,im,jm,km

im = i
jm = j
km = k

if (f.eq.1) then
   im = 1-i
else if (f.eq.2) then
   im = 2*ni+1-i
else if (f.eq.3) then
   jm = 1-j
else if (f.eq.4) then
   jm = 2*nj+1-j
else if (f.eq.5) then
   km = 1-k
else if (f.eq.6) then
   km = 2*nk+1-k
else
   stop 'errore in subroutine mirror: numero di faccia sbagliato'
end if

return
end
