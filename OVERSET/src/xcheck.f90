!===============================================================================
!  Verifica che tutte le celle chimera abbiano una valida dipendenza
!===============================================================================
subroutine xcheck(igr)
use prec
use moddef
use modpar
implicit none
character (len=2)  :: fi
character (len=80) :: str
integer (kind=I4P) :: i,j,k,igr,ibl
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: icc,tcc,cou,m,p
integer (kind=I4P) :: imin,jmin,kmin,imax,jmax,kmax
integer (kind=I4P) :: bmin,iimn,jjmn,kkmn,bmax,iimx,jjmx,kkmx
real (kind=R8P) :: rat,tot,dmin,dmax,mmin,mmax,dloc
type (point) :: bc,cc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  stampa di controllo
write(*,'(a,i3,a)') 'inizio xcheck  (igr =',igr,')'

!...............................................................................
! Controlla che non siano rimaste celle standard sulle cornici
!
cou = 0
do ibl=1,nbl
   nii = blogrd(ibl,igr)%ni
   njj = blogrd(ibl,igr)%nj
   nkk = blogrd(ibl,igr)%nk

   m = 0
   i = 0
   do k=0,nkk+1
   do j=0,njj+1
      icc = blochi(ibl,igr)%cella(i,j,k)%t
      if (icc.eq.stdflag) then
         call warn(ndeb,m,ibl,nii,njj,nkk,stdflag,i,j,k)
         cou = cou+1
      end if
   end do
   end do
   i = nii+1
   do k=0,nkk+1
   do j=0,njj+1
      icc = blochi(ibl,igr)%cella(i,j,k)%t
      if (icc.eq.stdflag) then
         call warn(ndeb,m,ibl,nii,njj,nkk,stdflag,i,j,k)
         cou = cou+1
      end if
   end do
   end do

   j = 0
   do k=0,nkk+1
   do i=1,nii
      icc = blochi(ibl,igr)%cella(i,j,k)%t
      if (icc.eq.stdflag) then
         call warn(ndeb,m,ibl,nii,njj,nkk,stdflag,i,j,k)
         cou = cou+1
      end if
   end do
   end do
   j = njj+1
   do k=0,nkk+1
   do i=1,nii
      icc = blochi(ibl,igr)%cella(i,j,k)%t
      if (icc.eq.stdflag) then
         call warn(ndeb,m,ibl,nii,njj,nkk,stdflag,i,j,k)
         cou = cou+1
      end if
   end do
   end do

   k = 0
   do j=1,njj
   do i=1,nii
      icc = blochi(ibl,igr)%cella(i,j,k)%t
      if (icc.eq.stdflag) then
         call warn(ndeb,m,ibl,nii,njj,nkk,stdflag,i,j,k)
         cou = cou+1
      end if
   end do
   end do
   k = nkk+1
   do j=1,njj
   do i=1,nii
      icc = blochi(ibl,igr)%cella(i,j,k)%t
      if (icc.eq.stdflag) then
         call warn(ndeb,m,ibl,nii,njj,nkk,stdflag,i,j,k)
         cou = cou+1
      end if
   end do
   end do

end do

if (cou.gt.0) then
   write(*,'(/,a,i4,a,/)') 'Trovate ',cou,' celle di cornice con icc = 0'
   write(ndeb,'(/,a,i4,a,/)') 'Trovate ',cou,' celle di cornice con icc = 0'
!  stop
end if

!...............................................................................
! Controlla che le celle interne siano solo celle standard, chimera (tcc=27),
! pareti interne (tcc=28), spiaggie (tcc=29)
! calcola la percentuale di celle chimera interne sul totale di celle interne
!
cou = 0
rat = 0d0
tot = 0d0
do ibl=1,nbl
   nii = blogrd(ibl,igr)%ni
   njj = blogrd(ibl,igr)%nj
   nkk = blogrd(ibl,igr)%nk
   do k=1,nkk
      do j=1,njj
         do i=1,nii
            tot = tot+1d0
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if (tcc.gt.0) then
               rat = rat+1d0
               if (tcc.lt.innerx .or. tcc.gt.innerb) then
                  write(*,'(a,i3,a,4(i5))') 'tcc = ',tcc,'  per ',ibl,i,j,k
                  cou = cou+1
               end if
            end if
         end do
      end do
   end do
end do

if (cou.gt.0) then
   write(*,*) 'Trovate ',cou,' celle con tcc <> ',innerx,' o ', innerw
   stop
end if

write(*,'(a,f6.3)') '   rapporto chimera/totali (senza cornici) = ', rat/tot
write(ndeb,'(a,i3.3,a,$)') 'reticolo #',igr,'#  rapporto chimera/totali '
write(ndeb,'(a,f6.3)') '(senza cornici) = ', rat/tot

!...............................................................................
! calcola la percentuale di celle chimera includendo le cornici
!
rat = 0d0
tot = 0d0
do ibl=1,nbl
   do k= -blogrd(ibl,igr)%gc(5)+1 , blogrd(ibl,igr)%nk+blogrd(ibl,igr)%gc(6)
   do j= -blogrd(ibl,igr)%gc(3)+1 , blogrd(ibl,igr)%nj+blogrd(ibl,igr)%gc(4)
   do i= -blogrd(ibl,igr)%gc(1)+1 , blogrd(ibl,igr)%ni+blogrd(ibl,igr)%gc(2)
            tot = tot+1d0
            tcc = blochi(ibl,igr)%cella(i,j,k)%t
            if (tcc.gt.0) rat = rat+1d0
   end do
   end do
   end do
end do

write(*,'(a,f6.3)') '   rapporto chimera/totali (incluse cornici) = ', rat/tot
write(ndeb,'(a,i3.3,a,$)') 'reticolo #',igr,'#  rapporto chimera/totali '
write(ndeb,'(a,f6.3)') '(incluse cornici) = ', rat/tot

!...............................................................................
! Calcola la funzione distanza:
! distanza (centro cella)<-->(baricentro donatori) / (ingombro donatori)
! calcolata per tcc = 21..26,27,41..46
! uguale a 0 altrove

write(ndeb,'(/,a,i3.3,a)') 'reticolo #',igr,'#  distanza chimera-donatori {{{1'

! inizializzazione funzione distanza
do ibl=1,nbl
   blochi(ibl,igr)%cella(:,:,:)%dist = 0.0d0
end do

mmax = -1d30
do ibl=1,nbl
   nii = blogrd(ibl,igr)%ni
   njj = blogrd(ibl,igr)%nj
   nkk = blogrd(ibl,igr)%nk
   write(ndeb,'(a,i3.3)',advance='no') '   blocco n.',ibl
   dmin = 1d20
   dmax = -1d20
   do k=1,nkk
   do j=1,njj
   do i=1,nii
!  do k=0,nkk+1
!  do j=0,njj+1
!  do i=0,nii+1
      dloc = 0.0d0
      tcc = blochi(ibl,igr)%cella(i,j,k)%t
      if (tcc.gt.0) then
         ! punti interni e facce chimera con donatori per centro cella
         if ((tcc.eq.innerx).or.(tcc.ge.offgen+1 .and. tcc.le.offgen+6)) then
            call bar_ing_donatori(igr,ibl,i,j,k,tot,bc)
            cc = blomet(ibl,igr)%cella(i,j,k)%cen
            dloc = sqrt(dist(cc,bc))/tot
         ! facce chimera con donatori per centro faccia
         else if (tcc.ge.offchi+1 .and. tcc.le.offchi+6) then
            call bar_ing_donatori(igr,ibl,i,j,k,tot,bc)
            p = tcc-offchi
            call centro_faccia_cornice(igr,ibl,i,j,k,p,cc)
            dloc = sqrt(dist(cc,bc))/tot
         end if
      end if
      if (dloc.lt.dmin) then
         imin = i
         jmin = j
         kmin = k
         dmin = dloc
      end if
      if (dloc.gt.dmax) then
         imax = i
         jmax = j
         kmax = k
         dmax = dloc
      end if
   end do
   end do
   end do
!  write(ndeb,'(a,e9.3,a,f9.5,a,3i4,a)') &
!           '      min = ',dmin,'      max = ',dmax,' (',imax,jmax,kmax,')'
   write(fi,'(a,i1)') 'I',nfig(namax)+1
   str = '(a,f9.5,' // fi // ',' // fi // ',' // fi // ',' // fi // ',' &
         // fi // ',' // fi // ')'
   write(ndeb,str) '     max = ',dmax,imax,jmax,kmax,nii,njj,nkk
   if (dmax.gt.mmax) then
      mmax = dmax
      bmax = ibl
      iimx = imax
      jjmx = jmax
      kkmx = kmax
   end if
end do
write(ndeb,'(a,f9.5,6x,a,i3.3)') '   MASSIMO ASSOLUTO = ',mmax,'BLOCCO N.',bmax

write(*,'(a)', advance='no') '   distanza massima            = '
write(*,'(f9.5,5x,a,i3.3,a,3i4,a)') mmax,'blocco n.',bmax,' (',iimx,jjmx,kkmx,')'

!...............................................................................
! Calcola la funzione rapporto:
! (dimensione cella chimera) / (dimensione media donatori)
! calcolata per tcc = 21..26,27,41..46
! uguale a 0 altrove

write(ndeb,'(/,a,i3.3,a)') 'reticolo #',igr,&
                           '#  rapporto dimensioni chimera-donatori {{{1'

mmin =  1d30
mmax = -1d30
do ibl=1,nbl
   nii = blogrd(ibl,igr)%ni
   njj = blogrd(ibl,igr)%nj
   nkk = blogrd(ibl,igr)%nk
   write(ndeb,'(a,i3.3)',advance='no') '   blocco n.',ibl
   dmin = 1d20
   dmax = -1d20
   do k=0,nkk+1
      do j=0,njj+1
         do i=0,nii+1
            dloc = blochi(ibl,igr)%cella(i,j,k)%size
            tcc =  blochi(ibl,igr)%cella(i,j,k)%t
            ! punti interni e facce chimera con donatori per centro cella
            if ((tcc.eq.innerx) .or. &
                (tcc.ge.offgen+1 .and. tcc.le.offgen+6) .or. &
                (tcc.ge.offchi+1 .and. tcc.le.offchi+6)) then
               call dim_media_donatori(igr,ibl,i,j,k,tot)
               rat = dloc/tot
               if (rat.lt.dmin) then
                  imin = i
                  jmin = j
                  kmin = k
                  dmin = rat
               end if
               if (rat.gt.dmax) then
                  imax = i
                  jmax = j
                  kmax = k
                  dmax = rat
               end if
            end if
         end do
      end do
   end do
!  write(ndeb,'(a,e9.3,a,f9.5,a,3i4,a)') &
!           '      min = ',dmin,'      max = ',dmax,' (',imax,jmax,kmax,')'
   write(fi,'(a,i1)') 'I',nfig(namax)+1
   str = '(a,e9.3,a,' // fi // ',",",' // fi // ',",",' // fi // ',' &
       // 'a,e9.3,a,' // fi // ',",",' // fi // ',",",' // fi // ',a)'
   write(ndeb,str) '   min = ',dmin,' (',imin,jmin,kmin,&
                   ')  max = ',dmax,' (',imax,jmax,kmax,')'
   if (dmax.gt.mmax) then
      mmax = dmax
      bmax = ibl
      iimx = imax
      jjmx = jmax
      kkmx = kmax
   end if
   if (dmin.lt.mmin) then
      mmin = dmin
      bmin = ibl
      iimn = imin
      jjmn = jmin
      kkmn = kmin
   end if
end do
write(ndeb,'(a,e10.3,6x,a,i3.3)') '   MINIMO  ASSOLUTO = ',mmin,'BLOCCO N.',bmin
write(ndeb,'(a,e10.3,6x,a,i3.3)') '   MASSIMO ASSOLUTO = ',mmax,'BLOCCO N.',bmax

write(*,'(a)', advance='no') '   rapporto minimo dimensioni  = '
write(*,'(e10.3,5x,a,i3.3,a,3i4,a)') mmin,'blocco n.',bmin,' (',iimn,jjmn,kkmn,')'
write(*,'(a)', advance='no') '   rapporto massimo dimensioni = '
write(*,'(e10.3,5x,a,i3.3,a,3i4,a)') mmax,'blocco n.',bmax,' (',iimx,jjmx,kkmx,')'

end subroutine

!===============================================================================
! calcolo baricentro ed ingombro complessivo dei donatori
!===============================================================================
subroutine bar_ing_donatori(igr,ibl,i,j,k,dm,bc)
use prec
use moddef
use modpar, only : blochi,blogrd,blomet
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
integer (kind=I4P) :: bb,ii,jj,kk,iii,jjj,kkk,nd,p
real (kind=R8P) :: dm, ww
type (point) :: bc, cc,mnc,mxc

bc  = 0d0
nd  = blochi(ibl,igr)%cella(i,j,k)%n

! inizializza dm tramite le dimensioni della prima donatrice
mnc = 1d20
mxc = -1d20
bb = blochi(ibl,igr)%cella(i,j,k)%q(1)%b
ii = blochi(ibl,igr)%cella(i,j,k)%q(1)%i
jj = blochi(ibl,igr)%cella(i,j,k)%q(1)%j
kk = blochi(ibl,igr)%cella(i,j,k)%q(1)%k
do kkk=kk-1,kk
   do jjj=jj-1,jj
      do iii=ii-1,ii
         cc = blogrd(bb,igr)%nodo(iii,jjj,kkk)
         mnc%x = min(mnc%x,cc%x)
         mnc%y = min(mnc%y,cc%y)
         mnc%z = min(mnc%z,cc%z)
         mxc%x = max(mxc%x,cc%x)
         mxc%y = max(mxc%y,cc%y)
         mxc%z = max(mxc%z,cc%z)
      end do
   end do
end do
cc = mxc-mnc
dm = 0.5d0*norma(cc)

mnc = 1d20
mxc = -1d20
do p=1,nd
   ! blocco e indici del donatore
   bb = blochi(ibl,igr)%cella(i,j,k)%q(p)%b
   ii = blochi(ibl,igr)%cella(i,j,k)%q(p)%i
   jj = blochi(ibl,igr)%cella(i,j,k)%q(p)%j
   kk = blochi(ibl,igr)%cella(i,j,k)%q(p)%k
   ww = blochi(ibl,igr)%cella(i,j,k)%q(p)%w
   ! centro cella del donatore
   cc = blomet(bb,igr)%cella(ii,jj,kk)%cen
   ! baricentro donatori
   bc = bc + ww*cc
   ! ingombro dell'insieme dei donatori
   mnc%x = min(mnc%x,cc%x)
   mnc%y = min(mnc%y,cc%y)
   mnc%z = min(mnc%z,cc%z)
   mxc%x = max(mxc%x,cc%x)
   mxc%y = max(mxc%y,cc%y)
   mxc%z = max(mxc%z,cc%z)
end do
cc = mxc-mnc
dm = max(dm,0.5d0*norma(cc))

end subroutine

!===============================================================================
! calcolo una dimensione media dei donatori
!===============================================================================
subroutine dim_media_donatori(igr,ibl,i,j,k,dm)
use prec
use modpar, only : blochi
implicit none
integer (kind=I4P), intent(in) :: igr,ibl,i,j,k
integer (kind=I4P) :: bb,ii,jj,kk,nn
real (kind=R8P) :: dm,dl,ww

dm = 0d0
do nn=1,blochi(ibl,igr)%cella(i,j,k)%n
   bb = blochi(ibl,igr)%cella(i,j,k)%q(nn)%b
   ii = blochi(ibl,igr)%cella(i,j,k)%q(nn)%i
   jj = blochi(ibl,igr)%cella(i,j,k)%q(nn)%j
   kk = blochi(ibl,igr)%cella(i,j,k)%q(nn)%k
   ww = blochi(ibl,igr)%cella(i,j,k)%q(nn)%w
   dl = blochi(bb,igr)%cella(ii,jj,kk)%size
   dm = dm + ww*dl
end do

end subroutine

!===============================================================================
! stampa warning cornici indefinite
!===============================================================================
subroutine warn(ndeb,m,ibl,imax,jmax,kmax,stdflag,i,j,k)
use prec
implicit none
integer (kind=I4P), intent(inout) :: m
integer (kind=I4P), intent(in) :: ndeb,ibl,imax,jmax,kmax,stdflag,i,j,k

if (m.eq.0) then
   write(*,'(a)') "****************************"
   write(*,'(a,4(i5))') "BLOCCO  ",ibl,imax,jmax,kmax
   write(*,'(a)') "****************************"
   write(*,'(a,i2,a,3(i5))') 'icc = ',stdflag,'  in ',i,j,k

   write(ndeb,'(a)') "****************************"
   write(ndeb,'(a,4(i5))') "BLOCCO  ",ibl,imax,jmax,kmax
   write(ndeb,'(a)') "****************************"
   write(ndeb,'(a,i2,a,3(i5))') 'icc = ',stdflag,'  in ',i,j,k

   m = 1
else
   write(*,'(a,i2,a,3(i5))') 'icc = ',stdflag,'  in ',i,j,k
   write(ndeb,'(a,i2,a,3(i5))') 'icc = ',stdflag,'  in ',i,j,k
end if

end subroutine
