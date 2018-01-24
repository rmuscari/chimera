!===============================================================================
! cerca donatori per punti interni di strato limite
!===============================================================================
subroutine strato_limite_main(igr,p)
!{{{1
use prec
use moddef, only : point, operator(*)
use modpar
implicit none
integer (kind=I4P), intent(in) :: igr,p
integer (kind=I4P) :: nii,njj,nkk
integer (kind=I4P) :: tgr,ibl,nf,nl(6),is,js,ks,i,j,k
integer (kind=I4P) :: args(8)  ! vettore di comodo per passare gli argomenti
integer (kind=I4P) :: famn,famp
type (point) :: n

! famiglia del patch
famp = mod(bp(4,p),100)

! memorizza in variabili locali i parenti di "p"
j = 1
famn = 0
if (famp.gt.0) then
   blt = boulaf(famp)
   do while (family2patch(j,famp).gt.0)
      if (family2patch(j,famp).ne.p) then
         famn = famn+1
         famw(famn) = family2patch(j,famp)
      end if
      j = j+1
   end do
else
   blt = boulap(p)
end if

tgr = 2**(igr-1)
ibl = bp(1,p)
nf  = bp(2,p)
nl  = bp(5:10,p)/tgr

nii = blogrd(ibl,igr)%ni
njj = blogrd(ibl,igr)%nj
nkk = blogrd(ibl,igr)%nk

if (nf.eq.1) then
   is = 1
   args = (/  igr,ibl,nii,njj,nkk,is,0,0 /)
   do k=nl(5)+1,nl(6)
   do j=nl(3)+1,nl(4)
      i  = 1
      n = real(is,kind(0d0)) * blomet(ibl,igr)%cella(i,j,k)%sni
      call strato_limite_subr(args,i,j,k,p,famn,n)
   end do
   end do
else if (nf.eq.2) then
   is = -1
   args = (/  igr,ibl,nii,njj,nkk,is,0,0 /)
   do k=nl(5)+1,nl(6)
   do j=nl(3)+1,nl(4)
      i  = nii
      n = real(is,kind(0d0)) * blomet(ibl,igr)%cella(i,j,k)%sni
      call strato_limite_subr(args,i,j,k,p,famn,n)
   end do
   end do

else if (nf.eq.3) then
   js = 1
   args = (/  igr,ibl,nii,njj,nkk,0,js,0 /)
   do k=nl(5)+1,nl(6)
   do i=nl(1)+1,nl(2)
      j  = 1
      n = real(js,kind(0d0)) * blomet(ibl,igr)%cella(i,j,k)%snj
      call strato_limite_subr(args,i,j,k,p,famn,n)
   end do
   end do
else if (nf.eq.4) then
   js = -1
   args = (/  igr,ibl,nii,njj,nkk,0,js,0 /)
   do k=nl(5)+1,nl(6)
   do i=nl(1)+1,nl(2)
      j  = njj
      n = real(js,kind(0d0)) * blomet(ibl,igr)%cella(i,j,k)%snj
      call strato_limite_subr(args,i,j,k,p,famn,n)
   end do
   end do

else if (nf.eq.5) then
   ks = 1
   args = (/  igr,ibl,nii,njj,nkk,0,0,ks /)
   do j=nl(3)+1,nl(4)
   do i=nl(1)+1,nl(2)
      k  = 1
      n = real(ks,kind(0d0)) * blomet(ibl,igr)%cella(i,j,k)%snk
      call strato_limite_subr(args,i,j,k,p,famn,n)
   end do
   end do
else if (nf.eq.6) then
   ks = -1
   args = (/  igr,ibl,nii,njj,nkk,0,0,ks /)
   do j=nl(3)+1,nl(4)
   do i=nl(1)+1,nl(2)
      k  = nkk
      n = real(ks,kind(0d0)) * blomet(ibl,igr)%cella(i,j,k)%snk
      call strato_limite_subr(args,i,j,k,p,famn,n)
   end do
   end do

end if

return
!1}}}
end subroutine strato_limite_main

!===============================================================================
! trova i donatori per le celle di strato limite a partire dal "nearest wall
! vertex"
!===============================================================================
subroutine strato_limite_subr(args,i,j,k,p,famn,n)
!{{{2
use prec
use moddef
use modpar, only : innerw,innerx,donmax,blochi,blomet,lvmax,blt,nbl,lv,pri
implicit none
integer (kind=I4P), intent(in) :: args(8)
integer (kind=I4P) :: i,j,k
integer (kind=I4P), intent(in) :: p,famn
type (point), intent(in) :: n
integer (kind=I4P) :: igr,ibl,nii,njj,nkk,is,js,ks
integer (kind=I4P) :: pbl,ip,jp,kp
integer (kind=I4P) :: td,nd,pp,nwv(9),niter,jbl
integer (kind=I4P) :: iq,jq,kq,dq
integer (kind=I4P) :: test_indici
logical (kind=I1P) :: found,test
real (kind=R8P)    :: dw
type (donor)     :: md(donmax)
type (point)     :: e1,e2,c
integer (kind=I4P) :: iii,jjj,kkk
integer (kind=I1P), allocatable, dimension(:) :: lvbad

! alloca lvbad (in vista della parellizzazione tramite openmp è meglio che sia
! una variabile locale e che non compaia in nessun modulo)
allocate(lvbad(0:lvmax))

igr = args(1)
ibl = args(2)
nii = args(3)
njj = args(4)
nkk = args(5)
is  = args(6)
js  = args(7)
ks  = args(8)

nwv(1) = 0
lvbad(0:lvmax) = 1 ! per la marcatura delle pareti interne
pbl = 0 ! solo per non far apparire il warning quando compilo con il gfortran

! se il patch non ha familiari non cerca donatori
if (famn.gt.0) then
   ! prima cella a parete
   call getcenter(igr,ibl,i,j,k,c)
   ! cerca donatore per la prima cella a parete
   call strato_limite_findonor(igr,p,famn,lvbad,nwv,e1,e2,c,n)
   pbl = nwv(1)
end if

! un loop su tutti i punti di strato limite va comunque fatto:
! o per assegnare i donatori
! o per cercare i punti che cascano nelle scatole
if (nwv(1).gt.0) then

   ! la prima cella a parete ha un donatore: cerca i donatori per tutte le celle
   ! di strato limite

   dw = blochi(ibl,igr)%cella(i,j,k)%dist
   pp = blochi(ibl,igr)%cella(i,j,k)%wall
   test_indici = (i-1)*(nii-i)*abs(is) &  ! se lo s.l. è spesso rispetto alle
               + (j-1)*(njj-j)*abs(js) &  ! dimensioni del blocco può succedere
               + (k-1)*(nkk-k)*abs(ks)    ! di sforare un indice

   do while ((dw.le.blt) .and. (pp.eq.p) .and. (test_indici.ge.0))

      ! se pp.ne.p vuol dire che la parete più vicina alla cella non è più
      ! quella per la quale sto calcolando lo strato limite (ad es. mi trovo in
      ! corrispondenza di un angolo)

      ! cerca iterativamente il vertice più vicino a "c" muovendosi prima lungo
      ! la coordinata uscente dalla parete e poi nel piano parallelo
      c = blomet(ibl,igr)%cella(i,j,k)%cen
      niter = 1
      found = .false.
      do while (.not.found .and. niter.lt.4)
         ! cerca il vertice con distanza da parete più prossima a partire da nwv
         ! lo memorizza in ip,jp,kp lasciando invariato nwv
         call ricerca_lungo_normale(igr,nwv,ip,jp,kp,dw)
         ! sul livello determinato da ip o jp o kp (dipende dalla faccia)
         ! determina la faccetta che contiene la proiezione di "c"
         ! memorizza il vertice della faccetta più vicino a "c" in ip,jp,kp
         ! aggiorna gli indici di nwv "paralleli" al piano (ad es. se la parete
         ! è sul piano J=0 aggiorna nwv(3) e nwv(5))
         iq = ip
         jq = jp
         kq = kp
         call ricerca_sul_piano(igr,nwv,ip,jp,kp,e1,e2,c,n)
         if (nwv(1).gt.0) then
            ! ha trovato un donatore; verifica se si è arrivato a convergenza
            dq = abs(ip-iq)+abs(jp-jq)+abs(kp-kq)
            if (dq.eq.0) then
               ! individua donatori (attorno al vertice più vicino)
               ! e calcola pesi
               call trilinear_strato_limite(igr,nwv,ip,jp,kp,nd,md,dw,e1,e2,c)
               kkk = 0
               do jjj=1,nd
                  iii = blochi(md(jjj)%b,igr)% &
                        cella(md(jjj)%i,md(jjj)%j,md(jjj)%k)%t
                  kkk = max(iii,kkk)
               end do
               if (kkk.eq.0) found = .true.
            end if
            niter = niter+1
         else
            ! non ha trovato donatori; forza l'uscita dal ciclo tramite niter
            niter = 10
         end if
      end do

      if (found) then
         ! crea la cella chimera
         td = innerx
         call newchi(igr,ibl,i,j,k,td,nd,md)
      else
         ! controlla che il punto non caschi nella scatola di un altro livello
         ! lvbad(0:lvmax) NON va reinizializzato
         ! il seguente loop serve per escludere i blocchi che contengono "c"
         ! anche se non hanno patch di parete nella stessa famiglia
         do jbl=1,nbl
            if (lvbad(lv(jbl)).gt.0 .and. jbl.ne.ibl) then
               call nelblocco(igr,jbl,test,c)
               if (test) lvbad(lv(jbl)) = -1
            end if
         end do
         call check_parete_interna(pp,lv(ibl),pri(ibl),lvbad,c)
         if (pp.ge.0) then
            md(1) = donor(ibl,i,j,k,real(pp,R8P))
            call newchi(igr,ibl,i,j,k,innerw,1,md)
         end if
         ! nwv(1) è stato azzerato in ricerca_sul_piano e va ripristinato
         nwv(1) = pbl
      end if

      ! passa al punto successivo di s.l.
      i = i+is
      j = j+js
      k = k+ks
      dw = blochi(ibl,igr)%cella(i,j,k)%dist
      pp = blochi(ibl,igr)%cella(i,j,k)%wall
      test_indici = (i-1)*(nii-i)*abs(is) &
                  + (j-1)*(njj-j)*abs(js) &
                  + (k-1)*(nkk-k)*abs(ks)
   end do

else

   ! la prima cella a parete NON ha un donatore: controlla se i punti di strato
   ! limite cascano dentro qualche scatola

   dw = blochi(ibl,igr)%cella(i,j,k)%dist
   pp = blochi(ibl,igr)%cella(i,j,k)%wall
   test_indici = (i-1)*(nii-i)*abs(is) &
               + (j-1)*(njj-j)*abs(js) &
               + (k-1)*(nkk-k)*abs(ks)

   do while ((dw.le.blt) .and. (pp.eq.p) .and. (test_indici.ge.0))

      ! controlla che il punto non caschi in una scatola
      c = blomet(ibl,igr)%cella(i,j,k)%cen
      do jbl=1,nbl
         if (lvbad(lv(jbl)).gt.0 .and. jbl.ne.ibl) then
            call nelblocco(igr,jbl,test,c)
            if (test) lvbad(lv(jbl)) = -1
         end if
      end do
      call check_parete_interna(pp,lv(ibl),pri(ibl),lvbad,c)
      if (pp.ge.0) then
         md(1) = donor(ibl,i,j,k,real(pp,R8P))
         call newchi(igr,ibl,i,j,k,innerw,1,md)
      end if

      ! passa al punto successivo di s.l.
      i = i+is
      j = j+js
      k = k+ks
      dw = blochi(ibl,igr)%cella(i,j,k)%dist
      pp = blochi(ibl,igr)%cella(i,j,k)%wall
      test_indici = (i-1)*(nii-i)*abs(is) &
                  + (j-1)*(njj-j)*abs(js) &
                  + (k-1)*(nkk-k)*abs(ks)
   end do

end if

deallocate(lvbad)

return
!2}}}
end subroutine strato_limite_subr

!===============================================================================
! cerca il vertice a distanza controvariante minima da "c"
!===============================================================================
subroutine strato_limite_findonor(igr,pat,famn,lvbad,nwv,e1,e2,c,n)
!{{{9
use prec
use moddef
use modpar, only : blogrd,blomet,lvmax,blt,lv,famw,ingopatch,bp
implicit none

integer (kind=I4P), intent(in) :: igr,pat,famn
integer (kind=I1P) :: lvbad(0:lvmax)
integer (kind=I4P) :: nwv(9)
type (point), intent(in) :: c,n
!
logical (kind=I1P) :: coocontr,controllo_normali,controllo_distanza
logical (kind=I1P) :: fuori_ingombro
!
integer (kind=I4P) :: p,q,nl(6),tgr
integer (kind=I4P) :: ip,jp,kp,pbl,nf
integer (kind=I4P) :: imin,jmin,kmin,pmin
integer (kind=I4P) :: is,js,ks
logical (kind=I1P) :: cs,ct
real (kind=R8P)    :: x1,x2,dis_norm,dis_plan
type (point)     :: q1,q2,e1,e2,cp,ry
type (point)     :: norm_face

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! cerca il vertice di parete più vicino
tgr = 2**(igr-1)
do q=1,famn
   p   = famw(q)
   if (fuori_ingombro(c,ingopatch(1,p))) go to 13
   pbl = bp(1,p)
   ! se il livello del blocco è già stato marcato come non-perforante (lvbad<0)
   ! e il punto "c" ha già trovato un valido donatore (nwv(1)>0) salta al patch
   ! successivo
   if (lvbad(lv(pbl)).lt.0 .and. nwv(1).gt.0) go to 13
   nf  = bp(2,p)
   nl  = bp(5:10,p)/tgr

   ! cerca sul patch più vicino la faccetta che contiene la proiezione
   if (nf.lt.3) then
      is = 3-2*nf
      js = 0
      ks = 0
      ip = nl(1)
      do kp=nl(5),nl(6)-1
      do jp=nl(3),nl(4)-1
         norm_face = real(is,kind(0d0)) * blomet(pbl,igr)%cella(ip,jp+1,kp+1)%sni
         if (controllo_normali(n,norm_face)) then
            cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
            q1 = blogrd(pbl,igr)%nodo(ip,jp+1,kp)
            q2 = blogrd(pbl,igr)%nodo(ip,jp,kp+1)
            e1 = q1-cp
            e2 = q2-cp
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip
               jmin = jp+nint(x1)
               kmin = kp+nint(x2)
               go to 19
            end if
            cp = blogrd(pbl,igr)%nodo(ip,jp+1,kp+1)
            e1 = q2-cp  !!!
            e2 = q1-cp  !!!
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip
               jmin = jp+1-nint(x1)
               kmin = kp+1-nint(x2)
               go to 19
            end if
         end if
      end do
      end do

   else if (nf.gt.4) then
      is = 0
      js = 0
      ks = 11-2*nf
      kp = nl(5)
      do jp=nl(3),nl(4)-1
      do ip=nl(1),nl(2)-1
         norm_face = real(ks,kind(0d0)) * blomet(pbl,igr)%cella(ip+1,jp+1,kp)%snk
         if (controllo_normali(n,norm_face)) then
            cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
            q1 = blogrd(pbl,igr)%nodo(ip+1,jp,kp)
            q2 = blogrd(pbl,igr)%nodo(ip,jp+1,kp)
            e1 = q1-cp
            e2 = q2-cp
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+nint(x1)
               jmin = jp+nint(x2)
               kmin = kp
               go to 19
            end if
            cp = blogrd(pbl,igr)%nodo(ip+1,jp+1,kp)
            e1 = q2-cp  !!!
            e2 = q1-cp  !!!
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+1-nint(x1)
               jmin = jp+1-nint(x2)
               kmin = kp
               go to 19
            end if
         end if
      end do
      end do

   else
      is = 0
      js = 7-2*nf
      ks = 0
      jp = nl(3)
      do kp=nl(5),nl(6)-1
      do ip=nl(1),nl(2)-1
         norm_face = real(js,kind(0d0)) * blomet(pbl,igr)%cella(ip+1,jp,kp+1)%snj
         if (controllo_normali(n,norm_face)) then
            cp = blogrd(pbl,igr)%nodo(ip,jp,kp)
            q1 = blogrd(pbl,igr)%nodo(ip+1,jp,kp)
            q2 = blogrd(pbl,igr)%nodo(ip,jp,kp+1)
            e1 = q1-cp
            e2 = q2-cp
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+nint(x1)
               jmin = jp
               kmin = kp+nint(x2)
               go to 19
            end if
            cp = blogrd(pbl,igr)%nodo(ip+1,jp,kp+1)
            e1 = q2-cp  !!!
            e2 = q1-cp  !!!
            ry = c-cp
            !
            dis_norm = ry.dot.norm_face
            dis_plan = normq(ry-dis_norm*norm_face)
            cs = controllo_distanza(dis_plan,max(normq(e1),normq(e2)))
            ct = controllo_distanza(dis_norm**2,4d0*blt)
            !
            if (cs .and. ct .and. coocontr(ry,e1,e2,x1,x2) ) then
               pmin = p
               imin = ip+1-nint(x1)
               jmin = jp
               kmin = kp+1-nint(x2)
               go to 19
            end if
         end if
      end do
      end do
   end if
   ! se arriva qui non ha trovato donatori
   go to 13

19 continue

   ! il livello del patch va escluso da quelli che bucano il punto
   lvbad(lv(pbl)) = -1

   ! ha trovato un possibile donatore (D) per la cella di strato limite (C)
   ! * se la priorità di D è più alta di C, accetta senza'altro il donatore
   ! * se le priorità sono uguali, accetta D solo se ha una dimensione in pianta
   !   più piccola di C
   !
   ! N.B.: per ciascuna famiglia il loop per la ricerca di possibili donatori
   ! viene eseguito a partire dai patch a priorità più alta. Quindi se si
   ! verifica che la priorità di D è uguale a quella di C, vuol dire che non ha
   ! trovato alcun possibile donatore sui patch a priorità più alta.
   if ((nwv(1).eq.0)) then
      if (  (bp(4,pmin).gt.bp(4,pat)) .or. &
           ((bp(4,pmin).eq.bp(4,pat)).and.(normq(n).gt.normq(norm_face))) ) then
         nwv(1) = pbl
         nwv(2) = nf
         nwv(3) = imin
         nwv(4) = jmin
         nwv(5) = kmin
         nwv(6) = is
         nwv(7) = js
         nwv(8) = ks
         nwv(9) = pmin
         e1 = vers0(e1)
         cp = e1.cross.e2
         e2 = cp.cross.e1
         e2 = vers0(e2)
      end if
   end if

! fine loop sui patch della stessa famiglia
13 continue
end do

!9}}}
end subroutine strato_limite_findonor
