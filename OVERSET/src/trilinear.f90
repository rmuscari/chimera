!===============================================================================
! Dati gli 8 vertici di un esaedro, calcola l'interpolazione trilineare in un
! punto dato
!===============================================================================
subroutine trilinear(p,q,u)
use prec
use moddef
implicit none

real (kind=R8P),    parameter :: smmax = 1.0d-8
real (kind=R8P),    parameter :: comax = 1.0d0+1.0d-5
integer (kind=I4P), parameter :: itmax = 21

type (point) :: p(0:1,0:1,0:1),q
real (kind=R8P) :: u(8)
real (kind=R8P) :: a,b,c,m(4,3)
real (kind=R8P) :: sm
integer (kind=I4P) :: iter
logical (kind=I1P) :: cont

call tri2or(p,q,a,b,c)
cont = .true.
iter = 0
do while (cont)

   call matrici(p,q,a,b,c,m)
   call sislin(m)

   a = a+m(1,1)
   b = b+m(1,2)
   c = c+m(1,3)
   iter = iter+1

   sm = m(1,1)**2 + m(1,2)**2 + m(1,3)**2
   cont = sm.gt.smmax .and. iter.lt.itmax

end do

cont = abs(a).gt.comax .or. &
       abs(b).gt.comax .or. &
       abs(c).gt.comax .or. &
       iter  .ge.itmax .or. &
       sm    .ge.smmax

if (cont) then
!  write(*,'(a,e13.5,i4)') 'sm,iter           :',sm,iter
!  write(*,'(a,3e13.5)') 'a,b,c (trilinear) :',a,b,c
   call tri2or(p,q,a,b,c)
!  write(*,'(a,3e13.5)') 'a,b,c (tri2or)    :',a,b,c
end if

u(1) = 0.125d0*(1-c)*(1-b)*(1-a) ! p(0,0,0)
u(2) = 0.125d0*(1-c)*(1-b)*(1+a) ! p(1,0,0)
u(3) = 0.125d0*(1-c)*(1+b)*(1-a) ! p(0,1,0)
u(4) = 0.125d0*(1-c)*(1+b)*(1+a) ! p(1,1,0)
u(5) = 0.125d0*(1+c)*(1-b)*(1-a) ! p(0,0,1)
u(6) = 0.125d0*(1+c)*(1-b)*(1+a) ! p(1,0,1)
u(7) = 0.125d0*(1+c)*(1+b)*(1-a) ! p(0,1,1)
u(8) = 0.125d0*(1+c)*(1+b)*(1+a) ! p(1,1,1)

return
end

!-------------------------------------------------------------------------------
subroutine matrici(p,q,a,b,c,m)
use prec
use moddef
implicit none
type (point) :: f
type (point) :: p(0:1,0:1,0:1),q
real (kind=R8P) :: a,b,c
real (kind=R8P) :: m(4,3)

! derivate cambiate di segno di x,y,z rispetto ad "a"
f = 0.125d0*( (1-c)*(1-b)*p(0,0,0) - (1-c)*(1-b)*p(1,0,0) +&
              (1-c)*(1+b)*p(0,1,0) - (1-c)*(1+b)*p(1,1,0) +&
              (1+c)*(1-b)*p(0,0,1) - (1+c)*(1-b)*p(1,0,1) +&
              (1+c)*(1+b)*p(0,1,1) - (1+c)*(1+b)*p(1,1,1)  )
m(1,1) = f%x
m(1,2) = f%y
m(1,3) = f%z

! derivate cambiate di segno di x,y,z rispetto ad "b"
f = 0.125d0*( (1-c)*(1+a)*p(1,0,0) - (1-c)*(1+a)*p(1,1,0) +&
              (1-c)*(1-a)*p(0,0,0) - (1-c)*(1-a)*p(0,1,0) +&
              (1+c)*(1+a)*p(1,0,1) - (1+c)*(1+a)*p(1,1,1) +&
              (1+c)*(1-a)*p(0,0,1) - (1+c)*(1-a)*p(0,1,1)  )
m(2,1) = f%x
m(2,2) = f%y
m(2,3) = f%z

! derivate cambiate di segno di x,y,z rispetto ad "c"
f = 0.125d0*( (1-b)*(1+a)*p(1,0,0) - (1+b)*(1+a)*p(1,1,1) +&
              (1+b)*(1+a)*p(1,1,0) - (1-b)*(1+a)*p(1,0,1) +&
              (1+b)*(1-a)*p(0,1,0) - (1-b)*(1-a)*p(0,0,1) +&
              (1-b)*(1-a)*p(0,0,0) - (1+b)*(1-a)*p(0,1,1)  )
m(3,1) = f%x
m(3,2) = f%y
m(3,3) = f%z

! funzione da annullare
f = 0.125d0*( (1-c)*(1-b)*(1-a)*p(0,0,0) + (1-c)*(1-b)*(1+a)*p(1,0,0) +&
              (1-c)*(1+b)*(1-a)*p(0,1,0) + (1-c)*(1+b)*(1+a)*p(1,1,0) +&
              (1+c)*(1-b)*(1-a)*p(0,0,1) + (1+c)*(1-b)*(1+a)*p(1,0,1) +&
              (1+c)*(1+b)*(1-a)*p(0,1,1) + (1+c)*(1+b)*(1+a)*p(1,1,1) ) - q
m(4,1) = f%x
m(4,2) = f%y
m(4,3) = f%z

return
end

!===============================================================================
! soluzione di un sistema di equazioni lineari   M(3,3)*X(3)=N(3)
! N.B.:  N è memorizzato in M(4,:)
!===============================================================================
subroutine sislin(m)
use prec
implicit none
real (kind=R8P) :: m(4,3)
real (kind=R8P) :: w(4)
integer (kind=I4P) :: i,j,k

! loop di triangolarizzazione
do j=1,2

   ! ricerca del massimo m(j,j:3)
   k = j
   do i=j+1,3
      if (abs(m(j,i)) .gt. abs(m(j,k)))  k = i
   enddo

   ! inversione delle righe
   if (k.ne.i) then
      w = m(:,k)
      m(:,k) = m(:,j)
      m(:,j) = w
   endif

   ! triangolarizzazione
   w = m(:,j)/m(j,j)
   do k=j+1,3
      m(:,k) = m(:,k) - m(j,k)*w
      m(j,k) = 0.0d0
   enddo

enddo

! sostituzione all'indietro
w(3) = m(4,3)/m(3,3)
do j=2,1,-1
   w(j) = ( m(4,j) - sum(m(j+1:3,j)*w(j+1:3)) ) / m(j,j)
enddo

! memorizza la soluzione in m(1,:)
m(1,:) = w(1:3)

return
end

!===============================================================================
subroutine tri2or(p,q,a,b,c)
use prec
use moddef
implicit none
type (point) :: p(0:1,0:1,0:1),q
real (kind=R8P) :: a,b,c
type (point) :: f,qo

qo = q - 0.125d0*( p(0,0,0) + p(1,0,0) +&
                   p(0,1,0) + p(1,1,0) +&
                   p(0,0,1) + p(1,0,1) +&
                   p(0,1,1) + p(1,1,1)  )

! derivate di x,y,z rispetto ad "a" in (0,0,0)
f = p(1,1,1) - p(0,1,1) +&
    p(1,1,0) - p(0,1,0) +&
    p(1,0,1) - p(0,0,1) +&
    p(1,0,0) - p(0,0,0)

f = 0.125d0*f
a = qo.dot.f/normq(f)

! derivate di x,y,z rispetto ad "b"
f = p(1,1,1) - p(1,0,1) +&
    p(1,1,0) - p(1,0,0) +&
    p(0,1,1) - p(0,0,1) +&
    p(0,1,0) - p(0,0,0)

f = 0.125d0*f
b = qo.dot.f/normq(f)

! derivate di x,y,z rispetto ad "c"
f = p(1,1,1) - p(1,1,0) +&
    p(1,0,1) - p(1,0,0) +&
    p(0,1,1) - p(0,1,0) +&
    p(0,0,1) - p(0,0,0)

f = 0.125d0*f
c = qo.dot.f/normq(f)

a = min(1.0d0,max(-1.0d0,a))
b = min(1.0d0,max(-1.0d0,b))
c = min(1.0d0,max(-1.0d0,c))

return
end
