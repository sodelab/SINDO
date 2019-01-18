!---------------------------------------------------------------------
!
!      Based on: 
!          // Press, Teukolsky, Vetterling, and Flannery
!          // Numerical Recipes, 2nd edition,
!          // Cambridge University Press (1992)
!
!      Extended to 3D version by KY
!
!---------------------------------------------------------------------
!
!       Program test
!
!       Implicit None
!       Integer :: n1,n2,n3
!       Integer :: i,j,k,nn
!       Real(8) :: x1,x2,x3,z, xl,xh,dx
!       Real(8), dimension(:), allocatable :: q1,q2,q3
!       Real(8), dimension(:,:,:), allocatable :: V0,Vs
!
!          Open(10,file='q3q2q1s.pot',status='old')
!          Read(10,*) 
!          Read(10,*) n3,n2,n1
!          Allocate(q1(n1),q2(n2),q3(n3),V0(n3,n2,n1),Vs(n3,n2,n1))
!          Read(10,*) 
!          Do i=1,n1
!          Do j=1,n2
!          Do k=1,n3
!             Read(10,*) q3(k),q2(j),q1(i),V0(k,j,i)
!          End do
!          End do
!          End do
!          Call DSplie3(q3,q2,q1,V0,n3,n2,n1,Vs)
!
!         !test1
!          x3=-0.52826035; x2=-0.34556191; x1=-0.34007271
!          Call DSplin3(q3,q2,q1,V0,Vs,n3,n2,n1,x3,x2,x1,z)
!          Write(6,'(3f8.3,f12.6)') x3,x2,x1,z
!
!         !test2
!          xl=-0.34D+00; xh=0.34D+00
!          nn=31
!          dx=(xh-xl)/(nn-1)
!          x3=-0.52826035D+00; x2=-0.26217983D+00
!          Do i=1,nn
!             x1=xl+dble(i-1)*dx 
!             Call DSplin3(q3,q2,q1,V0,Vs,n3,n2,n1,x3,x2,x1,z)
!             Write(6,'(3f8.3,f12.6)') x3,x2,x1,z
!          End do
!
!          Deallocate(q3,q2,q1,V0,Vs)
!
!       End
!
!---------------------------------------------------------------------
!
! ----     Three-dimensional spline interpolation routines
!
!---------------------------------------------------------------------

       Subroutine DSplie3(x3a,x2a,x1a,ya,m3,m2,m1,y2a)

       Implicit None

         Integer :: m1,m2,m3
         Real(8) :: x3a(m3),x2a(m2),x1a(m1),ya(m3,m2,m1),y2a(m3,m2,m1)

         Integer :: j

            Do j=1,m1
               Call DSplie2(x3a,x2a,ya(:,:,j),m3,m2,y2a(:,:,j))
            End do

       End subroutine

!---------------------------------------------------------------------

       Subroutine DSplin3(x3a,x2a,x1a,ya,y2a,m3,m2,m1,x3,x2,x1,y)

       Implicit None

         Integer :: m1,m2,m3
         Real(8) :: x3a(m3),x2a(m2),x1a(m1),ya(m3,m2,m1),y2a(m3,m2,m1)
         Real(8) :: x1,x2,x3,y

         Integer :: j,k
         Real(8) :: ytmp(m1),y2tmp(m1)

            Do j=1,m1
               Call DSplin2(x3a,x2a,ya(:,:,j),y2a(:,:,j), &
                            m3,m2,x3,x2,ytmp(j))
            End do 
            Call DSpline(x1a,ytmp,m1,1.D+32,1.D+32,y2tmp)
            Call DSplint(x1a,ytmp,y2tmp,m1,x1,y)

       End subroutine

!---------------------------------------------------------------------
