!---------------------------------------------------------------------
!
!      Taken from: 
!          // Press, Teukolsky, Vetterling, and Flannery
!          // Numerical Recipes, 2nd edition,
!          // Cambridge University Press (1992)
!
!
!---------------------------------------------------------------------
!
!       Program test
!
!       Implicit None
!       Integer :: nx,ny
!       Integer :: i,j,nn
!       Real(8) :: xx,yy,z, xl,xh,yl,yh,dx,dy
!       Real(8), dimension(:), allocatable :: qx,qy
!       Real(8), dimension(:,:), allocatable :: Vint,Vs
!
!          Open(10,file='q3q2s.pot',status='old')
!          Read(10,*) 
!          Read(10,*) nx,ny
!          Allocate(qx(nx),qy(ny),Vint(nx,ny),Vs(nx,ny))
!          Read(10,*) 
!          Do i=1,ny
!          Do j=1,nx
!             Read(10,*) qx(j),qy(i),Vint(j,i)
!          End do
!          End do
!          !Do i=1,ny
!          !Do j=1,nx
!          !   write(6,'(2f8.3,f12.6)') qx(j),qy(i),Vint(j,i)
!          !End do
!          !End do
!          Call DSplie2(qx,qy,Vint,nx,ny,Vs)
!
!          xl=-0.34D+00; xh=0.34D+00
!          nn=31
!          dx=(xh-xl)/(nn-1)
!          yy=-0.187808D+00
!
!          Do i=1,nn
!             xx=xl+dble(i-1)*dx 
!             Call DSplin2(qx,qy,Vint,Vs,nx,ny,xx,yy,z)
!             Write(6,'(2f8.3,f12.6)') xx,yy,z
!          End do
!
!          Deallocate(qx,qy,Vint,Vs)
!
!       End
!
!---------------------------------------------------------------------
!
! ----     Two-dimensional spline interpolation routines
!
!---------------------------------------------------------------------

       Subroutine DSplie2(x1a,x2a,ya,m,n,y2a)

       Implicit None

         Integer :: m,n
         Real(8) :: x1a(m),x2a(n),ya(m,n),y2a(m,n)

         Integer :: j

            Do j=1,n
               Call DSpline(x1a,ya(:,j),m,1.D+32,1.D+32,y2a(:,j))
            End do

!            Do j=1,m
!               Call DSpline(x2a,ya(j,:),n,1.D+32,1.D+32,y2a(j,:))
!            End do

       End subroutine

!---------------------------------------------------------------------

       Subroutine DSplin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)

       Implicit None

         Integer :: m,n
         Real(8) :: x1a(m),x2a(n),ya(m,n),y2a(m,n)
         Real(8) :: x1,x2,y

         Integer :: j,k
         Real(8) :: ytmp(n),y2tmp(n)

            Do j=1,n
               Call DSplint(x1a,ya(:,j),y2a(:,j),m,x1,ytmp(j))
            End do 
            Call DSpline(x2a,ytmp,n,1.D+32,1.D+32,y2tmp)
            Call DSplint(x2a,ytmp,y2tmp,n,x2,y)

!         Real(8) :: ytmp(m),y2tmp(m)
!
!            Do j=1,m
!               Call DSplint(x2a,ya(j,:),y2a(j,:),n,x2,ytmp(j))
!            End do 
!            Call DSpline(x1a,ytmp,n,1.D+32,1.D+32,y2tmp)
!            Call DSplint(x1a,ytmp,y2tmp,m,x1,y)

       End subroutine

!---------------------------------------------------------------------
