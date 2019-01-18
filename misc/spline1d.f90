!---+----1----+----2----+----3----+----4----+----5----+----6----+----72!
!
!      Taken from: 
!          // Press, Teukolsky, Vetterling, and Flannery
!          // Numerical Recipes, 2nd edition,
!          // Cambridge University Press (1992)
!
!     ------------------------------------------------------------------
!
!      Program test
!
!      Integer :: mm,nn
!      Real(8) :: qqi,dqi,Vint
!      Real(8),dimension(:), allocatable :: qi,Vi,Vs
!
!         Open(10,file='q1s.pot',status='old')
!         Read(10,*) 
!         Read(10,*) nn
!         Allocate(qi(nn),Vi(nn))
!         Read(10,*) 
!         Do i=1,nn
!            Read(10,*) qi(i),Vi(i) 
!         End do 
!         Close(10)
!
!         mm=31
!         dqi=(qi(nn)-qi(1))/dble(mm-1)
!
!         Allocate(Vs(nn))
!         Call DSpline(qi,Vi,nn,1.D+32,1.D+32,Vs)
!
!         Do i=1,nn
!            Call DSplint(qi,Vi,Vs,nn,qi(i),Vint)
!            Write(6,'(f8.5,f12.6)') qi(i),Vint
!         End do
!         Write(6,*)
!
!         Do i=1,mm
!            qqi=qi(1)+dble(i-1)*dqi
!            Call DSplint(qi,Vi,Vs,nn,qqi,Vint)
!            Write(6,'(f8.5,f12.6)') qqi,Vint
!         End do
!         Deallocate(qi,Vi,Vs)
!
!      End
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----72!
!
!     Given arrays x and y of length nsz containing a tabulated 
!     function, i.e., yi=f(xi), with x1 < x2 < ... < xnsz, and given 
!     values yp1 and ypn for the first derivative of the interpolating 
!     function at points 1 and nsz, respectively, this routine returns 
!     an array of y2 of length nsz that contains the second derivatives 
!     of the interpolating function at the tabulated points xi.  If yp1 
!     and/or ypn are equal to 1*10^30 or larger, the routine is signaled 
!     to set the corresponding boundary condition for a natural spline, 
!     with zero second derivative on that boundary.
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----72!

      Subroutine DSpline(x,y,nsz,yp1,ypn,y2)

      Implicit None

        Integer :: nsz
        Real(8) :: yp1,ypn,x(nsz),y(nsz),y2(nsz)

        Integer :: i,k
        Real(8) :: p,qn,sig,un,u(nsz)

          if(yp1 > 1.D+30) then
             ! natural cubic spline
             y2(1)=0.D+00 
             u(1)=0.D+00
          else
             ! a specified first derivative
             y2(1)=-0.5D+00
             u(1)=(3.D+00/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
          endif

          Do i=2,nsz-1
             sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
             p=sig*y2(i-1)+2.D+00
             y2(i)=(sig-1.D+00)/p
             u(i)=(6.D+00*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
                  /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
          End do

          if(ypn > 1.D+30) then
             ! natural cubic spline
             qn=0.D+00
             un=0.D+00
          else
             ! a specified first derivative
             qn=0.5D+00
             un=(3.D+00/(x(nsz)-x(nsz-1)))* &
                (ypn-(y(nsz)-y(nsz-1))/(x(nsz)-x(nsz-1)))
          endif

          y2(nsz)=(un-qn*u(nsz-1))/(qn*y2(nsz-1)+1.D+00)
          Do k=nsz-1,1,-1
             y2(k)=y2(k)*y2(k+1)+u(k)
          End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----72!
!
!     Given the arrays xa and ya, which tabulate a function (with the 
!     xai's in increasing or decreasing order), and given the array y2a,
!     which is the output from spline above, and given a value of x, 
!     this routine returns a cubic-spline interpolated value.  The 
!     arrays xa, ya and y2a are all of the same size.
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----72!

      Subroutine DSplint(xa,ya,y2a,nsz,x,y)

      Implicit None

        Integer :: nsz
        Real(8) :: x,y,xa(nsz),ya(nsz),y2a(nsz)

        Integer :: k,khi,klo
        Real(8) :: a,b,h

          klo=1
          khi=nsz

        1 Continue
          if(khi-klo>1) then
             k=(khi+klo)/2
             if(xa(k)>x) then 
                khi=k
             else
                klo=k
             endif
             goto 1
          endif

          h=xa(khi)-xa(klo)
          if(h == 0.D+00) stop 'bad xa input in splint'
          a=(xa(khi)-x)/h
          b=(x-xa(klo))/h
          y=a*ya(klo)+b*ya(khi)+ &
                   ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h*h)/6.D+00

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----72!
