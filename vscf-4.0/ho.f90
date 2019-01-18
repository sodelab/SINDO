!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine HO_xmat(n,omega,x)

      Implicit None

         Integer :: n,j,k
         Real(8) :: omega,x(0:n-1,0:n-1)

         Real(8) :: PI,const
         Real(8), parameter :: h=6.62606896D-34,    &  ! Planck's constant in J.s
                     c=2.99792458D+10, &  ! velocity of light in cm.s-1
                     Na=6.02214D23          ! Avogadoro's constant

        !-------------------------------------
        ! >>  Qmat in Angs(amu)1/2
        !-------------------------------------

         PI=Acos(-1.D+00)
         const=SQRT(h*Na*1.D+23/c/omega)/(2.D+00*PI)

         x=0.D+00
         Do j=1,n-1
            k=j-1
            x(k,j)=Sqrt(dble(j)*0.5D+00)*const
            x(j,k)=x(k,j)
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine HO_kinmat(nx,omega,T)

      Implicit None

         Integer :: i,j,k,l
         Real(8) :: omgh
         Real(8), dimension(0:nx,0:nx) :: T

         Real(8), parameter :: H2wvn=2.194746D+05

         Integer:: nx
         Real(8) :: omega

         omgh=omega/H2wvn*0.5D+00

         T=0.D+00
         Do i=0,nx

            T(i,i)=(dble(i)+0.5D+00)*omgh
            if(i<nx-1) then
               T(i+2,i)=-0.5D+00*SQRT(dble((i+1)*(i+2)))*omgh
               T(i,i+2)=T(i+2,i)
            endif

         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
