!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
      Module Hess_data

      Implicit None

        Integer  :: Nfree
        Character            :: title*80
        Double Precision :: E0,dx0,dy0,dz0
        Double Precision, allocatable :: q0(:),delq(:)
        Double Precision, allocatable :: Ene(:),Grad(:,:),Hess(:,:)

        Contains

        Double Precision Function pHcmp(im,jm,km)

          Implicit None
          Integer :: im,jm,km

          pHcmp=Hess(Nfree*(im-1)+jm,2*km-1)

        End Function
        
        Double Precision Function mHcmp(im,jm,km)

          Implicit None
          Integer :: im,jm,km

          mHcmp=Hess(Nfree*(im-1)+jm,2*km)

        End Function

        Double Precision Function H0cmp(im,jm)

          Implicit None
          Integer :: im,jm

          H0cmp=Hess(Nfree*(im-1)+jm,0)

        End Function


      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8

      Subroutine mkhs_main(iopt)

      Implicit None

      Integer :: iopt
      Integer :: i,j,k
      Double Precision :: De,De1,De2,De3,delx2


      !--------------------------------------------------------------
      !  Energy
      !
          if(iopt==0) then 
             Write(10,100)
          else
             Write(10,101)
          endif
          Write(10,'(3x,D18.12)') Ene(0)

      100 Format('# Energy / hartree ')
      101 Format('# Dipole / debye ')

      !--------------------------------------------------------------
      !  Geometry
      !
          Write(10,110)
          Write(10,'(3d18.7)') q0

      110 Format('# Geometry / Angs amu1/2')

      !--------------------------------------------------------------
      !  1MR
      !
          Write(10,200) title

          ! ==  gi  ==
          if(iopt==0) then 
             Write(10,210)
          else
             Write(10,211)
          endif
          Do i=1,Nfree
             Write(10,250) i,Grad(i,0)
          End do

          ! ==  hii  ==
          if(iopt==0) then 
             Write(10,220)
          else
             Write(10,221)
          endif
          Do i=1,Nfree
             Write(10,250) i,H0cmp(i,i)
          End do

          ! ==  tiii  ==
          if(iopt==0) then
             Write(10,230)
          else
             Write(10,231)
          endif
          Do i=1,Nfree
             delx2=delq(i)*2.D+00
             De=(pHcmp(i,i,i)-mHcmp(i,i,i))/delx2
             Write(10,250) i,De
          End do

          ! ==  uiiii  ==
          if(iopt==0) then
             Write(10,240)
          else
             Write(10,241)
          endif
          Do i=1,Nfree
             delx2=delq(i)*delq(i)
             De=(pHcmp(i,i,i)+mHcmp(i,i,i)-2.D+00*H0cmp(i,i))/delx2
             Write(10,250) i,De
          End do

      200 Format('# 1MR ',a)
      210 Format('# Gradient / hartree Angs^-1 amu^-1/2 ')
      211 Format('# Gradient / debye Angs^-1 amu^-1/2 ')
      220 Format('# Hessian(i,i) / hartree Angs^-2 amu^-1 ')
      221 Format('# Hessian(i,i) / debye Angs^-2 amu^-1 ')
      230 Format('# Cubic(i,i,i) / hartree Angs^-3 amu^-3/2 ')
      231 Format('# Cubic(i,i,i) / debye Angs^-3 amu^-3/2 ')
      240 Format('# Quartic(i,i,i,i) / hartree Angs^-4 amu^-2 ')
      241 Format('# Quartic(i,i,i,i) / debye Angs^-4 amu^-2 ')
      250 Format(i4,D20.10)


      !--------------------------------------------------------------
      !  2MR
      !
          Write(10,300) title 

          ! ==  hij  ==
          if(iopt==0) then
             Write(10,310)
          else
             Write(10,311)
          endif
          Do i=2,Nfree
             Do j=1,i-1
                Write(10,350) i,j,H0cmp(j,i)
             End do
          End do
 
          ! ==  uiijj  ==
          if(iopt==0) then
             Write(10,320)
          else
             Write(10,321)
          endif
          Do i=2,Nfree
             Do j=1,i-1
                delx2=delq(j)*delq(j)
                De1=(pHcmp(i,i,j)+mHcmp(i,i,j)-2.D+00*H0cmp(i,i))/delx2
                delx2=delq(i)*delq(i)
                De2=(pHcmp(j,j,i)+mHcmp(j,j,i)-2.D+00*H0cmp(j,j))/delx2
                De=(De1+De2)*0.5D+00
                Write(10,350) i,j,De
             End do
          End do

          if(iopt==0) then
             Write(10,330)
          else
             Write(10,331)
          endif
          Do i=2,Nfree
             Do j=1,i-1

                ! ==  tiij  ==
                delx2=delq(j)*2.D+00
                De1=(pHcmp(i,i,j)-mHcmp(i,i,j))/delx2
                delx2=delq(i)*2.D+00
                De2=(pHcmp(i,j,i)-mHcmp(i,j,i))/delx2
                De=(De1+De2*2.D+00)/3.D+00
                Write(10,350) i,j,De

                ! ==  tjji  ==
                delx2=delq(i)*2.D+00
                De1=(pHcmp(j,j,i)-mHcmp(j,j,i))/delx2
                delx2=delq(j)*2.D+00
                De2=(pHcmp(i,j,j)-mHcmp(i,j,j))/delx2
                De=(De1+De2*2.D+00)/3.D+00
                Write(10,350) j,i,De

             End do
          End do

          if(iopt==0) then
             Write(10,340)
          else
             Write(10,341)
          endif

          Do i=2,Nfree
             Do j=1,i-1

                ! ==  uiiij  ==
                delx2=delq(i)*delq(i)
                De=(pHcmp(i,j,i)+mHcmp(i,j,i)-H0cmp(i,j)*2.D+00)/delx2
                Write(10,350) i,j,De

                ! ==  ujjji  ==
                delx2=delq(j)*delq(j)
                De=(pHcmp(i,j,j)+mHcmp(i,j,j)-H0cmp(i,j)*2.D+00)/delx2
                Write(10,350) j,i,De

             End do
          End do
 
      300 Format('# 2MR ',a)
      310 Format('# Hessian(i,j) / hartree Angs^-2 amu^-1 ')
      311 Format('# Hessian(i,j) / debye Angs^-2 amu^-1 ')
      320 Format('# Quartic(i,i,j,j) / hartree Angs^-4 amu^-2 ')
      321 Format('# Quartic(i,i,j,j) / debye Angs^-4 amu^-2 ')
      330 Format('# Cubic(i,i,j) / hartree Angs^-3 amu^-3/2 ')
      331 Format('# Cubic(i,i,j) / debye Angs^-3 amu^-3/2 ')
      340 Format('# Quartic(i,i,i,j) / hartree Angs^-4 amu^-2 ')
      341 Format('# Quartic(i,i,i,j) / debye Angs^-4 amu^-2 ')
      350 Format(2i4,D20.10)


      !--------------------------------------------------------------
      !  3MR
      !
          Write(10,400) title

          ! ==  tijk  ==
          if(iopt==0) then
             Write(10,410)
          else
             Write(10,411)
          endif
          Do i=3,Nfree
          Do j=2,i-1
          Do k=1,j-1

             delx2=delq(i)*2.D+00
             De1=(pHcmp(j,k,i)-mHcmp(j,k,i))/delx2
             delx2=delq(j)*2.D+00
             De2=(pHcmp(k,i,j)-mHcmp(k,i,j))/delx2
             delx2=delq(k)*2.D+00
             De3=(pHcmp(i,j,k)-mHcmp(i,j,k))/delx2
             De=(De1+De2+De3)/3.D+00
             Write(10,430) i,j,k,De

          End do
          End do
          End do

          if(iopt==0) then
             Write(10,420)
          else
             Write(10,421)
          endif

          Do i=3,Nfree
          Do j=2,i-1
          Do k=1,j-1

             ! ==  uiijk  ==
             delx2=delq(i)*delq(i)
             De=(pHcmp(j,k,i)+mHcmp(j,k,i)-2.D+00*H0cmp(j,k))/delx2
             Write(10,430) i,j,k,De

             ! ==  uijjk  ==
             delx2=delq(j)*delq(j)
             De=(pHcmp(i,k,j)+mHcmp(i,k,j)-2.D+00*H0cmp(i,k))/delx2
             Write(10,430) j,k,i,De

             ! ==  uijkk  ==
             delx2=delq(k)*delq(k)
             De=(pHcmp(i,j,k)+mHcmp(i,j,k)-2.D+00*H0cmp(i,j))/delx2
             Write(10,430) k,i,j,De

          End do
          End do
          End do
 
      400 Format('# 3MR ',a)
      410 Format('# Cubic(i,j,k) / hartree Angs^-3 amu^-3/2 ')
      411 Format('# Cubic(i,j,k) / debye Angs^-3 amu^-3/2 ')
      420 Format('# Quartic(i,i,j,k) / hartree Angs^-4 amu^-2 ')
      421 Format('# Quartic(i,i,j,k) / debye Angs^-4 amu^-2 ')
      430 Format(3i4,D20.10)

 
      End Subroutine
 
      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8


      End module
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
!
      Subroutine read_hBase_data(N,del,qq,tl)

      Use Hess_data

      Integer :: N
      Real(8) :: del(N),qq(N)
      Character :: tl*80

      Integer :: i,j,k,l,m,ii,jj

!--------------------------------------------------------------
!
! --    Nfree,delq,q0
        Nfree=N

        Allocate(delq(Nfree))
        delq=del

        Allocate(q0(Nfree))
        q0=qq

        title=tl

        Write(6,'('' -> Ene, Grad, and Hess at each point'')')
!       Restore- Ene,Grad,Hess

        ii=2*N
        jj=N*N
        Allocate(Ene(0:ii),Grad(N,0:ii),Hess(jj,0:ii))

        Open(100,file='tmp1',status='old')

        ! Equilibrium Geometry
        Read(100,*,err=1000) i,Ene(0)
        Read(100,*) 
        Read(100,*) Grad(:,0)
        Read(100,*) 
        Read(100,*) Hess(:,0)

        Write(6,'(''  > Equilibrium geometry'')')
        Call writeData(N,Ene(0),Grad(:,0),Hess(:,0))
        Write(6,*)

        Do i=1,Nfree
           Read(100,*,err=1000)
           Do j=1,2
              l=2*(i-1)+j
              Read(100,*) k,Ene(l)
              Read(100,*) 
              Read(100,*) Grad(:,l)
              Read(100,*) 
              Read(100,*) Hess(:,l)

              Write(6,'(''  > MODE '',i3,'' ('',i1,'')'')') i,j
              Call writeData(N,Ene(l),Grad(:,l),Hess(:,l))

           End do
           Write(6,*)
        End do

        Close(100)

        Return

   1000 Continue
        Write(6,*) 'Error while reading tmp1'
        Stop

        Contains

        Subroutine writeData(N,E,G,H)

        Implicit None

        Integer :: i,N
        Real(8) :: E,G(N),H(N*N)

              Write(6,100) E
              Write(6,110) 
              Write(6,115) G
              Write(6,120) 
              if(N>3) then
                 Do i=1,N
                    Write(6,125) i,H((i-1)*N+1:(i-1)*N+3)
                    Write(6,126) H((i-1)*N+4:i*N)
                 End do
              else
                 Do i=1,N
                    Write(6,125) i,H((i-1)*N+1:i*N)
                 End do
              endif
          100 Format(4x,'Energy:  ',f18.12)
          110 Format(4x,'Gradient:')
          115 Format(8x,3f12.6)
          120 Format(4x,'Hessian:')
          125 Format(4x,i4,3f12.6)
          126 Format(8x,3f12.6)

        End subroutine

      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
      Subroutine hBase_hs

      Use Hess_data

      Implicit None

         Open(10,file='out.hs',status='unknown')
         Call mkhs_main(0)
         Close(10)

      End
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
