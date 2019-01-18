!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
      Module Ene_data

      Implicit None

        Integer              :: Nfree,MR
        Integer, allocatable :: Qmode(:)
        Character            :: title*80
        Double Precision     :: E0,dx0,dy0,dz0
        Double Precision, allocatable :: q0(:),delq(:)
        Double Precision, allocatable :: E1(:),E2(:),E3(:)
        Double Precision, allocatable :: dx1(:),dx2(:),dx3(:)
        Double Precision, allocatable :: dy1(:),dy2(:),dy3(:)
        Double Precision, allocatable :: dz1(:),dz2(:),dz3(:)

        Contains

        Double Precision Function E1MR(q,a)

          Implicit None
          Integer, intent(in) :: q,a

          E1MR=E1((q-1)*6+a)

        End Function
        
        Double Precision Function E2MR(qi,qj,a)

          Implicit None
          Integer, intent(in) :: qi,qj,a
          Integer :: qqi,qqj,tmp

          qqi=qi; qqj=qj
          if(qqi<qqj) then
             tmp=qqi; qqi=qqj; qqj=tmp
          endif 

          E2MR=E2(((qqi-1)*(qqi-2)/2+qqj-1)*12+a)

        End Function

        Double Precision Function E3MR(qi,qj,qk,a)

          Implicit None
          Integer, intent(in) :: qi,qj,qk,a
          Integer :: qqi,qqj,qqk,tmp

          qqi=qi; qqj=qj; qqk=qk
          if(qqj<qqk) then
             tmp=qqj; qqj=qqk; qqk=tmp
          endif 
          if(qqi<qqj) then
             tmp=qqi; qqi=qqj; qqj=tmp
          endif 

          E3MR=E3(((qqi-1)*(qqi-2)*(qqi-3)/6+(qqj-1)*(qqj-2)/2+qqk-1)*8+a)

        End Function


      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8


        Subroutine mkhs_main(iopt)

        Implicit None

        Integer :: iopt
        Integer :: i,j,k
        Double Precision :: De,delx2


      !--------------------------------------------------------------
      !  Energy
          if(iopt==0) then 
             Write(10,100)
          else
             Write(10,101)
          endif
          Write(10,'(3x,D18.12)') E0

      100 Format('# Energy / hartree ')
      101 Format('# Dipole / debye ')

      !--------------------------------------------------------------
      !  Geometry
          Write(10,110)
          Write(10,'(3d18.7)') q0

      110 Format('# Geometry / Angs amu1/2')

      !--------------------------------------------------------------
      !  1MR
          Write(10,200) title

          ! ==  gi  ==
          if(iopt==0) then 
             Write(10,210)
          else
             Write(10,211)
          endif
          Do i=1,Nfree
             delx2=delq(i)*2.D+00
             De=(E1MR(i,3)-E1MR(i,4))/delx2
             Write(10,250) Qmode(i),De
          End do

          ! ==  hii  ==
          if(iopt==0) then 
             Write(10,220)
          else
             Write(10,221)
          endif
          Do i=1,Nfree
             delx2=delq(i)*delq(i)
             De=(E1MR(i,3)-2.D+00*E0+E1MR(i,4))/delx2
             Write(10,250) Qmode(i),De
          End do

          ! ==  tiii  ==
          if(iopt==0) then
             Write(10,230)
          else
             Write(10,231)
          endif
          Do i=1,Nfree
             delx2=delq(i)*delq(i)*delq(i)*8.D+00
             De=(E1MR(i,1)-3.D+00*E1MR(i,3)+3.D+00*E1MR(i,4)-E1MR(i,6))/delx2
             Write(10,250) Qmode(i),De
          End do

          ! ==  uiiii  ==
          if(iopt==0) then
             Write(10,240)
          else
             Write(10,241)
          endif
          Do i=1,Nfree
             delx2=delq(i)*delq(i)*delq(i)*delq(i)
             De=(E1MR(i,2)-4D+00*E1MR(i,3)+6.D+00*E0-4.D+00*E1MR(i,4) &
                                          +E1MR(i,5))/delx2
             Write(10,250) Qmode(i),De
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

          if(MR==1) return

      !--------------------------------------------------------------
      !  2MR
      !
          Write(10,300) title

          ! ==  hii  ==
          if(iopt==0) then
             Write(10,310)
          else
             Write(10,311)
          endif
          Do i=2,Nfree
             Do j=1,i-1
                delx2=delq(i)*delq(j)*4.D+00
                De=(E2MR(i,j,2)-E2MR(i,j,5)-E2MR(i,j,8)+E2MR(i,j,11))/delx2
                Write(10,350) Qmode(i),Qmode(j),De
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
                delx2=delq(i)*delq(i)*delq(j)*delq(j)
                De=(E2MR(i,j,2)+E2MR(i,j,5)+E2MR(i,j,8)+E2MR(i,j,11) &
                   -2.D+00*E1MR(i,3)-2.D+00*E1MR(i,4) &
                   -2.D+00*E1MR(j,3)-2.D+00*E1MR(j,4) &
                   +4.D+00*E0)/delx2
                Write(10,350) Qmode(i),Qmode(j),De
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
                delx2=delq(i)*delq(i)*delq(j)*2.D+00
                De=(E2MR(i,j,2)-2.D+00*E1MR(j,3)+E2MR(i,j,8) &
                   -E2MR(i,j,5)+2.D+00*E1MR(j,4)-E2MR(i,j,11))/delx2
                Write(10,350) Qmode(i),Qmode(j),De

                ! ==  tjji  ==
                delx2=delq(i)*delq(j)*delq(j)*2.D+00
                De=(E2MR(i,j,2)-2.D+00*E1MR(i,3)+E2MR(i,j,5) &
                   -E2MR(i,j,8)+2.D+00*E1MR(i,4)-E2MR(i,j,11))/delx2
                Write(10,350) Qmode(j),Qmode(i),De

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
                delx2=delq(i)*delq(i)*delq(i)*delq(j)*16.D+00
                De=(E2MR(i,j,1)-3.D+00*E2MR(i,j,2) &
                   +3.D+00*E2MR(i,j,8)-E2MR(i,j,9) &
                   -E2MR(i,j,4)+3.D+00*E2MR(i,j,5) &
                   -3.D+00*E2MR(i,j,11)+E2MR(i,j,12))/delx2
                Write(10,350) Qmode(i),Qmode(j),De

                ! ==  ujjji  ==
                delx2=delq(i)*delq(j)*delq(j)*delq(j)*16.D+00
                De=(E2MR(i,j,3)-3.D+00*E2MR(i,j,2) &
                   +3.D+00*E2MR(i,j,5)-E2MR(i,j,6) &
                   -E2MR(i,j,7)+3.D+00*E2MR(i,j,8) &
                   -3.D+00*E2MR(i,j,11)+E2MR(i,j,10))/delx2
                Write(10,350) Qmode(j),Qmode(i),De

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

          if(MR==2) return

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

             delx2=delq(i)*delq(j)*delq(k)*8.D+00
             De=(E3MR(i,j,k,1)-E3MR(i,j,k,2)-E3MR(i,j,k,3) &
                +E3MR(i,j,k,4)-E3MR(i,j,k,5)+E3MR(i,j,k,6) &
                +E3MR(i,j,k,7)-E3MR(i,j,k,8))/delx2
             Write(10,430) Qmode(i),Qmode(j),Qmode(k),De

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
             delx2=delq(i)*delq(i)*delq(j)*delq(k)*4.D+00
             De=(E3MR(i,j,k,1)+E3MR(i,j,k,2)-E3MR(i,j,k,3) &
                -E3MR(i,j,k,4)-E3MR(i,j,k,5)-E3MR(i,j,k,6) &
                +E3MR(i,j,k,7)+E3MR(i,j,k,8) &
                -2.D+00*(E2MR(j,k,2)-E2MR(j,k,5) &
                -E2MR(j,k,8)+E2MR(j,k,11)))/delx2
             Write(10,430) Qmode(i),Qmode(j),Qmode(k),De

             ! ==  uijjk  ==
             delx2=delq(i)*delq(j)*delq(j)*delq(k)*4.D+00
             De=(E3MR(i,j,k,1)-E3MR(i,j,k,2)+E3MR(i,j,k,3) &
                -E3MR(i,j,k,4)-E3MR(i,j,k,5)+E3MR(i,j,k,6) &
                -E3MR(i,j,k,7)+E3MR(i,j,k,8) &
                -2.D+00*(E2MR(i,k,2)-E2MR(i,k,5) &
                -E2MR(i,k,8)+E2MR(i,k,11)))/delx2
             Write(10,430) Qmode(j),Qmode(k),Qmode(i),De

             ! ==  uijkk  ==
             delx2=delq(i)*delq(j)*delq(k)*delq(k)*4.D+00
             De=(E3MR(i,j,k,1)-E3MR(i,j,k,2)-E3MR(i,j,k,3) &
                +E3MR(i,j,k,4)+E3MR(i,j,k,5)-E3MR(i,j,k,6) &
                -E3MR(i,j,k,7)+E3MR(i,j,k,8) &
                -2.D+00*(E2MR(i,j,2)-E2MR(i,j,5) &
                -E2MR(i,j,8)+E2MR(i,j,11)))/delx2
             Write(10,430) Qmode(k),Qmode(i),Qmode(j),De

          End do
          End do
          End do
 
      400 Format('# 3MR ',a)
      410 Format('# Cubic(i,j,k) / hartree Angs^-3 amu^-3/2 ')
      411 Format('# Cubic(i,j,k) / debye Angs^-3 amu^-3/2 ')
      420 Format('# Quartic(i,i,j,k) / hartree Angs^-4 amu^-2 ')
      421 Format('# Quartic(i,i,j,k) / debye Angs^-4 amu^-2 ')
      430 Format(3i4,D20.10)

          if(MR==3) return

 
      End Subroutine
 
      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8


      End module
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
!
      Subroutine read_eBase_data(MR_in,N,del,qq,qm,tl)

      Use Ene_data

      Integer :: MR_in,N,Msym(N),qm(N)
      Real(8) :: del(N),qq(N)

      Integer :: i,j,k,l,m,ii,jj
      Character :: Sym*3,tl*80

!--------------------------------------------------------------
!
! --    Nfree,MR,delq,q0
        Nfree=N
        MR=MR_in

        Allocate(delq(Nfree))
        delq=del

        Allocate(q0(Nfree))
        q0=qq

        Allocate(Qmode(Nfree))
        Qmode=qm

        title=tl

!       Restore- E,dx,dy,dz
!        Allocate(E1(6*Nfree),E2(6*Nfree*(Nfree-1)),E3(4*Nfree*(Nfree-1)*(Nfree-2)/3))
!        Allocate(dx1(6*Nfree),dx2(6*Nfree*(Nfree-1)),dx3(4*Nfree*(Nfree-1)*(Nfree-2)/3))
!        Allocate(dy1(6*Nfree),dy2(6*Nfree*(Nfree-1)),dy3(4*Nfree*(Nfree-1)*(Nfree-2)/3))
!        Allocate(dz1(6*Nfree),dz2(6*Nfree*(Nfree-1)),dz3(4*Nfree*(Nfree-1)*(Nfree-2)/3))

        Open(100,file='tmp1',status='old')

        Read(100,100,err=1000) E0,dx0,dy0,dz0
  100   Format(4x,f22.12,3f14.8)

        Write(6,'('' -> Equilibrium geometry'')')
        Write(6,100) E0,dx0,dy0,dz0

        !  --    E(1MR)   --
        !  E1=  +3di,       E4=  -1di,
        !  E2=  +2di,       E5=  -2di,
        !  E3=  +1di,       E6=  -3di.
 
        Allocate(E1(6*Nfree),dx1(6*Nfree),dy1(6*Nfree),dz1(6*Nfree))
        E1=0.D+00
        Do i=1,Nfree
           Read(100,*,err=1000)
           Do j=1,6
              ii=(i-1)*6+j
              Read(100,100,err=1000) E1(ii),dx1(ii),dy1(ii),dz1(ii)
           End do
        End do

        !Write(6,'('' -> 1MR'')')
        Do i=1,Nfree
           Write(6,'(''  > MODE='',i4)') qm(i)
           Do j=1,6
              ii=(i-1)*6+j
              Write(6,100) E1(ii),dx1(ii),dy1(ii),dz1(ii)
           End do
        End do

        if(MR==1) goto 200


        !  --    E(2MR)   --
        !  E1=(+3di,+1dj), E4=(+3di,-1dj), E7=(-1di,+3dj), E10=(-1di,-3dj),
        !  E2=(+1di,+1dj), E5=(+1di,-1dj), E8=(-1di,+1dj), E11=(-1di,-1dj),
        !  E3=(+1di,+3dj), E6=(+1di,-3dj), E9=(-3di,+1dj), E12=(-3di,-1dj).

        Allocate(E2(6*Nfree*(Nfree-1)))
        Allocate(dx2(6*Nfree*(Nfree-1)))
        Allocate(dy2(6*Nfree*(Nfree-1)))
        Allocate(dz2(6*Nfree*(Nfree-1)))
        E2=0.D+00
        l=0
        Do i=1,Nfree
           Do j=1,i-1
              Read(100,*,err=1000)
              Do k=1,12
                 ii=l+k
                 Read(100,100,err=1000) E2(ii),dx2(ii),dy2(ii),dz2(ii)
              End do
              l=l+12
           End do
        End do

        !Write(6,'('' -> 2MR'')')
        l=0
        Do i=1,Nfree
           Do j=1,i-1
              Write(6,'(''  > MODE='',2i4)') qm(i),qm(j)
              Do k=1,12
                 ii=l+k
                 Write(6,100) E2(ii),dx2(ii),dy2(ii),dz2(ii)
              End do
              l=l+12
           End do
        End do

        if(MR==2) goto 200

        !   --    E(3MR)   --
        !   E1=(+di,+dj,+dk), E3=(+di,-dj,+dk), E5=(+di,+dj,-dk), E7=(+di,-dj,-dk),
        !   E2=(-di,+dj,+dk), E4=(-di,-dj,+dk), E6=(-di,+dj,-dk), E8=(-di,-dj,-dk).

        Allocate(E3(4*Nfree*(Nfree-1)*(Nfree-2)/3))
        Allocate(dx3(4*Nfree*(Nfree-1)*(Nfree-2)/3))
        Allocate(dy3(4*Nfree*(Nfree-1)*(Nfree-2)/3))
        Allocate(dz3(4*Nfree*(Nfree-1)*(Nfree-2)/3))
        E3=0.D+00
        m=0
        Do i=1,Nfree
           Do j=1,i-1
              Do k=1,j-1
                 Read(100,*,err=1000)
                 Do l=1,8
                    ii=m+l
                    Read(100,100,err=1000) E3(ii),dx3(ii),dy3(ii),dz3(ii)
                 End do
                 m=m+8
              End do
           End do
        End do

        !Write(6,'('' -> 3MR'')')
        m=0
        Do i=3,Nfree
           Do j=2,i-1
              Do k=1,j-1
                 Write(6,'(''  > MODE='',3i4)') qm(i),qm(j),qm(k)
                 Do l=1,8
                    ii=m+l
                    Write(6,100) E3(ii),dx3(ii),dy3(ii),dz3(ii)
                 End do
                 m=m+8
              End do
           End do
        End do
!
!--------------------------------------------------------------

    200 Continue
        Close(100)

!        Call System('rm tmp1')
        Return

   1000 Continue
        Write(6,*) 'Error while reading tmp1'
        Stop

      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
      Subroutine eBase_hs(Is,dpl)

      Use Ene_data

      Implicit None

      Integer :: Is(3)
      Logical :: dpl

         Open(10,file='out.hs',status='unknown')
         Call mkhs_main(0)
         Close(10)
         if(.not. dpl) return

         E0=dx0; E1=dx1; E2=dx2; E3=dx3
         Select case(Is(1))
            Case(1)
               Open(10,file='outx.hs',status='unknown')
            Case(2)
               Open(10,file='outy.hs',status='unknown')
            Case(3)
               Open(10,file='outz.hs',status='unknown')
         End select
         Call mkhs_main(1)
         Close(10)

         E0=dy0; E1=dy1; E2=dy2; E3=dy3
         Select case(Is(2))
            Case(1)
               Open(10,file='outx.hs',status='unknown')
            Case(2)
               Open(10,file='outy.hs',status='unknown')
            Case(3)
               Open(10,file='outz.hs',status='unknown')
         End select
         Call mkhs_main(1)
         Close(10)

         E0=dz0; E1=dz1; E2=dz2; E3=dz3
         Select case(Is(3))
            Case(1)
               Open(10,file='outx.hs',status='unknown')
            Case(2)
               Open(10,file='outy.hs',status='unknown')
            Case(3)
               Open(10,file='outz.hs',status='unknown')
         End select
         Call mkhs_main(1)
         Close(10)

      End
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
