!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Program main

        Implicit None

        Integer,parameter :: max_Nfree=50,max_nho=33,max_NP=9

        Integer :: i,j,k,l
        Integer :: ifl

        Integer :: Nfree,MR,NP
        Integer, dimension(max_Nfree) :: mode,nho
        Real(8) :: P0(max_NP),omega(max_Nfree),quad(max_nho,max_Nfree)
        Character :: title*30

          Call readParams(Nfree,MR,NP,mode,nho,omega)
          if(NP>max_NP) then
             write(6,*) 'Maximum NP is ',max_NP
             stop
          endif

          Call DVRGrids(Nfree,nho,max_nho,omega,quad)

          Write(6,*) 'Title?'
          Read(5,'(A)') title

          ifl=10
          Open(ifl,file='rpt2.dat',status='old')
          Read(ifl,*) i,P0(1:NP)

          Do i=1,Nfree
             Call mkq1pot(ifl,NP,mode(i),nho(i),quad(:,i),P0,title)
          End do
          if(MR==1) goto 100

          write(6,*) 
          write(6,*) 'ERROR: MR>1 is not ready'
          write(6,*) 'ERROR: Only the data for MR=1 is written'
          goto 100

          Do i=2,Nfree
          Do j=1,i-1
             Call mkq2pot(ifl,NP,mode(i),mode(j),nho(i),nho(j), & 
                          quad(:,i),quad(:,j),P0,title)
          End do
          End do
          if(MR==2) goto 100

          Do i=3,Nfree
          Do j=2,i-1
          Do k=1,j-1
             Call mkq3pot(ifl,NP,mode(i),mode(j),mode(k),nho(i),nho(j),nho(k), &
                       quad(:,i),quad(:,j),quad(:,k),P0,title)
          End do
          End do
          End do

      100 Continue
          Close(ifl)

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine mkq1pot(ifl,NP,mode,nho,qq,P0,title)

        Implicit None

        Integer :: i,j,k,l,ifl,NP,nho,mode,num
        Real(8) :: qq(nho),P0(NP),P1(NP,nho)
        Character :: title*30

        Character :: im,im2*2,fl*20

          Do i=1,nho
             Read(ifl,*) num,P1(:,i)
             !P1(:,i)=P1(:,i)-P0
          End do

          Do k=1,NP
             Call get_fname1(mode,fl)
             l=Len_Trim(fl)
             write(fl(l-3:),'(''-'',i1,''.pot'')') k
             Open(11,file=fl,status='new')
             Write(11,'(A)') title
             Write(11,'(''# Number of grids'')') 
             Write(11,*) nho
             Write(11,'(''# qi,P'',i2.2,'' (i='',i4,'')'')') k,mode
             Do i=1,nho
                j=nho-i+1
                Write(11,'(f12.8,3x,f15.10)') qq(j),P1(k,j)
             End do
             
             Close(11)

          End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     This routine is not ready
!
      Subroutine mkq2pot(ifl,NP,imode,jmode,nhoi,nhoj,qqi,qqj,P0,title)

        Implicit None

        Integer :: i,j,k,l,ifl,NP,nhoi,nhoj,imode,jmode,num
        Real(8) :: qqi(nhoi),qqj(nhoj),P0(NP),V2(nhoj,nhoi)
        Real(8) :: Vi(nhoi),Vj(nhoj)
        Character :: title*30

        Character :: im,jm,fl*20

          Do i=1,nhoi
          Do j=1,nhoj
             Read(ifl,*) num,V2(j,i)
          End do
          End do
          !V2=V2-P0

          i=(nhoi+1)/2
          Do j=1,nhoj
             Vj(j)=V2(j,i)
          End do
          j=(nhoj+1)/2
          Do i=1,nhoi
             Vi(i)=V2(j,i)
          End do

          !Write(im,'(i1)') imode
          !Write(jm,'(i1)') jmode
          !fl='q'//im//'q'//jm//'.pot'
          Call get_fname2(imode,jmode,fl)

          Open(11,file=fl,status='new')
          Write(11,'(A)') title
          Write(11,'(''# Number of grids (nj,ni) '')') 
          Write(11,*) nhoj,nhoi
          Write(11,'(''# qj,qi,V (i='',i4,'',j='',i4,'')'')') imode,jmode
          Do i=1,nhoi
             k=nhoi-i+1
          Do j=1,nhoj
             l=nhoj-j+1
             Write(11,'(2f12.8,3x,f15.10)') qqj(l),qqi(k),V2(l,k)-Vi(k)-Vj(l)
          End do
          End do
             
          Close(11)


      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     This routine is not ready
!
      Subroutine mkq3pot(ifl,NP,imode,jmode,kmode,nhoi,nhoj,nhok, &
                         qqi,qqj,qqk,P0,title)

        Implicit None

        Integer :: i,j,k,l,m,n,ifl,NP,nhoi,nhoj,nhok,imode,jmode,kmode,num
        Real(8) :: qqi(nhoi),qqj(nhoj),qqk(nhok),P0(NP),V3(nhok,nhoj,nhoi), &
                   Vi(nhoi),Vj(nhoj),Vk(nhok), &
                   Vij(nhoj,nhoi),Vjk(nhok,nhoj),Vik(nhok,nhoi), tmp
        Character :: title*30

        Character :: im,jm,km,fl*20

          Do i=1,nhoi
          Do j=1,nhoj
          Do k=1,nhok
             Read(ifl,*) num,V3(k,j,i)
          End do
          End do
          End do
          !V3=V3-P0

          i=(nhoi+1)/2
          Do j=1,nhoj
          Do k=1,nhok 
             Vjk(k,j)=V3(k,j,i)
          End do
          End do

          j=(nhoj+1)/2
          Do i=1,nhoi
          Do k=1,nhok 
             Vik(k,i)=V3(k,j,i)
          End do
          End do

          k=(nhok+1)/2
          Do i=1,nhoi
          Do j=1,nhoj 
             Vij(j,i)=V3(k,j,i)
          End do
          End do

          j=(nhoj+1)/2
          k=(nhok+1)/2
          Do i=1,nhoi
             Vi(i)=V3(k,j,i)
          End do

          i=(nhoi+1)/2
          k=(nhok+1)/2
          Do j=1,nhoj
             Vj(j)=V3(k,j,i)
          End do

          i=(nhoi+1)/2
          j=(nhoj+1)/2
          Do k=1,nhok
             Vk(k)=V3(k,j,i)
          End do

          !Write(im,'(i1)') imode
          !Write(jm,'(i1)') jmode
          !Write(km,'(i1)') kmode
          !fl='q'//im//'q'//jm//'q'//km//'.pot'
          Call get_fname3(imode,jmode,kmode,fl)

          Open(11,file=fl,status='new')
          Write(11,'(A)') title
          Write(11,'(''# Number of grids (nk,nj,ni) '')') 
          Write(11,*) nhok,nhoj,nhoi
          Write(11,'(''# qk,qj,qi,V (i='',i4,'',j='',i4,'',k='',i4,'')'')')  &
                      imode,jmode,kmode
          Do i=1,nhoi
             l=nhoi-i+1
          Do j=1,nhoj
             m=nhoj-j+1
          Do k=1,nhok
             n=nhok-k+1
             tmp=V3(n,m,l)-Vij(m,l)-Vik(n,l)-Vjk(n,m)+Vi(l)+Vj(m)+Vk(n)
             Write(11,'(3f12.8,f15.10)') qqk(n),qqj(m),qqi(l),tmp
          End do
          End do
          End do
             
          Close(11)


      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Subroutine DVRGrids(Nfree,nho,max_nho,omega,quad)

        Implicit None

        Integer :: i,j,k,l,max_nho,iout

        Integer :: Nfree,nho(Nfree)
        Real(8) :: Omega(Nfree),quad(max_nho,Nfree)
        Real(8), dimension(:,:), allocatable :: xx
        Real(8), dimension(:,:), allocatable :: xdvr,wdvr
        Real(8), dimension(:), allocatable :: qq,xl,HO
        Real(8) :: const,PI
        Real(8), parameter :: h=6.62606896D-34,    &  ! Planck's constant in J.s
                              c=2.99792458D+10, &  ! velocity of light in cm.s-1
                              Na=6.02214D+23       ! Avogadoro's constant

          !-------------------------------------
          ! >>  Quadrature points / Angs(amu)1/2
          !-------------------------------------

          PI=Acos(-1.D+00)
          Do i=1,Nfree
             l=nho(i)
             const=SQRT(h*Na*1.D+23/c/Omega(i))/(2.D+00*PI)

             Allocate(xx(0:l-1,0:l-1),xl(l*l),qq(l),HO(l))

             xx=0.D+00
             Do j=1,l-1
                k=j-1
                xx(k,j)=Sqrt(dble(j)*0.5D+00)*const
                xx(j,k)=xx(k,j)
             End do
             !Do j=0,l-1
             !   write(6,'(20f8.4)') xx(:,j)
             !End do

             Call huckeler(l,l,xx,qq,xl)
             !write(6,*)
             !write(6,'(20f8.4)') qq
             !write(6,*)
             !Do j=1,l
             !   write(6,'(20f8.4)') xl((j-1)*l+1:j*l)
             !End do

             Do j=1,l
                quad(j,i)=qq(j) 
             End do

             Call HOwfu(0,l,Omega(i),qq,HO)
             !write(iout,'(11f12.6)') (HO(j),j=1,l)

             Deallocate(xx,xl,qq,HO)

          End do


          !-------------------------------------

      End subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

!         n : vibrational quantum number
!     omega : harmonic frequeny (cm-1)
!         q : normal coordinate (Angs amu^1/2)

      Subroutine HOwfu(n,int,omega,qq,HO)

      Implicit None

         Integer, intent(in) :: n,int
         Real(8), intent(in) :: omega
         Real(8), dimension(:)  :: qq(int)
         Real(8), dimension(:,:):: HO(int,0:n)

         Integer :: i,j,k
         Real(8) :: PI,y,r,Nmz,tmp1,hh,hh1,hh2
         Real(8), parameter :: h=6.62606896D-34,    &  ! Planck's constant in J.s
                               c=2.99792458D+10, &  ! velocity of light in cm.s-1
                               Na=6.02214D+23       ! Avogadoro's constant

            PI=Acos(-1.D+00)
            r=2.D+00*PI * SQRT(c*omega/h/Na/1.0D+23)
            Do i=1,int
               y=qq(i)*r
               HO(i,0)=SQRT(r/SQRT(PI))*exp(-0.5D+00*y*y)

               tmp1=1.0D+00
               Do j=1,n
                  tmp1=tmp1*dble(j)*2.D+00
                  Nmz=SQRT(r/SQRT(PI)/tmp1)

                  ! Hermite polynomial
                  Select case(j)
                     case(1)
                       hh=2.D+00*y
                       hh2=1.D+00
                       hh1=hh
                     case default
                       hh=2.D+00*y*hh1-2.D+00*dble(j-1)*hh2
                       hh2=hh1
                       hh1=hh
                  End select

                  HO(i,j)=Nmz*exp(-0.5D+00*y*y)*hh
               End do
            End do

      End Subroutine HOwfu
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine readParams(Nf,MR,NP,mode,nho,omega)

      Implicit None

        Integer :: Nf,MR,NP
        Integer, dimension(*) :: mode,nho
        Real(8), dimension(*) :: omega

        Integer :: i,j,k,l
        Integer :: Nfree
        Real(8), dimension(:), allocatable :: omg

          ! Interactive Input
          Write(6,*) '=== Interactive Setup ==='
          Write(6,*) 'Enter the number of modes:'
          Read(5,*) Nf
          Write(6,*) 'Enter the number of mode couplings:'
          Read(5,*) MR

          Write(6,*) 'Enter the harmonic frequencies:'
          Read(5,*) omega(1:Nf)

          Write(6,*) 'Enter the number of grids:'
          Read(5,*) nho(1:Nf)

          Write(6,*) 'Enter the number of properties:'
          Read(5,*) NP

          Do i=1,Nf
             mode(i)=i
          End do
          
          !dbg Write(6,*) Nf,MR
          !dbg Write(6,'(100i3)') mode(1:Nf)
          !dbg Write(6,'(100i3)')  nho(1:Nf)
          !dbg Write(6,'(3f12.4)') omega(1:Nf)

      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
