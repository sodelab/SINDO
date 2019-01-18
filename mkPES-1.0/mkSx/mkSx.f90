!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Module mkSx_mod

        Integer :: Nat,Nat3,Nfree

        !     nt1,nt2,nt3                    ::  Num. of 1, 2, and 3MR terms 
        !     qq1(nt1),qq2(2,nt2),qq3(3,nt3) ::  Label of mode coupling terms
        Integer :: nt1,nt2,nt3
        Integer, allocatable :: qq1(:),qq2(:,:),qq3(:,:)

        !     nho(Nfree)        :: The number of basis functions for each mode
        !     int               :: Max number of Hermitte-Gauss grid
        !     quad(int,nt1)     :: Quadrature points / Angs(amu)1/2
        !     CHODVR            :: Contracted-HO DVR
        Integer, dimension(:), allocatable :: nho,nhalf
        Integer :: int
        Real(8), dimension(:,:), allocatable :: quad
        Logical :: CHODVR

        !     Potential energy
        Real(8) :: V0
        Integer, dimension(:), allocatable :: idx1,idx2,idx3,idx4

        !     Irst,Irst2 :: File indicators
        Integer :: Irst,Irst2

        !     Title of the PEF
        Character :: title*30

      End Module
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      PROGRAM main

      Implicit None

          Call v2AI_Const()
          Call Sx_Inp()
          Call Sx_nMR_PES()
          Call Sx_finalz
          Call v2AI_Dest()

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Sx_Inp()

      USE mkSx_mod

      Implicit None

        Integer, parameter :: nmax=1000

        Integer :: i,j,k,l
        Integer :: Get_Nat,Get_Nfree
        Integer :: vmax(nmax)
        Integer, dimension(nmax) :: mr1,mr2,mr3
        Integer, dimension(:), allocatable :: mode
        Real(8), dimension(:), allocatable :: qd
        Logical :: op

        Namelist /mkSx/vmax,mr1,mr2,mr3,title,CHODVR

!----------------------------------------------------------------------

        Write(6,100)
  100   Format(/,'--------------------(    MKSX MODULE    )--------------------',/)

! --    Read parameters
        Call spr_SetMaxmem(100)
        Nat=Get_Nat()
        Nfree=Get_Nfree()
        vmax=-1
        mr1=0
        mr2=0
        mr3=0
        title=''

        Rewind(5)
        Read(5,mkSx)

        write(6,110) 
    110 Format(' -> Options')

        if(title=='') then
           write(6,*) 'title is empty!'
           Stop
        endif

! --    3MR Terms
        nt3=0
        Do i=1,nmax,3
           if(mr3(i)==0) exit
           nt3=nt3+1
        End do
        Allocate(qq3(3,nt3))
        Do i=1,nt3
           qq3(1,i)=mr3(i*3-2)
           qq3(2,i)=mr3(i*3-1)
           qq3(3,i)=mr3(i*3)
           if(qq3(2,i)<qq3(3,i)) then
              j=qq3(2,i)
              qq3(2,i)=qq3(3,i)
              qq3(3,i)=j
           endif
           if(qq3(1,i)<qq3(3,i)) then
              j=qq3(1,i)
              qq3(1,i)=qq3(2,i)
              qq3(2,i)=qq3(3,i)
              qq3(3,i)=j
           elseif(qq3(1,i)<qq3(2,i)) then
              j=qq3(1,i)
              qq3(1,i)=qq3(2,i)
              qq3(2,i)=j
           endif
        End do

! --    2MR Terms
        nt2=0
        Do i=1,nmax,2
           if(mr2(i)==0) exit
           if(mr2(i)<mr2(i+1)) then
              j=mr2(i)
              mr2(i)=mr2(i+1)
              mr2(i+1)=j
           endif
           nt2=nt2+1
        End do
        Do i=1,nt3
           k=1
           Do j=1,nt2
              if(mr2(k)==qq3(1,i) .and. mr2(k+1)==qq3(2,i)) goto 10
              k=k+2
           End do
           mr2(nt2*2+1)=qq3(1,i)
           mr2(nt2*2+2)=qq3(2,i)
           nt2=nt2+1
        10 Continue

           k=1
           Do j=1,nt2
              if(mr2(k)==qq3(1,i) .and. mr2(k+1)==qq3(3,i)) goto 11
              k=k+2
           End do
           mr2(nt2*2+1)=qq3(1,i)
           mr2(nt2*2+2)=qq3(3,i)
           nt2=nt2+1
        11 Continue

           k=1
           Do j=1,nt2
              if(mr2(k)==qq3(2,i) .and. mr2(k+1)==qq3(3,i)) goto 12
              k=k+2
           End do
           mr2(nt2*2+1)=qq3(2,i)
           mr2(nt2*2+2)=qq3(3,i)
           nt2=nt2+1
        12 Continue

        End do

        Allocate(qq2(2,nt2))
        Do i=1,nt2
           qq2(1,i)=mr2(i*2-1)
           qq2(2,i)=mr2(i*2)
        End do

! --    1MR Terms
        nt1=0
        Do i=1,nmax
           if(mr1(i)==0) exit
           nt1=nt1+1
        End do
        Do i=1,nt2
           Do j=1,nt1
              if(mr1(j)==qq2(1,i)) goto 20
           End do
           mr1(nt1+1)=qq2(1,i)
           nt1=nt1+1
        20 Continue

           Do j=1,nt1
              if(mr1(j)==qq2(2,i)) goto 21
           End do
           mr1(nt1+1)=qq2(2,i)
           nt1=nt1+1
        21 Continue

        End do

        Allocate(qq1(nt1))
        qq1=mr1(:nt1)

!--     Print
!--     1MR
        if(nt1>0) then
           write(6,120) nt1
           write(6,121) qq1
    120    Format(2x,'> 1MR Terms :',i4)
    121    Format(5x,10(i2.2,2x))
        endif

!--     2MR
        if(nt2>0) then
           write(6,130) nt2
    130    Format(2x,'> 2MR Terms :',i4)
!f90           write(6,130) qq2
!f90    131    Format(5x,5('(',i2.2,',',i2.2,')',1x))
           j=1
           Do while(.true.)
              write(6,131) 
              Do i=1,5
                 write(6,132) qq2(:,j)
                 j=j+1
                 if(j>nt2) exit
              End do
              write(6,'(1x)')
              if(j>nt2) exit
           End do
    131    Format('     ',$)
    132    Format('(',i2.2,',',i2.2,') ',$)
        endif

!--     3MR
        if(nt3>0) then
           write(6,140) nt3
    140    Format(2x,'> 3MR Terms :',i4)
!f90           write(6,141) qq3
!f90       141 Format(5x,3('(',i2.2,',',i2.2,',',i2.2,')',1x))
           j=1
           Do while(.true.)
              write(6,131) 
              Do i=1,5
                 write(6,142) qq3(:,j)
                 j=j+1
                 if(j>nt3) exit
              End do
              write(6,'(1x)')
              if(j>nt3) exit
           End do
    142    Format('(',i2.2,',',i2.2,',',i2.2,') ',$)
        endif

!--     Restart Files
        Irst=210
        Do while(.true.)
           Inquire(Irst,opened=op)
           if(op) then
             Irst=Irst+1
           else
             exit
           endif
        End do
        Open(Irst,file='restart-r',status='unknown')

        Irst2=220
        Do while(.true.)
           Inquire(Irst2,opened=op)
           if(op) then
             Irst2=Irst2+1
           else
             exit
           endif
        End do
        Open(Irst2,file='restart-w',status='unknown')

        write(6,150)
    150 Format(/,' -> Restart(read)  FILE :  [ restart-r ]',/, &
                 ' -> Restart(write) FILE :  [ restart-w ]',/)

!--     Quadrature Grid Points
        Allocate(nho(Nfree),nhalf(Nfree),mode(Nfree))
        int=0
        Do i=1,Nfree
           if(mod(vmax(i),2)/=0) vmax(i)=vmax(i)+1
           nho(i)=vmax(i)+1
           if(int<nho(i)) int=nho(i)
           nhalf(i)=(nho(i)+1)/2
        End do
        Do i=1,Nfree
           mode(i)=i
        End do
        Allocate(quad(int,Nfree),qd(int))
        if(CHODVR) then
           Call CHODVRGrid(Nfree,mode,int,nho,quad)
        else
           Call HODVRGrid(Nfree,mode,int,nho,quad)
        endif
        Do i=1,Nfree
           qd=quad(1:nho(i),i)
           Do j=1,nho(i)
              quad(j,i)=qd(nho(i)+1-j)
           End do
        End do
        write(6,200)
        Do i=1,Nfree
           write(6,210) mode(i)
           write(6,220) quad(1:nho(i),i)
        End do
        Deallocate(mode,qd)

    200 Format(2x,'> DVR Grids')
    210 Format(6x,'Mode :',i3)
    220 Format(6x,6f9.4)

        End
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Subroutine Sx_finalz

      USE mkSx_mod

      Implicit None

      Integer :: in,io

          Call spr_Getio(in,io)
          Write(io,100)
      100 Format(/,'(  FINALIZE MKSX MODULE  )',/)

          Deallocate(quad,nho)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Sx_nMR_PES()

      USE mkSx_mod

      Implicit None

      ! ------------------------------------------------------------------------

      Integer :: io,in
      !Integer :: GetKey
      Logical :: dpl,get_dpl

      ! ------------------------------------------------------------------------

          dpl=Get_dpl()
          Call spr_GetIO(in,io)

          !----------------------
          ! >>  Q=0
          !----------------------
          write(io,100) 
      100 Format(' -> Energies on grids ',/, &
                2x,'> Q=0')
          Select case(GetKey())
            case(0)
              Call Sx_e_0(io)
            case(1)
              !Call Sx_ed_0(io)
          End select

          !----------------------
          ! >>  1MR PES 
          !----------------------
          Select case(GetKey())
            case(0)
              Call Sx_e_1MR(io)
            case(1)
              !Call Sx_ed_1MR(io)
          End select

          if(nt2==0) return

          !----------------------
          ! >>  2MR PES 
          !----------------------
          Select case(GetKey())
            case(0)
              Call Sx_e_2MR(io)
            case(1)
              !Call Sx_ed_2MR(io)
          End select

          if(nt3==0) return

          !----------------------
          ! >>  3MR PES 
          !----------------------
          Select case(GetKey())
            case(0)
              Call Sx_e_3MR(io)
            case(1)
              !Call Sx_ed_3MR(io)
          End select

          ! --------------------------------------------------------------------

      Contains

      Function GetKey()

      Implicit None

      Integer :: GetKey, key

!           if(.not. dpl .and. .not. nmrcc) then
           if(.not. dpl) then
             ! Energy
             key=0
           elseif(dpl) then
             ! Energy and dipole
             key=1
!           elseif(nmrcc) then
!             ! Energy and NMR-CC
!             key=2
           endif
           GetKey=key

         End function

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine Sx_e_0(io)

      USE mkSx_mod

      Implicit None

      Integer :: io
      Real(8), dimension(:), allocatable   :: qq

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree))
          qq=0.D+00
          Call Sx_airun_e(0,qq,V0)
          Write(io,2000) 0,V0
          Deallocate(qq)

     2000 Format(i8,f18.10)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Sx_e_1MR(io)

      USE mkSx_mod

      Implicit None

      Integer :: io

      Integer :: i,im,ist,nh,nh0,iq1,idx,Mi,si(3)
      Real(8), dimension(:), allocatable :: qq,qu,V1

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree))
          qq=0.D+00

          Do i=1,nt1

             im=qq1(i)
             Write(io,1000) im
             Call chkq1pot(im,ist)
             if(ist<0) cycle

             idx=0
             nh=nhalf(im)
             nh0=nho(im)
             Allocate(V1(nh0),qu(nh0))
             qu=quad(1:nh0,im)

             Call SymOp(im,Mi,si)

             Call PEF1D(Nfree,qq,CHODVR,idx,im,Mi,nh0,qu,V0,V1)
             qq(im)=0.D+00

             Call intpl1(nh0,1,qu,V1)
             V1=V1-V0
             Call mkq1pot(im,nh0,qu,V1,title)
             Call reset_restartw

             Deallocate(V1,qu)

          End do

          Deallocate(qq)

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,f8.3)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine PEF1D(Nfree,qq,CHODVR,idx,ii,Mi,ni,qqi,V0,V1)

      Implicit None

      Integer :: Nfree,ii,idx,Mi,ni
      Real(8) :: qq(Nfree),qqi(ni),V1(ni),V0
      Logical :: CHODVR

      Integer :: iq1,nhi
      Real(8) :: qiorg

         nhi=(ni+1)/2
         qiorg=qq(ii)

         Do iq1=1,nhi-1
            qq(ii)=qqi(iq1)
            idx=idx+1
            Call Sx_airun_e(idx,qq,V1(iq1))
            Write(6,2000) idx,V1(iq1),qq(ii)
         End do

         if(Mi/=0) then
            V1(nhi)=V0
            Write(6,2000) 0,V1(nhi),0.D+00

            Do iq1=nhi+1,ni
               V1(iq1)=V1(ni-iq1+1)
               Write(6,2000) 0,V1(iq1),qqi(iq1)
            End do

         else
            if(.not. CHODVR) then
               V1(nhi)=V0
               Write(6,2000) 0,V1(nhi),0.D+00

            else
               qq(ii)=qqi(nhi)
               idx=idx+1
               Call Sx_airun_e(idx,qq,V1(nhi))
               Write(6,2000) idx,V1(nhi),qq(ii)

            endif

            Do iq1=nhi+1,ni
               qq(ii)=qqi(iq1)
               idx=idx+1
               Call Sx_airun_e(idx,qq,V1(iq1))
               Write(6,2000) idx,V1(iq1),qq(ii)
            End do

         endif

         qq(ii)=qiorg

    2000 Format(i8,f18.10,100f8.3)

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Sx_e_2MR(io)

      USE mkSx_mod

      Implicit None

      Integer :: io

      Integer :: i,j,k,im,jm,iq1,iq2,idx,Mi,Mj,si(3),sj(3)
      Integer :: nhi,nhj,nhoi,nhoj,ist
      Real(8), dimension(:), allocatable :: qui,quj,V1i,V1j,qq
      Real(8), dimension(:,:), allocatable :: V2

      ! ------------------------------------------------------------------------

          Do k=1,nt2

             im=qq2(1,k)
             jm=qq2(2,k)
             Write(io,1000) im,jm
             Call chkq2pot(im,jm,ist)
             if(ist<0) cycle

             idx=0
             nhi=nhalf(im)
             nhj=nhalf(jm)
             nhoi=nho(im)
             nhoj=nho(jm)
             Call SymOp(im,Mi,si)
             Call SymOp(jm,Mj,sj)
             Allocate(qui(nhoi),quj(nhoj),qq(Nfree),V2(nhoj,nhoi),V1i(nhoi),V1j(nhoj))
             qui=quad(1:nhoi,im)
             quj=quad(1:nhoj,jm)
             qq=0.D+00

             Call getV1(im,nhoi,qui,V1i)
             Call getV1(jm,nhoj,quj,V1j)

             Call PEF2D(Nfree,qq,CHODVR,idx,im,jm,Mi,Mj,nhoi,nhoj,qui,quj, &
                        V0,V1i,V1j,V2)
             Call intpl2(nhoi,nhoj,1,quj,V2)
             V2=V2-V0
             Call mkq2pot(im,jm,nhoi,nhoj,qui,quj,V2,title,V1i,V1j)
             Call reset_restartw

             Deallocate(qui,quj,qq,V2,V1i,V1j)

          End do

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,2f8.3)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine PEF2D(Nfree,qq,CHODVR,idx,ii,jj,Mi,Mj,ni,nj,qqi,qqj,&
                       V0,V1i,V1j,V2)

      Implicit None

      Integer :: Nfree,idx,ii,jj,Mi,Mj,ni,nj
      Real(8) :: qq(Nfree),qqi(ni),qqj(nj),V0,V1i(ni),V1j(nj),V2(nj,ni)
      Logical :: CHODVR

      Integer :: iq1,iq2,nhi,nhj,type
      Real(8) :: qiorg,qjorg

         nhi=(ni+1)/2
         nhj=(nj+1)/2
         qiorg=qq(ii)
         qjorg=qq(jj)
 
         type=-1
         if(Mi/=0) then
            if(Mj==Mi) then
               type=1
            elseif(Mj/=0) then
               type=4
            else
               type=3
            endif
         else
            type=2
         endif
         !dbg write(6,'(2i3,'', type='',i3)') ii,jj,type
 
         Select Case(type)

            Case(1)
               Do iq2=1,nj
                  V2(iq2,nhi)=V0+V1j(iq2)
                  Write(6,2000) 0,V2(iq2,nhi),qqi(nhi),qqj(iq2)
               End do
               Do iq1=1,ni
                  V2(nhj,iq1)=V0+V1i(iq1)
                  Write(6,2000) 0,V2(nhj,iq1),qqi(iq1),qqj(nhj)
               End do
 
               ! qi<0 region
               Do iq1=1,nhi-1
                  qq(ii)=qqi(iq1)
 
                  Do iq2=1,nj
                     if(iq2==nhj) cycle
 
                     qq(jj)=qqj(iq2)
                     idx=idx+1
                     Call Sx_airun_e(idx,qq,V2(iq2,iq1))
                     Write(6,2000) idx,V2(iq2,iq1),qq(ii),qq(jj)
 
                  End do
               End do
 
               ! qi>0 region : generate by symmetry
               Do iq1=nhi+1,ni
                  Do iq2=1,nj
                     if(iq2==nhj) cycle
 
                     V2(iq2,iq1)=V2(nj-iq2+1,ni-iq1+1)
                     Write(6,2000) 0,V2(iq2,iq1),qqi(iq1),qqj(iq2)
                  End do
               End do
 
            Case(2)
               Do iq1=1,ni
                  if(iq1==nhi .and. abs(qqi(nhi)) < 1.D-06) then
                     Do iq2=1,nj
                        V2(iq2,nhi)=V0+V1j(iq2)
                        Write(6,2000) 0,V2(iq2,nhi),qqi(nhi),qqj(iq2)
                     End do
                     cycle
                  endif

                  qq(ii)=qqi(iq1)
                  Write(6,2001) ii,qq(ii)
                  Call PEF1D(Nfree,qq,CHODVR,idx,jj,Mj,nj,qqj, &
                             V0+V1i(iq1),V2(:,iq1))
               End do
 
            Case(3)
               Do iq2=1,nj
                  if(iq2==nhj .and. abs(qqj(nhj)) < 1.D-06) then
                     Do iq1=1,ni
                        V2(nhj,iq1)=V0+V1i(iq1)
                        Write(6,2000) 0,V2(nhj,iq1),qqi(iq1),qqj(nhj)
                     End do
                     cycle
                  endif
                  qq(jj)=qqj(iq2)
                  Write(6,2001) jj,qq(jj)
                  Call PEF1D(Nfree,qq,CHODVR,idx,ii,Mi,ni,qqi, &
                             V0+V1j(iq2),V2(iq2,:))
               End do

            Case(4)
               Do iq1=1,nhi-1
                  qq(ii)=qqi(iq1)
                  Write(6,2001) ii,qq(ii)
                  Call PEF1D(Nfree,qq,CHODVR,idx,jj,Mj,nj,qqj, &
                             V0+V1i(iq1),V2(:,iq1))
               End do
               Do iq2=1,nj
                  V2(iq2,nhi)=V0+V1j(iq2)
                  Write(6,2000) 0,V2(iq2,nhi),qqi(nhi),qqj(iq2)
               End do

               ! qi>0 region : generate by symmetry
               Do iq1=nhi+1,ni
                  Do iq2=1,nj
                     V2(iq2,iq1)=V2(iq2,ni-iq1+1)
                     Write(6,2000) 0,V2(iq2,iq1),qqi(iq1),qqj(iq2)
                  End do
               End do
 
         End select

         qq(ii)=qiorg
         qq(jj)=qjorg
 
    2000 Format(i8,f18.10,100f8.3)
    2001 Format(4x,'> MODE',i3,'= ',f8.3)
 
      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Sx_e_3MR(io)

      USE mkSx_mod

      Implicit None

      Integer :: io

      Integer :: i,j,k,ij,jk,ik,ijk,iq1,iq2,iq3,idx,jdx, &
                 Mi,Mj,Mk,si(3),sj(3),sk(3)
      Integer :: im,jm,km,nhi,nhj,nhk,nhoi,nhoj,nhok,ist
      Real(8), dimension(:), allocatable :: qq,qui,quj,quk,V1i,V1j,V1k
      Real(8), dimension(:,:), allocatable :: V2ij,V2ik,V2jk
      Real(8), dimension(:,:,:), allocatable :: V3

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree))

          ! --------------------------------------------------------------------
          ! >>>>>>>
          Do i=1,nt3

             im=qq3(1,i)
             jm=qq3(2,i)
             km=qq3(3,i)
             Write(io,1000) im,jm,km
             Call chkq3pot(im,jm,km,ist)
             if(ist<0) cycle

             idx=0
             nhi=nhalf(im)
             nhj=nhalf(jm)
             nhk=nhalf(km)
             nhoi=nho(im)
             nhoj=nho(jm)
             nhok=nho(km)
             Call SymOp(im,Mi,si)
             Call SymOp(jm,Mj,sj)
             Call SymOp(km,Mk,sk)
             Allocate(qui(nhoi),quj(nhoj),quk(nhok))
             Allocate(V1i(nhoi),V1j(nhoj),V1k(nhok))
             Allocate(V2ij(nhoj,nhoi),V2ik(nhok,nhoi),V2jk(nhok,nhoj))
             Allocate(V3(nhok,nhoj,nhoi))
             qui=quad(1:nhoi,im)
             quj=quad(1:nhoj,jm)
             quk=quad(1:nhok,km)

             Call getV1(im,nhoi,qui,V1i)
             Call getV1(jm,nhoj,quj,V1j)
             Call getV1(km,nhok,quk,V1k)
             Call getV2(im,jm,nhoi,nhoj,qui,quj,V2ij)
             Call getV2(im,km,nhoi,nhok,qui,quk,V2ik)
             Call getV2(jm,km,nhoj,nhok,quj,quk,V2jk)

             qq=0.D+00
             Call PEF3D(Nfree,qq,CHODVR,idx,im,jm,km,Mi,Mj,Mk,nhoi,nhoj,nhok,&
                        qui,quj,quk,V0,V1i,V1j,V1k,V2ij,V2ik,V2jk,V3)
             Call intpl3(nhoi,nhoj,nhok,1,quk,V3)
             V3=V3-V0
             Call mkq3pot(im,jm,km,nhoi,nhoj,nhok,qui,quj,quk,V3,title, &
                          V1i,V1j,V1k,V2ij,V2ik,V2jk)
             Call reset_restartw

             Deallocate(qui,quj,quk,V1i,V1j,V1k,V2ij,V2ik,V2jk,V3)

          End do

          Deallocate(qq)

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,3f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine PEF3D(Nfree,qq,CHODVR,idx,ii,jj,kk,Mi,Mj,Mk,ni,nj,nk, &
                       qqi,qqj,qqk,V0,V1i,V1j,V1k,V2ij,V2ik,V2jk,V3)

      Implicit None

      Integer :: Nfree,idx,ii,jj,kk,Mi,Mj,Mk,ni,nj,nk
      Real(8) :: qq(Nfree),qqi(ni),qqj(nj),qqk(nk),V0,V1i(ni),V1j(nj),V1k(nk), &
                 V2ij(nj,ni),V2ik(nk,ni),V2jk(nk,nj),V3(nk,nj,ni)
      Logical :: CHODVR

      Integer :: iq1,iq2,iq3,nhi,nhj,nhk,type,getChr
      Real(8) :: qiorg,qjorg,qkorg
      Character :: Sym*3,getSym*3 

         nhi=(ni+1)/2
         nhj=(nj+1)/2
         nhk=(nk+1)/2
         qiorg=qq(ii)
         qjorg=qq(jj)
         qkorg=qq(kk)

         Sym=getSym()

         type=-1
         if(Mi==0) then
            type=2
         elseif(Mj==0) then
            type=3
         elseif(Mk==0) then
            type=4
         else
            if(Mi==Mj .and. Mj==Mk) then
               type=1
            endif
            if(Mj==Mk .and. Mi/=Mj) then
               type=5
            endif
            if(Mi==Mk .and. Mi/=Mj) then
               type=6
            endif
            if(Mi==Mj .and. Mj/=Mk) then
               type=7
            endif
            if(Mi/=Mj .and. Mj/=Mk .and. Mi/=Mk) then
               if(Sym=='D2H' .and. getChr(Mi,Mj,Mk)==1) then
                  type=5
               else
                  type=8
               endif
            endif
         endif

         !dbg write(6,'(3i3,'', type='',i3)') ii,jj,kk,type

         Select Case(type)
            Case(1)

               Call copy_qi0()
               Call copy_qj0()
               Call copy_qk0()

               Do iq1=1,ni
                  if(iq1==nhi) cycle
                  qq(ii)=qqi(iq1)
               Do iq2=1,nj
                  if(iq2==nhj) cycle
                  qq(jj)=qqj(iq2)
               Do iq3=1,nhk-1
                  qq(kk)=qqk(iq3)

                  idx=idx+1
                  Call Sx_airun_e(idx,qq,V3(iq3,iq2,iq1))
                  Write(6,2000) idx,V3(iq3,iq2,iq1),qq(ii),qq(jj),qq(kk)

               End do
               End do
               End do

               ! qk>0 region
               Do iq1=1,ni
                  if(iq1==nhi) cycle
               Do iq2=1,nj
                  if(iq2==nhj) cycle
               Do iq3=nhk+1,nk
                  V3(iq3,iq2,iq1)=V3(nk-iq3+1,nj-iq2+1,ni-iq1+1)
                  Write(6,2000) 0,V3(iq3,iq2,iq1),qqi(iq1),qqj(iq2),qqk(iq3)

               End do
               End do
               End do

            Case(2)
               Do iq1=1,ni
                  if(iq1==nhi .and. abs(qqi(nhi)) < 1.D-06) then
                     Call copy_qi0()
                     cycle
                  endif

                  qq(ii)=qqi(iq1)
                  Write(6,2001) ii,qq(ii)
                  Call PEF2D(Nfree,qq,CHODVR,idx,jj,kk,Mj,Mk,nj,nk,qqj,qqk, &
                             V0+V1i(iq1),V1j+V2ij(:,iq1),V1k+V2ik(:,iq1),V3(:,:,iq1))
               End do

            Case(3)
               Do iq2=1,nj
                  if(iq2==nhj .and. abs(qqj(nhj)) < 1.D-06) then
                     Call copy_qj0()
                     cycle
                  endif

                  qq(jj)=qqj(iq2)
                  Write(6,2001) jj,qq(jj)
                  Call PEF2D(Nfree,qq,CHODVR,idx,ii,kk,Mi,Mk,ni,nk,qqi,qqk, &
                             V0+V1j(iq2),V1i+V2ij(iq2,:),V1k+V2jk(:,iq2),V3(:,iq2,:))
               End do

            Case(4)
               Do iq3=1,nk
                  if(iq3==nhk .and. abs(qqk(nhk)) < 1.D-06) then
                     Call copy_qk0()
                     cycle
                  endif

                  qq(kk)=qqk(iq3)
                  Write(6,2001) kk,qq(kk)
                  Call PEF2D(Nfree,qq,CHODVR,idx,ii,jj,Mi,Mj,ni,nj,qqi,qqj, &
                             V0+V1k(iq3),V1i+V2ik(iq3,:),V1j+V2jk(iq3,:),V3(iq3,:,:))
               End do

            Case(5)
               Do iq1=1,nhi-1
                  qq(ii)=qqi(iq1)
                  Write(6,2001) ii,qq(ii)
                  Call PEF2D(Nfree,qq,CHODVR,idx,jj,kk,Mj,Mk,nj,nk,qqj,qqk, &
                             V0+V1i(iq1),V1j+V2ij(:,iq1),V1k+V2ik(:,iq1),V3(:,:,iq1))
               End do
               Call copy_qi0()
               Do iq1=nhi+1,ni
               Do iq2=1,nj
               Do iq3=1,nk
                  V3(iq3,iq2,iq1)=V3(iq3,iq2,ni+1-iq1)
                  Write(6,2000) 0,V3(iq3,iq2,iq1),qqi(iq1),qqj(iq2),qqk(iq3)
               End do
               End do
               End do

            Case(6)
               Do iq2=1,nhj-1
                  qq(jj)=qqj(iq2)
                  Write(6,2001) jj,qq(jj)
                  Call PEF2D(Nfree,qq,CHODVR,idx,ii,kk,Mi,Mk,ni,nk,qqi,qqk, &
                             V0+V1j(iq2),V1i+V2ij(iq2,:),V1k+V2jk(:,iq2),V3(:,iq2,:))
               End do
               Call copy_qj0()
               Do iq1=1,ni
               Do iq2=nhj+1,nj
               DO iq3=1,nk
                  V3(iq3,iq2,iq1)=V3(iq3,nj+1-iq2,iq1)
                  Write(6,2000) 0,V3(iq3,iq2,iq1),qqi(iq1),qqj(iq2),qqk(iq3)
               End do
               End do
               End do

            Case(7)
               Do iq3=1,nhk-1
                  qq(kk)=qqk(iq3)
                  Write(6,2001) kk,qq(kk)
                  Call PEF2D(Nfree,qq,CHODVR,idx,ii,jj,Mi,Mj,ni,nj,qqi,qqj, &
                             V0+V1k(iq3),V1i+V2ik(iq3,:),V1j+V2jk(iq3,:),V3(iq3,:,:))
               End do
               Call copy_qk0()
               Do iq1=1,ni
               Do iq2=1,nj
               Do iq3=nhk+1,nk
                  V3(iq3,iq2,iq1)=V3(nk+1-iq3,iq2,iq1)
                  Write(6,2000) 0,V3(iq3,iq2,iq1),qqi(iq1),qqj(iq2),qqk(iq3)
               End do
               End do
               End do

            Case(8)
               Call copy_qi0()
               Call copy_qj0()
               Call copy_qk0()

               Do iq1=1,nhi-1
                  qq(ii)=qqi(iq1)
               Do iq2=1,nhj-1
                  qq(jj)=qqj(iq2)
               Do iq3=1,nk
                  if(iq3==nhk) cycle
                  qq(kk)=qqk(iq3)

                  idx=idx+1
                  Call Sx_airun_e(idx,qq,V3(iq3,iq2,iq1))
                  Write(6,2000) idx,V3(iq3,iq2,iq1),qq(ii),qq(jj),qq(kk)

               End do
               End do
               End do

               ! qi<0,qj>0 region
               Do iq1=1,nhi-1
               Do iq2=nhj+1,nj
               Do iq3=1,nk
                  if(iq3==nhk) cycle
                  V3(iq3,iq2,iq1)=V3(nk-iq3+1,nj-iq2+1,iq1)
                  Write(6,2000) 0,V3(iq3,iq2,iq1),qqi(iq1),qqj(iq2),qqk(iq3)

               End do
               End do
               End do

               ! qi>0 region
               Do iq1=nhi+1,ni
               Do iq2=1,nj
                  if(iq2==nhj) cycle
               Do iq3=1,nk
                  if(iq3==nhk) cycle
                  V3(iq3,iq2,iq1)=V3(nk-iq3+1,iq2,ni-iq1+1)
                  Write(6,2000) 0,V3(iq3,iq2,iq1),qqi(iq1),qqj(iq2),qqk(iq3)

               End do
               End do
               End do

         End select 

         qq(ii)=qiorg
         qq(jj)=qjorg
         qq(kk)=qkorg

    2000 Format(i8,f18.10,3f8.3)
    2001 Format(4x,'> MODE',i3,'= ',f8.3)

         Contains 

         Subroutine copy_qi0()

            ! qi=0
            Do iq2=1,nj
            Do iq3=1,nk
               V3(iq3,iq2,nhi)=V0+V1j(iq2)+V1k(iq3)+V2jk(iq3,iq2)
               Write(6,2000) 0,V3(iq3,iq2,nhi),qqi(nhi),qqj(iq2),qqk(iq3)
            End do
            End do
       2000 Format(i8,f18.10,3f8.3)

         End subroutine

         Subroutine copy_qj0()

            ! qj=0
            Do iq1=1,ni
            Do iq3=1,nk
               V3(iq3,nhj,iq1)=V0+V1i(iq1)+V1k(iq3)+V2ik(iq3,iq1)
               Write(6,2000) 0,V3(iq3,nhj,iq1),qqi(iq1),0.D+00,qqk(iq3)
            End do
            End do
       2000 Format(i8,f18.10,3f8.3)

         End subroutine

         Subroutine copy_qk0()

            ! qk=0
            Do iq1=1,ni
            Do iq2=1,nj
               V3(nhk,iq2,iq1)=V0+V1i(iq1)+V1j(iq2)+V2ij(iq2,iq1)
               Write(6,2000) 0,V3(nhk,iq2,iq1),qqi(iq1),qqj(iq2),0.D+00
            End do
            End do
       2000 Format(i8,f18.10,3f8.3)

         End Subroutine

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine Sx_airun_e(indx0,qq,E)

      USE mkSx_mod

      Implicit None

      Integer :: indx0,indx1
      Real(8) :: qq(Nfree),E
      Integer :: i

         Read(Irst,*,end=100,err=100) indx1,E
         if(indx0 .ne. indx1) then 
            Write(6,*) 'Invalid line is detected in the restart file:'
            !Write(6,'(i8,f20.10)') indx1,E
            Call write_e(6,indx1,E)
            Write(6,'('' It must have an index '',i8)') indx0
            Stop
         endif
         Call write_e(Irst2,indx0,E)
         Call myFlsh(Irst2)
         return

  100    Continue
         Call Run_e(qq,E)
         Call write_e(Irst2,indx0,E)
         Call myFlsh(Irst2)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine write_e(ifl,indx,E)

      Implicit None

          Integer :: ifl,indx
          Real(8) :: E

          Write(ifl,500) indx,E
      500 Format(i8,f20.10)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine chkq1pot(mode,istat)

        Implicit None

        Integer :: i,j,mode,istat
        Character :: fl*30,com*40

          istat=0
          Call get_fname1(mode,fl)
          Open(11,file=fl,status='new',err=10)
          Close(11)
          com='rm '//fl
          Call System(com)

          return

       10 continue
          Close(11)
          istat=-1

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mkq1pot(mode,nho,qq,V1,title)

        Implicit None

        Integer :: i,j,nho,mode
        Real(8) :: qq(nho),V1(nho)
        Character :: title*30,fl*30

          Call get_fname1(mode,fl)
          Open(11,file=fl,status='new')
          Write(11,'(A)') title
          Write(11,'(''# Number of grids'')')
          Write(11,*) nho
          Write(11,'(''# qi,V (i='',i4,'')'')') mode
          Do i=1,nho
             Write(11,'(f12.8,f15.10)') qq(i),V1(i)
          End do

          Close(11)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine chkq2pot(imode,jmode,istat)

        Implicit None

        Integer :: i,j,imode,jmode,istat
        Character :: fl*20,com*40

          istat=0
          Call get_fname2(imode,jmode,fl)
          Open(11,file=fl,status='new',err=10)
          Close(11)
          com='rm '//fl
          Call System(com)

          return

       10 continue
          Close(11)
          istat=-1

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkq2pot(imode,jmode,nhoi,nhoj,qqi,qqj,V2,title,Vi,Vj)

        Implicit None

        Integer :: i,j,k,l,nhoi,nhoj,imode,jmode,num
        Real(8) :: qqi(nhoi),qqj(nhoj),V2(nhoj,nhoi)
        Real(8) :: Vi(nhoi),Vj(nhoj)
        Character :: title*30

        Character :: im,jm,fl*20

          Call get_fname2(imode,jmode,fl)

          Open(11,file=fl,status='new')
          Write(11,'(A)') title
          Write(11,'(''# Number of grids (nj,ni) '')') 
          Write(11,*) nhoj,nhoi
          Write(11,'(''# qj,qi,V (i='',i4,'',j='',i4,'')'')') imode,jmode
          Do i=1,nhoi
          Do j=1,nhoj
             Write(11,'(2f12.8,f15.10)') qqj(j),qqi(i),V2(j,i)-Vi(i)-Vj(j)
          End do
          End do
             
          Close(11)


      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine chkq3pot(imode,jmode,kmode,istat)

        Implicit None

        Integer :: i,j,imode,jmode,kmode,istat
        Character :: fl*20,com*40

          istat=0
          Call get_fname3(imode,jmode,kmode,fl)
          Open(11,file=fl,status='new',err=10)
          Close(11)
          com='rm '//fl
          Call System(com)

          return

       10 continue
          Close(11)
          istat=-1

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkq3pot(imode,jmode,kmode,nhoi,nhoj,nhok,qqi,qqj,qqk,V3, &
                         title,V1i,V1j,V1k,V2ij,V2ik,V2jk)

        Implicit None

        Integer :: i,j,k,l,m,n,ifl,nhoi,nhoj,nhok,imode,jmode,kmode,num
        Real(8) :: qqi(nhoi),qqj(nhoj),qqk(nhok),E0,V3(nhok,nhoj,nhoi), &
                   V1i(nhoi),V1j(nhoj),V1k(nhok), &
                   V2ij(nhoj,nhoi),V2jk(nhok,nhoj),V2ik(nhok,nhoi), tmp
        Character :: title*30

        Character :: im,jm,km,fl*20

          Call get_fname3(imode,jmode,kmode,fl)

          Open(11,file=fl,status='new')
          Write(11,'(A)') title
          Write(11,'(''# Number of grids (nk,nj,ni) '')') 
          Write(11,*) nhok,nhoj,nhoi
          Write(11,'(''# qk,qj,qi,V (i='',i4,'',j='',i4,'',k='',i4,'')'')')  &
                      imode,jmode,kmode
          Do i=1,nhoi
          Do j=1,nhoj
          Do k=1,nhok
             tmp=V3(k,j,i)-V2ij(j,i)-V2ik(k,i)-V2jk(k,j)-V1i(i)-V1j(j)-V1k(k)
             Write(11,'(3f12.8,f15.10)') qqk(k),qqj(j),qqi(i),tmp
          End do
          End do
          End do
             
          Close(11)


      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine getV1(mode,nho,q1,V1)

        Implicit None

        Integer :: i,j,nho,mode
        Real(8) :: q1(nho),V1(nho),spfn1
        Real(8), allocatable :: qq(:),VV(:)
        Character :: fl*30

!          Call get_fname1(mode,fl)
!          Open(11,file=fl,status='old')
!          Read(11,*)
!          Read(11,*)
!          Read(11,*) i
!          Read(11,*)
!          Allocate(qq(i),VV(i))
!          Do i=1,nho
!             Read(11,*) qq(i),VV(i)
!          End do
!          Close(11)
          Call Spfn1_gen2(1,mode)
          Do i=1,nho
             V1(i)=Spfn1(q1(i))
          End do
          Call Spfn1_Dispose()

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine getV2(imode,jmode,nhoi,nhoj,qqi,qqj,V2)

        Implicit None

        Integer :: i,j,k,imode,jmode,nhoi,nhoj
        Real(8) :: qqj(nhoj),qqi(nhoi),V2(nhoj*nhoi),Spfn2
        Character :: fl*30

!          Call get_fname2(imode,jmode,fl)
!          Open(11,file=fl,status='old')
!          Do i=1,4
!             Read(11,*)
!          End do
!          Do i=1,nhoi*nhoj
!             Read(11,*) qqj,qqi,V2(i)
!          End do
!          Close(11)
          Call Spfn2_gen2(1,imode,jmode)
          k=1
          Do i=1,nhoi
          Do j=1,nhoj
             V2(k)=Spfn2(qqi(i),qqj(j))
             k=k+1
          End do
          End do
          Call Spfn2_Dispose()

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine reset_restartw

      USE mkSx_mod

      Implicit None

        Close(Irst2)
        Call System('rm restart-w')

        Open(Irst2,file='restart-w',status='unknown')
        Call write_e(Irst2,0,V0)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Subroutine CHODVRGrid(Nfree,qmod,maxnho,nho,quad)

        Implicit None

        Integer :: Nfree,qmod(Nfree),maxnho,nho(Nfree)
        Real(8) :: quad(maxnho,Nfree)

        Integer :: i,j,k,l,n
        Integer :: rCHO,nHO1,nCHO1,Mn,sn(3)

        Real(8) :: const,PI
        Real(8), dimension(:,:), allocatable :: xx,xho,Cwfn
        Real(8), dimension(:), allocatable :: qq,xl,Omega

          !-------------------------------------
          ! >>  Quadrature points / Angs(amu)1/2
          !-------------------------------------

          Allocate(Omega(Nfree))
          Call Get_omega(Omega)

          Call file_indicator(50,rCHO)
          Open(rCHO,file='cho-r.wfn',status='old',form='unformatted')

          PI=Acos(-1.D+00)
          Do n=1,Nfree

             Read(rCHO) nHO1,nCHO1
             if(nCHO1<nho(n)) then
                write(6,*) 'ERROR: Error occured while setting up CHO-DVR grid points.'
                write(6,*) 'ERROR: vmax is too large. It must be smaller than :',nCHO1
                Stop
             endif
             Allocate(Cwfn(nHO1,nCHO1))
             Read(rCHO) Cwfn
             !write(6,*) nHO1,nCHO1

             Allocate(xho(nHO1,nHO1))
             Call HO_xmat(nHO1,Omega(n),xho)
             !Do j=0,nHO1-1
             !   write(6,'(20f8.4)') xho(:,j)
             !End do

             ! xx = CHO*xho*CHO
             l=nho(n)
             Allocate(xx(l,l),xl(l*l),qq(l))
             Do i=1,nho(n)
                xx(i,i)=0.D+00
                Do k=1,nHO1
                Do l=1,nHO1
                   xx(i,i)=xx(i,i) + Cwfn(l,i)*xho(l,k)*Cwfn(k,i)
                End do
                End do
             End do

             Do i=1,nho(n)
             Do j=1,i-1
                xx(j,i)=0.D+00
                Do k=1,nHO1
                Do l=1,nHO1
                   xx(j,i)=xx(j,i) + Cwfn(l,j)*xho(l,k)*Cwfn(k,i)
                End do
                End do
                xx(i,j)=xx(j,i)
             End do
             End do
             !Do j=1,nho(n)
             !   write(6,'(20f8.4)') xx(:,j)
             !End do

             Call huckeler(nho(n),nho(n),xx,qq,xl)
             !Call diag2(nho(n),nho(n),xx,xl,qq)
             !write(6,*)
             !write(6,'(20f8.4)') qq
             !write(6,*)
             !Do j=1,l
             !   write(6,'(20f8.4)') xl((j-1)*l+1:j*l)
             !End do

             Call SymOp(n,Mn,sn)
             if(Mn/=0) then
                qq((nho(n)+1)/2)=0.D+00
             endif

             Do j=1,nho(n)
                quad(j,n)=qq(j) 
             End do

             Deallocate(xx,xl,qq,Cwfn,xho)

          End do
          Close(rCHO)

          Deallocate(Omega)

          !-------------------------------------

      End subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Subroutine HODVRGrid(Nfree,qmod,maxnho,nho,quad)

        Implicit None

        Integer :: Nfree,qmod(Nfree),maxnho,nho(Nfree)
        Real(8) :: quad(maxnho,Nfree)

        Integer :: i,j,k,l

        Real(8) :: const,PI
        Real(8), dimension(:,:), allocatable :: xx
        Real(8), dimension(:), allocatable :: qq,xl,HO,Omega

          !-------------------------------------
          ! >>  Quadrature points / Angs(amu)1/2
          !-------------------------------------

          Allocate(Omega(Nfree))
          Call Get_omega(Omega)

          PI=Acos(-1.D+00)
          Do i=1,Nfree
             l=nho(i)
             Allocate(xx(0:l-1,0:l-1),xl(l*l),qq(l),HO(l))

             Call HO_xmat(nho(i),Omega(i),xx)
             !Do j=0,l-1
             !   write(6,'(20f8.4)') xx(:,j)
             !End do

             Call huckeler(l,l,xx,qq,xl)
             qq((nho(i)+1)/2)=0.D+00
             !write(6,*)
             !write(6,'(20f8.4)') qq
             !write(6,*)
             !Do j=1,l
             !   write(6,'(20f8.4)') xl((j-1)*l+1:j*l)
             !End do

             Do j=1,l
                quad(j,i)=qq(j) 
             End do

             Deallocate(xx,xl,qq,HO)

          End do
          Deallocate(Omega)


          !-------------------------------------

      End subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
