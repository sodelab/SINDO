!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Module mkrpt_mod


        Integer :: Nat,Nat3,Nfree
        Integer, dimension(:), allocatable :: qmod

        !     nho(Nfree)          :: The number of basis functions for each mode
        !     int=maxval(nho)     :: Max number of Hermitte-Gauss grids

        Integer, dimension(:), allocatable :: nho
        Integer :: int

        !     MR                  :: Mode coupling representations
        !     quad(int,Nfree)     :: Quadrature points / Angs(amu)1/2

        Integer :: MR
        Real(8), dimension(:,:), allocatable :: quad

        !     Potential energy
        Real(8) :: V0
        Real(8), dimension(:), allocatable :: V1,V2,V3

        !     Dipole moment
        Real(8) :: d0(3)
        Real(8), dimension(:,:), allocatable :: d1,d2,d3

        !     NMR-CC
        Real(8), dimension(:), allocatable :: ncc0
        Real(8), dimension(:,:), allocatable :: ncc1

        !     Irpt,Irpt2,Irst,Irst2 :: File indicators
        !

        Integer :: Irpt,Irpt2,Irst,Irst2
        Integer, dimension(:), allocatable :: idx1,idx2,idx3,idx4

        Logical :: dryrun
        Logical :: CHODVR

      End Module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      PROGRAM main

      Implicit None

          Call v2AI_Const()
          Call rpt_Inp()
          Call mkrpt_nMR_PES()
          Call mkrpt_finalz
          Call v2AI_Dest()

      End 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine rpt_Inp()

      USE mkrpt_mod

      Implicit None

        Integer, parameter :: nmax=100

        Integer :: i,j,k,l
        Integer :: Get_Nat,Get_Nfree
        Integer :: vmax_base
        Integer, dimension(nmax) :: mode,vmax
        Logical :: op,Get_dpl,Get_nmrcc,Get_dryrun

        Namelist /mkrpt/MR,Nfree,mode,vmax,vmax_base,CHODVR

!----------------------------------------------------------------------
        Write(6,100)
  100   Format(/,'--------------------(    MKRPT MODULE   )--------------------',/)

! --    Read parameters
        Nat=Get_Nat()
        Nfree=Get_Nfree()
        Do i=1,Nfree
           mode(i)=i
        End do
        MR=2
        vmax=-1

        Rewind(5)
        Read(5,mkrpt)

! --    Mode
        Allocate(qmod(Nfree))
        qmod=mode(1:Nfree)
        Do i=1,Nfree
           k=qmod(i)
           Do j=i+1,Nfree
              if(k>qmod(j)) then
                 l=k
                 k=qmod(j)
                 qmod(j)=l
              endif
           End do
           qmod(i)=k
        End do 

        write(6,120) Nfree,qmod
        write(6,121) MR
    120 Format(' -> Options',/, &
              2x,'> Nfree:',i8,/, &
              2x,'> Mode :',100i3)
    121 Format(2x,'> MR   :',i8)
  
! --    # of HO function
        Allocate(nho(Nfree))
        nho=vmax_base
        Do i=1,Nfree
           if(vmax(i)/=-1) nho(i)=vmax(i)
           if(mod(nho(i),2)/=0) nho(i)=nho(i)+1
        End do
        nho=nho+1
        !  # of grids
        int=maxval(nho)
        Allocate(quad(int,Nfree))
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!MK        
        if(CHODVR) then
           Call CHODVRGrid(Nfree,mode,int,nho,quad)
        else
           Call HODVRGrid(Nfree,mode,int,nho,quad)
        endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
        write(*,*) "Except!" 
        write(6,200)
        Do i=1,Nfree
           write(6,210) qmod(i)
           write(6,220) quad(:,i)
        End do
    200 Format(2x,'> DVR Grids')
    210 Format(6x,'Mode :',i3)
    220 Format(6x,6f9.4)
        
        Irpt=200
        Do while(.true.)
           Inquire(Irpt,opened=op)
           if(op) then
             Irpt=Irpt+1
           else
             exit
           endif
        End do
        Open(Irpt,file='rpt.dat',status='unknown')
        Open(155,file='PES',status='unknown')
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

        if(Get_dpl() .or. Get_nmrcc()) then
           Irpt2=230
           Do while(.true.)
              Inquire(Irpt2,opened=op)
              if(op) then
                Irpt2=Irpt2+1
              else
                exit
              endif
           End do
           Open(Irpt2,file='rpt2.dat',status='unknown')

           write(6,151)
        else
           write(6,150)
        endif

    150 Format(/,' -> Data           FILE :  [ rpt.dat ]',/, &
                 ' -> Restart(read)  FILE :  [ restart-r ]',/, &
                 ' -> Restart(write) FILE :  [ restart-w ]',/)
    151 Format(/,' -> Data (energy)  FILE :  [ rpt.dat ]',/, &
                 '         (others)  FILE :  [ rpt2.dat ]',/, &
                 ' -> Restart(read)  FILE :  [ restart-r ]',/, &
                 ' -> Restart(write) FILE :  [ restart-w ]',/)

        dryrun=Get_dryrun()

        return


      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Subroutine mkrpt_finalz

      USE mkrpt_mod

      Implicit None 

      Integer :: in,io

          Call spr_Getio(in,io)
          Write(io,100)
      100 Format(/,'(  FINALIZE MKRPT MODULE  )',/)

          Deallocate(nho,quad)

          if(allocated(idx1)) Deallocate(idx1)
          if(allocated(idx2)) Deallocate(idx2)
          if(allocated(idx3)) Deallocate(idx3)

          if(allocated(V1)) Deallocate(V1)
          if(allocated(V2)) Deallocate(V2)
          if(allocated(V3)) Deallocate(V3)

          if(allocated(d1)) Deallocate(d1)
          if(allocated(d2)) Deallocate(d2)
          if(allocated(d3)) Deallocate(d3)

          if(allocated(ncc0)) Deallocate(ncc0)
          if(allocated(ncc1)) Deallocate(ncc1)


      End subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mkrpt_nMR_PES()

      USE mkrpt_mod

      Implicit None

      ! ------------------------------------------------------------------------

      Integer :: io,in
      !Integer :: GetKey
      Logical :: dpl,get_dpl,nmrcc,get_nmrcc

      ! ------------------------------------------------------------------------

          dpl=Get_dpl()
          nmrcc=Get_nmrcc()
          Call spr_GetIO(in,io)

          !----------------------
          ! >>  Q=0
          !----------------------
          write(io,100) 
      100 Format(' -> Energies on grids ',/, &
                2x,'> Q=0')
          Select case(GetKey())
            case(0)
              Call mkrpt_e_0(io)
            case(1)
              Call mkrpt_ed_0(io)
            case(2)
              Call mkrpt_en_0(io)
          End select

          !----------------------
          ! >>  1MR PES 
          !----------------------
          Select case(GetKey())
            case(0)
              Call mkrpt_e_1MR(io)
            case(1)
              Call mkrpt_ed_1MR(io)
            case(2)
              Call mkrpt_en_1MR(io)
          End select

          if(MR==1) return

          !----------------------
          ! >>  2MR PES 
          !----------------------
          Select case(GetKey())
            case(0)
              Call mkrpt_e_2MR(io)
            case(1)
              !Call mkrpt_ed_2MR(io)
          End select

          if(MR==2) return

          !----------------------
          ! >>  3MR PES 
          !----------------------
          Select case(GetKey())
            case(0)
              Call mkrpt_e_3MR(io)
            case(1)
              !Call mkrpt_ed_3MR(io)
          End select

          ! --------------------------------------------------------------------

      Contains

      Function GetKey()

      Implicit None

      Integer :: GetKey, key
!      Integer :: key

           if(.not. dpl .and. .not. nmrcc) then
             ! Energy
             key=0
           elseif(dpl) then
             ! Energy and dipole
             key=1
           elseif(nmrcc) then
             ! Energy and NMR-CC
             key=2
           endif
           GetKey=key

         End function

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_e_0(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io
      Real(8), dimension(:), allocatable   :: qq

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree))
          qq=0.D+00
          Call mkrpt_airun_e(0,qq,V0)
          Write(io,2000) 0,V0
          Write(Irpt,'(i8,f18.10)') 0,V0
          Deallocate(qq)

     2000 Format(i8,f18.10)!MK 18>20 like mr1-2-3

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_e_1MR(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io

      Integer :: i,iq1,idx,Mi,si(3)
      Integer, dimension(:), allocatable :: nhalf
      Real(8), dimension(:), allocatable :: qq
      Real(8), dimension(:,:), allocatable :: VV

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree),nhalf(Nfree))
          qq=0.D+00
          Do i=1,Nfree
             nhalf(i)=(nho(i)+1)/2
          End do

          Allocate(idx1(Nfree+1))
          idx1(1)=0
          Do i=2,Nfree+1
             idx1(i)=idx1(i-1)+nho(i-1)
          End do
          Allocate(V1(idx1(Nfree)+nho(Nfree)))

          Do i=1,Nfree
             Write(io,1000) qmod(i)
             Call SymOp(qmod(i),Mi,si)

             if(Mi==0) then
                Do iq1=1,nhalf(i)-1
                   qq(i)=quad(iq1,i)
                   idx=idx1(i)+iq1
                   Call mkrpt_airun_e(idx,qq,V1(idx))
                   Write(io,2000) idx,V1(idx),qq(i)
                End do

                qq(i)=0.D+00
                idx=idx1(i)+nhalf(i)
                V1(idx)=V0
                Write(io,2000) idx,V1(idx),qq(i)

                Do iq1=nhalf(i)+1,nho(i)
                   qq(i)=quad(iq1,i)
                   idx=idx1(i)+iq1
                   Call mkrpt_airun_e(idx,qq,V1(idx))
                   Write(io,2000) idx,V1(idx),qq(i)
                End do

             else
                Do iq1=1,nhalf(i)-1
                   qq(i)=quad(iq1,i)
                   idx=idx1(i)+iq1
                   Call mkrpt_airun_e(idx,qq,V1(idx))
                   Write(io,2000) idx,V1(idx),qq(i)
                End do

                qq(i)=0.D+00
                idx=idx1(i)+nhalf(i)
                V1(idx)=V0
                Write(io,2000) idx,V1(idx),qq(i)

                Do iq1=1,nhalf(i)-1
                   V1(idx+iq1)=V1(idx-iq1)
                End do
                Do iq1=nhalf(i)+1,nho(i)
                   qq(i)=quad(iq1,i)
                   idx=idx1(i)+iq1
                   Write(io,2000) idx,V1(idx),qq(i)
                End do

             endif
             qq(i)=0.D+00

!             Allocate(VV(nho(i),1))
!             VV(:,1)=V1(idx1(i)+1:idx1(i+1))
!             Call intpl1(nho(i),1,quad(:,i),VV)
!             V1(idx1(i)+1:idx1(i+1))=VV(:,1)
!             Deallocate(VV)
             !MK if(.not. dryrun) &
              !MK  Call intpl1(nho(i),1,quad(:,i),V1(idx1(i)+1:idx1(i+1)))

          End do

          ! --------------------------------------------------------------------
          Do i=1,idx1(Nfree)+nho(Nfree)
             Write(Irpt,'(i8,f18.10)') i,V1(i)
          End do
          ! --------------------------------------------------------------------

          Deallocate(qq,nhalf)

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mkrpt_e_2MR(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io

      Integer :: i,j,ij,iq1,iq2,idx,jdx,Mi,Mj,si(3),sj(3)
      Integer, dimension(:), allocatable :: nhalf
      Real(8), dimension(:), allocatable :: qq
      Real(8), dimension(:,:), allocatable :: VV

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree),nhalf(Nfree))
          qq=0.D+00
          Do i=1,Nfree
             nhalf(i)=(nho(i)+1)/2
          End do

          ! --------------------------------------------------------------------
          Allocate(idx2(Nfree*(Nfree-1)/2+1))
          ij=0; idx2(1)=0
          Do i=2,Nfree
             Do j=1,i-1
                ij= ij+1
                idx2(ij+1)=idx2(ij)+nho(i)*nho(j)
             End do
          End do
          !dbg write(6,'(10i4)') idx2,m
          Allocate(V2(idx2(Nfree*(Nfree-1)/2+1)))
          ! -----------------------------------------------------------------

          Do i=2,Nfree
          Call SymOp(qmod(i),Mi,si)
          Do j=1,i-1
             Call SymOp(qmod(j),Mj,sj)
             ij= i*(i-1)/2 - (i-1) +j
             !dbg write(io,'(3i3)') i,j,ij
             Write(io,1000) qmod(i),qmod(j)

             ! qj=0 
             Do iq1=1,nho(i)
                qq(i)=quad(iq1,i)
                iq2=nhalf(j)
                qq(j)=quad(iq2,j)
                idx=idx2(ij)+(iq1-1)*nho(j)+iq2
                V2(idx)=V1(idx1(i)+iq1)
                Write(io,2000) idx,V2(idx),qq(i),qq(j)
             End do

             ! qi=0
             iq1=nhalf(i)
             qq(i)=quad(iq1,i)
             Do iq2=1,nho(j)
                qq(j)=quad(iq2,j)
                idx=idx2(ij)+(iq1-1)*nho(j)+iq2
                V2(idx)=V1(idx1(j)+iq2)
                Write(io,2000) idx,V2(idx),qq(i),qq(j)
             End do

             if(Mi==0) then

                if(Mj==0) then

                   Call PEF2D(0,0,Nfree,i,j,idx2(ij),nho(i),nho(j), &
                        quad(:,i),quad(:,j),V2(idx2(ij)+1:idx2(ij+1)))

                else

                   Call PEF2D(0,1,Nfree,i,j,idx2(ij),nho(i),nho(j), &
                        quad(:,i),quad(:,j),V2(idx2(ij)+1:idx2(ij+1)))

                   ! qj<0 region : generate by symmetry
                   Do iq1=1,nho(i)
                      if(iq1==nhalf(i)) cycle
                      qq(i)=quad(iq1,i)
                      Do iq2=nhalf(j)+1,nho(j)
                         qq(j)=quad(iq2,j)
                         idx=idx2(ij)+(iq1-1)*nho(j)+iq2
                         jdx=idx2(ij)+(iq1-1)*nho(j)+nho(j)-iq2+1
                         V2(idx)=V2(jdx)
                         Write(io,2000) idx,V2(idx),qq(i),qq(j)
                      End do
                   End do

                endif

             else
                if(Mj==0) then

                   Call PEF2D(1,0,Nfree,i,j,idx2(ij),nho(i),nho(j), &
                        quad(:,i),quad(:,j),V2(idx2(ij)+1:idx2(ij+1)))

                   ! qi<0 region : generate by symmetry
                   Do iq1=nhalf(i)+1,nho(i)
                      qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                         idx=idx2(ij)+(iq1-1)*nho(j)+iq2
                         jdx=idx2(ij)+(nho(i)-iq1)*nho(j)+iq2
                         V2(idx)=V2(jdx)
                         Write(io,2000) idx,V2(idx),qq(i),qq(j)
                      End do
                   End do

                elseif(Mi==Mj) then
                   Call PEF2D(1,0,Nfree,i,j,idx2(ij),nho(i),nho(j), &
                        quad(:,i),quad(:,j),V2(idx2(ij)+1:idx2(ij+1)))

                   ! qi<0 region : generate by symmetry
                   Do iq1=nhalf(i)+1,nho(i)
                      qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                         idx=idx2(ij)+(iq1-1)*nho(j)+iq2
                         jdx=idx2(ij)+(nho(i)-iq1)*nho(j)+nho(j)-iq2+1
                         V2(idx)=V2(jdx)
                         Write(io,2000) idx,V2(idx),qq(i),qq(j)
                      End do
                   End do

                else
                   Call PEF2D(1,1,Nfree,i,j,idx2(ij),nho(i),nho(j), &
                        quad(:,i),quad(:,j),V2(idx2(ij)+1:idx2(ij+1)))

                   ! qi>0, qj<0 region
                   Do iq1=1,nhalf(i)-1
                      qq(i)=quad(iq1,i)
                      Do iq2=nhalf(j)+1,nho(j)
                         qq(j)=quad(iq2,j)
                         idx=idx2(ij)+(iq1-1)*nho(j)+iq2
                         jdx=idx2(ij)+(iq1-1)*nho(j)+nho(j)-iq2+1
                         V2(idx)=V2(jdx)
                         Write(io,2000) idx,V2(idx),qq(i),qq(j)
                      End do
                   End do

                   ! qi<0 region : generate by symmetry
                   Do iq1=nhalf(i)+1,nho(i)
                      qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                         idx=idx2(ij)+(iq1-1)*nho(j)+iq2
                         jdx=idx2(ij)+(nho(i)-iq1)*nho(j)+iq2
                         V2(idx)=V2(jdx)
                         Write(io,2000) idx,V2(idx),qq(i),qq(j)
                      End do
                   End do

                endif
             endif

             qq(i)=0.D+00
             qq(j)=0.D+00

!             Allocate(VV(nho(i)*nho(j),1))
!             VV(:,1)=V2(idx2(ij)+1:idx2(ij+1))
!             Call intpl2(nho(i),nho(j),1,quad(:,j),VV)
!             V2(idx2(ij)+1:idx2(ij+1))=VV(:,1)
!             Deallocate(VV)
!MK             if(.not. dryrun) &
!MK                Call intpl2(nho(i),nho(j),1,quad(:,j),V2(idx2(ij)+1:idx2(ij+1)))

          End do
          End do

          ! --------------------------------------------------------------------
          Do i=1,idx2(Nfree*(Nfree-1)/2+1)
             Write(Irpt,'(i8,f20.10)') i,V2(i)
          End do
          ! --------------------------------------------------------------------
          Deallocate(qq,nhalf)

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,2f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     ir>0  ... qi>0
!     ir<0  ... qi<0
!     ir=0  ... diabled

!     jr>0  ... qj>0
!     jr<0  ... qj<0
!     jr=0  ... diabled 

      Subroutine PEF2D(ir,jr,Nfree,ii,jj,idx0,ni,nj,qqi,qqj,V2)

      Implicit None

      Integer :: ir,jr,Nfree,ii,jj,idx0,ni,nj
      Real(8) :: qqi(ni),qqj(nj),V2(ni*nj)

      Integer :: iq1,iq2,nhi,nhj,idx
      Real(8) :: qq(Nfree),dp(3)

         nhi=(ni+1)/2
         nhj=(nj+1)/2
         qq=0.D+00

         ! qi>0
         if(ir>=0) then

         ! qi>0, qj>0 region
         if(jr>=0) then
         Do iq1=1,nhi-1
            qq(ii)=qqi(iq1)
            Do iq2=1,nhj-1
               qq(jj)=qqj(iq2)
               idx=(iq1-1)*nj+iq2
               Call mkrpt_airun_e(idx+idx0,qq,V2(idx))
               Write(6,2000) idx+idx0,V2(idx),qq(ii),qq(jj)
            End do
         End do
         endif

         ! qi>0, qj<0 region
         if(jr<=0) then
         Do iq1=1,nhi-1
            qq(ii)=qqi(iq1)
            Do iq2=nhj+1,nj
               qq(jj)=qqj(iq2)
               idx=(iq1-1)*nj+iq2
               Call mkrpt_airun_e(idx+idx0,qq,V2(idx))
               Write(6,2000) idx+idx0,V2(idx),qq(ii),qq(jj)
            End do
         End do
         endif
         endif

         if(ir<=0) then

         ! qi<0, qj>0 region
         if(jr>=0) then
         Do iq1=nhi+1,ni
            qq(ii)=qqi(iq1)
            Do iq2=1,nhj-1
               qq(jj)=qqj(iq2)
               idx=(iq1-1)*nj+iq2
               Call mkrpt_airun_e(idx+idx0,qq,V2(idx))
               Write(6,2000) idx+idx0,V2(idx),qq(ii),qq(jj)
            End do
         End do
         endif

         ! qi<0, qj<0 region
         if(jr<=0) then
         Do iq1=nhi+1,ni
            qq(ii)=qqi(iq1)
            Do iq2=nhj+1,nj
               qq(jj)=qqj(iq2)
               idx=(iq1-1)*nj+iq2
               Call mkrpt_airun_e(idx+idx0,qq,V2(idx))
               Write(6,2000) idx+idx0,V2(idx),qq(ii),qq(jj)
            End do
         End do
         endif

         endif

    2000 Format(i8,f18.10,100f8.3)

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mkrpt_e_3MR(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io

      Integer :: i,j,k,ij,jk,ik,ijk,iq1,iq2,iq3,idx,jdx, &
                 Mi,Mj,Mk,si(3),sj(3),sk(3)
      Integer, dimension(:), allocatable :: nhalf
      Real(8), dimension(:), allocatable :: qq
      Real(8), dimension(:,:), allocatable :: VV

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree),nhalf(Nfree))
          qq=0.D+00
          Do i=1,Nfree
             nhalf(i)=(nho(i)+1)/2
          End do

          ! --------------------------------------------------------------------
          Allocate(idx3(Nfree*(Nfree-1)*(Nfree-2)/6+1))

          ijk=0; idx3(1)=0
          Do i=3,Nfree
          Do j=2,i-1
          Do k=1,j-1
             ijk=((i-1)*i*(2*i-1)/6-3*(i-1)*i/2)/2 + (i-1) + (j-2)*(j-1)/2 + k
             idx3(ijk+1)=idx3(ijk)+nho(i)*nho(j)*nho(k)
             !dbg write(6,'(4i3)') i,j,k,ijk
          End do
          End do
          End do
          !dbg write(6,'(10i6)') idx3,m

          Allocate(V3(idx3(Nfree*(Nfree-1)*(Nfree-2)/6+1)))
          V3=0.D+00

          ! --------------------------------------------------------------------
          ! >>>>>>>
          ijk=0
          Do i=3,Nfree
          Call SymOp(qmod(i),Mi,si)
          Do j=2,i-1
          Call SymOp(qmod(j),Mj,sj)
          Do k=1,j-1
          Call SymOp(qmod(k),Mk,sk)
             Write(io,1000) qmod(i),qmod(j),qmod(k)
             ijk=ijk+1
             ij= i*(i-1)/2 -(i-1)+j
             jk= j*(j-1)/2 -(j-1)+k
             ik= i*(i-1)/2 -(i-1)+k
             !dbg write(6,'(3i3,2x,3i3)') i,j,k,ij,jk,ik

             ! qi=0
             iq1=nhalf(i)
                qq(i)=quad(iq1,i)
             Do iq2=1,nho(j)
                qq(j)=quad(iq2,j)
             Do iq3=1,nho(k)
                qq(k)=quad(iq3,k)
                idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                jdx=idx2(jk)+(iq2-1)*nho(k)+iq3
                V3(idx)=V2(jdx)
                Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
             End do
             End do

             ! qj=0
             Do iq1=1,nho(i)
                qq(i)=quad(iq1,i)
             iq2=nhalf(j)
                qq(j)=quad(iq2,j)
             Do iq3=1,nho(k)
                qq(k)=quad(iq3,k)
                idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                jdx=idx2(ik)+(iq1-1)*nho(k)+iq3
                V3(idx)=V2(jdx)
                Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
             End do
             End do

             ! qk=0
             Do iq1=1,nho(i)
                qq(i)=quad(iq1,i)
             Do iq2=1,nho(j)
                qq(j)=quad(iq2,j)
             iq3=nhalf(k)
                qq(k)=quad(iq3,k)
                idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                jdx=idx2(ij)+(iq1-1)*nho(j)+iq2
                V3(idx)=V2(jdx)
                Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
             End do
             End do

             if(Mi==0) then
                if(Mj==0) then
                   if(Mk==0) then
                      ! OK
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         Call mkrpt_airun_e(idx,qq,V3(idx))
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                   else
                      ! OK
                      ! qk>0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nhalf(k)-1
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         Call mkrpt_airun_e(idx,qq,V3(idx))
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                      ! qk<0 region : generate by symmetry
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                      Do iq3=nhalf(k)+1,nho(k)
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                      +(iq2-1)*nho(k)+nho(k)-iq3+1
                         V3(idx)=V3(jdx)
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
                
                      End do
                      End do
                      End do

                   endif

                else
                   if(Mk==0) then
                      ! qj>0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nhalf(j)-1
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         Call mkrpt_airun_e(idx,qq,V3(idx))
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                      ! qj<0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=nhalf(j)+1,nho(j)
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) & 
                                      +(nho(j)-iq2)*nho(k)+iq3
                         V3(idx)=V3(jdx)
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                   elseif(Mk==Mj) then
                      ! qj>0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nhalf(j)-1
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         Call mkrpt_airun_e(idx,qq,V3(idx))
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                      ! qj<0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=nhalf(j)+1,nho(j)
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                      +(nho(j)-iq2)*nho(k)+nho(k)-iq3+1
                         V3(idx)=V3(jdx)
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                   else !(Mk/=Mj)
                      ! OK
                      ! qj>0,qk>0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nhalf(j)-1
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nhalf(k)-1
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         Call mkrpt_airun_e(idx,qq,V3(idx))
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                      ! qj>0,qk<0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nhalf(j)-1
                         qq(j)=quad(iq2,j)
                      Do iq3=nhalf(k)+1,nho(k)
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k) &
                                      +nho(k)-iq3+1
                         V3(idx)=V3(jdx)
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                      ! qj<0 region
                      Do iq1=1,nho(i)
                         if(iq1==nhalf(i)) cycle
                         qq(i)=quad(iq1,i)
                      Do iq2=nhalf(j)+1,nho(j)
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                      +(nho(j)-iq2)*nho(k)+iq3
                         V3(idx)=V3(jdx)
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
                
                      End do
                      End do
                      End do


                   endif

                endif

             else

                if(Mj==0) then
                   if(Mk==0 .or. Mk==Mi) then
                      ! qi>0 region
                      Do iq1=1,nhalf(i)-1
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)

                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         Call mkrpt_airun_e(idx,qq,V3(idx))
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do

                      End do
                      End do

                      if(Mk==0) then
                         ! OK
                         ! qi<0 region
                         Do iq1=nhalf(i)+1,nho(i)
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                         End do
                         End do
                         End do

                      else
                         ! Mk==Mi
                         ! OK
                         ! qi<0 region
                         Do iq1=nhalf(i)+1,nho(i)
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+nho(k)-iq3+1
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                         End do
                         End do
                         End do

                      endif

                   else
                      ! OK
                      ! Mk /=0 .and. Mi /= Mk
                      ! qi>0,qk>0 region
                      Do iq1=1,nhalf(i)-1
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nhalf(k)-1
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         Call mkrpt_airun_e(idx,qq,V3(idx))
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                      ! qi>0,qk<0 region
                      Do iq1=1,nhalf(i)-1
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                      Do iq3=nhalf(k)+1,nho(k)
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k) &
                                      + nho(k)-iq3+1
                         V3(idx)=V3(jdx)
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do

                      ! qi<0 region
                      Do iq1=nhalf(i)+1,nho(i)
                         qq(i)=quad(iq1,i)
                      Do iq2=1,nho(j)
                         if(iq2==nhalf(j)) cycle
                         qq(j)=quad(iq2,j)
                      Do iq3=1,nho(k)
                         if(iq3==nhalf(k)) cycle
                         qq(k)=quad(iq3,k)
                         idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k)+(iq2-1)*nho(k)+iq3
                         jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                      +(iq2-1)*nho(k)+iq3
                         V3(idx)=V3(jdx)
                         Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)

                      End do
                      End do
                      End do
                   endif
                else
                   if(Mk==0) then
                      if(Mi==Mj) then
                         ! OK
                         ! qi>0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            Call mkrpt_airun_e(idx,qq,V3(idx))
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi<0 region
                         Do iq1=nhalf(i)+1,nho(i)
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(nho(j)-iq2)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                      else
                         ! CHECKED
                         ! qi>0,qj>0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nhalf(j)-1
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            Call mkrpt_airun_e(idx,qq,V3(idx))
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi>0,qj<0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=nhalf(j)+1,nho(j)
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(nho(j)-iq2)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi<0 region
                         Do iq1=nhalf(i)+1,nho(i)
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                      endif
                   else
                      if(Mi==Mj .and. Mj/=Mk) then
                         ! OK
                         ! qi>0,qk>0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nhalf(k)-1
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            Call mkrpt_airun_e(idx,qq,V3(idx))
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi<0,qk>0 region
                         Do iq1=nhalf(i)+1,nho(i)
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nhalf(k)-1
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(nho(j)-iq2)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qk<0 region
                         Do iq1=1,nho(i)
                            if(iq1==nhalf(i)) cycle
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=nhalf(k)+1,nho(k)
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+nho(k)-iq3+1
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                      elseif(Mi==Mk .and. Mk/=Mj) then
                         ! qj>0,qk>0 region
                         Do iq1=1,nho(i)
                            if(iq1==nhalf(i)) cycle
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nhalf(j)-1
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nhalf(k)-1
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            Call mkrpt_airun_e(idx,qq,V3(idx))
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qj>0,qk<0 region
                         Do iq1=1,nho(i)
                            if(iq1==nhalf(i)) cycle
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nhalf(j)-1
                            qq(j)=quad(iq2,j)
                         Do iq3=nhalf(k)+1,nho(k)
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+nho(k)-iq3+1
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qj<0 region
                         Do iq1=1,nho(i)
                            if(iq1==nhalf(i)) cycle
                            qq(i)=quad(iq1,i)
                         Do iq2=nhalf(j)+1,nho(j)
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(nho(j)-iq2)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                      elseif(Mi/=Mj .and. Mj==Mk) then
                         ! qi>0,qj>0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nhalf(j)-1
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            Call mkrpt_airun_e(idx,qq,V3(idx))
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi>0,qj<0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=nhalf(j)+1,nho(j)
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(nho(j)-iq2)*nho(k)+nho(k)-iq3+1
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi<0 region
                         Do iq1=nhalf(i)+1,nho(i)
                            qq(i)=quad(iq1,i)
                         Do iq2=nhalf(j)+1,nho(j)
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                      elseif(Mi==Mj .and. Mj==Mk) then
                         ! qk>0 region
                         Do iq1=1,nho(i)
                            if(iq1==nhalf(i)) cycle
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nhalf(k)-1
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            Call mkrpt_airun_e(idx,qq,V3(idx))
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qk<0 region
                         Do iq1=1,nho(i)
                            if(iq1==nhalf(i)) cycle
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=nhalf(k)+1,nho(k)
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(nho(j)-iq2)*nho(k)+nho(k)-iq3+1
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                      else !(Mi/=Mj/=Mk)
                         ! qi>0,qj>0,qk>0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nhalf(j)-1
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nhalf(k)-1
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            Call mkrpt_airun_e(idx,qq,V3(idx))
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi>0,qj>0,qk<0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nhalf(j)-1
                            qq(j)=quad(iq2,j)
                         Do iq3=nhalf(k)+1,nho(k)
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+nho(k)-iq3+1
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi>0,qj<0 region
                         Do iq1=1,nhalf(i)-1
                            qq(i)=quad(iq1,i)
                         Do iq2=nhalf(j)+1,nho(j)
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(nho(j)-iq2)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                         ! qi<0 region
                         Do iq1=nhalf(i)+1,nho(i)
                            qq(i)=quad(iq1,i)
                         Do iq2=1,nho(j)
                            if(iq2==nhalf(j)) cycle
                            qq(j)=quad(iq2,j)
                         Do iq3=1,nho(k)
                            if(iq3==nhalf(k)) cycle
                            qq(k)=quad(iq3,k)
                            idx=idx3(ijk)+(iq1-1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            jdx=idx3(ijk)+(nho(i)-iq1)*nho(j)*nho(k) &
                                         +(iq2-1)*nho(k)+iq3
                            V3(idx)=V3(jdx)
                            Write(io,2000) idx,V3(idx),qq(i),qq(j),qq(k)
      
                         End do
                         End do
                         End do

                      endif
                   endif
                endif

             endif
             qq(i)=0.D+00
             qq(j)=0.D+00
             qq(k)=0.D+00

!MK             if(.not. dryrun) &
 !MK               Call intpl3(nho(i),nho(j),nho(k),1, & 
  !MK                          quad(:,k),V3(idx3(ijk)+1:idx3(ijk+1)))
          End do
          End do
          End do

          ! --------------------------------------------------------------------
          Do i=1,idx3(Nfree*(Nfree-1)*(Nfree-2)/6+1)
             Write(Irpt,'(i8,f20.10)') i,V3(i)
          End do
          ! --------------------------------------------------------------------
          Deallocate(qq,nhalf)

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,3f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_airun_e(indx0,qq,E)

      USE mkrpt_mod

      Implicit None

      Integer :: indx0,indx1
      Real(8) :: qq(Nfree),E
      Integer :: i,Nf,Get_Nfree
      Real(8), dimension(:), allocatable :: qfull

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
         Nf=Get_Nfree()
         Allocate(qfull(Nf))
         qfull=0.D+00
         Do i=1,Nfree
            qfull(qmod(i))=qq(i)
         End do
         Call Run_e(qfull,E)
         Call write_e(Irst2,indx0,E)
         Call myFlsh(Irst2)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_ed_0(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io
      Real(8), dimension(:), allocatable   :: qq

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree))
          qq=0.D+00
          Call mkrpt_airun_ed(0,qq,V0,d0)
          Write(io,2000) 0,V0
          Write(Irpt,'(i8,f20.10)') 0,V0
          Write(Irpt2,'(i8,3f12.8)') 0,d0
          Deallocate(qq)

     2000 Format(i8,f18.10)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_ed_1MR(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io

      Integer :: i,j,iq1,idx
      Integer, dimension(:), allocatable :: nhalf
      Real(8), dimension(:), allocatable :: qq
      Real(8), dimension(:,:), allocatable :: VV

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree),nhalf(Nfree))
          qq=0.D+00
          Do i=1,Nfree
             nhalf(i)=(nho(i)+1)/2
          End do

          Allocate(idx1(Nfree+1))
          idx1(1)=0
          Do i=2,Nfree+1
             idx1(i)=idx1(i-1)+nho(i-1)
          End do
          Allocate(V1(idx1(Nfree)+nho(Nfree)))
          Allocate(d1(3,idx1(Nfree)+nho(Nfree)))

          Do i=1,Nfree
             Write(io,1000) qmod(i)

             Do iq1=1,nhalf(i)-1
                qq(i)=quad(iq1,i)
                idx=idx1(i)+iq1
                Call mkrpt_airun_ed(idx,qq,V1(idx),d1(:,idx))
                Write(io,2000) idx,V1(idx),qq(i)
             End do

             qq(i)=0.D+00
             idx=idx1(i)+nhalf(i)
             V1(idx)=V0
             d1(:,idx)=d0
             Write(io,2000) idx,V1(idx),qq(i)

             Do iq1=nhalf(i)+1,nho(i)
                qq(i)=quad(iq1,i)
                idx=idx1(i)+iq1
                Call mkrpt_airun_en(idx,qq,V1(idx),d1(:,idx))
                Write(io,2000) idx,V1(idx),qq(i)
             End do

             qq(i)=0.D+00

             Allocate(VV(nho(i),4))
             VV(:,1)=V1(idx1(i)+1:idx1(i+1))
             VV(:,2)=d1(1,idx1(i)+1:idx1(i+1))
             VV(:,3)=d1(2,idx1(i)+1:idx1(i+1))
             VV(:,4)=d1(3,idx1(i)+1:idx1(i+1))
   !MK          if(.not. dryrun) &
   !MK          Call intpl1(nho(i),4,quad(:,i),VV)
             V1(idx1(i)+1:idx1(i+1))=VV(:,1)
             d1(1,idx1(i)+1:idx1(i+1))=VV(:,2)
             d1(2,idx1(i)+1:idx1(i+1))=VV(:,3)
             d1(3,idx1(i)+1:idx1(i+1))=VV(:,4)
             Deallocate(VV)

          End do

          ! --------------------------------------------------------------------
          Do i=1,idx1(Nfree)+nho(Nfree)
             Write(Irpt,'(i8,f20.10)') i,V1(i)
             Write(Irpt2,'(i8,3f12.4)') i,d1(:,i)
          End do
          ! --------------------------------------------------------------------

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_airun_ed(indx0,qq,E,d)

      USE mkrpt_mod

      Implicit None

      Integer :: indx0,indx1
      Real(8) :: qq(Nfree),E,d(3)
      Integer :: i,Nf,Get_Nfree
      Real(8), dimension(:), allocatable :: qfull

         Read(Irst,*,end=100,err=100) indx1,E,d
         if(indx0 .ne. indx1) then 
            Write(6,*) 'Invalid line is detected in the restart file:'
            Call Write_ed(6,indx1,E,d)
            Write(6,'('' It must have an index '',i8)') indx0
            Stop
         endif
         Call write_ed(Irst2,indx0,E,d)
         Call myFlsh(Irst2)
         return

  100    Continue
         Nf=Get_Nfree()
         Allocate(qfull(Nf))
         qfull=0.D+00
         Do i=1,Nfree
            qfull(qmod(i))=qq(i)
         End do
         Call Run_ed(qfull,E,d)
         Call write_ed(Irst2,indx0,E,d)
         Call myFlsh(Irst2)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_en_0(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io
      Real(8), dimension(:), allocatable   :: qq

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree),ncc0(Nat*(Nat-1)/2))
          qq=0.D+00
          Call mkrpt_airun_en(0,qq,V0,ncc0)
          Write(io,2000) 0,V0
          Write(Irpt,'(i8,f20.10)') 0,V0
          Write(Irpt2,'(i8,100f12.4)') 0,ncc0
          Deallocate(qq)

     2000 Format(i8,f18.10,100f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_en_1MR(io)

      USE mkrpt_mod

      Implicit None

      Integer :: io

      Integer :: i,j,iq1,idx
      Integer, dimension(:), allocatable :: nhalf
      Real(8), dimension(:), allocatable :: qq
      Real(8), dimension(:,:), allocatable :: VV

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree),nhalf(Nfree))
          qq=0.D+00
          Do i=1,Nfree
             nhalf(i)=(nho(i)+1)/2
          End do

          Allocate(idx1(Nfree+1))
          idx1(1)=0
          Do i=2,Nfree+1
             idx1(i)=idx1(i-1)+nho(i-1)
          End do
          Allocate(V1(idx1(Nfree)+nho(Nfree)))
          Allocate(ncc1(Nat*(Nat-1)/2,idx1(Nfree)+nho(Nfree)))

          Do i=1,Nfree
             Write(io,1000) qmod(i)

             Do iq1=1,nhalf(i)-1
                qq(i)=quad(iq1,i)
                idx=idx1(i)+iq1
                Call mkrpt_airun_en(idx,qq,V1(idx),ncc1(:,idx))
                Write(io,2000) idx,V1(idx),qq(i)
             End do

             qq(i)=0.D+00
             idx=idx1(i)+nhalf(i)
             V1(idx)=V0
             ncc1(:,idx)=ncc0
             Write(io,2000) idx,V1(idx),qq(i)

             Do iq1=nhalf(i)+1,nho(i)
                qq(i)=quad(iq1,i)
                idx=idx1(i)+iq1
                Call mkrpt_airun_en(idx,qq,V1(idx),ncc1(:,idx))
                Write(io,2000) idx,V1(idx),qq(i)
             End do

             qq(i)=0.D+00

             Allocate(VV(nho(i),1+Nat*(Nat-1)/2))
             VV(:,1)=V1(idx1(i)+1:idx1(i+1))
             Do j=1,Nat*(Nat-1)/2
                VV(:,j+1)=ncc1(j,idx1(i)+1:idx1(i+1))
             End do
             if(.not. dryrun) &
             Call intpl1(nho(i),1+Nat*(Nat-1)/2,quad(:,i),VV)
             V1(idx1(i)+1:idx1(i+1))=VV(:,1)
             Do j=1,Nat*(Nat-1)/2
                ncc1(j,idx1(i)+1:idx1(i+1))=VV(:,j+1)
             End do
             Deallocate(VV)

          End do

          ! --------------------------------------------------------------------
          Do i=1,idx1(Nfree)+nho(Nfree)
             Write(Irpt,'(i8,f20.10)') i,V1(i)
             Write(Irpt2,'(i8,100f12.4)') i,ncc1(:,i)
          End do
          ! --------------------------------------------------------------------

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mkrpt_airun_en(indx0,qq,E,ncc)

      USE mkrpt_mod

      Implicit None

      Integer :: indx0,indx1
      Real(8) :: qq(Nfree),E,ncc(Nat*(Nat-1)/2)
      Integer :: i,Nf,Get_Nfree
      Real(8), dimension(:), allocatable :: qfull

         Read(Irst,*,end=100,err=100) indx1,E,ncc
         if(indx0 .ne. indx1) then 
            Write(6,*) 'Invalid line is detected in the restart file:'
            !Write(6,500) indx1,E,ncc
            Call write_en(6,indx1,E,(Nat*(Nat-1)/2),ncc)
            Write(6,'('' It must have an index '',i8)') indx0
            Stop
         endif
         Call write_en(Irst2,indx0,E,(Nat*(Nat-1)/2),ncc)
         Call myFlsh(Irst2)
         return

  100    Continue
         Nf=Get_Nfree()
         Allocate(qfull(Nf))
         qfull=0.D+00
         Do i=1,Nfree
            qfull(qmod(i))=qq(i)
         End do
         Call Run_en(qfull,E,ncc)
         Call write_en(Irst2,indx0,E,(Nat*(Nat-1)/2),ncc)
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

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine write_ed(ifl,indx,E,d)

      Implicit None

          Integer :: ifl,indx,ii
          Real(8) :: E,d(3)

          Write(ifl,500) indx,E,d
      500 Format(i8,f20.10,3f12.8)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine write_en(ifl,indx,E,ii,ncc)

      Implicit None

          Integer :: ifl,indx,ii
          Real(8) :: E,ncc(ii)

          Write(ifl,500) indx,E,ncc
      500 Format(i8,f20.10,100f12.4)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!Added by MK
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
          write(*,*) "CHODVRGrid"
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

