!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/07/19
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

Module Vib_mod

        ! -- (SYSTEM parameters) ------------------------------------------------
        !     Nfree         :: Number of vibrational degrees of freedom
        !     omegaf(Nfree) :: Harmonic frequencies (cm-1)

        Integer :: Nfree
        Real(8), dimension(:), allocatable :: omegaf

        ! -- (RUN OPTION paramters) --------------------------------------------- 
        !     lvscf :: Run VSCF if true. (default=.t.)
        !     lvci  :: Run VCI if true. (default=.f.)
        !     lvcc  :: Run VCC if true. (default=.f.)
        !     lvpt  :: Run VPT if true. (default=.f.)
        !     ltdh  :: Run TDH if true. (default=.f.)
        !     ltdv  :: Run TDV if true. (default=.f.)
        !     lvav  :: Run Vibrational Average if true. (default=.f.)

        Logical :: lvscf,lvci,lvcc,lvpt,ltdh,ltdv,lvav
        Logical :: summary !MK to dump info on fort.666

        ! -- (POTENTIAL parameters) --------------------------------------------- 
        !     MR         :: Mode Coupling Representation
        !     PotDir     :: Location of PEF data files
        !     Len_PotDir :: Length of PotDir
        
        Integer :: MR
        Character (len=80) :: PotDir
        Integer            :: Len_PotDir

        !     o  QFF
        !        nQx     :: Num. of xMR terms
        !        mQx     :: Mode of each terms
        !        coeffx  :: Coefficients
        !     o  Spline
        !        nSx     :: Num. of xMR terms
        !        mSx     :: Mode of each terms

        Integer :: nQ1,nQ2,nQ3, nS1,nS2,nS3
        Character :: tQ1*80,tQ2*80,tQ3*80
        Integer, dimension(:), allocatable :: mQ1, mS1
        Integer, dimension(:,:), allocatable :: mQ2,mQ3, mS2,mS3
        Real(8), dimension(:,:), allocatable :: coeff1,coeff2,coeff3

        ! -- (BASIS FUNCTION parameters) ----------------------------------------
        !     CHODVR      :: Switch to contracted HODVR 
        !     nHO(Nfree)  :: Num. of HO wfn 
        !     nCHO(Nfree) :: Num. of contracted HO wfn 
        !     maxCHO      :: Maximum of nCHO

        Logical :: CHODVR
        Integer :: maxCHO
        Integer, dimension(:), allocatable :: nHO,nCHO

        ! -- (POINTER) ---------------------------------------------------------
        !     idx1(Nfree) :: Pointer to A(nCHO(i),Nfree)
        !     idx2(Nfree) :: Pointer to A(nCHO(i),nCHO(i),Nfree)

        Integer, dimension(:), allocatable :: idx1,idx2

      CONTAINS

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine get_fname1(imod,fname)

      Implicit None

         Integer :: i,imod
         Character(*) :: fname

         fname=PotDir(:Len_PotDir)
         i=Len_PotDir
         Call get_header(imod,i,fname(i+1:))
         fname(i+1:)='.pot'

      End Subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      ! imod>jmod
      Subroutine get_fname2(imod,jmod,fname)

      Implicit None

         Integer :: imod,jmod,i,j
         Character(*) :: fname

         fname=PotDir(:Len_PotDir)
         i=Len_PotDir
         Call get_header(imod,i,fname(i+1:))
         Call get_header(jmod,i,fname(i+1:))
         fname(i+1:)='.pot'

      End Subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      ! i>j>k
      Subroutine get_fname3(imod,jmod,kmod,fname)

      Implicit None

         Integer :: i,j,k,imod,jmod,kmod
         Character(*) :: fname

         fname=PotDir(:Len_PotDir)
         i=Len_PotDir
         Call get_header(imod,i,fname(i+1:))
         Call get_header(jmod,i,fname(i+1:))
         Call get_header(kmod,i,fname(i+1:))
         fname(i+1:)='.pot'

      End Subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine get_header(imod,ipt,fp)

      Implicit None

         Integer :: imod,ipt
         Character(*) :: fp
         Character :: cm1*1,cm2*2,cm3*3

         if(imod<10) then
            write(cm1,'(i1)') imod
            fp='q'//cm1
            ipt=ipt+2
         elseif(imod<100) then
            write(cm2,'(i2)') imod
            fp='q'//cm2
            ipt=ipt+3
         else
            write(cm3,'(i3)') imod
            fp='q'//cm3
            ipt=ipt+4
         endif

      End Subroutine


      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine setupPEF(ierr)

      Implicit None

         Integer :: ierr

         Integer :: i,j,k,l,m,n

         Integer   :: ifl
         Character :: fp*80

         Integer :: spr_memalloc
         Integer, parameter :: nSmax=1000000,nQmax=1000000
         Integer :: modeS(3,nSmax),modeQ(3,nQmax)
         Real(8) :: coeff(6), &
                    qff_GetGi, qff_GetHii, qff_GetTiii, qff_GetUiiii, &
                    qff_GetHij, qff_GetTiij, qff_GetUiijj, qff_GetUiiij, &
                    qff_GetTijk, qff_GetUiijk

         Logical :: qff

       ! Initialize
         nS1=0; nQ1=0
         nS2=0; nQ2=0
         nS3=0; nQ3=0

       ! Location of DATA files
         Call GetEnv('POTDIR',PotDir)
         Len_PotDir=Len_Trim(PotDir)
         if(Len_PotDir /=0) then
            PotDir(Len_PotDir+1:)='/'
            Len_PotDir=Len_PotDir+1
         else
            PotDir='./'
            Len_PotDir=2
         endif

       ! Check if 001.hs exists
         qff=.false.
         Open(unit=7,file=PotDir(:Len_PotDir)//'001.hs',status='OLD',err=10)
         Close(7)

       ! Read QFF
         qff=.true.
         Call qff_Construct(ierr)
         if(ierr<0) return
         Call qff_GetTitle(tQ1,tQ2,tQ3)

      10 Continue

       ! 1MR-PEFs
         coeff=0.D+00
         Do i=1,Nfree
          ! Disable Mode i
            if(nCHO(i)==0) cycle

          ! SPLINE function
            Call get_fname1(i,fp)
            Call file_indicator(12,ifl)
            Open(unit=ifl,file=fp,status='OLD',err=100)
            nS1=nS1+1
            modeS(1,nS1)=i
            Close(ifl)
            cycle

            100 Continue

          ! QFF
            if(qff) then
               coeff(1)=qff_GetGi(i) 
               coeff(2)=qff_GetHii(i) 
               coeff(3)=qff_GetTiii(i) 
               coeff(4)=qff_GetUiiii(i) 
               Do j=1,4
                  if(abs(coeff(j))>1.D-06) then
                     nQ1=nQ1+1
                     modeQ(1,nQ1)=i
                     exit
                  endif
               End do
               cycle
             endif

          ! Disable current Mode
            nCHO(i)=0

         End do

         if(nS1/=0) then
            ierr=spr_memalloc(1,dble(nS1*4))
            allocate(mS1(nS1))
            mS1=modeS(1,1:nS1)
         endif

         if(nQ1/=0) then
            ierr=spr_memalloc(1,dble(nQ1)*4.D+00+dble(nQ1)*4.D+00*8.D+00)
            allocate(mQ1(nQ1),coeff1(4,nQ1))
            mQ1=modeQ(1,1:nQ1)
            Do i=1,nQ1
               coeff1(1,i)=qff_GetGi(mQ1(i))
               coeff1(2,i)=qff_GetHii(mQ1(i))
               coeff1(3,i)=qff_GetTiii(mQ1(i))
               coeff1(4,i)=qff_GetUiiii(mQ1(i))
            End do
         endif
         if(MR==1) goto 1000


       ! 2MR-PEFs
         coeff=0.D+00
         Do i=1,Nfree
          ! Disable Mode i
            if(nCHO(i)==0) cycle

         Do j=1,i-1
          ! Disable Mode j
            if(nCHO(j)==0) cycle

          ! SPLINE function
            Call get_fname2(i,j,fp)
            Call file_indicator(12,ifl)
            Open(unit=ifl,file=fp,status='OLD',err=200)
            nS2=nS2+1
            modeS(1,nS2)=i
            modeS(2,nS2)=j
            Close(ifl)
            cycle

            200 Continue

          ! QFF
            if(qff) then
               k=(i-1)*(i-2)/2 + j
               coeff(1)=qff_GetHij(k)
               coeff(2)=qff_GetTiij(k,0)
               coeff(3)=qff_GetTiij(k,1)
               coeff(4)=qff_GetUiijj(k)
               coeff(5)=qff_GetUiiij(k,0)
               coeff(6)=qff_GetUiiij(k,1)
               Do k=1,6
                  if(abs(coeff(k))>1.D-06) then
                     nQ2=nQ2+1
                     modeQ(1,nQ2)=i
                     modeQ(2,nQ2)=j
                     exit
                  endif
               End do
               cycle
            endif

         End do
         End do

         if(nS2/=0) then
            ierr=spr_memalloc(1,dble(nS2)*2.D+00*4.D+00)
            allocate(mS2(2,nS2))
            mS2(1,:)=modeS(1,1:nS2)
            mS2(2,:)=modeS(2,1:nS2)
         endif

         if(nQ2/=0) then
            ierr=spr_memalloc(1,dble(nQ2)*2.D+00*4.D+00+dble(nQ2)*6.D+00*8.D+00)
            allocate(mQ2(2,nQ2),coeff2(6,nQ2))
            mQ2(1,:)=modeQ(1,1:nQ2)
            mQ2(2,:)=modeQ(2,1:nQ2)
            Do n=1,nQ2
               i=mQ2(1,n)
               j=mQ2(2,n)
               k=(i-1)*(i-2)/2 + j
               coeff2(1,n)=qff_GetHij(k)
               coeff2(2,n)=qff_GetTiij(k,0)
               coeff2(3,n)=qff_GetTiij(k,1)
               coeff2(4,n)=qff_GetUiijj(k)
               coeff2(5,n)=qff_GetUiiij(k,0)
               coeff2(6,n)=qff_GetUiiij(k,1)
            End do
         endif
         if(MR==2) goto 1000


       ! 3MR-PEFs
         coeff=0.D+00
         Do i=3,Nfree
          ! Disable Mode i
            if(nCHO(i)==0) cycle

         Do j=2,i-1
          ! Disable Mode j
            if(nCHO(j)==0) cycle

         Do k=1,j-1
          ! Disable Mode k
            if(nCHO(k)==0) cycle


          ! SPLINE function
            Call get_fname3(i,j,k,fp)
            Call file_indicator(12,ifl)
            Open(unit=ifl,file=fp,status='OLD',err=300)
            nS3=nS3+1
            modeS(1,nS3)=i
            modeS(2,nS3)=j
            modeS(3,nS3)=k
            Close(ifl)
            cycle

            300 Continue

          ! QFF
            if(qff) then
               l=(i-1)*(i-2)*(i-3)/6 + (j-1)*(j-2)/2 + k
               coeff(1)=qff_GetTijk(l)
               coeff(2)=qff_GetUiijk(l,0)
               coeff(3)=qff_GetUiijk(l,1)
               coeff(4)=qff_GetUiijk(l,2)
               Do l=1,4
                  if(abs(coeff(l))>1.D-06) then
                     nQ3=nQ3+1
                     modeQ(1,nQ3)=i
                     modeQ(2,nQ3)=j
                     modeQ(3,nQ3)=k
                     exit
                  endif
               End do
               cycle
            endif

         End do
         End do
         End do

         if(nS3/=0) then
            ierr=spr_memalloc(1,dble(nS3)*3.D+00*4.D+00)
            allocate(mS3(3,nS3))
            mS3(1,:)=modeS(1,1:nS3)
            mS3(2,:)=modeS(2,1:nS3)
            mS3(3,:)=modeS(3,1:nS3)
         endif

         if(nQ3/=0) then
            ierr=spr_memalloc(1,dble(nQ3)*3.D+00*4.D+00+dble(nQ3)*4.D+00*8.D+00)
            allocate(mQ3(3,nQ3),coeff3(4,nQ3))
            mQ3(1,:)=modeQ(1,1:nQ3)
            mQ3(2,:)=modeQ(2,1:nQ3)
            mQ3(3,:)=modeQ(3,1:nQ3)
            Do n=1,nQ3
               i=mQ3(1,n)
               j=mQ3(2,n)
               k=mQ3(3,n)
               l=(i-1)*(i-2)*(i-3)/6 + (j-1)*(j-2)/2 + k
               coeff3(1,n)=qff_GetTijk(l)
               coeff3(2,n)=qff_GetUiijk(l,0)
               coeff3(3,n)=qff_GetUiijk(l,1)
               coeff3(4,n)=qff_GetUiijk(l,2)
            End do
         endif
         if(MR==3) goto 1000

 1000    Continue
         if(qff) Call qff_Destruct


      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine genGrid(ierr)

      Implicit None

         Integer :: ierr
         Integer :: spr_memalloc

         Integer :: i,j,k,l, n
         Integer :: nGi,nGj,nGk
         Real(8) :: b1si,b2si,b1sj,b2sj,b1sk,b2sk
         Real(8), dimension(:), allocatable :: qsi,qsj,qsk

         ! ONLY for debug
         Logical, parameter :: uniform_grid=.false.

         ! setup CHO
         Call Modal_Construct(ierr)
         if(ierr<0) return

         ! setup nGridx and Vx
         Call PEGrid_Construct(ierr)
         if(ierr<0) return

         ! Determine optimal nCHO
         Do n=1,nS1
            Call Spfn1_Gen1(0,n)
            Call Spfn1_getnGrid(nGi)
            Allocate(qsi(nGi))
            Call Spfn1_getGrid(qsi)
            b1si=abs(qsi(1))
            b2si=abs(qsi(nGi))
            Deallocate(qsi)
            Call Spfn1_Dispose()

            nGi=nCHO(mS1(n))
            Call optGrid(mS1(n),nGi,b1si,b2si)
            if(nGi /= nCHO(mS1(n))) then
               write(6,100) mS1(n),nCHO(mS1(n))-1,nGi-1
               write(6,*)
               nCHO(mS1(n)) = nGi
            endif
            Call PEGrid_setnGrid1(nGi,n)

         End do
     100 Format(9x,'WARNING:  MISMATCH OF PEF-GRID AND DVR-GRID IS DETECTED',/, &
                9x,'WARNING:  RUNNING WITH SMALLER BASIS SETS FOR',/, &
                9x,'WARNING:        MODE=',i4,/, &
                9x,'WARNING:        VMAX=',i4,' ->',i4)

         Do n=1,nS2
            Call Spfn2_Gen1(0,n)
            Call Spfn2_getnGrid(nGi,nGj)
            Allocate(qsi(nGi),qsj(nGj))
            Call Spfn2_getGrid(qsi,qsj)
            b1si=abs(qsi(1))
            b2si=abs(qsi(nGi))
            b1sj=abs(qsj(1))
            b2sj=abs(qsj(nGj))
            Deallocate(qsi,qsj)
            Call Spfn2_Dispose()

            if(nGi > nCHO(mS2(1,n))) nGi=nCHO(mS2(1,n))
            if(nGj > nCHO(mS2(2,n))) nGj=nCHO(mS2(2,n))

            Call optGrid(mS2(1,n),nGi,b1si,b2si)
            Call optGrid(mS2(2,n),nGj,b1sj,b2sj)
            Call PEGrid_setnGrid2(nGi,nGj,n)

         End do

         Do n=1,nS3
            Call Spfn3_Gen1(0,n)
            Call Spfn3_getnGrid(nGi,nGj,nGk)
            Allocate(qsi(nGi),qsj(nGj),qsk(nGk))
            Call Spfn3_getGrid(qsi,qsj,qsk)
            b1si=abs(qsi(1))
            b2si=abs(qsi(nGi))
            b1sj=abs(qsj(1))
            b2sj=abs(qsj(nGj))
            b1sk=abs(qsk(1))
            b2sk=abs(qsk(nGk))
            Deallocate(qsi,qsj,qsk)
            Call Spfn3_Dispose()

            if(nGi > nCHO(mS3(1,n))) nGi=nCHO(mS3(1,n))
            if(nGj > nCHO(mS3(2,n))) nGj=nCHO(mS3(2,n))
            if(nGk > nCHO(mS3(3,n))) nGk=nCHO(mS3(3,n))

            Call optGrid(mS3(1,n),nGi,b1si,b2si)
            Call optGrid(mS3(2,n),nGj,b1sj,b2sj)
            Call optGrid(mS3(3,n),nGk,b1sk,b2sk)
            Call PEGrid_setnGrid3(nGi,nGj,nGk,n)

         End do

         if(uniform_grid) then
            Write(6,'(9x,''RUN WITH THE UNIFORM_GRID OPTION'')')
            Do n=1,nS1
               i=mS1(n)
               Call PEGrid_setnGrid1(nCHO(i),n)
            End do
            Do n=1,nS2
               i=mS2(1,n); j=mS2(2,n)
               Call PEGrid_setnGrid2(nCHO(i),nCHO(j),n)
            End do
            Do n=1,nS3
               i=mS3(1,n); j=mS3(2,n); k=mS3(3,n)
               Call PEGrid_setnGrid3(nCHO(i),nCHO(j),nCHO(k),n)
            End do
         endif

       ! Setup idx1,idx2
         ierr=spr_memalloc(1,dble(Nfree)*8.D+00)
         if(ierr<0) return
         allocate(idx1(Nfree),idx2(Nfree))
         idx1(1)=0
         idx2(1)=0
         Do i=2,Nfree
            idx1(i)=idx1(i-1)+nCHO(i-1)
            idx2(i)=idx2(i-1)+nCHO(i-1)*nCHO(i-1)
         End do

         Call Modal_open(ierr)
         if(ierr<0) return

         Do i=1,Nfree
            Call Modal_add(i,nCHO(i))
         End do
         Do n=1,nS2
            i=mS2(1,n)
            j=mS2(2,n)
            Call PEGrid_getnGrid2(nGi,nGj,n)
            Call Modal_add(i,nGi)
            Call Modal_add(j,nGj)
         End do
         Do n=1,nS3
            i=mS3(1,n)
            j=mS3(2,n)
            k=mS3(3,n)
            Call PEGrid_getnGrid3(nGi,nGj,nGk,n)
            Call Modal_add(i,nGi)
            Call Modal_add(j,nGj)
            Call Modal_add(k,nGk)
         End do
         Call Modal_Close(ierr)
         if(ierr<0) return

         Call Modal_genDVR_ALL()
         Call PEGrid_onGrid(ierr)
         if(ierr<0) return

      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine optGrid(mode,nG,b1s,b2s)

      Implicit None 

      Integer :: mode,nG,nGm
      Real(8) :: b1s,b2s,b1m,b2m
      Real(8), parameter :: allow_extrapolation=1.10D+00
      Real(8) :: qm(nG),xdvr(nG,nG)

         b1s=b1s*allow_extrapolation
         b2s=b2s*allow_extrapolation
         nGm=nG

      10 Continue
         Call Modal_genDVR(mode,nGm,xdvr,qm)
         b1m=abs(qm(1))
         b2m=abs(qm(nGm))
         !dbg write(6,'(2i4)') mode,nGm
         !dbg write(6,'(2f8.4)') b1s,b1m
         !dbg write(6,'(2f8.4)') b2s,b2m
         if(b1s < b1m .or. b2s < b2m) then
            nGm=nGm-1
            goto 10
         endif

         nG=nGm

      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      End Module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
