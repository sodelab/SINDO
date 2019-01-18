!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/12/09
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Module Vav_mod

        Implicit None

        ! -----------------------------------------------------------------------
        !     r_v  :: Vibrational averaged structure
        Logical :: r_v

        ! -- (Average of molecular properties) ---------------------------------- 
        !     MR2                 :: Mode representation
        !     NP                  :: Number of properties
        !
        !     np1                 :: Number of 1MR terms
        !     mp1(2,np1)          :: Index of 1MR terms
        !                              1: NP, 2: Mode
        !
        !     P0(NP)              :: Value at Q=0
        !     P1(maxCHO,NP,Nfree) :: 1MR terms

        Integer :: MR2,NP
        Real(8), dimension(:), allocatable   :: P0
        Real(8), dimension(:,:,:), allocatable :: P1

        !     PrptDir     :: Location of property DATA files
        !     Len_PrptDir :: Length of PrptDir
        
        Character (len=80) :: PrptDir
        Integer            :: Len_PrptDir

      CONTAINS

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine get_fname1(imod,iNP,fname)

      Implicit None

         Integer :: i,imod,iNP
         Character(*) :: fname

         fname=PrptDir(:Len_PrptDir)
         i=Len_PrptDir
         Call get_header(imod,i,fname(i+1:))
         if(iNP<10) then
            write(fname(i+1:),'(''-'',i1,''.pot'')') iNP
         elseif(iNP<100) then
            write(fname(i+1:),'(''-'',i2,''.pot'')') iNP
         endif

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

      Subroutine setupPRPT(io)

      USE Vib_mod, only: Nfree,nCHO,maxCHO

      Implicit None

         Integer :: ierr,spr_memalloc,io

         Integer   :: ifl
         Character :: fp*80,title*40

         Integer :: i,j,k,nGi,nsf
         Real(8), allocatable :: qi(:),qsf(:),psf(:),psfs(:)

         Real(8), parameter :: allow_ext=1.1D+00

       ! Initialize
         ierr=spr_memalloc(-1,dble(maxCHO*NP*Nfree*8+NP*8))
         Allocate(P0(NP),P1(maxCHO,NP,Nfree))
         P0=0.D+00
         P1=0.D+00

       ! Location of DATA files
         Call GetEnv('POTDIR',PrptDir)
         Len_PrptDir=Len_Trim(PrptDir)
         if(Len_PrptDir /=0) then
            PrptDir(Len_PrptDir+1:)='/'
            Len_PrptDir=Len_PrptDir+1
         else
            PrptDir='./'
            Len_PrptDir=2
         endif

         Call file_indicator(12,ifl)

       ! 1MR Data
         Do i=1,Nfree
          ! Disable Mode i
            if(nCHO(i)==0) cycle

            nGi=nCHO(i)
            Allocate(qi(nGi))
            Call Modal_getQ(i,nGi,qi)

            Do j=1,NP
            ! SPLINE function
               Call get_fname1(i,j,fp)
               Open(unit=ifl,file=fp,status='OLD')
               Read(ifl,'(a)') title
               Read(ifl,*)
               Read(ifl,*) nsf
               Allocate(qsf(nsf),psf(nsf),psfs(nsf))
               Read(ifl,*)
               Do k=1,nsf
                  Read(ifl,*) qsf(k),psf(k)
               End do
               Close(ifl)
               Call DSpline(qsf,psf,nsf,1.D+32,1.D+32,psfs)

               if(abs(qi(1)) < abs(qsf(1))*allow_ext .and. &
                  abs(qi(nGi)) < abs(qsf(nsf))*allow_ext ) then

                  Call DSplint(qsf,psf,psfs,nsf,0.D+00,P0(j)) 
                  Do k=1,nGi
                     Call DSplint(qsf,psf,psfs,nsf,qi(k),P1(k,j,i)) 
                     P1(k,j,i)=P1(k,j,i) - P0(j)
                  End do

               else
                  write(6,*) 'ERROR: TOO LARGE DVR GRID EXTENDING OVER THE PROPERTY SURFACE'
                  write(6,*) 'ERROR: RERUN VMAX WITH A SMALLER VALUE'
                  Stop

               endif

               Deallocate(qsf,psf,psfs)

               write(io,151) i,j,nGi,title
            End do

            Deallocate(qi)

         End do

         if(MR2==1) goto 1000

        !=========================
        !  MR=2,3,.. is not ready
        !=========================

    1000 Continue

        !dbg  write(6,'(3f8.2)') P0
        !dbg  Do i=1,Nfree
        !dbg     write(6,'(i4)') i
        !dbg     Do j=1,NP
        !dbg        write(6,'(i4,11f8.2)') j,P1(:,j,i)
        !dbg     End do
        !dbg  End do

    151 Format(8x,'o MODE=',i3,', P[',i2.2,'], GRID=',i4,3x,a40)

      End Subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      End Module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
