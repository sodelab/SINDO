!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/06/08
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Module Vcor_mod

        USE Vib_mod

        ! -- (FILE description) -------------------------------------------------
        !      o rwf  :: File indicator for 'vscf-w.wfn'

        Integer :: rwf

        ! -- (Reference VSCF & Target states) -----------------------------------
        !     ss         :: State specific VSCF reference
        !     ref(Nfree) :: Reference VSCF configuration
        !     tar(Nfree) :: Target vibrational state
        !     Nstate              :: The number of states to be calculated
        !     Label(Nstate*Nfree) :: Vibrational states to be calculated

        Logical :: ss
        Integer :: Nstate
        Integer, allocatable :: label0(:),tar(:),ref(:)

        ! -- (Virtual VSCF configurations) --------------------------------------

        !      o nLvl       :: Level of excitation
        Integer :: nLvl

        !      > P-space selection
        !      o nGen   :: Num. of interaction to generate configurations
        !      o Thresholds
        !           pth0 : Select if (E0-En) < pth0         (default 500 cm-1) 
        !           pth1 : Select if <0|H|n>/(E0-En) > pth1 (default 0.1)
        !           pth2 : Discard CI weight < pth2         (default 0.05) 
        Integer :: nGen
        Real(8) :: pth0,pth1,pth2

        !      > P-space (Resonant states)
        !      o maxpCUP                  :: Max order of mode coupling
        !      o npCUP                    :: Order of mode coupling
        !      o npCnf                    :: Number of VSCF configurations
        !      o pCnf(maxpCUP*2+1,npCnf)  :: Label  of VSCF configurations
        !
        !      > Q-space
        !      o maxqCUP                  :: Max order of mode coupling
        !      o nqCUP                    :: Order of mode coupling
        !      o nqCnf                    :: Number of VSCF configurations
        !      o qCnf(maxqCUP*2+1,nqCnf)  :: Label  of VSCF configurations
        !      o nqCnfi(maxqCUP)          :: Number of n-mode configurations

        Integer, allocatable :: pCnf(:,:),qCnf(:,:),nqCnfi(:)
        Integer :: maxpCUP,maxqCUP,npCUP,nqCUP
        Integer :: npCnf,nqCnf

        !     zp         :: true if the target is the zero-point
        Logical :: zp

        ! -- (Modals in DVR) ----------------------------------------------------
        !
        !      o type_Table1(nCHO(i),Nfree) :: Pointer to block
        !      o block :: DVR coefficients and weights (density)

        Integer, allocatable :: type_Table1(:)
        Real(8), allocatable :: block(:)

        ! -- (Kinetic energy) ---------------------------------------------------

        Real(8), allocatable :: Tmat(:)

        ! -- (Q-Matrix) ---------------------------------------------------------
        !     Qmat(4,maxCHO,maxCHO,Nfree) :: Matrix elements of Q,Q^2,Q^3,Q^4
        !                                    in terms of one-mode functions

        Real(8), dimension(:,:,:,:), allocatable :: Qmat

      Contains

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine genKinmat(mm,Tm)

      Implicit None 

         Integer :: mm,i,j,k,l,nHOm,nCHOm
         Real(8) :: Tm(nCHO(mm),nCHO(mm))
         Real(8), allocatable :: CHOm(:,:),CHOn(:,:),Thm(:,:),Cwfn(:,:)


         nCHOm=nCHO(mm)
         if(.not. CHODVR) then
            nHOm=nCHO(mm)
         else
            nHOm=nHO(mm)
         endif
         !dbg write(6,*) mm,nCHOm,nHOm

         Allocate(CHOm(nHOm,nCHOm),Thm(nHOm,nHOm))
         Call HO_kinmat(nHOm-1,omegaf(mm),Thm)
         !dbg Do i=1,nHOm
         !dbg    write(6,'(11f8.4)') Thm(:,i)
         !dbg End do
         !dbg write(6,*)

         if(.not. CHODVR) then 
            Call Modal_getCwfn(mm,CHOm)

         else
            Allocate(CHOn(nHOm,nCHOm))
            Call Modal_getCHO(mm,CHOn)
            Allocate(Cwfn(nCHOm,nCHOm))
            Call Modal_getCwfn(mm,Cwfn)

            Do i=1,nCHOm
            Do j=1,nHOm
               CHOm(j,i)=0.D+00
               Do k=1,nCHOm
                  CHOm(j,i)=CHOm(j,i) + CHOn(j,k)*Cwfn(k,i)
               End do
            End do
            End do
            !dbg Do i=1,nCHOm
            !dbg    write(6,'(11f8.4)') CHOn(:,i)
            !dbg End do
            !dbg write(6,*)
            !dbg Do i=1,nCHOm
            !dbg    write(6,'(11f8.4)') Cwfn(:,i)
            !dbg End do
            !dbg write(6,*)

            Deallocate(CHOn,Cwfn)
            
         endif

         !  Tm = CHOm*Thm*CHOm
         Do i=1,nCHOm
            Tm(i,i)=0.D+00
            Do k=1,nHOm
            Do l=1,nHOm
               Tm(i,i)=Tm(i,i) + CHOm(l,i)*Thm(l,k)*CHOm(k,i)
            End do
            End do
         End do

         Do i=1,nCHOm
         Do j=1,i-1
            Tm(j,i)=0.D+00
            Do k=1,nHOm
            Do l=1,nHOm
               Tm(j,i)=Tm(j,i) + CHOm(l,j)*Thm(l,k)*CHOm(k,i)
            End do
            End do
            Tm(i,j)=Tm(j,i)
         End do
         End do

         !dbg Do i=1,nCHOm
         !dbg    write(6,'(11f8.4)') Tm(:,i)
         !dbg End do
         !dbg write(6,*)

         Deallocate(CHOm,Thm)

      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      End module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
