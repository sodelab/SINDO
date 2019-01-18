!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/15
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Module Vscf_mod

        USE Vib_mod

        ! -- (FILE description) -------------------------------------------------
        !     rwf :: File indicator of 'vscf-r.wfn'
        !     wwf :: File indicator of 'vscf-w.wfn'

        Integer :: rwf,wwf

        ! -- (STATE parameters) ------------------------------------------------- 
        !     Nstate              :: The number of states to be calculated
        !     Nst                 :: Index of current vibrational state
        !     Label(Nstate*Nfree) :: Vibrational states to be calculated

        Integer :: Nstate,Nst
        Integer, dimension(:), allocatable:: label

        ! -- (VSCF parameters) -------------------------------------------------- 
        !     restart :: Restart option
        !     Maxitr  :: Max cycle of VSCF iteration
        !     Ethresh :: Threshold energy for convergence criteria / in cm-1
        !     Etot    :: Total energy

        Logical :: restart
        Integer :: Maxitr
        Real(8) :: Ethresh, Etot
        Double Precision, parameter :: pi=3.14159265358979323846
        ! -- (VSCF one-mode variables) ------------------------------------------
        !     Tmat(nCHO(i),nCHO(i),Nfree):: Kinetic energy matrix
        !     type_Table1(nCHO(i),Nfree) :: Pointer to block
        !     block                      :: Mean-Field potential and DVR weights

        Integer, dimension(:), allocatable :: type_Table1
        Real(8), dimension(:), allocatable :: Tmat,block


      CONTAINS 

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine setKinmat()

      Implicit None 

         Integer :: i
         Real(8), allocatable :: Tm(:)

         Do i=1,Nfree
            Allocate(Tm(nCHO(i)*nCHO(i)))
            Call genKinmat(i,Tm)
            Tmat(idx2(i)+1:idx2(i)+nCHO(i)*nCHO(i))=Tm

            Deallocate(Tm)
         End do

      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine genKinmat(mm,Tm)

      Implicit None 

         Integer :: mm,i,j,k,l
         Real(8) :: Tm(nCHO(mm),nCHO(mm))
         Real(8), allocatable :: CHOm(:,:),Thm(:,:)

         if(.not. CHODVR) then 
            Call HO_kinmat(nCHO(mm)-1,omegaf(mm),Tm)

         else
            Allocate(CHOm(nHO(mm),nCHO(mm)),Thm(nHO(mm),nHO(mm)))
            Call HO_kinmat(nHO(mm)-1,omegaf(mm),Thm)
            Call modal_getCHO(mm,CHOm)

            !  Tm = CHOm*Thm*CHOm
            Do i=1,nCHO(mm)
               Tm(i,i)=0.D+00
               Do k=1,nHO(mm)
               Do l=1,nHO(mm)
                  Tm(i,i)=Tm(i,i) + CHOm(l,i)*Thm(l,k)*CHOm(k,i)
               End do
               End do
            End do

            Do i=1,nCHO(mm)
            Do j=1,i-1
               Tm(j,i)=0.D+00
               Do k=1,nHO(mm)
               Do l=1,nHO(mm)
                  Tm(j,i)=Tm(j,i) + CHOm(l,j)*Thm(l,k)*CHOm(k,i)
               End do
               End do
               Tm(i,j)=Tm(j,i)
            End do
            End do

            Deallocate(CHOm,Thm)

         endif
         !write(6,*) mm
         !Do i=1,nCHO(mm)
         !   write(6,'(11f8.4)') Tm(:,i)
         !End do

      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      End Module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
