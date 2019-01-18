!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/08
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Module PEGrid_mod

        USE Vib_mod, only: nS1,nS2,nS3

        ! 1-mode terms
        !  nGrid1(Nfree) :: Number of grid points
        !    ptV1(Nfree) :: Pointer of V1
        !      V1(Nfree) :: Potential energy values on Grid
        Integer, dimension(:), allocatable :: nGrid1,ptV1
        Real(8), dimension(:), allocatable :: V1

        ! 2-mode terms
        !  nGrid2(2,nS2) :: Number of grid points
        !    ptV2(nS2)   :: Pointer of V2
        !      V2(nS2)   :: Potential energy values on Grid
        Integer, dimension(:,:), allocatable :: nGrid2
        Integer, dimension(:), allocatable :: ptV2
        Real(8), dimension(:), allocatable :: V2

        ! 3-mode terms
        !  nGrid3(3,nS3) :: Number of grid points
        !    ptV3(nS3) :: Pointer of V3
        !      V3(nS3) :: Potential energy values on Grid
        Integer, dimension(:,:), allocatable :: nGrid3
        Integer, dimension(:), allocatable :: ptV3
        Real(8), dimension(:), allocatable :: V3

      CONTAINS

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine setV1(n,nGi,Vi)

      USE Vib_mod, only: mS1

      Implicit None

         Integer :: n,i,nGi
         Real(8) :: Vi(nGi),qi(nGi),Spfn1

         Call Modal_getQ(mS1(n),nGi,qi)
         Call Spfn1_Gen1(1,n)
         Do i=1,nGi
            Vi(i)=Spfn1(qi(i))
         End do
         Call Spfn1_Dispose()

         !dbg Write(6,*) mS1(n)
         !dbg Write(6,'(2f8.4)') (qi(i),Vi(i),i=1,nGi)

      End subroutine
 
      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine setV2(n,nGi,nGj,Vij)

      USE Vib_mod, only: mS2

      Implicit None

         Integer :: n,i,j,nGi,nGj
         Real(8) :: Vij(nGj,nGi),qi(nGi),qj(nGj),Spfn2

         Call Modal_getQ(mS2(1,n),nGi,qi)
         Call Modal_getQ(mS2(2,n),nGj,qj)
         Call Spfn2_Gen1(1,n)
         Do i=1,nGi
         Do j=1,nGj
            Vij(j,i)=Spfn2(qi(i),qj(j))
         End do
         End do
         Call Spfn2_Dispose()

         !dbg Write(6,*) mS2(:,n)
         !dbg Write(6,'(3f8.4)') ((qj(j),qi(i),Vij(j,i),j=1,nGi),i=1,nGi)

      End subroutine

 
      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine setV3(n,nGi,nGj,nGk,Vijk)

      USE Vib_mod, only: mS3

      Implicit None

         Integer :: n,i,j,k,nGi,nGj,nGk
         Real(8) :: Vijk(nGk,nGj,nGi),qi(nGi),qj(nGj),qk(nGk),Spfn3

         Call Modal_getQ(mS3(1,n),nGi,qi)
         Call Modal_getQ(mS3(2,n),nGj,qj)
         Call Modal_getQ(mS3(3,n),nGk,qk)
         Call Spfn3_Gen1(1,n)
         Do i=1,nGi
         Do j=1,nGj
         Do k=1,nGk
            Vijk(k,j,i)=Spfn3(qi(i),qj(j),qk(k))
         End do
         End do
         End do
         Call Spfn3_Dispose()

         !dbg Write(6,*) mS3(:,n)
         !dbg Write(6,'(4f8.4)') (((qk(k),qj(j),qi(i),Vijk(k,j,i),k=1,nGk),j=1,nGi),i=1,nGi)

      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
 
      End module
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
