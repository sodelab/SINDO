!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/15
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Module Modal_mod 

        USE Vib_mod, ONLY: Nfree,nHO,nCHO,CHODVR,idx1,idx2,omegaf

        ! -- (Type of integration) ---------------------------------------------
        !  Ntype  :: Number of types of modal
        !  type_Table1(nCHO,Nfree) :: Pointer to block of the modal  
        !  type_Table2(3,Ntype)    :: 1: mode, 2: nGrid, 3: pointer
        Integer :: Ntype
        Integer, allocatable :: type_Table1(:), type_Table2(:,:)

        ! -- (DVR basis and Grids) ---------------------------------------------
        !  memsz  :: size of block
        !  block  :: xdvr and qq 
        Integer :: memsz
        Real(8), dimension(:), allocatable :: block

        ! -- (Modal coefficients and energies) ---------------------------------
        !  Ene    :: modal energy
        !  Cwfn0  :: Coefficients
        Real(8), dimension(:), allocatable :: Cwfn0,Ene

        ! -- (Contraction coefficients) -----------------------------------------
        !  CHO(nHO,nCHO,Nfree) :: Contracted HO wfn
        !  ptCHO(Nfree):: Pointer of CHO
        Integer, dimension(:), allocatable :: ptCHO
        Real(8), dimension(:), allocatable :: CHO

      End module
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
