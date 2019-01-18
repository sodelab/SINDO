!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Module v2ai_mod

        ! -- (Information of the molecule) -------------------------------------
        Integer :: Nat,Mfree
        Double Precision, dimension(:,:), allocatable :: xeq,CL
        Double Precision, dimension(:), allocatable :: Zmass,Omega
        Logical :: Readmol
        Character :: path*40,Typ*3
        Character(2), dimension(:), allocatable :: Lbl

        ! -- (Ab initio run options) -------------------------------------------
        Integer, parameter :: nline=100
        Character :: Method*10,Com(nline)*120,Com2(nline)*120,RunTyp*3
        Logical :: dpl,nmrcc,dryrun
        Real(8) :: minlimit,maxlimit  !MK interpolate if energy exceeds limit (inpl())

        ! -- (Information of normal modes) -------------------------------------
        Integer :: Nfree,Isym(3)
        Character :: Sym*3
        Integer :: Nirp
        Integer, dimension(:), allocatable :: Msym,Nsym
        Integer, dimension(:,:), allocatable :: Lsym
        Integer, dimension(:), allocatable :: ActvMode,inActvMode

        ! -- (Character table of D2h) ------------------------------------------
        Integer :: ich(8,7)

        !               E C2z C2y C2x i xy xz yz
        !    B1g
        Data ich(:,1) / 1, 1, -1, -1, 1, 1,-1,-1 /
        !    B2g
        Data ich(:,2) / 1,-1,  1, -1, 1,-1, 1,-1 /
        !    B3g
        Data ich(:,3) / 1,-1, -1,  1, 1,-1,-1, 1 /
        !    Au
        Data ich(:,4) / 1, 1,  1,  1,-1,-1,-1,-1 /
        !    B1u
        Data ich(:,5) / 1, 1, -1, -1,-1,-1, 1, 1 /
        !    B2u
        Data ich(:,6) / 1,-1,  1, -1,-1, 1,-1, 1 /
        !    B3u
        Data ich(:,7) / 1,-1, -1,  1,-1, 1, 1,-1 /

      End Module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
