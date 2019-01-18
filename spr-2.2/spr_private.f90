!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!   Last modified  2007/05/09
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! SUPER MODULE : Fundamental paramters are given by this module.
!
! Information of the molecule 
!        Nat  : The number of atoms
!        Nat3 : The number of atoms * 3
!        Nfree: The number of internal degrees of freedom, i.e. Nat3-6(or Nat3-5)
!        Zmass: The mass of each atom
!        xin  : Input structure (Angs)
!        Freq : Harmonic Frequencies at xin (cm-1)
!        CL   : Normal displacement vector at xin
!
! File assignment
!        Inp  : Input file
!        Iout : Output file
!
! System parameter
!       rMaxmem : Maximum size of memory (byte)
!          rmem : Current size of memory (byte)
!          smem : Max memory used so far (byte)
!
MODULE spr_private
!
!
  Integer :: Nat,Nat3,Nfree
  Double Precision, dimension(:), allocatable :: Zmass,xin,Freq
  Double Precision, dimension(:,:), allocatable :: CL
!
  Integer :: Inp=5
  Integer :: Iout=6
!
  Real(8) :: rMaxmem,rmem,smem
!
END MODULE

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
