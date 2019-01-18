!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/07/19
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
MODULE qff_mod
!
!  Private variables for QFF
!
!    Nfree: The number of internal degrees of freedom, i.e. 3N-6
!
  Integer :: Nfree,Iout,MR
!
!    E0                : Energy
!    EC(Nfree)         : Geometry
!    Gi,Hii,Tiii,Uiiii    : 1MR
!    Hij,Tiij,Uiijj,Uiiij : 2MR
!    Tijk,Uiijk           : 3MR
!
  Double Precision:: E0
  Double Precision, dimension(:), Allocatable:: EC
  Double Precision, dimension(:), Allocatable:: Gi,Hii,Tiii,Uiiii
  Double Precision, dimension(:), Allocatable:: Hij,Tiij,Uiijj,Uiiij
  Double Precision, dimension(:), Allocatable:: Tijk,Uiijk
!
!    tl1,tl2,tl3  : Title of 1MR, 2MR, and 3MR QFF data
!
  Character :: tl1*80,tl2*80,tl3*80
!
END MODULE 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
