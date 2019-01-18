!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
!   Code description by K.Yagi:  yagi@qcl.t.u-tokyo.ac.jp
!
!   Last modified : 2007/11/15
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  Module RPOLYMER

    Integer :: Nat, Nat3, Nfree
    Integer :: Inp=100
    Logical :: linear
    Logical :: polymer

  End Module RPOLYMER
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RPOLYMER_Const(ch,N)

    USE RPOLYMER

    Implicit None

    Integer, intent(in):: N
    Character(len=*),intent(in)  :: ch

    Logical op

!chk      Write(6,'(A)') ch
      Do while(0==0) 
         Inquire(Inp,opened=op)
         if(op) then
            Inp=Inp+1
         else
            exit
         endif
      End do
      Open(Inp,file=ch,status='OLD')

      Nat=N
      Nat3=Nat*3
      linear=.false.
      polymer=.false.
      Nfree=Nat3-6

    return

  END SUBROUTINE RPOLYMER_Const
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RPOLYMER_Dest()

    USE RPOLYMER

    Implicit None

      Close(Inp)

    return

  END SUBROUTINE RPOLYMER_Dest
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RPOLYMER_Label(x)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i
    Character(2), dimension(Nat) :: x
    Character :: ch*120

      Rewind(Inp)

      Read(Inp,*)
      Read(Inp,*)

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(1:2),*) x(i)
      End do

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RPOLYMER_Label
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RPOLYMER_Freq(Omega)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Double Precision, dimension(Nfree), intent(out) :: Omega

    Integer :: i,j,ii
    Logical :: op
    Character :: ch*120

      ii=200
      Do while(0==0) 
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do

      Open(ii,file='ncaout',status='old',err=10)

      Rewind(ii)

   20 Read(ii,'(A)') ch
      i=Index(ch,'vibrational frequencies')
      if(i.EQ.0) GOTO 20

      Do i=1,Nfree
         Read(ii,*) Omega(i)
         Do j=1,Nat
            Read(ii,*)
         End do
      End do
      Close(ii)

   return

! -----------------------------------------------------------------------------

 10 Continue
    Write(6,*)
    Write(6,*)'ERROR> ncaout is not found!'
    Write(6,*)'ERROR> Program ends with an error.'
    Stop


  END SUBROUTINE RPOLYMER_Freq
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RPOLYMER_Geom(x)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j
    Double Precision :: B2A
    Double Precision, dimension(3,Nat):: x
    Character :: ch*120

      Rewind(Inp)

      Read(Inp,'(a)') ch
      Read(Inp,'(a)') ch

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(3:62),'(3f20.15)') (x(j,i),j=1,3)
      End do

      B2A = 0.52917724924d+00
      x=x*B2A

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RPOLYMER_Geom
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RPOLYMER_L(CL)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j,k,ii
    Logical :: op
    Double Precision, dimension(Nat3,Nfree):: CL
    Character :: ch*120

      ii=200
      Do while(0==0) 
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do

      Open(ii,file='ncaout',status='old',err=10)

      Rewind(ii)

   20 Read(ii,'(A)') ch
      i=Index(ch,'vibrational frequencies')
      if(i.EQ.0) GOTO 20

      Do i=1,Nfree
         Read(ii,*)
         Do j=1,Nat
            Read(ii,*) k,CL(3*(j-1)+1:3*j,i)
         End do
      End do
      Close(ii)

    return

! -----------------------------------------------------------------------------

 10 Continue
    Write(6,*)
    Write(6,*)'ERROR> ncaout is not found!'
    Write(6,*)'ERROR> Program ends with an error.'
    Stop

  END SUBROUTINE RPOLYMER_L
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RPOLYMER_Mass(x)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i
    Double Precision, dimension(Nat):: x
    Character(2), dimension(Nat):: labl
    Character :: ch*120
    Double Precision, dimension(36):: Ms_DATA

    Data Ms_DATA &
    !H(1)        He(2)      Li(3)     Be(4)     B(5)       C(6)
    /1.00782504, 4.0026032, 7.016003, 9.012182, 11.009305, 12.0, &
    !N(7)        O(8)        F(9)        Ne(10)    Na(11)
     14.0030740, 15.9949146, 18.9984032, 19.99244, 22.989768, & 
    !Mg(12)     Al(13)     Si(14)     P(15)      S(16)
     23.985042, 26.981539, 27.976927, 30.973762, 31.972071, &
    !Cl(17)      Ar(18)    
     34.9688527, 39.962384 , &
    !K(19)      Ca(20)     Sc(21)
     38.963707, 39.962591, 44.955910, &
    !Ti(22)     V(23)     Cr(24)    Mn(25)    Fe(26)    Co(27)
     47.947947, 50.94396, 51.94051, 54.93805, 55.93494, 58.93320, &
    !Ni(28)    Cu(29)    Zn(30)    Ga(31)    Ge(32)    As(33)    
     57.93535, 62.92960, 63.92914, 68.92558, 73.92118, 74.92159, &
    !Se(34)    Br(35)    Kr(36)    
     79.91652, 78.91834, 83.91151 /
 
!     Rb        Sr        Y
!     84.91179, 87.90562, 88.90585, &
!    !Zr        Nb        Mo        Ru         Rh         Pd
!     89.90470, 92.90638, 97.90541, 101.90435, 102.90550, 105.90348, &
!    !Ag         Cd         In         Sn         Sb         Te
!     106.90509, 113.90336, 114.90388, 119.90220, 120.90382, 129.90623, &
!    !I          Xe         Cs         Ba         La         Ce
!     126.90447, 131.90414, 132.90543, 137.90523, 138.90635, 139.90543, &
!    !Pr         Nd         Sm         Eu         Gd         Tb
!     140.90765, 141.90772, 151.91973, 152.92123, 157.92401, 158.92534,
!    !D()      C-13()    N-15()   O-17()   O-18()
!     2.0141018, 13.0033548, 15.001090, 16.999131, 17.999160
     

      Call RPOLYMER_Label(labl)

      Do i=1,Nat

         Select Case(trim(labl(i)))
            Case('H')
               x(i)=Ms_DATA(1)
            Case('HE')
               x(i)=Ms_DATA(2)
            Case('Li')
               x(i)=Ms_DATA(3)
            Case('BE')
               x(i)=Ms_DATA(4)
            Case('B')
               x(i)=Ms_DATA(5)
            Case('C')
               x(i)=Ms_DATA(6)
            Case('N')
               x(i)=Ms_DATA(7)
            Case('O')
               x(i)=Ms_DATA(8)
            Case('F')
               x(i)=Ms_DATA(9)
            Case('NE')
               x(i)=Ms_DATA(10)
            Case('NA')
               x(i)=Ms_DATA(11)
            Case('MA')
               x(i)=Ms_DATA(12)
            Case('AL')
               x(i)=Ms_DATA(13)
            Case('SI')
               x(i)=Ms_DATA(14)
            Case('P')
               x(i)=Ms_DATA(15)
            Case('S')
               x(i)=Ms_DATA(16)
            Case('CL')
               x(i)=Ms_DATA(17)
            Case('AR')
               x(i)=Ms_DATA(18)
            Case('K')
               x(i)=Ms_DATA(19)
            Case('CA')
               x(i)=Ms_DATA(20)
            Case('SC')
               x(i)=Ms_DATA(21)
            Case('TI')
               x(i)=Ms_DATA(22)
            Case('V')
               x(i)=Ms_DATA(23)
            Case('CR')
               x(i)=Ms_DATA(24)
            Case('MN')
               x(i)=Ms_DATA(25)
            Case('FE')
               x(i)=Ms_DATA(26)
            Case('CO')
               x(i)=Ms_DATA(27)
            Case('NI')
               x(i)=Ms_DATA(28)
            Case('CU')
               x(i)=Ms_DATA(29)
            Case('ZN')
               x(i)=Ms_DATA(30)
            Case('GA')
               x(i)=Ms_DATA(31)
            Case('GE')
               x(i)=Ms_DATA(32)
            Case('AS')
               x(i)=Ms_DATA(33)
            Case('SE')
               x(i)=Ms_DATA(34)
            Case('BR')
               x(i)=Ms_DATA(35)
            Case('KR')
               x(i)=Ms_DATA(36)
            Case default
              Write(6,*) 
              Write(6,*) &
              '>ERROR:: Ms_DATA for atom ',trim(labl(i)),' is not ready.'
              Write(6,*) '>ERROR:: Sorry, the program stops.'
              Stop
         End select

      End do

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RPOLYMER_Mass
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  FUNCTION RPOLYMER_Ene(opt)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: I
    Character(*) :: opt
    Double Precision :: RPOLYMER_Ene
    Character :: ch*120

      Rewind(Inp)
 
      Select Case(trim(opt))

          Case('HF')
 10          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'TOTAL SCF ENERGY')
             if(I==0) goto 10 
             Read(ch(31:50),*) RPOLYMER_Ene

          case('MP2')
 20          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'EMP2')
             if(I==0) goto 20
             Read(ch(8:27),*) RPOLYMER_Ene
          case('MK')   
 30          Continue
             Read(Inp,*,end=2000)RPOLYMER_Ene	
          
          Case default
             goto 1000

      End select

    Return

      1000 Continue
      Write(*,100) trim(opt)
      Stop

100 Format('RPOLYMER_Ene: ERROR OPT=',a10,/, &
           'avaliable options are: HF and MP2 and')

      2000 Continue
      Write(*,200) trim(opt)
      RPOLYMER_Ene=0.0
      return

200 Format('RPOLYMER_Ene: ERROR OPT=',a10,/, &
           'Energy not found in the output!')

! -----------------------------------------------------------------------------

  End FUNCTION RPOLYMER_Ene
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8

  Subroutine RPOLYMER_setLinear(option)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Logical :: option

      linear = option
      if(option) Nfree=Nat3-5

      return

! -----------------------------------------------------------------------------


  End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8

  Subroutine RPOLYMER_setPolymer(option)

    USE RPOLYMER

    Implicit None

! -----------------------------------------------------------------------------

    Logical :: option

      polymer = option
      if(option) Nfree=Nat3-4

      return

! -----------------------------------------------------------------------------


  End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8

