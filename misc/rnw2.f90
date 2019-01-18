!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
!   Code description by K.Yagi:  yagi@qcl.t.u-tokyo.ac.jp
!
!   Last modified : 2007/12/11
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  Module RNW2

    Integer :: Nat, Nat3, Nfree
    Integer :: Inp=100
    Logical :: linear

  End Module RNW2
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Const(ch,N)

    USE RNW2

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
      Nfree=Nat3-6

    return

  END SUBROUTINE RNW2_Const
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Dest()

    USE RNW2

    Implicit None

      Close(Inp)

    return

  END SUBROUTINE RNW2_Dest
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Label(x)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i
    Character(2), dimension(Nat) :: x
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=INDEX(ch,'Atom information')
      if(i.EQ.0) GOTO 10
      Read(Inp,*)
      Read(Inp,*)

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(3:5),*) x(i)
      End do

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RNW2_Label
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Freq(Omega)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Double Precision, dimension(Nfree), intent(out) :: Omega

    Integer :: i,j,ii
    Logical :: op
    Double Precision :: dummy

      ii=200
      Do while(0==0) 
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do

      Open(ii,file='nmode',status='old',err=10)

      if(.not. linear) then
         j=6
      else
         j=5
      endif

      Do i=1,j
         Read(ii,*) dummy
      End do

      Do i=1,Nfree
         Read(ii,*) Omega(i)
      End do
      Close(ii)

   return

! -----------------------------------------------------------------------------

 10 Continue
    Write(6,*)
    Write(6,*)'ERROR> nmode is not found!'
    Write(6,*)'ERROR> Program ends with an error.'
    Stop


  END SUBROUTINE RNW2_Freq
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Geom(x)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j
    Double Precision :: B2A
    Double Precision, dimension(3,Nat):: x
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=Index(ch,'Atom information')
      if(i.EQ.0) GOTO 10
      Read(Inp,'(a)') ch
      Read(Inp,'(a)') ch

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(15:59),'(3f15.7)') (x(j,i),j=1,3)
      End do

      B2A = 0.52917724924d+00
      x=x*B2A

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RNW2_Geom
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Grad(x)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j
    Double Precision, dimension(3,Nat):: x
    Character :: cc1*120

      Rewind(Inp)

   10 READ(Inp,'(A)') cc1
      i=Index(cc1,'X              Y              Z')
      IF(i.EQ.0) GOTO 10

      Read(Inp,'(a)') cc1
      Do i=1,Nat
         Read(Inp,'(a)') cc1
         Read(cc1(20:80),*) (x(j,i),j=1,3)
      End do

      Do i=1,Nat
         Do j=1,3
            x(j,i)=-x(j,i)
         End do
      End do
    
! -----------------------------------------------------------------------------

   return

  END SUBROUTINE RNW2_Grad
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Hess(H)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j,k,l,ll,itr,line
    Double Precision, dimension(Nat3,Nat3):: H
    Character :: cc1*120

      Rewind(Inp)

   10 Read(Inp,'(A)') cc1
      i=Index(cc1,'Force constants')
      if(i.EQ.0) goto 10

      if(Mod(Nat3,5).EQ.0) then
        itr=Nat3/5
      else
        itr=(Nat3-MOD(Nat3,5))/5+1
      endif

      Do l=1,itr
         ll=5*(l-1)
         Read(Inp,'(A)') cc1
         Do i=1,Nat3-ll
            line=5
            if(i.LT.5) line=i
            Read(Inp,'(A)') cc1
            Read(cc1,*) j,(H(i+ll,k+ll),k=1,line)
         End do
      End do

      Do j=1,nat3
         Do i=1,j
            H(i,j)=H(j,i)
         End do
      End do

! -----------------------------------------------------------------------------

   return

  END SUBROUTINE RNW2_Hess
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
 
!
  SUBROUTINE RNW2_L(CL)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j,ii
    Logical :: op
    Double Precision, dimension(Nat3,Nfree):: CL
    Double Precision :: dummy

      ii=200
      Do while(0==0) 
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do

      Open(ii,file='nmode',status='old',err=10)

      Do i=1,Nat3
         Read(ii,*) dummy
      End do

      if(.not. linear) then
         j=6
      else
         j=5
      endif

      Do i=1,j
         Do j=1,Nat3
            Read(ii,*) dummy
         End do
      End do

      Do i=1,Nfree
         Do j=1,Nat3
            Read(ii,*) CL(j,i)
         End do
      End do
      Close(ii)

    return

! -----------------------------------------------------------------------------

 10 Continue
    Write(6,*)
    Write(6,*)'ERROR> nmode is not found!'
    Write(6,*)'ERROR> Program ends with an error.'
    Stop

  END SUBROUTINE RNW2_L
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
! Dipole moment vector in Debye
!
  SUBROUTINE RNW2_Dpl(x)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i
    Real(8), dimension(3) :: x
    Integer :: ii
    Logical :: op

      ii=300
      Do while(0==0)
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do

      Open(ii,file='dipole',status='old',err=1000)

      Read(ii,*,end=1000) x

      Close(ii)

    Return

 1000 Continue

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RNW2_Dpl
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RNW2_Mass(x)

    USE RNW2

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
     

      Call RNW2_Label(labl)

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

  END SUBROUTINE RNW2_Mass
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
! iopt =0 :: RHF type WFN
!      =1 :: UHF type WFH, put alpha orbitals
!      =2 :: UHF type WFH, put bet orbitals
!
  SUBROUTINE RNW2_MO(iopt,n,eig,vec)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j,n,iopt
    Double Precision, dimension(n):: eig
    Double Precision, dimension(n,n):: vec
    Character :: ch*120


      ! Check if POP=FULL flag is set.
      Rewind(Inp)

   10 Read(Inp,'(A)',end=100) ch
      i=INDEX(ch,'POP=FULL')
      if(i.EQ.0) GOTO 10

      ! Read MO
      Rewind(Inp)
      eig=0.D+00
      vec=0.D+00

      if(iopt.eq.2) then
      15 Read(Inp,'(A)') ch
         i=INDEX(ch,'Beta Molecular Orbital')
         if(i.EQ.0) GOTO 15
      endif

      j=1
      Do while(0==0)

      20 Read(Inp,'(A)') ch
         i=INDEX(ch,'EIGEN')
         if(i.EQ.0) GOTO 20

         if(j+4 <= n) then 
            Read(ch(20:72),*) eig(j:j+4)

            Do i=1,n
               Read(Inp,'(A)') ch
               Read(ch(20:72),*) vec(i,j:j+4)
            End do
         else
            Read(ch(20:72),*) eig(j:n)
            Do i=1,n
               Read(Inp,'(A)') ch
               Read(ch(20:72),*) vec(i,j:n)
            End do
            exit
         endif
         j=j+5

      End do

    return

! ----------------------------------------------!

  100 Continue
      Write(6,*) 'ERROR in RNW2_MO'
      Write(6,*) 'Set POP=FULL flag in the input file to punch out MO'
      n=-1

    return

! -----------------------------------------------------------------------------

  END SUBROUTINE RNW2_MO
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  FUNCTION RNW2_Ene(option,track)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: ii,option
    Logical :: op
    Character(*) :: opt
    Double Precision :: RNW2_Ene

      ii=300
      Do while(0==0)
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do
      write(filename,'('nw.out',i3)')track
      Open(ii,file='eneragy',status='old',err=1000)

      Read(ii,*) RNW2_Ene

      Close(ii)
 
    Return

      1000 Continue
      Write(6,*)
      Write(6,*)'ERROR> energy is not found!'
      Write(6,*)'ERROR> Program ends with an error.'
      Stop

! -----------------------------------------------------------------------------

  End FUNCTION RNW2_Ene
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8


  Subroutine RNW2_setLinear(option)

    USE RNW2

    Implicit None

! -----------------------------------------------------------------------------

    Logical :: option

      linear = option
      if(option) Nfree=Nat3-5

      return

! -----------------------------------------------------------------------------

  End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
