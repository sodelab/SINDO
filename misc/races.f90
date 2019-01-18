!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
!   Code description by K.Yagi:  yagi@qcl.t.u-tokyo.ac.jp
!
!   Last modified : 2007/11/15
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  Module RACES

    Integer :: Nat, Nat3, Nfree
    Integer :: Inp=100
    Logical :: linear

  End Module RACES
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Const(ch,N)

    USE RACES

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
      Nfree=Nat3-6

    return

  END SUBROUTINE RACES_Const
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Dest()

    USE RACES

    Implicit None

      Close(Inp)

    return

  END SUBROUTINE RACES_Dest
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Label(x)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i
    Character(2), dimension(Nat) :: x
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=INDEX(ch,'C o o')
      if(i.EQ.0) GOTO 10
      Read(Inp,*)
      Read(Inp,*)

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(6:8),*) x(i)
      End do

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RACES_Label
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Freq(Omega)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Double Precision, dimension(Nfree), intent(out) :: Omega

    Integer :: i,j,ii
    Logical :: op

      ii=200
      Do while(0==0) 
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do

      Open(ii,file='NORMCO',status='old',err=10)

      Do i=1,Nat+1
         Read(ii,*)
      End do

      if(.not. linear) then 
         j=6*(3+Nat)
      else
         j=5*(3+Nat)
      endif
      Do i=1,j
         Read(ii,*)
      End do

      Do i=1,Nfree
         Read(ii,*)
         Read(ii,*) Omega(i)
         Read(ii,*)
         Do j=1,Nat
            Read(ii,*)
         End do
      End do
      Close(ii)

   return

! -----------------------------------------------------------------------------

 10 Continue
    Write(6,*)
    Write(6,*)'ERROR> NORMCO is not found!'
    Write(6,*)'ERROR> Program ends with an error.'
    Stop


  END SUBROUTINE RACES_Freq
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Geom(x)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j
    Double Precision :: B2A
    Double Precision, dimension(3,Nat):: x
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=Index(ch,'C o o r d i n a t e s')
      if(i.EQ.0) GOTO 10
      Read(Inp,'(a)') ch
      Read(Inp,'(a)') ch

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(21:65),'(3f15.8)') (x(j,i),j=1,3)
      End do

      B2A = 0.52917724924d+00
  
    ! x=x*B2A

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RACES_Geom
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Grad(x)

    USE RACES

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

  END SUBROUTINE RACES_Grad
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Hess(H)

    USE RACES

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

  END SUBROUTINE RACES_Hess
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
 
!
  SUBROUTINE RACES_L(CL)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j,ii
    Logical :: op
    Double Precision, dimension(Nat3,Nfree):: CL

      ii=200
      Do while(0==0) 
         Inquire(ii,opened=op)
         if(op) then
            ii=ii+1
         else
            exit
         endif
      End do

      Open(ii,file='NORMCO',status='old',err=10)

      Do i=1,Nat+1
         Read(ii,*)
      End do

      if(.not. linear) then 
         j=6*(3+Nat)
      else
         j=5*(3+Nat)
      endif
      Do i=1,j
         Read(ii,*)
      End do

      Do i=1,Nfree
         Read(ii,*)
         Read(ii,*)
         Read(ii,*)
         Do j=1,Nat
            Read(ii,*) CL(3*(j-1)+1:3*j,i)
         End do
      End do
      Close(ii)

    return

! -----------------------------------------------------------------------------

 10 Continue
    Write(6,*)
    Write(6,*)'ERROR> NORMCO is not found!'
    Write(6,*)'ERROR> Program ends with an error.'
    Stop

  END SUBROUTINE RACES_L
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
! Dipole moment vector in Debye
!
  SUBROUTINE RACES_Dpl(x)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i
    Real(8), dimension(3) :: x
    Character :: ch*120

      Rewind(Inp)

      ! SCF Dipole
   10 Read(Inp,'(A)') ch
      i=INDEX(ch,'dipole moment')
      if(i.EQ.0) GOTO 10

      Read(Inp,100)x

      ! Correlated Dipole
   20 Read(Inp,'(A)',end=1000) ch
      i=INDEX(ch,'dipole moment')
      if(i.EQ.0) GOTO 20

      Read(Inp,100)x
  100 Format(5x,3(6x,f16.10))

 1000 Continue

! -----------------------------------------------------------------------------

    return

  END SUBROUTINE RACES_Dpl
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  SUBROUTINE RACES_Mass(x)

    USE RACES

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
     

      Call RACES_Label(labl)

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

  END SUBROUTINE RACES_Mass
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
! iopt =0 :: RHF type WFN
!      =1 :: UHF type WFH, put alpha orbitals
!      =2 :: UHF type WFH, put bet orbitals
!
  SUBROUTINE RACES_MO(iopt,n,eig,vec)

    USE RACES

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
      Write(6,*) 'ERROR in RACES_MO'
      Write(6,*) 'Set POP=FULL flag in the input file to punch out MO'
      n=-1

    return

! -----------------------------------------------------------------------------

  END SUBROUTINE RACES_MO
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
!
  FUNCTION RACES_Ene(opt)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: I
    Character(*) :: opt
    Double Precision :: RACES_Ene
    Character :: ch*120

      Rewind(Inp)
 
      Select Case(trim(opt))

          Case('HF')
 10          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'E(SCF)')
             if(I==0) goto 10 
             Read(ch(13:33),*) RACES_Ene

          case('MP2')
 20          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'Total MBPT(2)')
             if(I==0) goto 20
             Read(ch(38:55),*) RACES_Ene

          Case('CCSD')
 30          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'CCSD        energy')
             if(I==0) goto 30
             Read(ch(37:57),*) RACES_Ene

          Case('EOM-CC')
 40          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'Total EOM-CCSD')
             if(I==0) goto 40
             Read(ch(36:56),*) RACES_Ene

 42          Continue
             Read(Inp,'(A)',end=44) ch
             I=Index(ch,'Total EOM-CCSD')
             if(I==0) goto 42
             Read(ch(36:56),*) RACES_Ene
             goto 42

 44          Continue

          Case('CCSD(T)')
 50          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'CCSD(T)        =')
             if(I==0) goto 50
             Read(ch(27:47),*) RACES_Ene

          Case('CCSDT')
 60          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'CCSDT       energy is')
             if(I==0) goto 60
             Read(ch(37:57),*) RACES_Ene

          Case default
             goto 1000

      End select

    Return

      1000 Continue
      Write(*,100) trim(opt)
      Stop

100 Format('RACES_Ene: ERROR OPT=',a10,/, &
           'avaliable options are: HF,MP2,CCSD,CCSD(T),CCSDT,EOM-CC')

      2000 Continue
      Write(*,200) trim(opt)
      RACES_Ene=0.0
      return

200 Format('RACES_Ene: ERROR OPT=',a10,/, &
           'Energy not found in the output!')

! -----------------------------------------------------------------------------

  End FUNCTION RACES_Ene
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8

  Subroutine RACES_NMRCC(PSO,DSO,FC,SD,JJ,AtmPairs)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Integer :: i,j,k
    Real(8), dimension(Nat*(Nat-1)/2) :: PSO,DSO,FC,SD,JJ
    Character(5), dimension(Nat*(Nat-1)/2) :: AtmPairs
    Character(2), dimension(Nat) :: Atms
    Character :: ch*120

      Rewind(Inp)
   10 Continue
      Read(Inp,'(A)',end=11) ch
      I=index(ch,'PSO contribution')
      if(I==0) goto 10

      Do i=1,7
         Read(Inp,'(A)') ch
      End do

      Read(Inp,'(A)') ch
      Read(ch(7:8),'(A)') Atms(1)
      j=0
      Do i=2,Nat
         Read(Inp,'(A)') ch
         Read(ch(7:80),*) Atms(i),(PSO(j+k),k=1,i-1)
         j=j+i-1
      End do

      k=1
      Do i=1,Nat
      Do j=1,i-1
         AtmPairs(k)=trim(Atms(i))//trim(Atms(j))
         k=k+1
      End do
      End do
      !dbg Write(6,'(A)') Atms
      !dbg Write(6,'(A)') AtmPairs
      goto 20

      !----------------------------------------------------------
      11 Continue
         Write(*,201)
         PSO=0.D+00
         Do i=1,Nat*(Nat-1)/2
            AtmPairs(i)='--'
         End do

     201 Format('RACES_NMRCC: ERROR ',/, &
                'PSO contributions are not found in the output!')
         Rewind(Inp)

      !----------------------------------------------------------

   20 Continue
      Read(Inp,'(A)',end=22) ch
      I=index(ch,'Fermi-contact contribution')
      if(I==0) goto 20

      Do i=1,8
         Read(Inp,'(A)') ch
      End do

      j=0
      Do i=2,Nat
         Read(Inp,'(A)') ch
         Read(ch(10:80),*) (FC(j+k),k=1,i-1)
         j=j+i-1
      End do
      goto 30

      !----------------------------------------------------------
      22 Continue
         Write(*,202)
         FC=0.D+00

     202 Format('RACES_NMRCC: ERROR ',/, &
                'Fermi-constact contributions are not found in the output!')
         Rewind(Inp)

      !----------------------------------------------------------

   30 Continue
      Read(Inp,'(A)',end=33) ch
      I=index(ch,'Spin-dipole contribution')
      if(I==0) goto 30

      Do i=1,8
         Read(Inp,'(A)') ch
      End do

      j=0
      Do i=2,Nat
         Read(Inp,'(A)') ch
         Read(ch(10:80),*) (SD(j+k),k=1,i-1)
         j=j+i-1
      End do
      goto 40

      !----------------------------------------------------------
      33 Continue
         Write(*,203)
         SD=0.D+00

     203 Format('RACES_NMRCC: ERROR ',/, &
                'Spin-dipole contributions are not found in the output!')
      !----------------------------------------------------------

      Rewind(Inp)
      J=0
   40 Continue
      Read(Inp,'(A)',end=44) ch
      I=index(ch,'DSO contribution')
      if(I==0) goto 40
      if(J==0) then
         ! The first output is for HF density matrix
         J=1
         goto 40
      endif

      Do i=1,8
         Read(Inp,'(A)') ch
      End do

      j=0
      Do i=2,Nat
         Read(Inp,'(A)') ch
         Read(ch(10:80),*) (DSO(j+k),k=1,i-1)
         j=j+i-1
      End do
      goto 50

      !----------------------------------------------------------

      44 Continue
         Write(*,204)
         DSO=0.D+00

     204 Format('RACES_NMRCC: ERROR ',/, &
                'DSO contributions are not found in the output!')
      !----------------------------------------------------------

   50 Continue
      JJ=PSO+DSO+FC+SD

      return

! -----------------------------------------------------------------------------


  End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8


  Subroutine RACES_setLinear(option)

    USE RACES

    Implicit None

! -----------------------------------------------------------------------------

    Logical :: option

      linear = option
      if(option) Nfree=Nat3-5

      return

! -----------------------------------------------------------------------------


  End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
