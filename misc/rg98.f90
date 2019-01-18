!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!   Code description by K.Yagi:  yagi@qcl.t.u-tokyo.ac.jp
!
!   Last modified : 2007/04/20
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  Module RG98

    Integer :: Nat, Nat3, Nfree
    Integer :: Inp=100
    Logical :: linear

  End Module RG98
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Const(ch,N)

    USE RG98

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

  END SUBROUTINE RG98_Const
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Dest()

    USE RG98

    Implicit None

      Close(Inp)

    return

  END SUBROUTINE RG98_Dest
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Force(force)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i,j
    Double Precision, dimension(3,Nat):: force
    Character :: ch*120

      Rewind(Inp)

      j=0
   10 Read(Inp,'(A)') ch
      i=Index(ch,'Forces (Hartrees/Bohr)')
      if(i.EQ.0) GOTO 10

      Read(Inp,'(a)') ch 
      Read(Inp,'(a)') ch 
      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(27:38),*) force(1,i)
         Read(ch(42:53),*) force(2,i)
         Read(ch(57:68),*) force(3,i)
      End do

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_Force

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Label(x)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i
    Integer, dimension(Nat) :: ANum
    Character(*), dimension(Nat) :: x
    Character :: ch*120

    Character(2) :: ALabel(18)
    Data ALabel / 'H', 'He', &
       'Li', 'Be',  'B',  'C', 'N', 'O',  'F', 'Ne', &
       'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'  /

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=INDEX(ch,'Atomic')
      if(i.EQ.0) GOTO 10
      Read(Inp,*)
      Read(Inp,*)

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(13:18),*) ANum(i)
         x(i)=ALabel(ANum(i))
      End do

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_Label
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Freq(ii,Omega)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer, intent(in)::ii
    Double Precision, dimension(ii), intent(out) :: Omega

    Integer :: a,b,c,i,j,k
    Character :: ch*120

      a=Mod(ii-1,3)
      b=(ii-a)/3+1

      Rewind(Inp)

      Do i=1,b

 10      Continue
         Read(Inp,'(A)') ch
         j=Index(ch,'Frequencies -- ')
         if(j==0) goto 10 

         Read(ch(16:72),'(f11.4,2(f23.4))') (Omega(k),k=3*(i-1)+1,3*i)

      End do
 

  END SUBROUTINE RG98_Freq
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Geom(x)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i,j
    Double Precision, dimension(3,Nat):: x
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=Index(ch,'              X           Y           Z')
      if(i.EQ.0) GOTO 10
      Read(Inp,'(a)') ch

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(35:70),'(3f12.6)') (x(j,i),j=1,3)
      End do

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_Geom
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Geom2(x)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i,j
    Double Precision, dimension(3,Nat):: x
    Character :: ch*120

      Rewind(Inp)

      j=0
   10 Read(Inp,'(A)') ch
      i=Index(ch,'Standard orientation')
      if(i.EQ.0) GOTO 10

      Read(Inp,'(a)') ch
      Read(Inp,'(a)') ch
      Read(Inp,'(a)') ch
      Read(Inp,'(a)') ch
      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(35:70),'(3f12.6)') (x(j,i),j=1,3)
      End do

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_Geom2
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Grad(x)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i,j,k,ipes,nl
    Real(8) :: x(3,Nat),xx(Nat3)
    Character :: cc1*120
    Logical :: op

      ! Try reading PES.out
      ipes=20
      Do while(0==0) 
         Inquire(ipes,opened=op)
         if(op) then
            ipes=ipes+1
         else
            exit
         endif
      End do
      Open(ipes,file='PES.out',status='old',err=10)
      Read(ipes,*)
      i=mod(3*Nat,6)
      if(i==0) then
         nl=3*Nat/6
      else
         nl=(3*Nat-i)/6+1
      endif
      Do i=1,nl
         Read(ipes,*)
      End do
      Do i=1,nl-1
         Read(ipes,*) xx(6*(i-1)+1:6*i)
      End do
      Read(ipes,*) xx(6*(nl-1)+1:Nat3)

      Close(ipes)

      k=1
      Do i=1,Nat
      Do j=1,3
         x(j,i)=xx(k)
         k=k+1
      End do
      End do

      return

   10 Continue

      Rewind(Inp)

   20 READ(Inp,'(A)') cc1
      i=Index(cc1,'X              Y              Z')
      IF(i.EQ.0) GOTO 20

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
    
! --------------------------------------------------------------------------------

   return

  END SUBROUTINE RG98_Grad
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Hess(H)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: ipes,nl
    Integer :: i,j,k,l,ll,itr,line
    Real(8) :: H(Nat3,Nat3),xx(Nat3*(Nat3+1)/2)
    Character :: cc1*120
    Logical :: op

      ! Try reading PES.out
      ipes=20
      Do while(0==0) 
         Inquire(ipes,opened=op)
         if(op) then
            ipes=ipes+1
         else
            exit
         endif
      End do
      Open(ipes,file='PES.out',status='old',err=10)
      Read(ipes,*)
      i=mod(3*Nat,6)
      if(i==0) then
         nl=3*Nat/6
      else
         nl=(3*Nat-i)/6+1
      endif
      Do i=1,nl*2
         Read(ipes,*)
      End do
      i=mod(Nat3*(Nat3+1)/2,6)
      if(i==0) then
         nl=Nat3*(Nat3+1)/2/6
      else
         nl=(Nat3*(Nat3+1)/2-i)/6+1
      endif
      Do i=1,nl-1
         Read(ipes,*) xx(6*(i-1)+1:6*i)
      End do
      Read(ipes,*) xx(6*(nl-1)+1:Nat3*(Nat3+1)/2)

      Close(ipes)

      k=1
      Do i=1,Nat3
      Do j=1,i
         H(i,j)=xx(k)
         H(j,i)=xx(k)
         k=k+1
      End do
      End do

      return

   10 Continue

      Rewind(Inp)

   20 Read(Inp,'(A)') cc1
      i=Index(cc1,'Force constants in')
      if(i.EQ.0) goto 20

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

! --------------------------------------------------------------------------------

   return

  END SUBROUTINE RG98_Hess
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
 
!      In the standard output of G98, the unit vectors of Normal coordinate 
! are represented in Cartesian coordinate and normalized.  In this routine,
! they are first transformed into mass weighted coordinate, and then are
! retransformed into Cartesian coordinate.  In this form, the representation 
! of geometry in Cartesian coordinate and in Normal coordinate is related by,
!
!    Delta(C_i) =  Sum(j=1,Nfree){X(i,j)*Q(j)}
!
! where C and Q denotes displacement in Cartesian coordinate and in Normal
! coordinate, respectivly.
 
!
  SUBROUTINE RG98_L(iopt,CL)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: iopt,itrn,i,j,k,l,ii,jj,nmod
    Double Precision:: xm,tmp
    Double Precision, dimension(Nat3):: Zmass
    Double Precision, dimension(Nat3,Nfree):: CL,CL2
    Character :: ch*130

      Rewind(Inp)

      nmod=mod(Nfree,5)
      if(nmod==0) then
         itrn=Nfree/5
      else
         itrn=(Nfree-nmod)/5
      endif
      Do i=1,itrn
         ii=5*(i-1)
   10    Read(Inp,'(A)') ch
         j=INDEX(ch,'Coord Atom Element:')
         if(j.EQ.0) GOTO 10

         Do j=1,Nat*3
            Read(Inp,'(A)') ch
            Read(ch(20:80),*) (CL2(j,ii+k),k=1,5)
         End do
      End do
      if(nmod/=0) then
         ii=5*itrn
   12    Read(Inp,'(A)') ch
         j=INDEX(ch,'Coord Atom Element:')
         if(j.EQ.0) GOTO 12

         Do j=1,Nat*3
            Read(Inp,'(A)') ch
            Read(ch(20:80),*) (CL2(j,ii+k),k=1,nmod)
         End do
      endif
 
!old      itrn=Nfree/3
!old      Do i=1,itrn
!old         ii=3*(i-1)
!old   10    Read(Inp,'(A)') ch
!old         j=INDEX(ch,'Atom AN')
!old         if(j.EQ.0) GOTO 10
!old
!old         jj=0
!old         Do j=1,Nat
!old            Read(Inp,'(A)') ch
!old            Read(ch(10:130),*) ((CL2(jj+l,ii+k),l=1,3),k=1,3)
!old            jj=jj+3
!old         End do
!old     
!old      End do
 
!dbg      Do i=1,Nfree
!dbg         Write(6,'(3f12.6)') CL2(:,i)
!dbg         Write(6,*)
!dbg      End do
!dbg      Stop

      Call RG98_Mass(Zmass)

! Transform into mass weighted coordinate
      Do i=1,Nfree
         l=1
         Do j=1,Nat
            xm=SQRT(Zmass(j))
            Do k=1,3
               CL(l,i)=xm*CL2(l,i)
               l=l+1
            End do
         End do
      End do

! Normalize
      Do i=1,Nfree
         tmp=0
         Do j=1,Nat3
            tmp = tmp + CL(j,i)**2
         End do
         tmp=SQRT(tmp)
         Do j=1,Nat3
            CL(j,i)=CL(j,i)/tmp
         End do
      End do

! Gamess output
      if(iopt==1) then
         Do i=1,Nfree
            l=1
            Do j=1,Nat
               xm=SQRT(Zmass(j))
               Do k=1,3
                  CL(l,i)=CL(l,i)/xm
                  l=l+1
               End do
            End do
         End do
      endif

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_L
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
! Dipole moment vector in Debye
!
  SUBROUTINE RG98_Dpl(x)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i
    Real(8), dimension(3) :: x
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=INDEX(ch,'Dipole')
      if(i.EQ.0) GOTO 10

      Read(Inp,100)x
  100 Format(3(6x,f14.8))

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_Dpl
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_MKPop(x,chr)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i
    Real(8), dimension(3)   :: x
    Real(8), dimension(Nat) :: chr
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=INDEX(ch,'Charges from ESP fit')
      if(i.EQ.0) GOTO 10

      Read(Inp,'(A)')ch
      Read(ch(27:60),*) x
      Read(Inp,*)
      Do i=1,Nat
         Read(Inp,'(A)')ch
         Read(ch(10:25),*) chr(i)
      End do

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_MKPop
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE RG98_Mass(x)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: i
    Double Precision, dimension(Nat):: x
    Character :: ch*120

      Rewind(Inp)

   10 Read(Inp,'(A)') ch
      i=INDEX(ch,'Kelvin.')
      if(i.EQ.0) GOTO 10

      Do i=1,Nat
         Read(Inp,'(a)') ch
         Read(ch(40:50),*) x(i)
      End do

! --------------------------------------------------------------------------------

    return

  END SUBROUTINE RG98_Mass
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
! iopt =0 :: RHF type WFN
!      =1 :: UHF type WFH, put alpha orbitals
!      =2 :: UHF type WFH, put bet orbitals
!
  SUBROUTINE RG98_MO(iopt,n,eig,vec)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

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
      Write(6,*) 'ERROR in RG98_MO'
      Write(6,*) 'Set POP=FULL flag in the input file to punch out MO'
      n=-1

    return

! --------------------------------------------------------------------------------

  END SUBROUTINE RG98_MO
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION RG98_Zpe(Iopt)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: I,Iopt
    Double Precision :: RG98_Zpe
    Character :: ch*120

      Rewind(Inp)

 10   Continue
      Read(Inp,'(A)') ch
      I=Index(ch,'Zero')
      if(I==0) goto 10 

      Select Case(Iopt)

          Case(0)
             Read(Inp,'(A)') ch
             Read(ch(32:43),'(f12.5)') RG98_Zpe

          Case(1)
             Read(ch(32:43),'(f12.1)') RG98_Zpe
             RG98_Zpe=RG98_Zpe*1.D-03

          Case(2)
 15          Continue
             Read(Inp,'(A)') ch
             I=Index(ch,'Zero')
             if(I==0) goto 15 

             Read(ch(47:59),'(f12.6)') RG98_Zpe

          Case default
             Write(*,*) 'RG98_Zpe: ERROR IOPT=',Iopt

      End Select

    Return

! --------------------------------------------------------------------------------

  End FUNCTION RG98_Zpe
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION RG98_Ene(opt)

    USE RG98

    Implicit None

! --------------------------------------------------------------------------------

    Integer :: I
    Character(*) :: opt
    Double Precision :: RG98_Ene
    Character :: ch*120

      Rewind(Inp)
 
      Select Case(trim(opt))

          Case('RHF')
 10          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'E(RHF)')
             if(I==0) goto 10 
             Read(ch(21:36),'(d25.14)') RG98_Ene

          Case('UHF')
 15          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'E(UHF)')
             if(I==0) goto 15 
             Read(ch(21:36),'(d25.14)') RG98_Ene

          case('MP2')
 20          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'EUMP2')
             if(I==0) goto 20
             Read(ch(36:63),*) RG98_Ene

          Case('B3LYP')
 30          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'E(RB+HF-LYP)')
             if(I==0) goto 30
             Read(ch(27:42),'(d25.14)') RG98_Ene

          Case('BLYP')
 35          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'E(RB-LYP)')
             if(I==0) goto 35
             Read(ch(24:39),'(d25.14)') RG98_Ene

          Case('QCISD')
 40          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'Largest amplitude')
             if(I==0) goto 40
             Backspace(Inp)
             Backspace(Inp)
             Backspace(Inp)
             Read(Inp,'(A)') ch
             Read(ch(44:64),*) RG98_Ene

          Case('CCSD(T)')
 50          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'CCSD(T)=')
             if(I==0) goto 50
             Read(ch(10:28),*) RG98_Ene

          Case('MP4(SDQ)')
 60          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'UMP4(SDQ)')
             if(I==0) goto 60
             Read(ch(45:64),*) RG98_Ene

          Case('INDO')
 70          Continue
             Read(Inp,'(A)',end=2000) ch
             I=Index(ch,'Energy=')
             if(I==0) goto 70
             Read(ch(9:26),*) RG98_Ene

          Case default
             goto 1000

      End select

    Return

      1000 Continue
      Write(*,100) trim(opt)
      Stop

100 Format('RG98_Ene: ERROR OPT=',a10,/, &
           'avaliable options are: UHF,MP2,B3LYP,QCISD,CCSD(T),MP4(SDQ)')

      2000 Continue
      Write(*,200) trim(opt)
      RG98_Ene=0.D+00
      return

200 Format('RG98_Ene: ERROR OPT=',a10,/, &
           'Energy not found in the output!')

! --------------------------------------------------------------------------------

  End FUNCTION RG98_Ene
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8


  Subroutine RG98_setLinear(option)

    USE RG98

    Implicit None

! -----------------------------------------------------------------------------

    Logical :: option

      linear = option
      if(option) Nfree=Nat3-5

      return

! -----------------------------------------------------------------------------


  End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
