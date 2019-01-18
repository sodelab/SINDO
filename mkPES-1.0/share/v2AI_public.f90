!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!      Program test
!
!          Call Read_inp(i)
!
!      End
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine v2ai_Const()

      USE v2ai_mod

         Call Read_Inp()
         Call nma_Construct(Nat,Nfree,Zmass,xeq,CL)

      End subroutine 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine v2ai_Dest()

      USE v2ai_mod

         Call nma_Destruct
         Deallocate(xeq,CL,Zmass,Lbl,Omega,Msym)
         if(Sym /= 'C1') Deallocate(Nsym,Lsym)

      End subroutine 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Read_Inp()

      USE v2ai_mod

      Implicit None

        Integer, parameter :: nmax=100

        Integer, dimension(nmax) :: ModeSym
        Integer :: i,j,k,l,m1,m2,maxNsym(1)
        Integer :: vmax_base
        Integer, dimension(nmax) :: vmax
        Character(3) :: ML(nmax)
!        Character :: path*40,Typ*3
        Logical :: op,linear,polymer

        Namelist /mode/Nat,Nfree,path,Typ,Sym,ModeSym,linear,ReadMol,polymer,minlimit,maxlimit
        Namelist /airun/ RunTyp,Method,Com,Com2,dpl,nmrcc,dryrun

!----------------------------------------------------------------------

!       Set default values and read parameters
        Nat=0
        Nfree=0
!        Nfree=4
        Sym='C1'
        ModeSym=0
!        linear=.true.
!        PRINT*,"Nfree inside = ",Nfree
        linear=.false.
        polymer=.false.
        path=''
        Typ=''; RunTyp='';Method='';Com='';Com2=''
        ReadMol=.false.
        dryrun=.false.
        minlimit=-1D+10
        maxlimit=1D+10
        Rewind(5)
        Read(5,mode)
        Rewind(5)
        Read(5,airun,end=10)
     10 Continue

! --    Nat,Nfree
!        PRINT*,"Nfree inside = ",Nfree
        if (linear .and. polymer) then
            write(6,*) 'no linear polymer please'
            stop
        endif
        if(.not. linear) then 
!=            Mfree=Nat*3-6
        else
            Mfree=Nat*3-5
        endif
        if(.not. polymer) then 
            Mfree=Nat*3-6
        else
            Mfree=Nat*3-4
        endif
        !Mfree = 8 != oos
        if(Nfree==0) Nfree=Mfree
        if(Nfree>nmax) then
           Write(6,'(''  ERROR: MAXIMUM NUMBER OF MODE IS '',i4)') nmax
           Write(6,'(''  ERROR: TERMINATED IN READ_INP '')')
           Stop
        endif
        !PRINT*,"Nfree inside = ",Nfree
!MK
        Allocate(xeq(3,Nat),CL(Nat*3,Nfree),Zmass(Nat),Lbl(Nat), &
                 Omega(Nfree),Msym(Nfree))
!MK
    !old Allocate(xeq(3,Nat),CL(Nat*3,Mfree),Zmass(Nat),Lbl(Nat), &
    !old             Omega(Mfree),Msym(Nfree))             

        Write(6,100) 
        Write(6,105) Nat,Nfree
  100   Format(/,'--------------------(  HESSIAN OUTPUT  )--------------------',/)
  105   Format(' -> Number of atoms =',i4,', degrees of freedom =',i4,'.')

! --    Read data from an output of Hessian calculation
! --    xeq,Zmass,Lbl,CL,Omega

        if(.not. ReadMol) then
           Write(6,*) '-> Normal Coordinates read from: ',path
           Write(6,*) '-> Type is (',Typ,')'
           if(Typ=='G98') Typ='GAU'
           if(Typ=='GAU') then
              Call RG98_Const(path,Nat)
              if(linear) Call RG98_setLinear(linear)
              Call RG98_Geom2(xeq)
              Call RG98_Mass(Zmass)
              Call RG98_Label(Lbl)
              Call RG98_L(0,CL)
              Call RG98_Freq(Mfree,Omega)
              Call RG98_Dest
           elseif(Typ=='ACE') then
              Call RACES_Const(path,Nat)
              if(linear) Call RACES_setLinear(linear)
              Call RACES_Geom(xeq)
              Call RACES_Mass(Zmass)
              Call RACES_Label(Lbl)
              Call RACES_L(CL)
              Call RACES_Freq(Omega)
              Call RACES_Dest
           elseif(Typ=='NWC') then
              Call RNWCHEM_Const(path,Nat)
              if(linear) Call RNWCHEM_setLinear(linear)
              Call RNWCHEM_Geom(xeq)
              Call RNWCHEM_Mass(Zmass)
              Call RNWCHEM_Label(Lbl)
              Call RNWCHEM_L(CL)
              Call RNWCHEM_Freq(Omega)
              Call RNWCHEM_Dest
              Method='Automatic'
           elseif(Typ=='PLM') then
              Call RPOLYMER_Const(path,Nat)
              if(linear) Call RPOLYMER_setLinear(linear)
              if(polymer) Call RPOLYMER_setPolymer(polymer)
              Call RPOLYMER_Geom(xeq)
              Call RPOLYMER_Mass(Zmass)
              Call RPOLYMER_Label(Lbl)
              Call RPOLYMER_L(CL)
              Call RPOLYMER_Freq(Omega)
              Call RPOLYMER_Dest
           elseif(Typ=='DTB') then
              Call RDFTB2(Nat,xeq,Zmass,Lbl,Omega,CL)
           endif
        else
           Call Read_Data(Nat,Nfree,xeq,Zmass,Lbl,CL,Omega)
        endif

        Write(6,*) '-> Mass of atoms is:'
        !Do i=1,Nat
        Write(6,'(''    '',3f14.5)') Zmass
        !End do
        Write(6,*) '-> The equilibrium geometry is:'
        Do i=1,Nat
           Write(6,110) Lbl(i),xeq(:,i)
        End do
  110   Format('      ',a2,f10.6,2f14.6)

        Write(6,*) '-> Harmonic Frequency is:'
        Write(6,'(''    '',3f14.4)') Omega

        Write(6,*) '-> Normal displacement vector is:'
        Do i=1,Nfree
           Write(6,'(''  > Mode='',i4)') i
           Write(6,'(''    '',3f14.8)') CL(:,i)
        End do

! --    Information of normal modes
! --    Nfree,Sym,ModeSym

        Write(6,120)
  120   Format(/,'--------------------(  MODE PARAMETER  )--------------------',/)
! --    Sym,Msym
        Do i=1,Nfree
           Msym(i)=ModeSym(i)
        End do
        Sym=trim(AdjustL(Sym))
        Call lw2up(Sym)
        Call getML(Nfree,Sym,Msym,ML,Nirp)
        if(Sym /= 'C1') then
           Allocate(Nsym(Nirp))
           Nsym=0
           Do i=1,Nfree
              Nsym(Msym(i)+1)=Nsym(Msym(i)+1)+1
           End do
           maxNsym=maxval(Nsym)
           Allocate(Lsym(maxNsym(1),Nirp))
           Lsym=0
           Nsym=0
           Do i=1,Nfree
              j=Msym(i)+1
              Nsym(j)=Nsym(j)+1
              Lsym(Nsym(j),j)=i
           End do
           !dbg Write(6,'(4i3)') Nsym
           !dbg Do i=1,Nirp
           !dbg    Write(6,'(5i3)') Lsym(:,i)
           !dbg End do
        endif
        Write(6,125) Nfree,Sym
        if(Nfree > 10) then
           Write(6,126) ML(1:10)
           Write(6,127) ML(11:Nfree)
        else
           Write(6,126) ML(1:Nfree)
        endif
  125   Format(' -> # of Mode    : ',i3,/, &
               '  > Irred. rep.  :  ',a3)
  126   Format('  > Mode Symmetry: ',10a4)
  127   Format('                   ',10a4)

! --    check xyz orientation
        if(dpl) then
           Select case(Sym)
              Case('C1')
                 Continue
              Case('C2v')
                 Call Coord_C2v(Nat,Nfree,xeq)
              Case('Cs')
                 Call Coord_Cs(Nat,Nfree,xeq)
              Case default
                 Write(6,*) '-- Error --'
                 Write(6,*) Sym,' is not supported'
                 Stop
           end Select
        endif
 
! --    Set options for ab initio runs
! --    Method,Com,Com2,RunTyp,dpl,nmrcc

        Write(6,130)
  130   Format(/,'--------------------(  AB INITIO RUNS  )--------------------',/)
        if(RunTyp=='') RunTyp=Typ
        if(RunTyp=='G98') RunTyp='GAU'
        Write(6,'('' -> RunTyp is ('',a3,'')'')') RunTyp
        if(Method=='' .and. RunTyp /= 'DTB' .and. RunTyp /= 'WRT' .and. RunTyp /= 'PES') then
           Write(6,*) 
           Write(6,*) 'ERROR> Method is not defined'
           Write(6,*) 'ERROR> Program will end'
           Stop
        endif
        Write(6,135) trim(Method)
        if(dpl) Write(6,136)
        if(nmrcc) Write(6,137)
  135   Format(' -> Method       : ',a10)
  136   Format(' -> Dipole       :  [  TRUE  ]')
  137   Format(' -> NMR-CC       :  [  TRUE  ]')
        Write(6,'(''  > Sample input for ab initio program:'',/)')
        Write(6,'(25(''/*''))')
        if(Runtyp=='PES') then
           write(*,*)'PES is used. MK'
        elseif(Runtyp=='WRT') then
          write(*,*)'General method "WRT" is used. MK'
          Call Generalinp(6,Nat,xeq)
        else
           Write(6,*) 'ERROR!'
           Write(6,*) 'Runtyp is not set in airun group'
           Write(6,*) 'ERROR!'
           Stop
        endif
        Write(6,'(25(''/*''))')

        return

      Contains

      Subroutine getML(Nf,Sym,Msym,ML,Nirp)

         Integer :: Nf,Msym(Nf),Nirp,maxsym(1),i
         Character :: Sym*3
         Character(3) :: ML(Nf)
         Character(3) :: ci(2),cs(2),c2(2),d2(4),c2v(4),c2h(4),d2h(8)

         Data ci  / 'Ag ','Au ' /
         Data cs  / 'A` ','A" ' /
         Data c2  / 'A  ','B  ' /
         Data d2  / 'A1 ','B1 ','B2 ','B3 ' /
         Data c2v / 'A1 ','A2 ','B1 ','B2 ' /
         Data c2h / 'Ag ','Au ','Bg ','Bu ' /
         Data d2h / 'Ag ','Au ','B1g','B1u','B2g','B2u','B3g','B3u'/

         maxsym=maxval(Msym)

         Select case(Sym)

            case('C1')
               ML='A  '
               Nirp=1

            case('CI')
               if(maxsym(1) > 2) then 
                  Write(6,*) '-- Error --'
                  Write(6,*) '  illegal number specified for mode symmetry '
                  Stop
               endif
               Do i=1,Nf
                  ML(i)=ci(Msym(i)+1)
               End do
               Nirp=2

            case('CS')
               if(maxsym(1) > 2) then 
                  Write(6,*) '-- Error --'
                  Write(6,*) '  illegal number specified for mode symmetry '
                  Stop
               endif
               Do i=1,Nf
                  ML(i)=cs(Msym(i)+1)
               End do
               Nirp=2

            case('C2')
               if(maxsym(1) > 2) then 
                  Write(6,*) '-- Error --'
                  Write(6,*) '  illegal number specified for mode symmetry '
                  Stop
               endif
               Do i=1,Nf
                  ML(i)=c2(Msym(i)+1)
               End do
               Nirp=2

            !case('D2')
            !   if(maxsym(1) > 4) then 
            !      Write(6,*) '-- Error --'
            !      Write(6,*) '  illegal number specified for mode symmetry '
            !      Stop
            !   endif
            !   Do i=1,Nf
            !      ML(i)=d2(Msym(i)+1)
            !   End do
            !   Nirp=4

            case('C2V')
               if(maxsym(1) > 4) then 
                  Write(6,*) '-- Error --'
                  Write(6,*) '  illegal number specified for mode symmetry '
                  Stop
               endif
               Do i=1,Nf
                  ML(i)=c2v(Msym(i)+1)
               End do
               Nirp=4

            case('C2H')
               if(maxsym(1) > 4) then 
                  Write(6,*) '-- Error --'
                  Write(6,*) '  illegal number specified for mode symmetry '
                  Stop
               endif
               Do i=1,Nf
                  ML(i)=c2h(Msym(i)+1)
               End do
               Nirp=4

            case('D2H')
               if(maxsym(1) > 8) then 
                  Write(6,*) '-- Error --'
                  Write(6,*) '  illegal number specified for mode symmetry '
                  Stop
               endif
               Do i=1,Nf
                  ML(i)=d2h(Msym(i)+1)
               End do
               Nirp=8

            case default
               Write(6,*) '-- Error --'
               Write(6,*) Sym,' is not supported'
               Stop

         End select

      End subroutine

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine SymOp(i,Ms,sign)

      USE v2ai_mod

      Implicit None

        Integer :: i,Ms,sign(3)

          Ms=Msym(i)

          sign=1
          if(Sym=='C2v') then

             Select case(Msym(i))

               case(1)
                 sign(2)=-1
               case(2)
                 sign(1)=-1

             End select

          elseif(Sym=='Cs') then

             if(Msym(i)==1) then
                sign(3)=-1
             endif

          endif

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function getSym()

      USE v2ai_mod

      Implicit None

        Character :: getSym*3

          getSym=Sym

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     getChr =0  (Ag)
!            =1  (not Ag)
!

      Function getChr(mi,mj,mk)

      USE v2ai_mod

      Implicit None

        Integer :: getChr,mi,mj,mk,jch,i

        getChr=0
        Do i=1,8
           jch=ich(i,mi)*ich(i,mj)*ich(i,mk)
           if(jch/=1) then 
              getChr=1 
              return
           endif
        End do

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!      Function get_Num_of_Irrp(Nir)
!
!      USE v2ai_mod
!
!      Implicit None
!
!      Integer :: get_Num_of_Irrp,Nir
!
!         if(Nir > Nirp) then 
!             Write(6,*) 'ERROR in Function Get_Num_of_Irrp'
!             Write(6,*) 'Nir is larger than Nirp'
!         endif
!
!         get_Num_of_Irrp = Nsym(Nir+1)
!
!      End Function
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!      Subroutine get_Modes_of_Irrp(Nir,Lm)
!
!      USE v2ai_mod
!
!      Implicit None
!
!      Integer :: i,Nir,Lm(*)
!
!         if(Nir > Nirp) then 
!             Write(6,*) 'ERROR in Function Get_Num_of_Irrp'
!             Write(6,*) 'Nir is larger than Nirp'
!         endif
!
!         Do i=1,Nsym(Nir+1)
!            Lm(i)=Lsym(i,Nir+1)
!         End do
!
!      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function getIrrep1(mode)

      USE v2ai_mod

      Implicit None

         Integer :: getIrrep1,mode

         getIrrep1=Msym(mode)

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function getIrrep2(mode1,mode2)

      USE v2ai_mod

      Implicit None

         Integer :: getIrrep2,mode1,mode2,Ms1,Ms2,Ms
         Integer :: mul4(3,3),mul8(7,7)

         Data mul4 / 0,3,2, 3,0,1, 2,1,0 /
         Data mul8 / 0,3,2,5,4,7,6, &
                     3,0,1,6,7,4,5, &
                     2,1,0,7,6,5,4, &
                     5,6,7,0,1,2,3, &
                     4,7,6,1,0,3,2, &
                     7,4,5,2,3,0,1, &
                     6,5,4,3,2,1,0  /

         Ms1=Msym(mode1)
         Ms2=Msym(mode2)

         if(Ms1==0) then
            getIrrep2=Ms2
            return
         endif

         if(Ms2==0) then
            getIrrep2=Ms1
            return
         endif

         if(Ms1==Ms2) then
            getIrrep2=0
            return
         endif

         if(Sym /= 'D2H') then
            getIrrep2=mul4(Ms1,Ms2)
         else
            getIrrep2=mul8(Ms1,Ms2)
         endif

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function get_Nat()

      USE v2ai_mod

      Integer :: Get_Nat

          get_Nat=Nat

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function get_Nfree()

      USE v2ai_mod

      Integer :: get_Nfree

          get_Nfree=Nfree

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function get_dpl()

      USE v2ai_mod

      Logical :: get_dpl

          get_dpl=dpl

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function get_nmrcc()

      USE v2ai_mod

      Logical :: get_nmrcc

          get_nmrcc=nmrcc

      End Function

!---+----1----+--!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function Get_dryrun()

      USE v2ai_mod

      Logical :: Get_dryrun

          Get_dryrun=dryrun

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

!
      Subroutine get_Omega(w)

      USE v2ai_mod

      Real(8) :: w(Nfree)

          w=Omega

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     Z-axis is the C2 axis.
!
      Subroutine Coord_C2v(Nat,Nfree,xeq)

      Implicit None

        Integer :: Nat,Nfree
        Real(8) :: xeq(3,Nat)

        Integer :: i,j,k,ix
        Real(8) :: xx(3)

          ix=0
          Do j=1,Nat
             xx=xeq(:,j)
             Do k=1,2
                xx(k)=-xx(k)
             End do
             !dbg Write(6,'(3f10.4)') xx
             Do k=1,Nat
                if(abs(xx(1)-xeq(1,k))<1.0D-7 .and. &
                   abs(xx(2)-xeq(2,k))<1.0D-7 .and. &
                   abs(xx(3)-xeq(3,k))<1.0D-7 ) then 
                   ix=1
                   exit
                endif
             End do
             !dbg Write(6,'(i3)') ix
             if(ix==0) then
                Write(6,*) 
                Write(6,*) '!--(ERROR)--!'
                Write(6,*) '   The input Z-axis must be the C2 axis'
                Stop
             endif
             ix=0
          End do

          Do i=1,2
             ix=0
             Do j=1,Nat
                if(abs(xeq(i,j))>1.0D-7) then
                   xx=xeq(:,j)
                   xx(i)=-xx(i)
                   Do k=1,Nat
                      if(abs(xx(1)-xeq(1,k))<1.0D-7 .and. &
                         abs(xx(2)-xeq(2,k))<1.0D-7 .and. &
                         abs(xx(3)-xeq(3,k))<1.0D-7 ) then 
                         ix=1
                         exit
                      endif
                   End do
                   if(ix==0) then
                      Write(6,*) 
                      Write(6,*) '!--(ERROR)--!'
                      Write(6,*) '   The input geometry is not C2v'
                      Stop
                   endif
                   ix=0
                endif
             End do
          End do

      End Subroutine 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     XY plane.
!
      Subroutine Coord_Cs(Nat,Nfree,xeq)

      Implicit None

        Integer :: Nat,Nfree,Isym(3)
        Real(8) :: xeq(3,Nat),CL(3,Nat,Nfree)

        Integer :: i,j,k,ix
        Real(8) :: xx(3)

          ix=0
          Do j=1,Nat
             xx=xeq(:,j)
             xx(3)=-xx(3)
             !dbg Write(6,'(3f10.4)') xx
             Do k=1,Nat
                if(abs(xx(1)-xeq(1,k))<1.0D-7 .and. &
                   abs(xx(2)-xeq(2,k))<1.0D-7 .and. &
                   abs(xx(3)-xeq(3,k))<1.0D-7 ) then 
                   ix=1
                   exit
                endif
             End do
             !dbg Write(6,'(i3)') ix
             if(ix==0) then
                Write(6,*) 
                Write(6,*) '!--(ERROR)--!'
                Write(6,*) '   XY plane must be the Cs plane'
                Stop
             endif
             ix=0
          End do

      End Subroutine 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine Read_Data(Nat,Nf,xeq,Zmass,Lbl,CL,Omg)

      Implicit None

        Integer :: Nat,Nf
        Real(8) :: xeq(Nat*3),Zmass(Nat),CL(3*Nat*Nf),Omg(Nf)
        Character(2) :: Lbl(Nat)

        Integer :: i,j,k
        Integer, parameter :: n=200,m=3*n-6
        Double Precision :: mass(n), x(3*n)
        Double Precision :: omega(n)!MK omega(m)
        Double Precision, dimension(3*n*m) :: L
        Character(2) :: Label(n)

        Namelist /mol/ mass,x,omega,L,Label

           Rewind(5)
           Read(5,mol)

           Do i=1,Nat
              Zmass(i)=mass(i)
              Lbl(i)=Label(i)
           End do
           Do i=1,Nf
              Omg(i)=omega(i)
           End do
           Do i=1,Nat*3
              xeq(i)=x(i)
           End do
           Do i=1,3*Nat*Nf
              CL(i)=L(i)
           End do


      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine intpl1(ni,nd,qi,VV)
      USE v2ai_mod !MK

      Implicit None

         Integer :: ni,nd
         Real(8) :: qi(ni),VV(ni,nd),q0(ni),V0(ni,nd),sV0(ni)

         Integer :: i,j,jj,npt,ipt(ni),i0,i1,i2
            if(qi(1)<qi(ni)) then
               i0=1;i1=ni;i2=1
            else
               i0=ni;i1=1;i2=-1
            endif
            jj=0
            npt=0;ipt=0
            q0=0.D+00
            V0=0.D+00
            !Do j=ni,1,-1
            Do j=i0,i1,i2
!MK
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
               if(VV(j,1)>minlimit .and. VV(j,1)<maxlimit) then !MK
                  jj=jj+1
                  q0(jj)=qi(j)
                  Do i=1,nd
                     V0(jj,i)=VV(j,i)
                  End do
               else
                  npt=npt+1
                  ipt(npt)=j
               endif
            End do

            if(npt==0) return
            !Write(6,'(11f10.4)') qi
            !Write(6,'(11f10.4)') VV(:,1)
            Write(123,*)npt," points interpolated."!MK
            Do i=1,nd
               Call DSpline(q0(1:jj),V0(1:jj,i),jj,1D+32,1D+32,sV0(1:jj))
               Do j=1,npt
                  Call DSplint(q0(1:jj),V0(1:jj,i),sV0(1:jj),jj, &
                               qi(ipt(j)),VV(ipt(j),i))
               End do
            End do

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine intpl2(ni,nj,nd,qj,VV)
      USE v2ai_mod !MK
      Implicit None

         Integer :: ni,nj,nd
         Real(8) :: qj(nj),VV(nj,ni,nd),q0(nj),V0(nj,nd),sV0(nj)

         Integer :: i,j,k,jj,npt,ipt(nj),j0,j1,j2

         if(qj(1)<qj(nj)) then
            j0=1;j1=nj;j2=1
         else
            j0=nj;j1=1;j2=-1
         endif

         Do i=1,ni

            jj=0
            npt=0;ipt=0
            q0=0.D+00
            V0=0.D+00
            !Do j=nj,1,-1
            Do j=j0,j1,j2
               if(VV(j,i,1)>minlimit .and. VV(j,i,1)<maxlimit) then 
                  jj=jj+1
                  q0(jj)=qj(j)
                  Do k=1,nd
                     V0(jj,k)=VV(j,i,k)
                  End do
               else
                  npt=npt+1
                  ipt(npt)=j
               endif
            End do
            if(npt==0) cycle
            !Write(6,'(11f10.4)') qj
            !Write(6,'(11f10.4)') VV(j,i)
            Write(123,*)npt," points interpolated." !MK  
            Do k=1,nd
               Call DSpline(q0(1:jj),V0(1:jj,k),jj,1D+32,1D+32,sV0(1:jj))
               Do j=1,npt
                  Call DSplint(q0(1:jj),V0(1:jj,k),sV0(1:jj),jj, & 
                               qj(ipt(j)),VV(ipt(j),i,k))
               End do
            End do

         End do

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine intpl3(ni,nj,nk,nd,qk,VV)
      USE v2ai_mod !MK
      Implicit None

         Integer :: ni,nj,nk,nd
         Real(8) :: qk(nk),VV(nk,nj,ni,nd),q0(nk),V0(nk,nd),sV0(nk)

         Integer :: d,i,j,k,jj,npt,ipt(nk),k0,k1,k2

         if(qk(1)<qk(nk)) then
            k0=1;k1=nk;k2=1
         else
            k0=nk;k1=1;k2=-1
         endif

         Do i=1,ni
         Do j=1,nj

            jj=0
            npt=0;ipt=0
            q0=0.D+00
            V0=0.D+00
            !Do k=nk,1,-1
            Do k=k0,k1,k2
               if(VV(k,j,i,1)> minlimit .and. VV(k,j,i,1)< maxlimit) then 
                  jj=jj+1
                  q0(jj)=qk(k)
                  Do d=1,nd
                     V0(jj,d)=VV(k,j,i,d)
                  End do
               else
                  npt=npt+1
                  ipt(npt)=k
               endif
            End do
            if(npt==0) cycle
            !Write(6,'(11f10.4)') qk
            !Write(6,'(11f10.4)') VV(:,j,i)
           Write(123,*)npt," points interpolated." !MK
            Do d=1,nd
               Call DSpline(q0(1:jj),V0(1:jj,d),jj,1D+32,1D+32,sV0(1:jj))
               Do k=1,npt
                  Call DSplint(q0(1:jj),V0(1:jj,d),sV0(1:jj),jj, &
                               qk(ipt(k)),VV(ipt(k),j,i,d))
               End do
            End do

         End do
         End do

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
