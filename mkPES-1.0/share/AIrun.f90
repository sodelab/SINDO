!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!MK Ab initio run

      Subroutine Run_e(qq,E)

      USE v2ai_mod, ONLY : Nat,Nfree,Method,RunTyp,dryrun

      Implicit None

      Double Precision :: qq(Nfree),x(3,Nat),E
      Double Precision :: RACES_Ene,RG98_Ene,RNWChem_Ene,RPOLYMER_Ene
      Character :: fi*12,fo*10,com*80,charfilenumber*5,chartrack*8
      Integer::filenumber,track
        !Initialize
         E=0.D+00        
         !Don't do AIrun and return 
         if(dryrun) return
 
        !Transform normal coordinates (qq) into Cartesian (x)
         Call nma_q2x(x,qq)
         Open(700,file='PESinfo',status='unknown')
         Open(701,file='WRTinfo',status='unknown')
         Select case(RunTyp)

!MK          
           Case('PES')
           Read(155,*)E
           Write(700,'(f15.8,20f7.3)')E,qq
!MK       
           Case('WRT')
            Open(70,file='track',status='old')
            Read(70,*)track
            E=track
            write(chartrack,'(i8)')track
            write(*,*) chartrack
            Write(701,'(a,20f7.3)')chartrack,qq
            fi='inp.'//adjustl(trim(chartrack))
            Open(7,file=adjustl(trim(fi)),status='new')
            Call Generalinp(7,Nat,x)
            Close(7)
            write(chartrack,'(i8)')track+1
            Call System('echo '//adjustl(trim(chartrack))//' > track')
            Close(70)
          End select

        !Remove files
       ! com='rm '//fi//fo
       ! if (Runtyp/='PES' .and. Runtyp/='WRT')Call System(com)

      End 
!MK
      Subroutine Generalinp(io,Nat,x)

      USE v2ai_mod, ONLY : Lbl,Com,nline

      Implicit None

      Integer :: io,Nat
      Double Precision :: x(3,Nat)

      Integer :: i,j,k
      i=1
        !general input
         Do while(Com(i)/='sindo:geometry') 
            Write(io,'(a)')trim(adjustl(Com(i)))
            i=i+1         
         enddo   
         Do j=1,Nat
            Write(io,'(a2,3f15.8)') Lbl(j),(x(k,j),k=1,3)
         End do
        Do j=i+1,nline
            if(Com(j)/='sindo:end') then
            	Write(io,'(a)')trim(adjustl(Com(j)))
            else 
                exit
            endif
        Enddo
      End
!MK
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine Run_ed(qq,E,d)

      USE v2ai_mod, ONLY : Nat,Nfree,Method,RunTyp,dpl,dryrun

      Implicit None

      Real(8) :: qq(Nfree),x(3,Nat),E,d(3)
      Real(8) :: RACES_Ene,RG98_Ene,RNWChem_Ene
      Character :: fi*6,fo*6,com*80

        !Initialize
         E=0.D+00; d=0.D+00
        !Don't do AIrun and return 
         if(dryrun) return

        !Transform normal coordinates (qq) into Cartesian (x)
         Call nma_q2x(x,qq)

         Select case(RunTyp)

           Case('ACE') 
              fi='ZMAT'
              fo='AC.out'
              !write(6,'(a,3f8.3)') fo,qq

             !ACESII input
              Open(7,file=fi,status='new')
              Call ACES2inp(7,Nat,x)
              Close(7)

             !Run ACESII
              com='./runaces2b '//fo
              Call System(com)

             !Read energy and dipole
              Call RACES_Const(fo,Nat)
              E=RACES_Ene(Method)
              Call RACES_Dpl(d)
              Call RACES_Dest

           Case('GAU') 
              fi='rn.com'
              fo='rn.out'
              !write(6,'(a,3f8.3)') fo,qq

             !GAUSSIAN input
              Open(7,file=fi,status='new')
              Call GAUinp(7,Nat,x)
              Close(7)

             !Run GAUSSIAN
              com='./runGaub '//fi
              Call System(com)

             !Read energy and dipole
              Call RG98_Const(fo,Nat)
              E=RG98_Ene(Method)
              Call RG98_Dpl(d)
              Call RG98_Dest

           Case('NWC')
            fi='NW.nw'
            fo='NW.out'

           !NWCHEM input
            Open(7,file=fi,status='new')
            Call NWCHEMinp(7,Nat,x)
            Close(7)

           !Run NWCHEM
            com='./runnwchemb '//fi//' '//fo
            Call System(com)

           !Read energy and dipole
            Call RNWCHEM_Const(fo,Nat)
            E=RNWCHEM_Ene(trim(Method))
            Call RNWCHEM_Dpl(d)
            Call RNWCHEM_Dest

         End select

        !Remove files
         com='rm '//fi//' '//fo
         Call System(com)

      End 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine Run_en(qq,E,ncc)

      USE v2ai_mod, ONLY : Nat,Nfree,Method,RunTyp,dryrun

      Implicit None

      Real(8) :: qq(Nfree),x(3,Nat),E
      Real(8), dimension(Nat*(Nat-1)/2) :: PSO,DSO,FC,SD,ncc
      Character(5), dimension(Nat*(Nat-1)/2) :: AtmPairs
      Real(8) :: RACES_Ene,RG98_Ene,RNWChem_Ene
      Character :: fi*6,fo*6,com*30

        !Initialize
         E=0.D+00; ncc=0.D+00
        !Don't do AIrun and return 
         if(dryrun) return
        !Transform normal coordinates (qq) into Cartesian (x)
         Call nma_q2x(x,qq)

         Select case(RunTyp)

           Case('ACE') 
              fi='ZMAT'
              fo='AC.out'
              !write(6,'(a,3f8.3)') fo,qq

             !ACESII input
              Open(7,file=fi,status='new')
              Call ACES2inp(7,Nat,x)
              Close(7)

             !Run ACESII
              com='./runaces2b '//fo
              Call System(com)

             !Read energy and dipole
              Call RACES_Const(fo,Nat)
              E=RACES_Ene(Method)
              Call RACES_NMRCC(PSO,DSO,FC,SD,ncc,AtmPairs)
              Call RACES_Dest

           Case default
              Write(6,*) 'ERROR!'
              Write(6,*) 'ERROR! NMR-CC is not ready for Runtyp=',RunTyp
              Stop

         End select

        !Remove files
         com='rm '//fi//' '//fo
         Call System(com)

        !debug
        ! Open(300,file='debug',status='unknown',access='append')
        ! Write(300,'(100a)') AtmPairs
        ! Close(300)

      End 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine RunHess(qq,Ene,gn,hn,idx)

      USE v2ai_mod, ONLY : Nat,Nfree,Method,Typ,RunTyp,path,dryrun

      Implicit None

      Integer :: idx
      Real(8), parameter :: B2A = 0.52917724924d+00
      Real(8) :: qq(Nfree),gn(Nfree),hn(Nfree,Nfree)
      Real(8) :: x(3,Nat),gx(3,Nat),hx(Nat*3,Nat*3)
      Real(8) :: Ene,RACES_Ene,RG98_Ene,RNWChem_Ene
      Character :: fi*6,fo*6,com*80

        !Initialize
         Ene=0.D+00; gn=0.D+00; hn=0.D+00
        !Don't do AIrun and return 
         if(dryrun) return
         if(idx==-1) then
           !Hessian is avaliable
            Select case(Typ)
              Case('GAU') 
                !Read energy, gradient, and Hessian
                 Call RG98_Const(path,Nat)
                 Ene=RG98_Ene(Method)
                 Call RG98_Grad(gx)
                 Call RG98_Hess(hx)
                 Call RG98_Dest
                 gx = gx/B2A
                 hx = hx/B2A/B2A
            End select

         else
           !Transform normal coordinates (qq) into Cartesian (x)
            Call nma_q2x(x,qq)

            Select case(RunTyp)

              Case('GAU') 
                 fi='rn.com'
                 fo='rn.out'
                 !write(6,'(a,3f8.3)') fo,qq

                !GAUSSIAN input
                 Open(7,file=fi,status='new')
                 Call GAUinp(7,Nat,x)
                 Close(7)

                !Run GAUSSIAN
                 com='./runGaub '//fi
                 Call System(com)

                !Read energy, gradient, and Hessian
                 Call RG98_Const(fo,Nat)
                 Ene=RG98_Ene(Method)
                 Call RG98_Grad(gx)
                 Call RG98_Hess(hx)
                 Call RG98_Dest
                 gx = gx/B2A
                 hx = hx/B2A/B2A

              Case('DTB')
                 fi='df.com'
                 fo='df.out'

                !DFTB input
                 Open(7,file=fi,status='new')
                 Call DFTBinp(7,Nat,x)
                 Close(7)

                !Run DFTB
                 com='./rundftb '//fi
                 Call System(com)

                !Read energy, gradient, and Hessian
                 Call RDFTB(Nat,Ene,gx,hx)

            End select

           !Remove files
            com='rm '//fi//' '//fo
            Call System(com)

         endif

         Call nma_g2qg(gx,gn)
         Call nma_h2qh(hx,hn)

      End 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine GAUinp(io,Nat,x)

      USE v2ai_mod, ONLY : Lbl,Com,Com2,nline

      Implicit None

      Integer :: io,Nat
      Real(8) :: x(3,Nat)

      Integer :: i,j,k

        !gaussian input
         Write(io,'(''$RunGauss'')')
         Do i=1,nline
            if(Com(i)=='/n') then
               Write(io,*)
            elseif(Com(i)=='/geom') then
               Do j=1,Nat
                  Write(io,'(a2,3f15.8)') Lbl(j),(x(k,j),k=1,3)
               End do
            elseif(Com(i)/='') then
               Write(io,'(a)') Com(i)
            endif
         End do
         Write(io,*)

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine GMSinp(io,Nat,x)

      Implicit None

      Integer :: io,Nat
      Real(8) :: x(3,Nat)

      Integer :: i,j

        !GAMESS input
         Write(io,'('' $CONTRL SCFTYP=RHF RUNTYP=ENERGY MULT=1 NZVAR=0'')')
         Write(io,'('' MPLEVL=2 ISPHER=1 UNITS=ANGS $END'')')
         Write(io,'('' $MP2 NACORE=1 MP2PRP=.TRUE. $END'')')
         Write(io,'('' $BASIS EXTFIL=.TRUE. GBASIS=CCPVTZ $END'')')
         Write(io,'('' $SYSTEM TIMLIM=10000 MEMORY=40000000 $END'')')
         Write(io,'('' $GUESS GUESS=huckel $END'')')
         Write(io,'('' $DATA'')')
         Write(io,'('' title'')')
         Write(io,'('' C1'')')
         Write(io,'('' C 6.0 '',3(2x,f15.8))') (X(i,1),i=1,3)
         Write(io,'('' H 1.0 '',3(2x,f15.8))') (X(i,2),i=1,3)
         Write(io,'('' H 1.0 '',3(2x,f15.8))') (X(i,3),i=1,3)
         Write(io,'('' H 1.0 '',3(2x,f15.8))') (X(i,4),i=1,3)
         Write(io,'('' H 1.0 '',3(2x,f15.8))') (X(i,5),i=1,3)
         Write(io,'('' $END'',/)')

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine ACES2inp(io,Nat,x)

      USE v2ai_mod, ONLY : Lbl,Com,nline

      Implicit None

      Integer :: io,Nat
      Real(8) :: x(3,Nat)

      Integer :: i,j,k

        !ACESII input
         Do i=1,nline
            if(Com(i)=='/n') then
               Write(io,*)
            elseif(Com(i)=='/geom') then
               Do j=1,Nat
                  Write(io,'(a2,3f15.8)') Lbl(j),(x(k,j),k=1,3)
               End do
               Write(io,*)
            elseif(Com(i)=='*ACES2') then
               Write(io,'(''*ACES2'')')
               Write(io,'(''COORD=CARTESIAN,UNITS=ANGSTROM'')')
            elseif(Com(i)/='') then
               Write(io,'(a)') Com(i)
            endif
         End do
         Write(io,*)

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!MK changed to 1.3

      Subroutine NWCHEMinp(io,Nat,x)

      USE v2ai_mod, ONLY : Lbl,Com

      Implicit None

      Integer :: io,Nat
      Real(8) :: x(3,Nat)

      Integer :: i,j

        !NWCHEM input
         Write(io,'(''start TITLE'')')
         Write(io,'(''echo'')')
         !Write(io,'(''geometry units angstrom abelian noautoz'')')
         Write(io,'(''geometry units angstrom noautoz'')')
         Do i=1,Nat
            Write(io,'(a2,3f15.8)') Lbl(i),(x(j,i),j=1,3)
         End do
         Write(io,'(''end'')')
         Do i=1,10
            if(Com(i)/='') Write(io,'(a)') Com(i)
         End do
         Write(io,'(''task tce energy'')')

      End


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine POLYMERinp(io,Nat,x)

      USE v2ai_mod, ONLY : Lbl,Com,nline

      Implicit None

      Integer :: io,Nat
      Real(8) :: x(3,Nat)

      Integer :: i,j,k

        !NWCHEM input
         Do i=1,nline
            if(Com(i)=='/n') then
               Write(io,*)
            elseif(Com(i)=='/geom') then
               Write(io,'(''GEOMETRY'')')
               Do j=1,Nat
                  Write(io,'(a2,3f15.8)') Lbl(j),(x(k,j),k=1,3)
               End do
               Write(io,'(''X 0.0 0.0 0.0'')')
            elseif(Com(i)/='') then
               Write(io,'(a)') Com(i)
            endif
         End do
         Write(io,'(''JOBTYPE'')')
         Write(io,'(''0'')')
         Write(io,'(''UNITS'')')
         Write(io,'(''ANGSTROM'')')

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine DFTBinp(io,Nat,x)

      USE v2ai_mod, ONLY : Lbl,Com,Com2,nline

      Implicit None

      Integer :: io,Nat
      Real(8) :: x(3,Nat)
      Character :: work*60

      Integer :: i,j,k
      i=1
        !DFTB input
         Write(io,'(i4)') Nat
         Write(io,'(''   2'')')
         Do j=1,Nat
            Write(io,'(a2,3f15.8)') Lbl(j),(x(k,j),k=1,3)
         End do

         Call GetEnv('DFTBWORK',work)
         Write(io,'(a)') trim(work) 

      End
