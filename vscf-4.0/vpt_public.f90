!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/11/02
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Vpt_Construct(ierr)

      USE Vpt_mod

      Implicit None

         Integer :: ierr
         Integer :: i,j,k,kk,l,in,iout,spr_memalloc

         Integer, parameter :: nmax=10000,m=1000
         Integer, dimension(m*nmax):: state
         Integer, dimension(m,nmax):: state_label

         Namelist /vpt/nLvl,nGen,pth0,pth1,pth2,dpt,maxpCUP,maxqCUP,prtth, &
                       ss,Nstate,state,state_label
    !nlvl=0 ss=.t. gives full MP2

      ! ----------------------------------------------------------------

         ierr=0
         Call spr_Getio(in,iout)

         ! ------------------------------------------------------------
         ! >> Read input 

         ! --- default ---

         nLvl=4
         nGen=3
         pth0=500
         pth1=1.D-01
         pth2=5.D-02
         dpt=0
         prtth=5.D-02
         maxpCUP=-1
         maxqCUP=-1

         ss=.false.
         Nstate=1 
         state=0
         state_label=0

         Rewind(in)
         Read(in,vpt,end=10)
      10 Continue

         write(iout,100)
     100 Format(//,'(  ENTER VPT MODULE  )',//)

         if(ss) then
            write(iout,110)
         else
            write(iout,112)
            if(Nstate>nmax) then
               ierr=-1
               Write(iout,'(''  ERROR: MAXIMUM NUMBER OF STATE IS LIMITIED TO '',i4)') nmax
               Write(iout,'(''  ERROR: TERMINATED IN VPT_CONSTRUCT '')')
               return
            endif

            Call Vcor_setTarget(state,state_label(1:Nfree,1:Nstate))

            if(Nfree<=30) then
               Do j=1,Nstate
                  write(iout,115) j,(label0(k+(j-1)*Nfree),k=1,Nfree)
               End do
            else
               kk=Nfree/30
               Do j=1,Nstate
                  write(iout,115) j,(label0(k+(j-1)*Nfree),k=1,30)
                  Do k=2,kk
                     write(iout,116) (label0(l+(j-1)*Nfree),l=30*(k-1)+1,30*k)
                  End do
                  write(iout,116) (label0(k+(j-1)*Nfree),k=30*kk+1,Nfree)
               End do
            endif
            write(iout,*)

         endif
     110 Format(3x,'>> VPT WITH STATE SPECIFIC VSCF REFERENCE',/)
     112 Format(3x,'>> VPT WITH ZERO-POINT VSCF REFERENCE',/)
     115 Format(7x,'TARGET STATE ',i3.3,':',3x,6(5i1,1x))
     116 Format(27x,6(5i1,1x))

         write(iout,140) 
         if(dpt>=1) then 
            write(iout,145) dpt
         elseif(dpt<0) then
            write(iout,144)
         endif
         if(nLvl/=0) then
            if(nLvl<5) then
               write(iout,150) nLvl
            else
               write(iout,'('' ERROR: MAXIMUM LEVEL IS 4 IN CURRENT VERSION'')')
               write(iout,'('' ERROR: RERUN WITH nLvl < 5 '')')
               write(iout,*)
               Stop
            endif
         else
            if(dpt/=0) then
               write(iout,*)
               write(iout,*)
               write(iout,'(''  WARNING: FULL-MP2 IS SPECIFIED'')')
               write(iout,'(''  WARNING: TURNING OFF DPT OPTION'')')
               write(iout,*)
               write(iout,*)
               dpt=0
            endif
            write(iout,151)
         endif
         if(dpt/=0) then
            write(iout,154) nGen
            write(iout,155) pth0
            write(iout,156) pth1
            write(iout,157) pth2
            write(iout,158) prtth
         else
            pth0=0.D+00
         endif
         if(dpt==0) then
            maxqCUP=MR
            write(iout,162) maxqCUP
         else
            if(maxqCUP>maxpCUP) maxpCUP=maxqCUP
            if(maxpCUP>0) then
               write(iout,160) maxpCUP
            else
               write(iout,161) 
            endif
            if(maxqCUP>0) then
               write(iout,162) maxqCUP
            else
               write(iout,163) 
            endif
         endif
         write(iout,*)

     140 Format(3x,'>> VPT OPTIONS',/)
     144 Format(7x,'o DPT      :      1+(2)')
     145 Format(7x,'o DPT      : ',i10)
     150 Format(7x,'o LEVEL    : ',i10)
     151 Format(7x,'o LEVEL    :   FULL-MP2')
     154 Format(7x,'o N_GEN    : ',i10)
     155 Format(7x,'o PTH0     : ',f10.1,' [CM-1]')
     156 Format(7x,'o PTH1     : ',f10.4)
     157 Format(7x,'o PTH2     : ',f10.4)
     158 Format(7x,'o PRINTTH  : ',f10.4)
     160 Format(7x,'o MAX-PCUP : ',i10)
     161 Format(7x,'o MAX-PCUP :  UNLIMITED')
     162 Format(7x,'o MAX-QCUP : ',i10)
     163 Format(7x,'o MAX-QCUP :  UNLIMITED')

         ! ------------------------------------------------------------

         if(lvscf) then
            Write(iout,170)
            Call file_indicator(10,rwf)
            Open(rwf,file='vscf-w.wfn',form='UNFORMATTED')

         else
             Write(iout,*)
             Write(iout,'(8x,''ERROR:  VPT REQUIRES VSCF REFERENCES '')')
             Write(iout,'(8x,''ERROR:  RE-RUN WITH VSCF=.TRUE. '')')
             Write(iout,*)
             Stop

         endif

     170 Format(7x,'o READ VSCF WFN  : vscf-w.wfn')

         ! ------------------------------------------------------------

         write(iout,200)
         Call spr_meminfo
         Call timer(1)
     200 Format(//,'(  SETUP OF VPT COMPLETED  )',//)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vpt_finalz()

      USE Vpt_mod

      Implicit None
      IF(lvscf) CLOSE(rwf)
      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vpt_main()

      USE Vpt_mod

      Implicit None

         Integer :: i,j,k,l,kk
         Integer :: in,io,ierr,spr_memalloc,Ist
         Integer :: lb(Nfree)

         Call spr_Getio(in,io)
         write(io,100)

         Call Vcor_Construct(ierr)
         if(ierr<0) then
            write(io,*) 'ERROR IN VPT_MAIN WHILE SETTING UP VCOR_MODULE'
            Stop
         endif

         write(io,120) 

         Ist=1
         Do while (.true.) 

            if(ss .or. Ist==1) then
               Call Vcor_rVSCF(ierr)
               if(ierr<0) exit
               if(Ist==1 .and. .not. zp) then
                  write(io,*) 'ERROR:  THE FIRST VSCF REFERENCE MUST BE THE ZERO POINT STATE'
                  write(io,*) 'ERROR:  RESET state IN &VSCF'
                  Stop
               endif
            else
               if(Ist==(Nstate+1)) exit
               tar=label0(Nfree*(Ist-1)+1:Nfree*Ist)
               zp=.false.
            endif

            Call Vcor_setKinmat()
            Call Vcor_setWfn()
            Call Vcor_setQmat()

            if(Nfree<31) then
               write(io,125) Ist,tar
            else
               write(io,125) Ist,(tar(i),i=1,30)
               kk=Nfree/30
               Do k=2,kk
                  write(io,126) (tar(i),i=30*(k-1)+1,30*k)
               End do
               write(io,126) (tar(i),i=30*kk+1,Nfree)
            endif
            write(io,*)

            if(dpt/=0) then
               Call Vcor_set_pCnf()

               if(npCnf>0) then
                  write(io,128) npCnf
                  if(Nfree<31) then
                     Do i=1,npCnf
                        Call Vcor_getLabel(i+1,lb)
                        write(io,130) i,lb
                     End do
                  else
                     kk=Nfree/30
                     Do i=1,npCnf
                        Call Vcor_getLabel(i+1,lb)
                        write(io,130) i,(lb(j),j=1,30)
                        Do k=2,kk
                           write(io,131) (lb(j),j=30*(k-1)+1,30*k)
                        End do
                        write(io,131) (lb(j),j=30*kk+1,Nfree)
                     End do
                  endif
                  write(io,*)

               endif
            endif

            if(dpt/=1) then
               if(nLvl>0) then
                  Call Vcor_set_qCnf()
               else
                  Call Vpt_set_qCnf()
               endif
               write(io,135)
               write(io,136) nqCnf
               Do i=1,nqCUP
                  if(nqCnfi(i)>0) write(io,137) i,nqCnfi(i)
               End do
            endif

            Call timer(1)
            Call Vpt_boost(io)

            Ist=Ist+1
            Call Vcor_ClearCnf

            Call timer(1)
            write(io,*)

         End do

         Call Vcor_finalz()

         write(io,200)
         Call spr_meminfo
         Call timer(1)

     100 Format(/'(  ENTER VPT MAIN MODULE  )',//)
     120 Format(/,3x,'[  LOOP OVER VPT STATES  ]',//)
     125 Format(5x,'>> STATE ',i3.3,':',3x,6(5i1,1x))
     126 Format(21x,6(5i1,1x))
     128 Format(9x,'o RESONANT STATES   (P-SPACE) : ',i6,/)
     130 Format(12x,'(',i3.3,'):  ',6(5i1,1x))
     131 Format(20x,6(5i1,1x))
     135 Format(9x,'o CORRELATED STATES (Q-SPACE)',/)
     136 Format(12x,'=== TOTAL ===  : ',i10)
     137 Format(19x,i2,'-MODE : ',i10)
     200 Format(/'(  EXIT VPT MAIN MODULE  )',//)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vpt_set_qCnf()

      USE Vpt_mod

      Implicit None

         Integer :: m2,m3,nq,i,j,k,n1,n2,n3,spr_memalloc
         Integer, allocatable :: cup(:),mode(:,:),vv(:,:)
         !dbg Integer :: ref0(Nfree)

         i=maxCHO-1
         m2=i*i
         m3=m2*i

         nq=Nfree*(Nfree-1)/2 * m2
         if(MR>2) then
            nq=nq + Nfree*(Nfree-1)*(Nfree-2)/6 * m3
         endif
         if(.not. ss) nq=nq + Nfree*i

         i=spr_memalloc(-1,dble(nq)*4.D+00+dble(nq)*dble(MR)*8.D+00)
         Allocate(cup(nq),mode(MR,nq),vv(MR,nq))
         !dbg write(6,'(i7)') nq

         nq=0
         if(.not. ss) then
            Do i=1,Nfree
               Do n1=0,nCHO(i)-1
                  if(n1==tar(i)) cycle
                  nq=nq+1
                  cup(nq)=1
                  mode(1,nq)=i
                  vv(1,nq)=n1-tar(i)
               End do
            End do
         endif
         !dbg write(6,'(i7)') nq

         Do i=1,Nfree
         Do j=1,i-1
            Do n1=0,nCHO(i)-1
               if(n1==tar(i)) cycle
               Do n2=0,nCHO(j)-1
                  if(n2==tar(j)) cycle
                  nq=nq+1
                  cup(nq)=2
                  mode(1,nq)=i
                  mode(2,nq)=j
                  vv(1,nq)=n1-tar(i)
                  vv(2,nq)=n2-tar(j)
               End do
            End do
         End do
         End do
         !dbg write(6,'(i7)') nq

         if(MR==2) goto 1000

         Do i=1,Nfree
         Do j=1,i-1
         Do k=1,j-1
            Do n1=0,nCHO(i)-1
               if(n1==tar(i)) cycle
               Do n2=0,nCHO(j)-1
                  if(n2==tar(j)) cycle
                  Do n3=0,nCHO(k)-1
                     if(n3==tar(k)) cycle
                     nq=nq+1
                     cup(nq)=3
                     mode(1,nq)=i
                     mode(2,nq)=j
                     mode(3,nq)=k
                     vv(1,nq)=n1-tar(i)
                     vv(2,nq)=n2-tar(j)
                     vv(3,nq)=n3-tar(k)
                  End do
               End do
            End do
         End do
         End do
         End do

    1000 Continue
         !dbg write(6,'(i5)') nq
         !dbg Do i=1,nq
         !dbg    ref0=tar
         !dbg    Do j=1,cup(i)
         !dbg       ref0(mode(j,i))=ref0(mode(j,i))+vv(j,i)
         !dbg    End do
         !dbg    write(6,'(i5,2x,10i3)') i,ref0
         !dbg End do

         Call Vcor_setup_qCnf(MR,nq)
         Do i=1,nq
            Call Vcor_add_qCnf(cup(i),mode(1:cup(i),i),vv(1:cup(i),i))
         End do

         Call spr_memdealloc(dble(size(cup))*4.D+00+dble(size(mode))*8.D+00)
         Deallocate(cup,mode,vv)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     VPT with state specific VSCF reference function
!

      Subroutine Vpt_boost(io)

      USE Vpt_mod

      Implicit None

         Integer :: io,spr_memalloc
         Integer :: i,j,k,ii,ij,nlb,lb(Nfree)
         Integer :: cp,Vcor_getCUP
         Real(8) :: cc,Eq,Hmat,Vcor_getHmat

         Integer :: nCI,npCI
         Real(8), allocatable :: Ciw(:,:),Cie(:),Ci(:),Hqp(:,:),Hp(:,:),Hpp(:),H0(:)

         Real(8) :: E0,E1,E2
         Real(8) :: Ezp,E,delEi,delEj
         Real(8), parameter :: H2wvn=2.194746D+05

         Save Ezp

         if(npCnf==0) then

            !  Vibrational MP2

            Call Vcor_getLabel(1,lb)
            E0=E0i(1)
            E1=Vcor_getHmat(0,1,1)-E0

            E2=0.D+00; nlb=1+npCnf
            Do i=1,nqCnf
               nlb=nlb+1
               cp=Vcor_getCUP(1,nlb)
               Hmat=Vcor_getHmat(cp,1,nlb)
               delEi=E0-E0i(nlb)
               if(abs(delEi)>1.D-04) E2=E2 + Hmat*Hmat/delEi
            End do

            if(.not. zp) then
               E=(E0+E1+E2)*H2wvn
               write(io,115) E,E-Ezp
               write(666,'(f10.3)',advance='no')E-Ezp!MK
            else
               Ezp=(E0+E1+E2)*H2wvn
               write(io,110) Ezp
               write(666,666)nLvl,nGen,dpt,maxpCUP,maxqCUP,prtth
               666 Format(/,'VPT(cm^-1)  for nLvl=',i2,' nGen=',i2,' dpt=',i2,' maxpCUP=',i2,' maxqCUP=',i2,' prtth=',e9.2 )         
               write(666,'(f10.3)',advance='no')Ezp !MK
            endif

         else

            ! Calculate H(PP)

            npCI=1+npCnf

            i=spr_memalloc(-1,(dble(npCI*(npCI+1)/2)+dble(npCI*npCI))*8.D+00 &
                              +dble(npCI)*16.D+00)
            Allocate(Hpp(npCI*(npCI+1)/2),Ciw(npCI,npCI),Cie(npCI),Ci(npCI))
  
            ! P-space CI
            k=1
            if(ss) then
               Do i=1,npCI
                  cp=Vcor_getCUP(1,i)
                  if(cp/=1 .and. cp<=MR) then
                     Hpp(k)=Vcor_getHmat(cp,1,i)
                  else
                     Hpp(k)=0.D+00
                  endif
                  k=k+1
   
                  Do j=2,i
                     cp=Vcor_getCUP(j,i)
                     if(cp<=MR) then
                        Hpp(k)=Vcor_getHmat(cp,j,i)
                     else
                        Hpp(k)=0.D+00
                     endif
                     k=k+1
                  End do
               End do

            else

               Do i=1,npCI
                  Do j=1,i
                     cp=Vcor_getCUP(j,i)
                     if(cp<=MR) then
                        Hpp(k)=Vcor_getHmat(cp,j,i)
                     else
                        Hpp(k)=0.D+00
                     endif
                     k=k+1
                  End do
               End do

            endif

            if(dpt==1) then
               Call diag(npCI,npCI,Hpp,Ciw,Cie)
               ! Do i=1,npCI
               !    write(6,'(12f8.4)') Ciw(:,i)
               ! End do
               ! write(6,*)

               Call Vpt_print(npCI,Ezp,Ciw,Cie)

            else

               ! Calculate H(QP)

               i=spr_memalloc(-1,dble(npCI*nqCnf)*8.D+00)
               Allocate(Hqp(nqCnf,npCI))

               if(ss) then
                  nlb=npCI
                  Do j=1,nqCnf
                     nlb=nlb+1
                     cp=Vcor_getCUP(1,nlb)
                     if(cp/=1 .and. cp<=MR) then
                        Hqp(j,1)=Vcor_getHmat(cp,1,nlb)
                     else
                        Hqp(j,1)=0.D+00
                     endif
                  End do

                  Do i=2,npCI
                     nlb=npCI
                     Do j=1,nqCnf
                        nlb=nlb+1
                        cp=Vcor_getCUP(i,nlb)
                        if(cp<=MR) then
                           Hqp(j,i)=Vcor_getHmat(cp,i,nlb)
                        else
                           Hqp(j,i)=0.D+00
                        endif
                     End do
                  End do

               else
                  Do i=1,npCI
                     nlb=npCI
                     Do j=1,nqCnf
                        nlb=nlb+1
                        cp=Vcor_getCUP(i,nlb)
                        if(cp<=MR) then
                           Hqp(j,i)=Vcor_getHmat(cp,i,nlb)
                        else
                           Hqp(j,i)=0.D+00
                        endif
                     End do
                  End do

               endif
            endif

            if(dpt==2) then

               ! QDPT2 (PT + CI)
               nCI=1+npCnf+nqCnf
               i=spr_memalloc(-1,dble(nCI)*8.D+00)
               Allocate(H0(nCI))

               Do i=1,nCI
                  H0(i)=E0i(i)
               End do

               ij=1
               Do i=1,npCI
                  Do j=1,i-1
                     Do k=1,nqCnf
                        delEi=H0(i)-H0(k+npCI)
                        delEj=H0(j)-H0(k+npCI)
                        !if(abs(delEi)>1.D-04 .and. abs(delEj)>1.D-04) &
                        Hpp(ij)=Hpp(ij) + Hqp(k,i)*Hqp(k,j) &
                               *0.5D+00*(1.D+00/delEi+1/delEj)
                     End do
                     ij=ij+1
                  End do

                  Do k=1,nqCnf
                     delEi=H0(i)-H0(k+npCI)
                     !if(abs(delEi)>1.D-04) &
                     Hpp(ij)=Hpp(ij) + Hqp(k,i)*Hqp(k,i)/delEi
                  End do
                  ij=ij+1

               End do

               Call diag(npCI,npCI,Hpp,Ciw,Cie)

               ! write(6,'(f12.4)') Cie*H2wvn-Ezp
               ! Do i=1,npCI
               !    write(6,'(20f8.4)') Ciw(:,i)
               ! End do

               Call Vpt_print(npCI,Ezp,Ciw,Cie)

               Call spr_memdealloc(dble(nCI)*8.D+00)
               Deallocate(H0)

            elseif(dpt==-2) then

               Call diag(npCI,npCI,Hpp,Ciw,Cie)

               Do i=1,npCI
                  Ci=Ciw(i,:)

                  cc=0.D+00; k=1 
                  Do j=1,npCI
                     if(cc<abs(Ci(j))) then
                        cc=abs(Ci(j))
                        k=j
                     endif
                  End do
                  !dbg write(6,*) i,k

                  if(k/=i) then
                     Ci=Ciw(:,i)
                     Ciw(:,i)=Ciw(:,k)
                     Ciw(:,k)=Ci
                     cc=Cie(i)
                     Cie(i)=Cie(k)
                     Cie(k)=cc
                  endif
               End do

               ! cc=Vcor_getHmat(0,1,1)
               ! write(6,'(f12.4)') cc*H2wvn-Ezp
               ! write(6,'(f12.4)') Cie*H2wvn-Ezp
               ! Do i=1,npCI
               !    write(6,'(12f8.4)') Ciw(:,i)
               ! End do
               ! write(6,*)

               Do k=1,npCI

                  Call Vcor_getLabel(k,lb)
                  E0=E0i(k)

                  E2=0.D+00
                  nlb=npCI
                  Do i=1,nqCnf
                     Hmat=0.D+00
                     Do j=1,npCI
                        Hmat=Hmat + Hqp(i,j)*Ciw(j,k)
                     End do

                     nlb=nlb+1
                     !Eq=E0i(nlb)
                     delEi=E0-E0i(nlb)
                     E2=E2 + Hmat*Hmat/delEi
                     !Eq=Vcor_getHmat(0,nlb,nlb)
                     !E2=E2 + Hmat*Hmat/(Cie(k)-Eq)
                  End do
                  Cie(k)=Cie(k)+E2

               End do

               Call Vpt_print(npCI,Ezp,Ciw,Cie)

            endif

            if(allocated(Hqp)) then
               Call spr_memdealloc(dble(npCI*nqCnf)*8.D+00)
               Deallocate(Hqp)
            endif

            Call spr_memdealloc(dble(npCI*(npCI+1)/2)*8.D+00+ &
                   dble(npCI*npCI)*8.D+00+dble(npCI)*16.D+00)
            Deallocate(Hpp,Ciw,Cie,Ci)

         endif

     110 Format(/,10x,'   E(VPT2)   =',f12.2)
     115 Format(/,10x,'   E(VPT2)   =',f12.2,/ &
                  10x,'   E(VPT2)-E0=',f12.2)

         Contains

         Function E0i(iSt)

         Implicit None

            Integer :: iSt,i,j,label(Nfree)
            Real(8) :: E0i,Ei

            E0i=0.D+00
            Call Vcor_getLabel(iSt,label)
            Do j=1,Nfree
               Call Modal_getEne_i(j,label(j),Ei)
               E0i=E0i + Ei
            End do

         End function

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vpt_print(npCI,Ezp,Ciw,Cie)

      USE Vpt_mod

      Implicit None

         Integer :: npCI,in,io,i,j,k,kk
         Integer :: Lbl(npCI),Label(Nfree)
         Real(8) :: E,Ezp,Ciw(npCI,npCI),Cie(npCI),Ci(npCI),Wi(npCI)
         Real(8), parameter :: H2wvn=2.194746D+05, &
                               thsh=1.D-03

         Call spr_Getio(in,io)

         Do i=1,npCI
            Ci=Ciw(:,i)
            Do j=1,npCI
               Wi(j)=Ci(j)*Ci(j)
            End do
            if(Wi(1)< prtth) cycle
            Call sort(npCI,Lbl,Wi)
            Call Vcor_getLabel(Lbl(1),Label)
            if(Nfree<31) then
               if(Lbl(1)/=1) then 
                  write(io,101) Label
               else
                  write(io,100) Label
               endif
            else
               kk=Nfree/30
               if(Lbl(1)/=1) then 
                  write(io,101) Label(1:30)
                  Do k=2,kk
                     write(io,102) Label(30*(k-1)+1:30*k)
                  End do
                  write(io,102) Label(30*kk+1:)
               else
                  write(io,100) Label(1:30)
                  Do k=2,kk
                     write(io,102) Label(30*(k-1)+1:30*k)
                  End do
                  write(io,102) Label(30*kk+1:)
               endif
            endif
            E=Cie(i)*H2wvn
            write(io,115) E,E-Ezp
            write(666,*)"??"
            write(666,'(f12.2)')E-Ezp!MK
            write(io,120)
            Do j=1,npCI 
               if(Wi(j) > thsh) then
                  Call Vcor_getLabel(Lbl(j),Label)
                  if(Nfree<31) then 
                     write(io,125) Ci(Lbl(j)),Wi(j),Label
                  else
                     write(io,125) Ci(Lbl(j)),Wi(j),Label(1:30)
                     kk=Nfree/30
                     Do k=2,kk
                        write(io,126) Label(30*(k-1)+1:30*k)
                     End do
                     write(io,126) Label(30*kk+1:)
                  endif
               else
                  exit
               endif
            End do

         End do

     100 Format(/,9x,'o REFERENCE STATE : ',6(5i1,1x))
     101 Format(/,9x,'o RESONANCE STATE : ',6(5i1,1x))
     102 Format(29x,6(5i1,1x))
     115 Format(/,10x,'   E(VPT)   =',f12.2,/ &
                  10x,'   E(VPT)-E0=',f12.2)
     120 Format(/,10x,'   COEFF.  WEIGHT      CONFIG.') 
     125 Format(13x,f6.3,2x,f6.3,6x,6(5i1,1x))
     126 Format(33x,6(5i1,1x))

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
