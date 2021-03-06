!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/12/19
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Vscf_Construct(i)

      USE Vscf_mod
     
      Implicit None

         Integer :: i,i1,i2,j,j1,j2,k,kk,l,in,iout,spr_memalloc

         Integer, parameter :: nmax=10000,m=1000
         Integer, dimension(m*nmax):: state
         Integer, dimension(m,nmax):: state_label

         Real(8) :: msz
 
         Logical :: op
 
         Namelist /vscf/Nstate,state,state_label,restart,Maxitr,Ethresh

         i=0
         Call spr_Getio(in,iout)

         if(Nfree>nmax) then
            i=-1
            Write(iout,'(''  ERROR: MAXIMUM NUMBER OF MODE IS '',i4)') nmax
            Write(iout,'(''  ERROR: TERMINATED IN VSCF_CONSTRUCT '')')
            return
         endif

         ! --+----2----+----3----+----4----+----5----+----6----+----7----+----80
         ! >> Read input 

         ! --- default ---
         ! - The number of state(s) to be calculated and quantum numbers of 
         ! - each mode.  Default is only the ground state. 
         Nstate=1
         state=0
         state_label=0
         restart=.false.
         Maxitr=10
         Ethresh=1.0D-03

         Rewind(in)
         Read(in,vscf,end=10)
      10 Continue

         if(Nstate>m) then 
            i=-1
            Write(iout,'(''  ERROR: MAXIMUM NUMBER OF STATE IS LIMITIED TO '',i4)') m
            Write(iout,'(''  ERROR: TERMINATED IN VSCF_CONSTRUCT '')')
            return
         endif

         i=spr_memalloc(1,dble(Nstate)*dble(Nfree)*4.D+00)
         if(i<0) return
         allocate(label(Nstate*Nfree))
         label=0

         Do i=1,Nstate*Nfree
            if(state(i)/=0) then 
               label=state(1:Nstate*Nfree)
               k=-1
               exit
            endif
         End do

         Do i1=1,Nstate
         Do i2=1,Nfree
            if(state_label(i2,i1)/=0) then

               if(k==-1) then
                  Write(iout,'(''  ERROR: BOTH state AND state_label ARE SPECIFIED.'')') 
                  Write(iout,'(''  ERROR: USE THE EITHER ONE TO SPECIFY THE STATES.'')') 
                  i=-1
                  return
               endif

               k=1
               Do j1=1,Nstate
               Do j2=1,Nfree
                  label(k)=state_label(j2,j1)
                  k=k+1
               End do
               End do
               exit

            endif
         End do
         End do

         k=1
         Do i1=1,Nstate
         Do i2=1,Nfree
            if(label(k)> nCHO(i2)-1) label(k)= nCHO(i2)-1
            if(label(k)< 0) label(k)=0
            k=k+1
         End do
         End do

         ! --+----2----+----3----+----4----+----5----+----6----+----7----+----80
         ! >> Output options

         write(iout,100)
         write(iout,105)
         write(iout,110)
         if(Nfree<=30) then
            Do j=1,Nstate
               write(iout,115) j,(label(k+(j-1)*Nfree),k=1,Nfree)
            End do
         else
            kk=Nfree/30
            Do j=1,Nstate
               write(iout,115) j,(label(k+(j-1)*Nfree),k=1,30)
               Do l=2,kk
                  write(iout,116) (label(k+(j-1)*Nfree),k=30*(l-1)+1,30*l)
               End do
               write(iout,116) (label(k+(j-1)*Nfree),k=30*kk+1,Nfree)
            End do
         endif
         if(lvpt .or. lvcc .or. lvci) then
            if(Ethresh > 1.D-06) Ethresh=1.D-06
         endif
         write(iout,130) restart,Maxitr, Ethresh

     100 Format(//,'(  ENTER VSCF MODULE  )',//)
     105 Format(2x,'---  VIBRATIONAL SELF-CONSISTENT FIELD CALCULATIONS   ---',//)
     110 Format(/,3x,'>> VIBRATIONAL STATES ')
     115 Format(7x,'STATE ',i3.3,':',3x,6(5i1,1x))
     116 Format(20x,6(5i1,1x))
     130 Format(/,3x,'>> SCF OPTION',/, &
         7x,'      RESTART :',l9,/, &
         7x,'MAX ITERATION :',i9,/, &
         7x,' CONV. THRESH :',e9.2)

         ! --+----2----+----3----+----4----+----5----+----6----+----7----+----80
         ! >> Data files

         if(restart) then
            rwf=20
            Do while(0==0)
               Inquire(rwf,opened=op)
               if(op) then
                  rwf=rwf+1
               else
                  exit
               endif
            End do
            Open(rwf,file='vscf-r.wfn',status='old',form='UNFORMATTED')
         endif

         wwf=30
         Do while(0==0)
            Inquire(wwf,opened=op)
            if(op) then
               wwf=wwf+1
            else
               exit
            endif
         End do
         Open(wwf,file='vscf-w.wfn',status='unknown',form='UNFORMATTED')

         write(iout,135)
     135 Format(/,3x,'>> FILES',/, &
                  7x,' WFN READ FILE : vscf-r.wfn',/, &
                  7x,'WFN WRITE FILE : vscf-w.wfn')

         ! --+----2----+----3----+----4----+----5----+----6----+----7----+----80
         ! >> Memory allocation

         j=idx2(Nfree) + nCHO(Nfree)*nCHO(Nfree)
         msz=dble(j)*8.D+00
         i=spr_memalloc(1,msz)
         if(i<0) return
         allocate(Tmat(j))

         Call setKinmat()

         j=idx1(Nfree) + nCHO(Nfree)
         msz=dble(j)*4.D+00
         i=spr_memalloc(1,msz)
         if(i<0) return
         allocate(type_Table1(j))
         Call Modal_gettypeTable1(type_Table1)
         i1=0
         Do i=1,Nfree
            Do j=1,nCHO(i)
               if(type_Table1(idx1(i)+j)<0) cycle
               !write(6,'(3i4)') i,j,i1
               type_Table1(idx1(i)+j)=i1
               i1=i1+j*2
            End do
         End do
         !write(6,'(3i4)') i1
         msz=dble(i1)*8.D+00
         i=spr_memalloc(1,msz)
         if(i<0) return
         allocate(block(i1))

         i=0
         write(iout,160)
         Call spr_meminfo
         Call timer(1)

     160 Format(//,'(  SETUP OF VSCF COMPLETED  )',//)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Vscf_finalz()

      USE Vscf_mod

      Implicit None

         Integer :: in,io
         Real(8) :: memsz

         Call spr_Getio(in,io)
         Write(io,100)
     100 Format(/,'(  FINALIZE VSCF MODULE  )',/)

         memsz=dble(size(label))*4.D+00  &
             + dble(size(Tmat))*8.D+00
         Call spr_memdealloc(memsz)
         Deallocate(label,Tmat)

         if(restart) Close(rwf)
         Close(wwf)

         Call spr_meminfo
         Call timer(1)

      End subroutine 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Vscf_main()
 
      USE Vscf_mod
 
      Implicit None
 
         Integer :: i,j,jj,k,kk,l,itr,nGj
         Integer :: in,io,maxn,ierr,spr_memalloc
 
         Real(8) :: Enew,Eold,delE,E0,Vscf_totE
         Real(8), parameter :: H2wvn=2.194746D+05
         Real(8), dimension(:), allocatable :: Ene
         Real(8), dimension(:,:), allocatable :: C,Hmat,Vmat
         Logical :: cnv
 
         Integer :: NP,vav_GetNP
         Logical :: r_v,vav_Getr_v
 
         Call spr_GetIO(in,io)
         maxn=maxval(nCHO)
 
         r_v=vav_Getr_v()
         NP=vav_GetNP()
 
         write(io,100)
         Do i=1,Nstate
 
            Nst=i
            if(Nfree<31) then
               Write(io,110)  i,label((i-1)*Nfree+1:i*Nfree)
            else
               Write(io,110)  i,label((i-1)*Nfree+1:(i-1)*Nfree+30)
               kk=Nfree/30
               Do k=2,kk
                  Write(io,111)  label((i-1)*Nfree+30*(k-1)+1:(i-1)*Nfree+30*k)
               End do
               Write(io,111)  label((i-1)*Nfree+30*kk+1:i*Nfree)
            endif
 
            if(restart) then
               Call Vscf_Guess(j)
               if(j==0) then
                  write(io,112)
               else
                  if(i==1) then 
                    write(io,113)
                  else
                    write(io,114)
                  endif
               endif
                  
            else
               if(i==1) then 
                 write(io,113)
               else
                 write(io,114)
               endif
 
            endif
 
            ! Start of VSCF LOOP
            cnv=.False.
            Eold=0.D+00
            itr=0
            Write(io,120)
            Do while(.not. cnv) 
 
               Call Vscf_update()
               Etot=Vscf_totE()
               Enew=Etot*H2wvn
               delE=Enew-Eold
               write(io,130) itr,Eold,Enew,delE
 
               if(abs(delE)>Ethresh) then
                  Eold=Enew
               else
                  cnv=.true.
                  exit
               endif
 
               itr=itr+1
               if(itr >= Maxitr+1) then 
                  write(io,1000) 
                  1000 Format(/,7x,' o MAX NUMBER OF ITERATION REACHED ...', &
                                   ' EXIT THE LOOP. ',/, &
                                7x,' o CURRENT VSCF STATE IS UNCONVERGED.')
                  exit
               endif
 
               Do j=1,Nfree
 
                  if(nCHO(j)==0) cycle

                  nGj=nCHO(j)
                  !dbg write(io,'(7x,''MODE='',i4,'', GRID='',i4)') j,nGj
                  ierr=spr_memalloc(0,dble(nGj*nGj)*24.D+00)
                  ierr=spr_memalloc(0,dble(nGj)*8.D+00)
                  if(ierr<0) return 
                  Allocate(C(nGj,nGj),Hmat(nGj,nGj),Vmat(nGj,nGj),Ene(nGj))
 
                  Ene=0.D+00
                  C=0.D+00
 
                  Call Vscf_getKinmat(j,Hmat)
                  Call Vscf_getVmat(j,Vmat)
                  Hmat=Hmat+Vmat
                  Call diag2(nGj,nGj,Hmat,C,Ene)
                  !Call huckeler(nGj,nGj,Hmat,Ene,C)
                  !Call diag0(nGj,nGj,Hmat,C,Ene)
                  Call Modal_setCwfn(j,C)
                  Call Modal_setEne(j,Ene)
 
                  !dbg write(io,'(9x,i4,f12.2)') 0,Ene(1)*H2wvn
                  !dbg write(io,'(9x,i4,f12.2)') (k-1,(Ene(k)-Ene(1))*H2wvn,k=2,11)
 
                  Call spr_memdealloc((dble(size(C))+dble(size(Hmat)) &
                                                    +dble(size(Vmat)))*8.D+00)
                  Call spr_memdealloc(dble(size(Ene))*8.D+00)
                  Deallocate(C,Hmat,Vmat,Ene)
 
               End do
 
            End do
!MK   writes fundamentals on fort.666 if summary=.t. @ &vscf >>>>>>>>>>>>         
            if (cnv .and. summary)then
            	if(i/=1) then 
               		write(666,'(f10.3)',advance='no')Enew-E0
            	else
          		write(666,666)ncho(1),MR,Maxitr, Ethresh
             666 Format('VSCF(cm^-1) for ncho(1)=',i2,' MR=',i1,' Maxitr=', i3, ' Ethresh=', e9.2)
          		write(666,'(f10.3)',advance='no')Enew
            	endif
            	
            else
             if (summary) write(666,*)"Not converged"
            endif
!MK <<<<<<<<<<<<<<<<<<<<<<<
            if(i/=1) then 
               		write(io,140) Enew,Enew-E0
            	else
               		write(io,145) Enew
             		E0=Enew
            endif	
 
            write(io,150)
            ierr=spr_memalloc(0,dble(maxn)*dble(Nfree)*8.D+00)
            if(ierr<0) return 
            Allocate(C(maxn,Nfree))
            Call Vscf_getConf(maxn,Nfree,C)
 
            jj=(Nfree-mod(Nfree,6))/6
            Do j=1,jj
               write(io,155) (k,k=(j-1)*6+1,j*6)
               Do k=1,maxn
                  write(io,156) k-1,(C(k,l),l=(j-1)*6+1,j*6)
               End do
            End do
            if(mod(Nfree,6)/=0) then
               write(io,155) (k,k=jj*6+1,Nfree)
               Do k=1,maxn
                  write(io,156) k-1,(C(k,l),l=jj*6+1,Nfree)
               End do
            endif
            Call spr_memdealloc(dble(size(C))*8.D+00)
            Deallocate(C)
 
            if(r_v) Call Vscf_ave()
            if(NP>0) Call Vscf_Prpt()
            Call Vscf_Dump()

         End do
         write(io,200)
   
     100 Format(/'(  ENTER VSCF MAIN MODULE  )',/ &
               /,3x,'[  LOOP OVER STATES  ]')
     110 Format(//,3x,'>> STATE ',i3.3,':',3x,6(5i1,1x))
     111 Format(19x,6(5i1,1x))
     112 Format(/,7x,'o RESTORE WFN FROM A PUNCH FILE')
     113 Format(/,7x,'o INITIAL GUESS FROM (CONTRACTED) HARMONIC OSCILLATOR ')
     114 Format(/,7x,'o USING THE PREVIOUS VSCF STATE AS INITAL WFN')
     120 Format(/,3x,'-- (ITERATION) ',8('-'),' (EOLD) ',8('-'), &
                    ' (ENEW) ',5('-'),' (DELTA E) --')
     130 Format(13x,i4,3f16.2)
     140 Format(3x,65('-'),/,7x,'E(VSCF)   =',f12.2,/,7x,'E(VSCF)-E0=',f12.2)
     145 Format(3x,65('-'),/,7x,'E(VSCF)=',f12.2)
     150 Format(/,7x,'o COEFFICIENTS')
     155 Format(/,9x,'-MODES-',1x,6i10)
     156 Format(10x,i5,2x,6f10.5)
     200 Format(//,3x,'[  END OF LOOP OVER STATES  ]',/ &
                //,'(  EXIT VSCF MAIN MODULE  )',/)

      End subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_getConf(maxn,Nf,C)
 
      USE Vscf_mod
 
      Implicit None
 
         Integer :: maxn,Nf
         Real(8) :: C(maxn,Nf)
 
         Integer :: i,j,n
         Integer :: lb(Nf)
         Real(8), dimension(:,:), allocatable :: Cwfn
 
         lb=label((Nst-1)*Nfree+1:Nst*Nfree)
         !dbg write(6,*) lb
 
         Do i=1,Nfree
            n=nCHO(i)
            Allocate(Cwfn(n,n))
            Call Modal_getCwfn(i,Cwfn)
            Do j=1,n
               C(j,i)=Cwfn(j,lb(i)+1)
            End do
            Do j=n+1,maxn
               C(j,i)=0.D+00
            End do
            Deallocate(Cwfn)
         End do
 
      End subroutine
 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Vscf_Dump()
 
      USE Vscf_mod
 
      Implicit None
 
         Integer :: i,j,k,l, nGi
         Real(8), allocatable :: Cwfn(:,:),Ene(:)
         Character :: wftyp*4
         Real(8), parameter :: H2wvn=2.194746D+05
 
         wftyp='VSCF'
         Write(wwf) wftyp
 
         Write(wwf) label((Nst-1)*Nfree+1:Nst*Nfree)
         Do i=1,Nfree
            nGi=nCHO(i)
            Write(wwf) nGi
            if(nCHO(i)==0) cycle
            Allocate(Cwfn(nGi,nGi),Ene(nGi))
            Call Modal_getCwfn(i,Cwfn)
            Call Modal_getEne(i,Ene)
            Write(wwf) Cwfn
            Write(wwf) Ene
            !dbg Write(6,'(3x,i4,100f8.1)') i,Ene*H2wvn
            Deallocate(Cwfn,Ene)
         End do
 
      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Function Vscf_totE()
 
      USE Vscf_mod
 
      Implicit None
 
         Real(8) :: Vscf_totE
 
         Real(8), parameter :: H2wvn=2.194746D+05
 
         Integer :: lb(Nfree)
 
         Integer :: i,j,k,ij,ijk,idx,jdx,int2
         Real(8) :: Ekin,Epot,omgh,tmp,tmp1,tmp2
 
         Integer :: n,nGi,nGj,nGk,mi,mj,mk
         Real(8), allocatable :: qav(:,:),qi(:,:),qj(:,:),qk(:,:)
         Real(8), parameter :: thresh=1.D-06
 
         Real(8), allocatable :: Cwfn(:,:),Cw(:),Tm(:,:)
         Real(8), allocatable :: Di(:),Dj(:),Dk(:)
         Real(8), allocatable :: Vi(:),Vj(:),Vk(:),Vij(:,:),Vijk(:,:,:)
 
         Integer :: spr_memalloc
 
         ! Label of the current state
         lb=label((Nst-1)*Nfree+1:Nst*Nfree)

         Ekin=0.D+00
         Do i=1,Nfree
            nGi=nCHO(i)
            Allocate(Cwfn(nGi,nGi),Cw(nGi),Tm(nGi,nGi))

            Call Modal_getCwfn(i,Cwfn)
            Cw=Cwfn(:,lb(i)+1)

            Call Vscf_getKinmat(i,Tm)

            Do j=1,nGi
               tmp=0.D+00
               Do k=1,nGi
                  tmp=tmp + Tm(j,k)*Cw(k)
               End do
               Ekin=Ekin + Cw(j)*tmp
            End do

            Deallocate(Cwfn,Cw,Tm)
         End do

         Epot=0.D+00
         Call Vscf_clearVmf()

         ! Potential energy and mean-field potential

         ! QFF
         if(nQ1 /=0 .or. nQ2 /=0 .or. nQ3 /=0) then
            i=spr_memalloc(-1,dble(Nfree)*4.D+00*8.D+00)
            Allocate(qav(4,Nfree))

            Do mi=1,Nfree
               nGi=nCHO(mi)
               Allocate(Di(nGi),qi(nGi,4))
               Call Vscf_getDwfn(mi,nGi,Di)
               Call Modal_getQ4(mi,qi)

               ! <Qi>
               tmp=0.D+00
               Do i=1,nGi
                  tmp=tmp + Di(i)*qi(i,1)
               End do
               qav(1,mi)=tmp
               ! <Qi^2>
               tmp=0.D+00
               Do i=1,nGi
                  tmp=tmp + Di(i)*qi(i,2)
               End do
               qav(2,mi)=tmp
               ! <Qi^3>
               tmp=0.D+00
               Do i=1,nGi
                  tmp=tmp + Di(i)*qi(i,3)
               End do
               qav(3,mi)=tmp
               ! <Qi^4>
               tmp=0.D+00
               Do i=1,nGi
                  tmp=tmp + Di(i)*qi(i,4)
               End do
               qav(4,mi)=tmp

               Deallocate(Di,qi)
            End do
         endif

         Do n=1,nQ1
            mi=mQ1(n)
            nGi=nCHO(mi)
            Allocate(Vi(nGi),qi(nGi,4))
            Call Modal_getQ4(mi,qi)

            Epot=Epot + qav(1,mi)*coeff1(1,n) &
                      + qav(2,mi)*coeff1(2,n) &
                      + qav(3,mi)*coeff1(3,n) &
                      + qav(4,mi)*coeff1(4,n) 
            Do i=1,nGi
               Vi(i)=   qi(i,1)*coeff1(1,n) &
                      + qi(i,2)*coeff1(2,n) &
                      + qi(i,3)*coeff1(3,n) &
                      + qi(i,4)*coeff1(4,n) 
            End do

            Call Vscf_addVmf(mi,nGi,Vi)
            Deallocate(Vi,qi)

         End do

         Do n=1,nS1
            mi=mS1(n)
            Call PEGrid_getnGrid1(nGi,n)
            Allocate(Vi(nGi),Di(nGi))

            Call PEGrid_getV1(n,Vi)
            Call Vscf_getDwfn(mi,nGi,Di)
            Call Vscf_addVmf(mi,nGi,Vi)
            tmp=0.D+00
            Do i=1,nGi
               tmp=tmp + Di(i)*Vi(i)
            End do
            Epot=Epot+tmp

            Deallocate(Vi,Di)

         End do

         if(MR==1) goto 999

         Do n=1,nQ2
            mi=mQ2(1,n)
            mj=mQ2(2,n)
            nGi=nCHO(mi)
            nGj=nCHO(mj)
            Allocate(Vi(nGi),Vj(nGj))
            Allocate(qi(nGi,4),qj(nGj,4))
            Call Modal_getQ4(mi,qi)
            Call Modal_getQ4(mj,qj)

            Epot=Epot + qav(1,mi)*qav(1,mj)*coeff2(1,n) &
                      + qav(2,mi)*qav(1,mj)*coeff2(2,n) &
                      + qav(1,mi)*qav(2,mj)*coeff2(3,n) &
                      + qav(2,mi)*qav(2,mj)*coeff2(4,n) &
                      + qav(3,mi)*qav(1,mj)*coeff2(5,n) &
                      + qav(1,mi)*qav(3,mj)*coeff2(6,n) 

            Do i=1,nGi
               Vi(i)= + qi(i,1)*qav(1,mj)*coeff2(1,n) &
                      + qi(i,2)*qav(1,mj)*coeff2(2,n) &
                      + qi(i,1)*qav(2,mj)*coeff2(3,n) &
                      + qi(i,2)*qav(2,mj)*coeff2(4,n) &
                      + qi(i,3)*qav(1,mj)*coeff2(5,n) &
                      + qi(i,1)*qav(3,mj)*coeff2(6,n) 
            End do

            Do j=1,nGj
               Vj(j)= + qav(1,mi)*qj(j,1)*coeff2(1,n) &
                      + qav(2,mi)*qj(j,1)*coeff2(2,n) &
                      + qav(1,mi)*qj(j,2)*coeff2(3,n) &
                      + qav(2,mi)*qj(j,2)*coeff2(4,n) &
                      + qav(3,mi)*qj(j,1)*coeff2(5,n) &
                      + qav(1,mi)*qj(j,3)*coeff2(6,n) 
            End do
            Call Vscf_addVmf(mi,nGi,Vi)
            Call Vscf_addVmf(mj,nGj,Vj)
            Deallocate(Vi,Vj,qi,qj)

         End do

         Do n=1,nS2
            mi=mS2(1,n)
            mj=mS2(2,n)
            Call PEGrid_getnGrid2(nGi,nGj,n)
            Allocate(Vij(nGj,nGi),Vi(nGi),Vj(nGj),Di(nGi),Dj(nGj))
            Call PEGrid_getV2(n,Vij)
            Call Vscf_getDwfn(mi,nGi,Di)
            Call Vscf_getDwfn(mj,nGj,Dj)

            tmp=0.D+00
            Do i=1,nGi
               Vi(i)=0.D+00
               Do j=1,nGj
                  Vi(i)=Vi(i) + Dj(j)*Vij(j,i)
               End do
               tmp=tmp + Di(i)*Vi(i)
            End do
            Epot=Epot+tmp

            Do j=1,nGj
               Vj(j)=0.D+00
               Do i=1,nGi
                  Vj(j)=Vj(j) + Di(i)*Vij(j,i)
               End do
            End do

            Call Vscf_addVmf(mi,nGi,Vi)
            Call Vscf_addVmf(mj,nGj,Vj)
            Deallocate(Vij,Vi,Vj,Di,Dj)

         End do

         if(MR==2) goto 999

         Do n=1,nQ3
            mi=mQ3(1,n)
            mj=mQ3(2,n)
            mk=mQ3(3,n)
            nGi=nCHO(mi)
            nGj=nCHO(mj)
            nGk=nCHO(mk)
            Allocate(Vi(nGi),Vj(nGj),Vk(nGk))
            Allocate(qi(nGi,4),qj(nGj,4),qk(nGk,4))
            Call Modal_getQ4(mi,qi)
            Call Modal_getQ4(mj,qj)
            Call Modal_getQ4(mk,qk)

            Epot=Epot + qav(1,mi)*qav(1,mj)*qav(1,mk)*coeff3(1,n) &
                      + qav(2,mi)*qav(1,mj)*qav(1,mk)*coeff3(2,n) &
                      + qav(1,mi)*qav(2,mj)*qav(1,mk)*coeff3(3,n) &
                      + qav(1,mi)*qav(1,mj)*qav(2,mk)*coeff3(4,n) 

            Do i=1,nGi
               Vi(i) =  qi(i,1)*qav(1,mj)*qav(1,mk)*coeff3(1,n) &
                      + qi(i,2)*qav(1,mj)*qav(1,mk)*coeff3(2,n) &
                      + qi(i,1)*qav(2,mj)*qav(1,mk)*coeff3(3,n) &
                      + qi(i,1)*qav(1,mj)*qav(2,mk)*coeff3(4,n) 
            End do

            Do j=1,nGj
               Vj(j)=   qav(1,mi)*qj(j,1)*qav(1,mk)*coeff3(1,n) &
                      + qav(2,mi)*qj(j,1)*qav(1,mk)*coeff3(2,n) &
                      + qav(1,mi)*qj(j,2)*qav(1,mk)*coeff3(3,n) &
                      + qav(1,mi)*qj(j,1)*qav(2,mk)*coeff3(4,n) 
            End do

            Do k=1,nGk
               Vk(k)=   qav(1,mi)*qav(1,mj)*qk(k,1)*coeff3(1,n) &
                      + qav(2,mi)*qav(1,mj)*qk(k,1)*coeff3(2,n) &
                      + qav(1,mi)*qav(2,mj)*qk(k,1)*coeff3(3,n) &
                      + qav(1,mi)*qav(1,mj)*qk(k,2)*coeff3(4,n) 
            End do

            Call Vscf_addVmf(mi,nGi,Vi)
            Call Vscf_addVmf(mj,nGj,Vj)
            Call Vscf_addVmf(mk,nGk,Vk)
            Deallocate(Vi,Vj,Vk,qi,qj,qk)

         End do

         Do n=1,nS3
            mi=mS3(1,n)
            mj=mS3(2,n)
            mk=mS3(3,n)
            Call PEGrid_getnGrid3(nGi,nGj,nGk,n)
            Allocate(Vijk(nGk,nGj,nGi),Vi(nGi),Vj(nGj),Vk(nGk))
            Allocate(Vij(nGj,nGi))
            Allocate(Di(nGi),Dj(nGj),Dk(nGk))
            Call PEGrid_getV3(n,Vijk)
            Call Vscf_getDwfn(mi,nGi,Di)
            Call Vscf_getDwfn(mj,nGj,Dj)
            Call Vscf_getDwfn(mk,nGk,Dk)

            Do i=1,nGi
            Do j=1,nGj
               Vij(j,i)=0.D+00
               Do k=1,nGk
                  Vij(j,i)=Vij(j,i) + Dk(k)*Vijk(k,j,i)
               End do
            End do
            End do

            Do i=1,nGi
               Vi(i)=0.D+00
               Do j=1,nGj
                  Vi(i)=Vi(i) + Dj(j)*Vij(j,i)
               End do
            End do

            Do j=1,nGj
               Vj(j)=0.D+00
               Do i=1,nGi
                  Vj(j)=Vj(j) + Di(i)*Vij(j,i)
               End do
            End do

            Do k=1,nGk
               Vk(k)=0.D+00
               Do i=1,nGi
                  tmp=0.D+00
                  Do j=1,nGj
                     tmp=tmp + Dj(j)*Vijk(k,j,i)
                  End do
                  Vk(k)=Vk(k) + tmp*Di(i)
            End do
            End do

            tmp=0.D+00
            Do i=1,nGi
               tmp=tmp + Di(i)*Vi(i)
            End do
            Epot=Epot + tmp

            Call Vscf_addVmf(mi,nGi,Vi)
            Call Vscf_addVmf(mj,nGj,Vj)
            Call Vscf_addVmf(mk,nGk,Vk)
            Deallocate(Vijk,Vij,Vi,Vj,Vk,Di,Dj,Dk)

         End do

         if(MR==3) goto 999

     999 Continue

         Vscf_totE=Ekin+Epot
         if(nQ1 /=0 .or. nQ2 /=0 .or. nQ3 /=0) then
            Call spr_memdealloc(dble(Nfree)*4.D+00*8.D+00)
            Deallocate(qav)
         endif
         return
 
      End Function
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!    --- Initial guess of one mode function ---
!
      Subroutine Vscf_Guess(istat)
 
      USE Vscf_mod
 
      Implicit None
 
         Integer :: inp,iout
 
         Integer :: i,istat,nGi,lb(Nfree)
         Real(8), allocatable :: Cwfn(:,:),Ene(:)
 
         Character :: wftyp*4
 
         istat=0
         Read(rwf,end=100) wftyp
         if(wftyp=='VSCF') then
            ! --> VSCF wfn <--
            Read(rwf) lb
 
            Do i=1,Nfree
               Read(rwf) nGi
               Allocate(Cwfn(nGi,nGi),Ene(nGi))
               Read(rwf) Cwfn
               Read(rwf) Ene
               if(nGi/=nCHO(i)) then
                  Call spr_Getio(inp,iout)
                  Write(iout,*)
                  Write(iout,*) 'FOR MODE          :',i
                  Write(iout,*) 'INPUT             :',nCHO(i)
                  Write(iout,*) 'FILE (vscf-r.wfn) :',nGi
                  Write(iout,*)
                  Write(iout,*) 'CHANGING THE NUM. OF BASIS SET IS NOT ALLOWED AT PRESENT'
                  Write(iout,*)
                  Stop
                  
               endif
               Call Modal_setCwfn(i,Cwfn)
               Call Modal_setEne(i,Ene)
               Deallocate(Cwfn,Ene)
            End do
 
         else
            Call spr_Getio(inp,iout)
 
            ! --> TDH wfn <--
            Write(iout,*) 'vscf-r.wfn CONTAINS TDH COMPLEX WFN'
            Write(iout,*) 'THIS OPTION IS NOT READY'
            Stop
         endif
         return
 
     100 Continue
         istat=-1
 
         return
 
      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_getKinmat(mode,Tm)
 
      USE Vscf_mod
 
      Implicit None
 
         Integer :: mode
         Real(8) :: Tm(nCHO(mode)*nCHO(mode))
 
         Tm=Tmat(idx2(mode)+1:idx2(mode)+nCHO(mode)*nCHO(mode))
 
      End subroutine
 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_setDwfn(mode,nGrid,Dwfn)

      USE Vscf_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: Dwfn(nGrid)

         pt=type_Table1(idx1(mode)+nGrid)+nGrid
         block(pt+1:pt+nGrid)=Dwfn

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_getDwfn(mode,nGrid,Dwfn)

      USE Vscf_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: Dwfn(nGrid)

         pt=type_Table1(idx1(mode)+nGrid)+nGrid
         Dwfn=block(pt+1:pt+nGrid)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_update()

      USE Vscf_mod

      Implicit None

         Integer :: i,j,k,l
         Integer :: nG,nGi,lb(Nfree)
         Real(8) :: xx
         Real(8), dimension(:,:), allocatable :: Cwfn,xdvr
         Real(8), dimension(:), allocatable :: Dwfn

         lb=label(Nfree*(Nst-1)+1:Nfree*Nst)+1
         Do i=1,Nfree

            nGi=nCHO(i)
            Allocate(Cwfn(nGi,nGi),xdvr(nGi,nGi),Dwfn(nGi))
            Call Modal_getCwfn(i,Cwfn)
            Call Modal_getxdvr(i,nGi,xdvr)

            Do j=1,nGi
               Dwfn(j)=0.D+00
               Do k=1,nGi
                  Dwfn(j)=Dwfn(j) + xdvr(k,j)*Cwfn(k,lb(i))
               End do
               Dwfn(j)=Dwfn(j)*Dwfn(j)
            End do
            Call Vscf_setDwfn(i,nGi,Dwfn)

            !dbg Call Vscf_getDwfn(i,nGi,Dwfn)
            !dbg write(6,'(2i4)') i,nGi
            !dbg write(6,'(11f8.4)') Dwfn

            Deallocate(Dwfn,xdvr)

            Do j=1,nGi-1
               if(type_Table1(idx1(i)+j)<0) cycle

               nG=j
               Allocate(xdvr(nG,nG),Dwfn(nG))
               Call Modal_getxdvr(i,nG,xdvr)

               Do k=1,nG
                  Dwfn(k)=0.D+00
                  Do l=1,nG
                     Dwfn(k)=Dwfn(k) + xdvr(l,k)*Cwfn(l,lb(i))
                  End do
                  Dwfn(k)=Dwfn(k)*Dwfn(k)
               End do
               Call Vscf_setDwfn(i,nG,Dwfn)

               !dbg Call Vscf_getDwfn(i,nG,Dwfn)
               !dbg write(6,'(2i4)') i,j
               !dbg write(6,'(11f8.4)') Dwfn

               Deallocate(Dwfn,xdvr)

            End do

            Deallocate(Cwfn)

         end do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_setVmf(mode,nGrid,Vmf)

      USE Vscf_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: Vmf(nGrid)

         pt=type_Table1(idx1(mode)+nGrid)
         block(pt+1:pt+nGrid)=Vmf

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_getVmf(mode,nGrid,Vmf)

      USE Vscf_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: Vmf(nGrid)

         pt=type_Table1(idx1(mode)+nGrid)
         Vmf=block(pt+1:pt+nGrid)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_addVmf(mode,nGrid,Vmf)

      USE Vscf_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: Vmf(nGrid)

         pt=type_Table1(idx1(mode)+nGrid)
         block(pt+1:pt+nGrid)=block(pt+1:pt+nGrid)+Vmf

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_clearVmf()

      USE Vscf_mod

      Implicit None

         Integer :: i,j,k,l, mm,nG
         Real(8), allocatable :: Vmf(:)

         Do i=1,Nfree
         Do j=1,nCHO(i)
            if(type_Table1(idx1(i)+j)<0) cycle
            mm=i
            nG=j

            Allocate(Vmf(nG))
            Vmf=0.D+00
            Call Vscf_setVmf(mm,nG,Vmf)
            Deallocate(Vmf)

         End do
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vscf_getVmat(mm,Vmat)

      USE Vscf_mod

      Implicit None

         Integer :: pt,mm,i,j,k,l,nG,nGi
         Real(8) :: Vmat(nCHO(mm),nCHO(mm))
         Real(8), allocatable :: Vmf(:),xdvr(:,:),Vms(:,:)

         Vmat=0.D+00
         nG=nCHO(mm)

         Allocate(Vmf(nG),xdvr(nG,nG))
         Call Vscf_getVmf(mm,nG,Vmf)
         Call Modal_getxdvr(mm,nG,xdvr)
         !dbg write(6,'(i4)') mm
         !dbg write(6,'(f8.4)') Vmf
         !dbg write(6,*)
         !dbg Do j=1,nG
         !dbg    write(6,'(11f8.4)') xdvr(:,j)
         !dbg End do
         !dbg write(6,*)

         Do j=1,nG
            Do k=1,j-1
               Vmat(j,k)=0.D+00
               Do l=1,nG
                  Vmat(j,k)=Vmat(j,k) + xdvr(j,l)*Vmf(l)*xdvr(k,l)
               End do
               Vmat(k,j)=Vmat(j,k)
            End do

            Vmat(j,j)=0.D+00
            Do l=1,nG
               Vmat(j,j)=Vmat(j,j) + xdvr(j,l)*Vmf(l)*xdvr(j,l)
            End do
         End do
         Deallocate(Vmf,xdvr)

         Do i=1,nG-1
            pt=idx1(mm)+i
            if(type_Table1(pt)<0) cycle

            nGi=i
            Allocate(Vmf(nGi),xdvr(nGi,nGi),Vms(nGi,nGi))
            Call Vscf_getVmf(mm,nGi,Vmf)
            Call Modal_getxdvr(mm,nGi,xdvr)

            Do j=1,nGi
               Do k=1,j-1
                  Vms(j,k)=0.D+00
                  Do l=1,nGi
                     Vms(j,k)=Vms(j,k) + xdvr(j,l)*Vmf(l)*xdvr(k,l)
                  End do
                  Vms(k,j)=Vms(j,k)
               End do

               Vms(j,j)=0.D+00
               Do l=1,nGi
                  Vms(j,j)=Vms(j,j) + xdvr(j,l)*Vmf(l)*xdvr(j,l)
               End do
            End do

            Do j=1,nGi
               Do k=1,j-1
                  Vmat(k,j)=Vmat(k,j)+Vms(k,j)
                  Vmat(j,k)=Vmat(k,j)
               End do
               Vmat(j,j)=Vmat(j,j)+Vms(j,j)
            End do
            Deallocate(Vmf,xdvr,Vms)

         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Vscf_ave

      USE Vscf_mod

      Implicit None

         Integer :: ie,spr_memalloc, Nat,spr_getNat
         Integer :: io,i,j,k,nGj
         Double precision :: A0
         Real(8), allocatable :: QQ(:),x0(:,:),Dwfn(:),qj(:),R0(:,:)
 
         Call spr_GetIO(j,io)

         ie=spr_memalloc(1,dble(Nfree)*8.D+00)
         Allocate(QQ(Nfree))

         write(io,100)
         QQ=0.D+00
         Do j=1,Nfree
            nGj=nCHO(j)
            Allocate(Dwfn(nGj),qj(nGj))
            Call Modal_getQ(j,nGj,qj)
            Call Vscf_getDwfn(j,nGj,Dwfn)
            Do k=1,nGj
               QQ(j)=QQ(j)+Dwfn(k)*qj(k)
            End do
            Deallocate(Dwfn,qj)
            !Write(io,'(i10)') nho(j)
            !Write(io,'(11f10.3)') wfn(:,j) 
            !Write(io,'(11f10.3)') quad(:,j)
         End do
         Write(io,200) QQ

         Write(io,110)
         Nat=spr_GetNat()
         ie=spr_memalloc(1,dble(Nat)*3.D+00*8.D+00)
         Allocate(x0(3,Nat),R0(Nat,Nat))
         Call vav_q2x(Nat,Nfree,x0,QQ)
         Write(io,200) x0
!MK
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!     
      R0=0.0;   
      Do j=1,Nat
       Do k=j+1,Nat
        R0(j,k)=sqrt((x0(1,j)-x0(1,k))**2+(x0(2,j)-x0(2,k))**2+(x0(3,j)-x0(3,k))**2)
        R0(k,j)=R0(j,k)
       End Do
      End Do
      Write(io,120)
      Write(io,200) R0
      Write(io,130)
!MK 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!     
      A0=0.0;   
      Do i=1,Nat
       Do j=1,Nat
         if(i==j)cycle
        Do k=1,Nat
         if(j==k .or. i==k .or. k>j)cycle 
        A0=acos((R0(i,j)**2.d0+R0(i,k)**2.d0-R0(k,j)**2.d0)/(2*R0(i,j)*R0(i,k)))*180.0D0/pi
        Write(6,300)j,i,k,A0
        End Do
       End Do
      End Do
!MK
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
         Call spr_memdealloc(dble(Nat)*3.D+00*8.D+00)
         Deallocate(x0,R0)

         Call spr_memdealloc(dble(Nfree)*8.D+00)
         Deallocate(QQ)

     100 Format(/,7x,'o VIBRATIONALLY AVERAGED STRUCTURE',//,7x,'  - Q0 -')
     110 Format(/,7x,'  - X0 -')
     120 Format(/,7x,'  - R0 -')
     130 Format(/,10x,'  - A0 -')
     200 Format(10x,3f12.6)
     300 Format(12x,3i2,' ^ ',f9.4)
     

      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Vscf_Prpt()

      USE vscf_mod
      USE vav_mod

      Implicit None

      Integer :: in,io
      Integer :: i,j,k,nGi

      Real(8), dimension(NP) :: P
      Real(8), dimension(maxCHO) :: Di

         P=0.D+00
         Do i=1,Nfree
            nGi=nCHO(i)
            if(nGi==0) cycle

            Call Vscf_getDwfn(i,nGi,Di)

            Do j=1,NP
               Do k=1,nGi
                   P(j)=P(j)+Di(k)*P1(k,j,i)
               End do
            End do
         End do

         if(MR2==1) goto 999

         !=========================
         !  MR=2,3,.. is not ready
         !=========================

     999 Continue

         Call spr_GetIO(in,io)
         Write(io,100)
         Do i=1,NP
            Write(io,200) i,P0(i),P(i),P0(i)+P(i)
         End do
     100 Format(/,7x,'o VIBRATIONAL AVERAGE OVER PROPERTY SURFACES',/, &
                  7x,'     [i]      P0             dP             <P>')
     200 Format(12x,i2,1x,3f15.6)


      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
