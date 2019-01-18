!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/11/16
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Vcor_Construct(ierr)

      USE Vcor_mod

      Implicit None


         Integer :: ierr,spr_memalloc
         Integer :: i,j,k,l,i1
         Real(8) :: msz

         ierr=0

         ! --+----2----+----3----+----4----+----5----+----6----+----7----+----80
         ! >> VSCF configurations

         npCUP=0
         nqCUP=0
         npCnf=0
         nqCnf=0

         ierr=spr_memalloc(1,dble(Nfree)*8.D+00)
         if(ierr<0) return
         Allocate(ref(Nfree),tar(Nfree))

         ! --+----2----+----3----+----4----+----5----+----6----+----7----+----80
         ! >> Kinetic energy matrix

         j=idx2(Nfree) + nCHO(Nfree)*nCHO(Nfree)
         msz=dble(j)*8.D+00
         ierr=spr_memalloc(1,msz)
         if(ierr<0) return
         Allocate(Tmat(j))

         ! --+----2----+----3----+----4----+----5----+----6----+----7----+----80
         ! >> DVR coefficients and weights

         j=idx1(Nfree) + nCHO(Nfree)
         msz=dble(j)*4.D+00
         ierr=spr_memalloc(1,msz)
         if(ierr<0) return
         Allocate(type_Table1(j))
         Call Modal_gettypeTable1(type_Table1)
         i1=0
         Do i=1,Nfree
            Do j=1,nCHO(i)
               if(type_Table1(idx1(i)+j)<0) cycle
               !dbg write(6,'(3i4)') i,j,i1
               type_Table1(idx1(i)+j)=i1
               i1=i1+j*nCHO(i)*2
            End do
         End do
         !dbg write(6,*) 'DEBUG ',i1
         msz=dble(i1)*8.D+00
         ierr=spr_memalloc(1,msz)
         if(ierr<0) return
         Allocate(block(i1))

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setTarget(state,state_label)

      USE Vcor_mod

      Implicit None

         Integer, dimension(Nstate*Nfree):: state
         Integer, dimension(Nfree,Nstate):: state_label

         Integer :: i,j,k,i1,i2,j1,j2,spr_memalloc

            i=spr_memalloc(1,dble(Nstate)*dble(Nfree)*4.D+00)
            if(i<0) return
            Allocate(label0(Nstate*Nfree))
            label0=0

            Do i=1,Nstate*Nfree
               if(state(i)/=0) then 
                  label0=state(1:Nstate*Nfree)
                  k=-1
                  exit
               endif
            End do

            Do i1=1,Nstate
            Do i2=1,Nfree
               if(state_label(i2,i1)/=0) then

                  if(k==-1) then
                     Write(6,'(''  ERROR: BOTH state AND state_label ARE SPECIFIED.'')') 
                     Write(6,'(''  ERROR: USE THE EITHER ONE TO SPECIFY THE STATES.'')') 
                     Stop
                  endif

                  k=1
                  Do j1=1,Nstate
                  Do j2=1,Nfree
                     label0(k)=state_label(j2,j1)
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
               if(label0(k)> nCHO(i2)-1) label0(k)= nCHO(i2)-1
               if(label0(k)< 0) label0(k)=0
               k=k+1
            End do
            End do


      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_finalz()

      USE Vcor_mod

      Implicit None

         Call spr_memdealloc(dble(size(ref))*4.D+00)
         Deallocate(ref)
         Call spr_memdealloc(dble(size(tar))*4.D+00)
         Deallocate(tar)
         Call spr_memdealloc(dble(size(Tmat))*8.D+00)
         Deallocate(Tmat)
         Call spr_memdealloc(dble(size(type_Table1))*4.D+00)
         Deallocate(type_Table1)
         Call spr_memdealloc(dble(size(block))*8.D+00)
         Deallocate(block)
         if(allocated(label0)) then
            Call spr_memdealloc(dble(size(label0))*4.D+00)
            Deallocate(label0)
         endif

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setKinmat()

      USE Vcor_mod

      Implicit None 

         Integer :: i
         Real(8), allocatable :: Tm(:)

         Do i=1,Nfree
            if(nCHO(i)==0) cycle
            Allocate(Tm(nCHO(i)*nCHO(i)))
            Call genKinmat(i,Tm)
            Tmat(idx2(i)+1:idx2(i)+nCHO(i)*nCHO(i))=Tm
            Deallocate(Tm)
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vcor_getKinmat(mode,mm,nn)
 
      USE Vcor_mod
 
      Implicit None
 
         Integer :: mode,mm,nn
         Real(8) :: Vcor_getKinmat

         if(nCHO(mode)/=0) then
            Vcor_getKinmat=Tmat(idx2(mode)+nCHO(mode)*mm+nn+1)
         else
            Vcor_getKinmat=0.D+00
         endif
 
      End Function
 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setWfn()

      USE Vcor_mod

      Implicit None

         Integer :: i,j,k,l,i1,nGi,nG
         Real(8) :: xx
         Real(8), dimension(:,:), allocatable :: Cwfn,Dwfn,Ewfn,xdvr

         Do i=1,Nfree
            if(nCHO(i)==0) cycle
            nGi=nCHO(i)
            Allocate(Cwfn(nGi,nGi),Dwfn(nGi,nGi),Ewfn(nGi,nGi),xdvr(nGi,nGi))
            Call Modal_getCwfn(i,Cwfn)
            Call Modal_getxdvr(i,nGi,xdvr)
            Do j=1,nGi
               Do k=1,nGi
                  Ewfn(k,j)=0.D+00
                  Do l=1,nGi
                     Ewfn(k,j)=Ewfn(k,j) + xdvr(l,k)*Cwfn(l,j)
                     Dwfn(k,j)=Ewfn(k,j)*Ewfn(k,j)
                  End do
               End do
            End do
            Call Vcor_setDwfn(i,nGi,Dwfn)
            Call Vcor_setEwfn(i,nGi,Ewfn)

            Deallocate(Dwfn)
            Deallocate(Ewfn)
            Deallocate(xdvr)

            Do i1=1,nGi-1
               if(type_Table1(idx1(i)+i1)<0) cycle

               nG=i1
               Allocate(Dwfn(nG,nGi),Ewfn(nG,nGi),xdvr(nG,nG))
               Call Modal_getxdvr(i,nG,xdvr)

               Do j=1,nGi
                  Do k=1,nG
                     Ewfn(k,j)=0.D+00
                     Do l=1,nG
                        Ewfn(k,j)=Ewfn(k,j) + xdvr(l,k)*Cwfn(l,j)
                        Dwfn(k,j)=Ewfn(k,j)*Ewfn(k,j)
                     End do
                  End do
               End do
               Call Vcor_setDwfn(i,nG,Dwfn)
               Call Vcor_setEwfn(i,nG,Ewfn)

               Deallocate(xdvr,Dwfn,Ewfn)

            End do

            Deallocate(Cwfn)

         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setQmat()

      USE Vcor_mod

      Implicit None

         Integer :: ier,spr_memalloc
         Integer :: i,j,q1,q2,qi,nGi

         Real(8) :: Di(maxCHO),tmp
         Real(8), allocatable :: qdi(:,:)

         if(nQ1==0 .and. nQ2==0 .and. nQ3==0) return

         ier=spr_memalloc(-1,dble(maxCHO*maxCHO*Nfree*4)*8.D+00)
         Allocate(Qmat(4,0:maxCHO-1,0:maxCHO-1,Nfree))

         Do i=1,Nfree
            if(nCHO(i)==0) cycle
            nGi=nCHO(i)
            Allocate(qdi(nGi,4))
            Call Modal_getQ4(i,qdi)

            Do q1=0,nCHO(i)-1
               Call Vcor_getDwfn(i,nGi,q1,Di)

               ! <Qi>,<Qi^2>,<Qi^3>,<Qi^4>
               Do j=1,4
                  tmp=0.D+00
                  Do qi=1,nGi
                     tmp=tmp + Di(qi)*qdi(qi,j)
                  End do
                  Qmat(j,q1,q1,i)=tmp
               End do

               Do q2=0,q1-1
                  Call Vcor_getXwfn(i,nGi,q1,q2,Di)
                  ! <Qi>,<Qi^2>,<Qi^3>,<Qi^4>
                  Do j=1,4
                     tmp=0.D+00
                     Do qi=1,nGi
                        tmp=tmp + Di(qi)*qdi(qi,j)
                     End do
                     Qmat(j,q2,q1,i)=tmp
                     Qmat(j,q1,q2,i)=tmp
                  End do

               End do
            End do

            Deallocate(qdi)

         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine Vcor_set_pCnf()

      USE Vcor_mod

      Implicit None

         Integer, parameter :: maxrs=1000,maxCUP0=20
         Integer ::maxCUP

         Integer :: i,j,k,l,irs,nRs,spr_memalloc
         Real(8) :: Er0,Er1,De0
         Integer :: tar0(Nfree),tar1(Nfree),nRs0(2)
         Integer, allocatable :: cup(:),mode(:,:),vv(:,:)
         Real(8), allocatable :: Ene1(:,:),x0(:)
         Integer :: npCI,jmax
         Real(8) :: Vcor_getHmat1,Vcor_getHmat0,aa
         Real(8), allocatable :: Hmat(:),Ciw(:,:),Cie(:)
         Integer, allocatable :: delCnf(:),iCnf(:),jCnf(:)
         Real(8), parameter :: H2wvn=2.194746D+05,dd=1.D-06

         Allocate(Ene1(0:maxCHO-1,MR))
         if(maxpCUP>0) then
            maxCUP=maxpCUP
         else
            maxCUP=maxCUP0
            if(maxCUP>Nfree) maxCUP=Nfree
         endif
         Allocate(cup(0:maxrs),mode(maxCUP,0:maxrs),vv(maxCUP,0:maxrs),x0(0:maxrs))

         nRs=0
         tar0=tar 
         De0=0.D+00
         !dbg write(6,'(''***'',50i1)') tar0
         Call Search()

         ! Skip if no resonant state
         if(nRs==0) goto 1000

         nRs0(2)=-1
         Do i=2,nGen
            nRs0(1)=nRs0(2)+1
            nRs0(2)=nRs-1

            Do irs=nRs0(1),nRs0(2)
               tar0=tar
               Do j=1,cup(irs)
                  tar0(mode(j,irs))=tar0(mode(j,irs))+vv(j,irs)
               End do
               !dbg write(6,'(''***'',f8.3,3x,50i1)') x0(irs),tar0

               ! De0 = E0-Ep1
               De0=0.D+00
               Do j=1,cup(irs)
                  Call Modal_getEne_i(mode(j,irs),tar0(mode(j,irs)),Er1)
                  Call Modal_getEne_i(mode(j,irs),tar(mode(j,irs)),Er0)
                  De0=De0+ Er0-Er1
               End do
               De0=De0*H2wvn

               Call Search()

            End do
         End do

         ! Setup P-space
         j=0
         Do i=0,nRs-1
            if(j<cup(i)) j=cup(i)
         End do

         Call Vcor_setup_pCnf(j,nRs)
         Do i=0,nRs-1
            Call Vcor_add_pCnf(cup(i),mode(1:cup(i),i),vv(1:cup(i),i))
         End do
         !dbg Do i=1,1+npCnf
         !dbg    Call Vcor_getLabel(i,tar1)
         !dbg    write(6,'(4x,i4,3x,50i1)') i,tar1
         !dbg End do
         !dbg write(6,*)

         ! P-space CI 
         npCI=1+npCnf
         Allocate(Hmat(npCI*(npCI+1)/2),Ciw(npCI,npCI),Cie(npCI),delCnf(0:nRs-1))

         k=1
         Do i=1,npCI
            Call Vcor_getLabel(i,tar0)
            Do j=1,i
               Call Vcor_getLabel(j,tar1)
               Hmat(k)=Vcor_getHmat1(tar0,tar1)
               k=k+1
            End do
         End do
         !dbg Do i=1,npCI
         !dbg    write(6,'(50f9.5)') (Hmat(i*(i-1)/2+j),j=1,i)
         !dbg End do
         !dbg write(6,*)

         Call diag(npCI,npCI,Hmat,Ciw,Cie)

         !dbg Do i=1,npCI
         !dbg    write(6,'(i4,f8.4)') i,Cie(i)
         !dbg    Do j=1,npCI
         !dbg       Call Vcor_getLabel(j,tar1)
         !dbg       write(6,'(4x,i4,2f8.4,3x,50i1)') j,Ciw(j,i),Ciw(j,i)*Ciw(j,i),tar1
         !dbg    End do
         !dbg End do
         !dbg write(6,*)

         Do i=1,npCI
         Do j=1,npCI
            Ciw(j,i)=Ciw(j,i)*Ciw(j,i)
         End do
         End do
         !dbg Do i=1,npCI
         !dbg    write(6,'(50f7.3)') (Ciw(i,j),j=1,npCI)
         !dbg End do
         !dbg write(6,*)

         ! Remove non-degenerate states
         Allocate(iCnf(npCI),jCnf(npCI))
         iCnf=0; jCnf=-1
         jCnf(1)=0
     100 Continue
         Do i=1,npCI
            if(jCnf(i)==0) then
               Do j=1,npCI
                  if(Ciw(i,j)>pth2) then 
                     Do k=1,npCI
                        if(jCnf(k)==-1 .and. Ciw(k,j)>pth2) jCnf(k)=0
                     End do
                  endif
               End do
               jCnf(i)=1
               goto 100
            endif
         End do
         !dbg write(6,'(i4)') jCnf
         !dbg write(6,*)

         delCnf=0
         k=0
         Do i=2,npCI
            if(jCnf(i)==-1) then
               k=k+1
               delCnf(i-2)=-1
            endif
         End do
         !Do i=1,npCI
         !   if(Ciw(1,i)>1.d+00-pth2) cycle

         !   aa=0.D+00
         !   Do j=2,npCI
         !      if(Ciw(j,i)>aa) then 
         !         aa=Ciw(j,i)
         !         jmax=j
         !      endif
         !   End do
         !   if(aa>pth2) then
         !      k=-1
         !      delCnf(jmax-2)=-1
         !      !dbg Call Vcor_getLabel(jmax,tar1)
         !      !dbg write(6,'(3x,f8.3,4x,10(5i1,1x))') Ciw(jmax,i),tar1
         !   endif

         !End do
         if(k==0) then 
            Deallocate(Hmat,Ciw,Cie,delCnf)
            goto 1000
         endif

         ! Reset P-space
         Call spr_memdealloc(dble(size(pCnf))*4.D+00)
         Deallocate(pCnf)
         npCnf=0
         npCUP=0

         k=0
         Do i=0,nRs-1
            if(delCnf(i)==0) k=k+1
         End do
         !dbg write(6,'(i4)') delCnf
         !dbg write(6,*) k

         j=0
         Do i=0,nRs-1
            if(j<cup(i) .and. delCnf(i)==0) j=cup(i)
         End do

         Call Vcor_setup_pCnf(j,k)
         Do i=0,nRs-1
            if(delCnf(i)==0) &
            Call Vcor_add_pCnf(cup(i),mode(1:cup(i),i),vv(1:cup(i),i))
         End do

         Deallocate(Hmat,Ciw,Cie,delCnf)

    1000 Continue
         Deallocate(Ene1,cup,mode,vv)

         Contains

         Subroutine Search()

         Implicit None

            Integer :: i,m1,m2,m3

            Do i=1,nS2
               m1=mS2(1,i)
               m2=mS2(2,i)
               Call Search2(m1,m2)
            End do

            Do i=1,nQ2
               m1=mQ2(1,i)
               m2=mQ2(2,i)
               Call Search2(m1,m2)
            End do

            if(MR==2) return

            Do i=1,nS3
               m1=mS3(1,i)
               m2=mS3(2,i)
               m3=mS3(3,i)
               Call Search3(m1,m2,m3)
            End do

            Do i=1,nQ3
               m1=mQ3(1,i)
               m2=mQ3(2,i)
               m3=mQ3(3,i)
               Call Search3(m1,m2,m3)
            End do

         End subroutine

         Subroutine Search2(m1,m2)

         Implicit None

            Integer :: m1,m2,mm(2),n1,n2

            if(tar0(m1)==0 .and. tar0(m2)==0) return
            tar1=tar0

            mm(1)=m1; mm(2)=m2
            Do j=1,2
               Call Modal_getEne(mm(j),Ene1(0:nCHO(mm(j))-1,j))
               Er1=Ene1(tar0(mm(j)),j)
               Do k=0,nCHO(mm(j))-1
                  Ene1(k,j)=(Er1-Ene1(k,j))*H2wvn
               End do
               !dbg Write(6,'(2i2)') mm(j),tar0(mm(j))
               !dbg Write(6,'(f12.3)') (Ene1(k,j),k=0,nCHO(mm(j))-1)
            End do

            n1=nCHO(m1)
            n2=nCHO(m2)

            Call Check2(m1,m2,-1,2)
            Call Check2(m1,m2,2,-1)
            Call Check2(m1,m2,-2,1)
            Call Check2(m1,m2,1,-2)

            if(nLvl==3) return

            Call Check2(m1,m2,-1,1)
            Call Check2(m1,m2,-1,3)
            Call Check2(m1,m2,1,-1)
            Call Check2(m1,m2,3,-1)
            Call Check2(m1,m2,-2,2)
            Call Check2(m1,m2,2,-2)
            Call Check2(m1,m2,-3,1)
            Call Check2(m1,m2,1,-3)

            if(nLvl==4) return

         End subroutine

         Subroutine Check2(mi,mj,dqi,dqj)

         Implicit None

            Integer :: mi,mj,dqi,dqj
            Integer :: i,j,k,l
            Real(8) :: De1,xx,yy,Vcor_getHmat0

            tar1(mi)=tar0(mi)+dqi
            tar1(mj)=tar0(mj)+dqj
            if(dqi<0) then
               if(tar1(mi)<0 .or. tar1(mj)>nCHO(mj)) return
            endif
            if(dqj<0) then
               if(tar1(mi)>nCHO(mi) .or. tar1(mj)<0) return
            endif

            De1=Ene1(tar1(mi),1)+Ene1(tar1(mj),2)+dd

            if(abs(De1+De0)<pth0) then 
               yy=Vcor_getHmat0(2,tar0,tar1)*H2wvn
               xx=yy/De1
               if(abs(xx)>pth1) then 
                  !dbg write(6,'(3x,3f8.1,f8.3,3x,50i1)') De1+De0,De1,yy,xx,tar1
                  x0(nRs)=abs(xx)
                  Call addRs2(mi,mj,dqi,dqj)
               endif
            endif

         End subroutine

         Subroutine addRs2(mi,mj,dqi,dqj)

         Implicit None

            Integer :: mi,mj,dqi,dqj
            Integer :: i,j,k,l,fl

            j=0
            Do i=1,Nfree
               if(tar1(i)/=tar(i)) then
                  j=j+1
                  mode(j,nRs)=i
                  vv(j,nRs)=tar1(i)-tar(i)
               endif
            End do
            cup(nRs)=j
            if(j>maxCUP) then
               if(maxpCUP>0) return
               Write(6,*) 'ERROR:  ERROR WHILE CONSTRUCTING P-SPACE'
               Write(6,*) 'ERROR:  EXCEEDED MAXIMUM ORDER OF COUPLING'
               Write(6,*) 'ERROR:  RESET maxCUP0'
               Write(6,*)
               Stop
            endif

            if(j==0) return

            if(nRs==0) then
               nRs=nRs+1
               return
            endif

            !write(6,'(i4)') nRs
            !write(6,'(3x,''*'',3x,10i3)') cup(nRs),mode(1:cup(nRs),nRs),vv(1:cup(nRs),nRs)
            Do i=0,nRs-1
               !write(6,'(3x,''-'',10i3)') i,cup(i),mode(1:cup(i),i),vv(1:cup(i),i)
               fl=1
               if(cup(i)==cup(nRs)) then
                  Do j=1,cup(i)
                     if(mode(j,i)/=mode(j,nRs)) then
                        fl=0
                        exit
                     endif
                  End do
                  if(fl==1) then
                     Do j=1,cup(i)
                        if(vv(j,i)/=vv(j,nRs)) then
                           fl=0
                           exit
                        endif
                     End do
                  endif
                  if(fl==1) exit
               else
                  fl=0
               endif
            End do
            !write(6,'(3x,''fl='',i3)') fl

            if(fl==0) nRs=nRs+1
            if(nRs>maxrs) then
               Write(6,*) 'ERROR:  EXCEEDED MAXIMUM NUMBER OF RESONANT STATES'
               Write(6,*) 'ERROR:  RESET maxrs'
               Stop
            endif
            !dbg write(6,'(i3,2x,3i2)') nRs,Rs0

         End subroutine

         Subroutine Search3(m1,m2,m3)

         Implicit None

         Integer :: m1,m2,m3,mm(3),n1,n2,n3

            if(tar0(m1)==0 .and. tar0(m2)==0 .and. tar0(m3)==0) return
            tar1=tar0

            mm(1)=m1; mm(2)=m2; mm(3)=m3
            Do j=1,3
               Call Modal_getEne(mm(j),Ene1(0:nCHO(mm(j))-1,j))
               Er1=Ene1(tar0(mm(j)),j)
               Do k=0,nCHO(mm(j))-1
                  Ene1(k,j)=(Er1-Ene1(k,j))*H2wvn
               End do
               !dbg Write(6,'(2i2)') mm(j),tar0(mm(j))
               !dbg Write(6,'(f12.3)') (Ene1(k,j),k=0,nCHO(mm(j))-1)
            End do

            n1=nCHO(m1)
            n2=nCHO(m2)
            n3=nCHO(m3)

            Call Check3(m1,m2,m3,-1,1,1)
            Call Check3(m1,m2,m3,1,-1,1)
            Call Check3(m1,m2,m3,1,1,-1)

            Call Check3(m1,m2,m3,-1,-1,1)
            Call Check3(m1,m2,m3,-1,1,-1)
            Call Check3(m1,m2,m3,1,-1,-1)

            if(nLvl==3) return

            Call Check3(m1,m2,m3,-1,2,1)
            Call Check3(m1,m2,m3,-1,1,2)
            Call Check3(m1,m2,m3,2,-1,1)
            Call Check3(m1,m2,m3,1,-1,2)
            Call Check3(m1,m2,m3,2,1,-1)
            Call Check3(m1,m2,m3,1,2,-1)

            Call Check3(m1,m2,m3,-1,-1,2)
            Call Check3(m1,m2,m3,-1,2,-1)
            Call Check3(m1,m2,m3,2,-1,-1)

            Call Check3(m1,m2,m3,-2,1,1)
            Call Check3(m1,m2,m3,1,-2,1)
            Call Check3(m1,m2,m3,1,1,-2)

            Call Check3(m1,m2,m3,-2,-1,1)
            Call Check3(m1,m2,m3,-1,-2,1)
            Call Check3(m1,m2,m3,-2,1,-1)
            Call Check3(m1,m2,m3,-1,1,-2)
            Call Check3(m1,m2,m3,1,-2,-1)
            Call Check3(m1,m2,m3,1,-1,-2)

            if(nLvl==4) return

         End subroutine

         Subroutine Check3(mi,mj,mk,dqi,dqj,dqk)

         Implicit None

            Integer :: mi,mj,mk,dqi,dqj,dqk
            Integer :: i,j,k,l
            Real(8) :: De1,xx,yy,Vcor_getHmat0

            !write(6,'(''check3- '',6i3)') tar1
            tar1(mi)=tar0(mi)+dqi
            if(tar1(mi)<0 .or. tar1(mi)>nCHO(mi)) return
            tar1(mj)=tar0(mj)+dqj
            if(tar1(mj)<0 .or. tar1(mj)>nCHO(mj)) return
            tar1(mk)=tar0(mk)+dqk
            if(tar1(mk)<0 .or. tar1(mk)>nCHO(mk)) return

            !write(6,'(''check3: '',6i3)') tar1
            De1=Ene1(tar1(mi),1)+Ene1(tar1(mj),2)+Ene1(tar1(mk),3)+dd

            if(abs(De1+De0)<pth0) then 
               yy=Vcor_getHmat0(3,tar0,tar1)*H2wvn
               xx=yy/De1
               if(abs(xx)>pth1) then 
                  !dbg write(6,'(3x,3f8.1,f8.3,3x,50i1)') De1+De0,De1,yy,xx,tar1
                  x0(nRs)=abs(xx)
                  Call addRs3(mi,mj,mk,dqi,dqj,dqk)
               endif
            endif

         End subroutine

         Subroutine addRs3(mi,mj,mk,dqi,dqj,dqk)

         Implicit None

            Integer :: mi,mj,mk,dqi,dqj,dqk
            Integer :: i,j,k,l,fl

            Integer :: Rs0(Nfree)

            Rs0=tar0
            Rs0(mi)=Rs0(mi)+dqi
            Rs0(mj)=Rs0(mj)+dqj
            Rs0(mk)=Rs0(mk)+dqk
            j=0
            Do i=1,Nfree
               if(Rs0(i)/=tar(i)) then
                  j=j+1
                  mode(j,nRs)=i
                  vv(j,nRs)=Rs0(i)-tar(i)
               endif
            End do
            cup(nRs)=j
            if(j>maxCUP) then
               if(maxpCUP>0) return
               Write(6,*) 'ERROR:  ERROR WHILE CONSTRUCTING P-SPACE'
               Write(6,*) 'ERROR:  EXCEEDED MAXIMUM ORDER OF COUPLING'
               Write(6,*) 'ERROR:  RESET maxCUP0'
               Stop
            endif

            if(j==0) return

            if(nRs==0) then
               nRs=nRs+1
               return
            endif

            !dbg write(6,'(10i3)') nRs,cup(nRs),mode(1:cup(nRs),nRs),vv(1:cup(nRs),nRs)
            Do i=0,nRs-1
               fl=1
               !dbg write(6,'(10i3)') i,cup(i),mode(1:cup(i),i),vv(1:cup(i),i)
               if(cup(i)==cup(nRs)) then
                  Do j=1,cup(i)
                     if(mode(j,i)/=mode(j,nRs)) then
                        fl=0
                        exit
                     endif
                  End do
                  if(fl==1) then
                     Do j=1,cup(i)
                        if(vv(j,i)/=vv(j,nRs)) then
                           fl=0
                           exit
                        endif
                     End do
                  endif
                  if(fl==1) exit
               else
                  fl=0
               endif
            End do

            if(fl==0) nRs=nRs+1
            if(nRs>maxrs) then
               Write(6,*) 'ERROR:  EXCEEDED MAXIMUM NUMBER OF RESONANT STATES'
               Write(6,*) 'ERROR:  RESET maxrs'
               Stop
            endif
            !dbg write(6,'(i3,2x,3i2)') nRs,Rs0

         End subroutine

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setup_pCnf(n,m)

      USE Vcor_mod

      Implicit None

         Integer :: m,n,i,spr_memalloc

         i=spr_memalloc(-1,dble(2*n+1)*dble(m)*4.D+00)
         Allocate(pCnf(2*n+1,m))
         npCnf=0
         npCUP=n

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_add_pCnf(CUP,mode,vv)

      USE Vcor_mod

      Implicit None

         Integer :: i,j,CUP,mode(CUP),vv(CUP)

         npCnf=npCnf+1
         pCnf(1,npCnf)=CUP

         j=2
         Do i=1,CUP
            pCnf(j,npCnf)=mode(i)
            pCnf(j+1,npCnf)=vv(i)
            j=j+2
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_get_pCnf(nn,CUP,mode,vv)

      USE Vcor_mod

      Implicit None

         Integer :: i,j,nn,CUP,mode(CUP),vv(CUP)

         CUP=pCnf(1,nn)

         j=2
         Do i=1,CUP
            mode(i)=pCnf(j,nn)
            vv(i)=pCnf(j+1,nn)
            j=j+2
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_set_qCnf()

      USE Vcor_mod

      Implicit None

         Integer :: spr_memalloc
         Integer :: maxq,maxCUP,nq,i,j,k,m1,m2,m3
         Integer, allocatable :: tar0(:),nq0(:),mode(:,:,:),vv(:,:,:),idx(:,:,:)

         Integer :: tcup
         Integer, allocatable :: tm(:),tv(:)

         i=spr_memalloc(-1,dble(Nfree)*4.D+00)
         Allocate(tar0(Nfree))

         if(maxqCUP>0) then
            maxCUP=maxqCUP
         else
            maxCUP=MR+npCUP
            if(maxCUP>Nfree) maxCUP=Nfree
         endif

         i=nS1+nQ1
         j=nS2+nQ2
         k=nS3+nQ3
         if(nLvl==3) then
            maxq=i*4 + j*8 + k*8
         elseif(nLvl==4) then
            maxq=i*8 + j*24 + k*32
         else
            maxq=i*10 + j*40 + k*80
         endif
         maxq=2.D+00*maxq*(npCnf+1)/maxCUP

         i=spr_memalloc(-1,dble(maxCUP)*4.D+00  &
                          +dble(maxq)*dble(maxCUP)*dble(maxCUP)*8.D+00 &
                          +dble(maxq)*dble(maxCUP)*8.D+00 )
         Allocate(nq0(maxCUP))
         Allocate(mode(maxCUP,0:maxq-1,maxCUP))
         Allocate(vv(maxCUP,0:maxq-1,maxCUP))
         Allocate(idx(0:1,0:maxq-1,maxCUP))
         idx=-1; nq0=0

         Do i=0,npCnf
            if(i/=0) then
               Call Vcor_getLabel(i+1,tar0)
            else
               tar0=tar
            endif

            if(.not. ss .or. npCnf>0 .or. lvci) then
               Do j=1,nS1
                  m1=mS1(j)
                  Call add1(m1)
               End do
               Do j=1,nQ1
                  m1=mQ1(j)
                  Call add1(m1)
               End do
            endif

            Do j=1,nS2
               m1=mS2(1,j); m2=mS2(2,j)
               Call add2(m1,m2)
            End do
            Do j=1,nQ2
               m1=mQ2(1,j); m2=mQ2(2,j)
               Call add2(m1,m2)
            End do

            Do j=1,nS3
               m1=mS3(1,j); m2=mS3(2,j); m3=mS3(3,j)
               Call add3(m1,m2,m3)
            End do
            Do j=1,nQ3
               m1=mQ3(1,j); m2=mQ3(2,j); m3=mQ3(3,j)
               Call add3(m1,m2,m3)
            End do

         End do

         ! write(6,'(20i5)') nq0
         ! Do k=1,maxCUP
         ! Do i=0,nq0(k)-1
         !    tar0=tar
         !    Do j=1,k
         !       tar0(mode(j,i,k))=tar0(mode(j,i,k))+vv(j,i,k)
         !    End do
         !    write(6,'(2i5,2x,30i2)') k,i+1,tar0
         ! End do
         ! End do

         ! Setup Q-space
         nq=0
         Do i=1,maxCUP
            nq=nq + nq0(i)
         End do
         Do i=maxCUP,1,-1
            if(nq0(i)>0) then
               j=i
               exit
            endif
         End do

         Call Vcor_setup_qCnf(j,nq)
         Do k=1,j
         Do i=0,nq0(k)-1
            Call Vcor_add_qCnf(k,mode(1:k,i,k),vv(1:k,i,k))
         End do
         End do
         !write(6,'(20i5)') nqCnfi

         Call spr_memdealloc(dble(Nfree)*4.D+00)
         Deallocate(tar0)
         Call spr_memdealloc(dble(size(mode))*4.D+00)
         Deallocate(mode)
         Call spr_memdealloc(dble(size(vv))*4.D+00)
         Deallocate(vv)
         Call spr_memdealloc(dble(size(nq0))*4.D+00)
         Deallocate(nq0)
         Call spr_memdealloc(dble(size(idx))*4.D+00)
         Deallocate(idx)

         Contains

         Subroutine add1(m1)

         Implicit None

            Integer :: ii,jj,m1,v1,n1

            v1=tar0(m1)
            n1=nCHO(m1)

            if(v1+3<n1) then
               tar0(m1)=v1+1
               Call add()
               tar0(m1)=v1+3
               Call add()
            elseif(v1+1<n1) then
               tar0(m1)=v1+1
               Call add()
            endif

            tar0(m1)=v1-1
            if(tar0(m1)>=0) Call add()
            tar0(m1)=v1-3
            if(tar0(m1)>=0) Call add()

            if(nLvl==3) then
               tar0(m1)=v1
               return
            endif

            if(v1+4<n1) then
               tar0(m1)=v1+2
               Call add()
               tar0(m1)=v1+4
               Call add()
            elseif(v1+2<n1) then
               tar0(m1)=v1+2
               Call add()
            endif

            tar0(m1)=v1-2
            if(tar0(m1)>=0) Call add()
            tar0(m1)=v1-4
            if(tar0(m1)>=0) Call add()

            if(nLvl==4) then
               tar0(m1)=v1
               return
            endif

            Do ii=5,nLvl
               tar0(m1)=v1+ii
               if(tar0(m1)<n1) Call add()

               tar0(m1)=v1-ii
               if(tar0(m1)>=0) Call add()
            End do
            tar0(m1)=v1

         End subroutine


         Subroutine add2(m1,m2)

         Implicit None

            Integer :: i1,i2,j1,j2,m1,m2,v1,v2,n1,n2

            v1=tar0(m1)
            v2=tar0(m2)
            n1=nCHO(m1)
            n2=nCHO(m2)

            if(nLvl==3) then
               Call add2_0(v1,v2,n1,n2,2,1)
               Call add2_0(v1,v2,n1,n2,1,2)

            else
               Do i1=2,nLvl
               Do i2=1,i1-1
                  j2=i1-i2
                  Call add2_0(v1,v2,n1,n2,i2,j2)
               End do
               End do

            endif
            tar0(m1)=v1
            tar0(m2)=v2

         End subroutine

         Subroutine add2_0(v1,v2,n1,n2,d1,d2)

         Implicit None

         Integer :: v1,v2,d1,d2,n1,n2

            tar0(m1)=v1+d1
            tar0(m2)=v2+d2
            if(tar0(m1)<n1 .and. tar0(m2)<n2) Call add()

            tar0(m1)=v1-d1
            tar0(m2)=v2+d2
            if(tar0(m1)>=0 .and. tar0(m2)<n2) Call add()

            tar0(m1)=v1+d1
            tar0(m2)=v2-d2
            if(tar0(m1)<n1 .and. tar0(m2)>=0) Call add()

            tar0(m1)=v1-d1
            tar0(m2)=v2-d2
            if(tar0(m1)>=0 .and. tar0(m2)>=0) Call add()

         End subroutine

         Subroutine add3(m1,m2,m3)

         Implicit None

            Integer :: m1,m2,m3,v1,v2,v3,n1,n2,n3
            Integer :: ii,j1,j2,j3,k1,k2,k3

            v1=tar0(m1)
            v2=tar0(m2)
            v3=tar0(m3)
            n1=nCHO(m1)
            n2=nCHO(m2)
            n3=nCHO(m3)

            Do ii=3,nLvl
               Do j1=1,ii-2
               Do j2=1,ii-2
                  j3=ii-j1-j2
                  if(j3<1) cycle
                  !dbg write(6,'(3i3)')j1,j2,j3
                  Call add3_0(v1,v2,v3,n1,n2,n3,j1,j2,j3)
               End do
               End do
            End do

            tar0(m1)=v1
            tar0(m2)=v2
            tar0(m3)=v3

         End subroutine

         Subroutine add3_0(v1,v2,v3,n1,n2,n3,d1,d2,d3)

         Implicit None

         Integer :: v1,v2,v3,d1,d2,d3,n1,n2,n3

            tar0(m1)=v1+d1
            tar0(m2)=v2+d2
            tar0(m3)=v3+d3
            if(tar0(m1)<n1 .and. tar0(m2)<n2 .and. tar0(m3)<n3) Call add()

            tar0(m1)=v1-d1
            tar0(m2)=v2+d2
            tar0(m3)=v3+d3
            if(tar0(m1)>=0 .and. tar0(m2)<n2 .and. tar0(m3)<n3) Call add()

            tar0(m1)=v1+d1
            tar0(m2)=v2-d2
            tar0(m3)=v3+d3
            if(tar0(m1)<n1 .and. tar0(m2)>=0 .and. tar0(m3)<n3) Call add()

            tar0(m1)=v1+d1
            tar0(m2)=v2+d2
            tar0(m3)=v3-d3
            if(tar0(m1)<n1 .and. tar0(m2)<n2 .and. tar0(m3)>=0) Call add()

            tar0(m1)=v1+d1
            tar0(m2)=v2-d2
            tar0(m3)=v3-d3
            if(tar0(m1)<n1 .and. tar0(m2)>=0 .and. tar0(m3)>=0) Call add()

            tar0(m1)=v1-d1
            tar0(m2)=v2+d2
            tar0(m3)=v3-d3
            if(tar0(m1)>=0 .and. tar0(m2)<n2 .and. tar0(m3)>=0) Call add()

            tar0(m1)=v1-d1
            tar0(m2)=v2-d2
            tar0(m3)=v3+d3
            if(tar0(m1)>=0 .and. tar0(m2)>=0 .and. tar0(m3)<n3) Call add()

            tar0(m1)=v1-d1
            tar0(m2)=v2-d2
            tar0(m3)=v3-d3
            if(tar0(m1)>=0 .and. tar0(m2)>=0 .and. tar0(m3)>=0) Call add()

         End subroutine

         Subroutine add()

         Implicit None

            Integer :: i,j,k,l,fl
            Integer :: compare,lr,nn
            Integer :: cup,mm1(maxCUP),vv1(maxCUP), &
                           mm2(maxCUP),vv2(maxCUP)

            cup=0
            Do i=1,Nfree
               if(tar0(i)/=tar(i)) then
                  cup=cup+1
                  mm1(cup)=i
                  vv1(cup)=tar0(i)-tar(i)
               endif
            End do
            if(cup>maxCUP) return
            if(cup==0) return

            if(npCnf>0) then 
               Do i=1,npCnf
                  fl=1
                  if(pCnf(1,i)==cup) then
                     Do j=1,cup
                        if(pCnf(2*j,i)/=mm1(j)) then
                           fl=0
                           exit
                        endif
                     End do
                     if(fl==1) then
                        Do j=1,cup
                           if(pCnf(2*j+1,i)/=vv1(j)) then
                              fl=0
                              exit
                           endif
                        End do 
                     endif
                     if(fl==1) exit
                  else
                     fl=0
                  endif
               End do
               if(fl==1) return
            endif

            if(nq0(cup)==0) then
               Do i=1,cup
                  mode(i,nq0(cup),cup)=mm1(i)
                  vv(i,nq0(cup),cup)=vv1(i)
               End do
               nq0(cup)=nq0(cup)+1 
               return
            endif

            nn=0
            Do j=1,cup
               mm2(j)=mode(j,nn,cup)
               vv2(j)=vv(j,nn,cup)
            End do
            Do while(.true.)
               lr=compare(cup,mm1,vv1,mm2,vv2)
               !write(6,'(5i3)') lr,cup
               !write(6,'(6i3)') mm1(1:cup),vv1(1:cup)
               !write(6,'(6i3)') mm2(1:cup),vv2(1:cup)
               !write(6,*)

               ! mm1=mm2
               if(lr<0) exit

               if(idx(lr,nn,cup)>=0) then
                  nn=idx(lr,nn,cup)
                  Do j=1,cup
                     mm2(j)=mode(j,nn,cup)
                     vv2(j)=vv(j,nn,cup)
                  End do

               else
                  idx(lr,nn,cup)=nq0(cup)
                  Do i=1,cup
                     mode(i,nq0(cup),cup)=mm1(i)
                     vv(i,nq0(cup),cup)=vv1(i)
                  End do
                  nq0(cup)=nq0(cup)+1
                  if(nq0(cup)>maxq) then
                     Write(6,*) 'ERROR:  EXCEEDED MAXIMUM NUMBER OF Q-SPACE'
                     Write(6,*) 'ERROR:  RESET maxq'
                     Stop
                  endif
                  exit

               endif

            End do

         End subroutine

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

         ! =-1: 1=2
         ! =0 : 1>2
         ! =1 : 1<2
         Function compare(cp,mm1,vv1,mm2,vv2)

         Implicit None

            Integer :: compare
            Integer :: cp,mm1(cp),vv1(cp),mm2(cp),vv2(cp)

            Integer :: i

            Do i=1,cp
               if(mm1(i)>mm2(i)) then
                  compare=0
                  return
               elseif(mm1(i)<mm2(i)) then
                  compare=1
                  return
               endif
            End do

            Do i=1,cp
               if(vv1(i)>vv2(i)) then
                  compare=0
                  return
               elseif(vv1(i)<vv2(i)) then
                  compare=1
                  return
               endif
            End do

            compare=-1

         End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setup_qCnf(n,m)

      USE Vcor_mod

      Implicit None

         Integer :: m,n,i,spr_memalloc

         i=spr_memalloc(-1,dble(2*n+1)*dble(m)*4.D+00+dble(n)*4.D+00)
         Allocate(qCnf(2*n+1,m),nqCnfi(n))
         nqCnf=0
         nqCUP=n
         nqCnfi=0

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_add_qCnf(CUP,mode,vv)

      USE Vcor_mod

      Implicit None

         Integer :: i,j,CUP,mode(CUP),vv(CUP)

         nqCnf=nqCnf+1
         nqCnfi(CUP)=nqCnfi(CUP)+1
         qCnf(1,nqCnf)=CUP

         j=2
         Do i=1,CUP
            qCnf(j,nqCnf)=mode(i)
            qCnf(j+1,nqCnf)=vv(i)
            j=j+2
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_getLabel(nn,label)

      USE Vcor_mod

      Implicit None

         Integer :: nn,label(Nfree)
         Integer :: i,j,nC,mm

         Do i=1,Nfree
            label(i)=tar(i)
         End do

         if(nn==1) return

         if(nn<npCnf+2) then
            mm=nn-1
            nC=pCnf(1,mm)
            j=2
            Do i=1,nC
               Label(pCnf(j,mm))=Label(pCnf(j,mm))+pCnf(j+1,mm)
               j=j+2
            End do
         else
            mm=nn-1-npCnf
            nC=qCnf(1,mm)
            j=2
            Do i=1,nC
               Label(qCnf(j,mm))=Label(qCnf(j,mm))+qCnf(j+1,mm)
               j=j+2
            End do
         endif

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_ClearCnf()

      USE Vcor_mod

      Implicit None

         if(allocated(pCnf)) then
            Call spr_memdealloc(dble(size(pCnf))*4.D+00)
            Deallocate(pCnf)
            npCnf=0
         endif
         if(allocated(qCnf)) then
            Call spr_memdealloc(dble(size(qCnf))*4.D+00)
            Deallocate(qCnf)
            nqCnf=0
         endif
         if(allocated(Qmat)) then
            Call spr_memdealloc(dble(size(Qmat))*8.D+00)
            Deallocate(Qmat)
         endif
         if(allocated(nqCnfi)) then
            Call spr_memdealloc(dble(size(nqCnfi))*4.D+00)
            Deallocate(nqCnfi)
         endif

         npCUP=0
         nqCUP=0
         npCnf=0
         nqCnf=0

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setEwfn(mode,nGrid,Ewfn)

      USE Vcor_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: Ewfn(nGrid*nCHO(mode))

         nGrid2=nGrid*nCHO(mode)
         pt=type_Table1(idx1(mode)+nGrid)
         block(pt+1:pt+nGrid2)=Ewfn

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_getEwfn(mode,nGrid,qm,Ewfn)

      USE Vcor_mod

      Implicit None

         Integer :: mode,nGrid,qm,pt
         Real(8) :: Ewfn(nGrid)

         pt=type_Table1(idx1(mode)+nGrid)+nGrid*qm
         Ewfn=block(pt+1:pt+nGrid)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_setDwfn(mode,nGrid,Dwfn)

      USE Vcor_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: Dwfn(nGrid*nCHO(mode))

         nGrid2=nGrid*nCHO(mode)
         pt=type_Table1(idx1(mode)+nGrid)+nGrid2
         block(pt+1:pt+nGrid2)=Dwfn

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_getDwfn(mode,nGrid,qm,Dwfn)

      USE Vcor_mod

      Implicit None

         Integer :: mode,nGrid,qm,nGrid2,pt
         Real(8) :: Dwfn(nGrid)

         nGrid2=nGrid*nCHO(mode)
         pt=type_Table1(idx1(mode)+nGrid)+nGrid2+nGrid*qm
         Dwfn=block(pt+1:pt+nGrid)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vcor_getXwfn(mode,nGrid,qm,qn,Xwfn)

      USE Vcor_mod

      Implicit None

         Integer :: mode,nGrid,qm,qn,i
         Real(8) :: Xwfn(nGrid),Em(nGrid),En(nGrid)

         Call Vcor_getEwfn(mode,nGrid,qm,Em)
         Call Vcor_getEwfn(mode,nGrid,qn,En)
         Do i=1,nGrid
            Xwfn(i)=Em(i)*En(i)
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vcor_getCUP(iin,jin)

      USE Vcor_mod

      Implicit None

         Integer :: iin,jin,Vcor_getCUP
         Integer :: ii,jj,i,j,k,l,cp
         Integer :: icnf,jcnf,icp,jcp
         Integer, allocatable ::   mi(:),mj(:),vi(:),vj(:)
                    !mi(nqCUP),mj(nqCUP),vi(nqCUP),vj(nqCUP)

         if(iin==jin) then
            Vcor_getCUP=0
            return
         endif

         ! ii > jj
         if(iin>jin) then
            ii=iin; jj=jin
         else
            ii=jin; jj=iin
         endif

         if(jj==1) then
           ! P0-Q
            if(ii> npCnf+1) then
               Vcor_getCUP=qCnf(1,ii-npCnf-1)
            else
           ! P0-P
               Vcor_getCUP=pCnf(1,ii-1)
            endif
            return
         endif

         ! Q-Q
         if(jj>npCnf+1) then
            icnf=ii-npCnf-1
            icp=qCnf(1,icnf)
            Allocate(mi(icp),vi(icp))
            Do i=1,icp
               mi(i)=qCnf(i*2,icnf)
               vi(i)=qCnf(i*2+1,icnf)
            End do

            jcnf=jj-npCnf-1
            jcp=qCnf(1,jcnf)
            Allocate(mj(jcp),vj(jcp))
            Do i=1,jcp
               mj(i)=qCnf(i*2,jcnf)
               vj(i)=qCnf(i*2+1,jcnf)
            End do

            cp=icp+jcp
            Do k=1,icp
            Do l=1,jcp
               if(mi(k)==mj(l)) then 
                  cp=cp-1
                  if(vi(k)==vj(l)) cp=cp-1
               endif
            End do
            End do

            Vcor_getCUP=cp

            Deallocate(mi,vi,mj,vj)
            return
         endif

         ! P-P
         if(ii<npCnf+2) then
            icnf=ii-1
            icp=pCnf(1,icnf)
            Allocate(mi(icp),vi(icp))
            Do i=1,icp
               mi(i)=pCnf(i*2,icnf)
               vi(i)=pCnf(i*2+1,icnf)
            End do

            jcnf=jj-1
            jcp=pCnf(1,jcnf)
            Allocate(mj(jcp),vj(jcp))
            Do i=1,jcp
               mj(i)=pCnf(i*2,jcnf)
               vj(i)=pCnf(i*2+1,jcnf)
            End do

            cp=icp+jcp
            Do k=1,icp
            Do l=1,jcp
               if(mi(k)==mj(l)) then 
                  cp=cp-1
                  if(vi(k)==vj(l)) cp=cp-1
               endif
            End do
            End do

            Vcor_getCUP=cp

            Deallocate(mi,vi,mj,vj)
            return
         endif

         ! P(j)-Q(i)
         icnf=ii-npCnf-1
         icp=qCnf(1,icnf)
         Allocate(mi(icp),vi(icp))
         Do i=1,icp
            mi(i)=qCnf(i*2,icnf)
            vi(i)=qCnf(i*2+1,icnf)
         End do

         jcnf=jj-1
         jcp=pCnf(1,jcnf)
         Allocate(mj(jcp),vj(jcp))
         Do i=1,jcp
            mj(i)=pCnf(i*2,jcnf)
            vj(i)=pCnf(i*2+1,jcnf)
         End do

         cp=icp+jcp
         Do k=1,icp
         Do l=1,jcp
            if(mi(k)==mj(l)) then 
               cp=cp-1
               if(vi(k)==vj(l)) cp=cp-1
            endif
         End do
         End do

         Vcor_getCUP=cp

         Deallocate(mi,vi,mj,vj)
         return

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vcor_getHmat(cp,mm,nn)

      USE Vcor_mod

      Implicit None

         Integer :: cp,mm,nn
         Integer :: nCnf(Nfree),mCnf(Nfree)
         Real(8) :: Vcor_getHmat, Vcor_getHmat0

         Call Vcor_getLabel(mm,mCnf)
         if(mm/=nn) then
            Call Vcor_getLabel(nn,nCnf)
         else
            nCnf=mCnf
         endif
         Vcor_getHmat=Vcor_getHmat0(cp,mCnf,nCnf)

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vcor_getHmat1(mCnf,nCnf)

      USE Vcor_mod

      Implicit None

         Integer :: i,cp,mm,nn
         Integer :: nCnf(Nfree),mCnf(Nfree)
         Real(8) :: Vcor_getHmat1, Vcor_getHmat0

         cp=0
         Do i=1,Nfree
            if(mCnf(i)/=nCnf(i)) cp=cp+1
         Enddo
         Vcor_getHmat1=Vcor_getHmat0(cp,mCnf,nCnf)

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vcor_getHmat0(cp,mCnf,nCnf)

      USE Vcor_mod

      Implicit None

         Integer :: cp
         Real(8) :: Vcor_getHmat0

         Integer :: nCnf(Nfree),mCnf(Nfree),mode(cp)

         Integer :: i,spr_memalloc
         Integer :: mi,mj,mk,nGi,nGj,nGk,qi,qj,qk,m0,m,n
         Real(8), allocatable :: Di(:),Dj(:),Dk(:)
         Real(8), allocatable :: Vi(:),Vij(:,:),Vijk(:,:,:)
         Real(8), allocatable :: qav(:,:),qdi(:,:)
         Real(8) :: Ekin,Epot,Vcor_getKinmat
         Real(8) :: tmp,tmp2,tmp3

         Vcor_getHmat0=0.D+00
         if(cp>0) then
            i=0
            Do mi=1,Nfree
               if(nCnf(mi)/=mCnf(mi)) then
                  i=i+1
                  mode(i)=mi
                  if(i==cp) exit
               endif
            End do
         endif

         ! QFF
         if(nQ1/=0 .or. nQ2/=0 .or. nQ3/=0) then
            i=spr_memalloc(-1,dble(Nfree)*4.D+00*8.D+00)
            Allocate(qav(4,Nfree))
            Do mi=1,Nfree
               if(nCHO(mi)==0) cycle
               Do i=1,4
                  qav(i,mi)=Qmat(i,mCnf(mi),nCnf(mi),mi)
               End do
            End do

         endif

         Ekin=0.D+00; Epot=0.D+00
         Select Case(cp)

        ! [ Diagonal ]----------------------------------------------------------- 
            Case(0) 
               Do mi=1,Nfree
                  Ekin=Ekin + Vcor_getKinmat(mi,mCnf(mi),mCnf(mi))
               End do

             ! 1-mode terms
               Do n=1,nQ1
                  mi=mQ1(n)
                  Epot=Epot + qav(1,mi)*coeff1(1,n) &
                            + qav(2,mi)*coeff1(2,n) &
                            + qav(3,mi)*coeff1(3,n) &
                            + qav(4,mi)*coeff1(4,n)
               End do

               Do n=1,nS1
                  mi=mS1(n)
                  Call PEGrid_getnGrid1(nGi,n)
                  Allocate(Vi(nGi),Di(nGi))

                  Call PEGrid_getV1(n,Vi)
                  Call Vcor_getDwfn(mi,nGi,mCnf(mi),Di)
                  tmp=0.D+00
                  Do qi=1,nGi
                     tmp=tmp + Di(qi)*Vi(qi)
                  End do
                  Epot=Epot+tmp

                  Deallocate(Vi,Di)
               End do

               if(MR==1) goto 100

             ! 2-mode terms
               Do n=1,nQ2
                  mi=mQ2(1,n)
                  mj=mQ2(2,n)
                  Epot=Epot + qav(1,mi)*qav(1,mj)*coeff2(1,n) &
                            + qav(2,mi)*qav(1,mj)*coeff2(2,n) &
                            + qav(1,mi)*qav(2,mj)*coeff2(3,n) &
                            + qav(2,mi)*qav(2,mj)*coeff2(4,n) &
                            + qav(3,mi)*qav(1,mj)*coeff2(5,n) &
                            + qav(1,mi)*qav(3,mj)*coeff2(6,n) 
               End do

               Do n=1,nS2
                  mi=mS2(1,n)
                  mj=mS2(2,n)
                  Call PEGrid_getnGrid2(nGi,nGj,n)
                  Allocate(Vij(nGj,nGi),Di(nGi),Dj(nGj))
                  Call PEGrid_getV2(n,Vij)
                  Call Vcor_getDwfn(mi,nGi,mCnf(mi),Di)
                  Call Vcor_getDwfn(mj,nGj,mCnf(mj),Dj)

                  tmp=0.D+00
                  Do qi=1,nGi
                     tmp2=0.D+00
                     Do qj=1,nGj
                        tmp2=tmp2 + Dj(qj)*Vij(qj,qi)
                     End do
                     tmp=tmp + Di(qi)*tmp2
                  End do
                  Epot=Epot+tmp

                  Deallocate(Vij,Di,Dj)
               End do

               if(MR==2) goto 100

             ! 3-mode terms
               Do n=1,nQ3
                  mi=mQ3(1,n)
                  mj=mQ3(2,n)
                  mk=mQ3(3,n)
                  Epot=Epot + qav(1,mi)*qav(1,mj)*qav(1,mk)*coeff3(1,n) &
                            + qav(2,mi)*qav(1,mj)*qav(1,mk)*coeff3(2,n) &
                            + qav(1,mi)*qav(2,mj)*qav(1,mk)*coeff3(3,n) &
                            + qav(1,mi)*qav(1,mj)*qav(2,mk)*coeff3(4,n) 
               End do

               Do n=1,nS3
                  mi=mS3(1,n)
                  mj=mS3(2,n)
                  mk=mS3(3,n)
                  Call PEGrid_getnGrid3(nGi,nGj,nGk,n)
                  Allocate(Vijk(nGk,nGj,nGi))
                  Allocate(Di(nGi),Dj(nGj),Dk(nGk))
                  Call PEGrid_getV3(n,Vijk)
                  Call Vcor_getDwfn(mi,nGi,mCnf(mi),Di)
                  Call Vcor_getDwfn(mj,nGj,mCnf(mj),Dj)
                  Call Vcor_getDwfn(mk,nGk,mCnf(mk),Dk)

                  tmp=0.D+00
                  Do qi=1,nGi
                     tmp2=0.D+00
                     Do qj=1,nGj
                        tmp3=0.D+00
                        Do qk=1,nGk
                           tmp3=tmp3 + Dk(qk)*Vijk(qk,qj,qi)
                        End do
                        tmp2=tmp2 + Dj(qj)*tmp3
                     End do
                     tmp=tmp + Di(qi)*tmp2
                  End do
                  Epot=Epot+tmp

                  Deallocate(Vijk,Di,Dj,Dk)

               End do

               if(MR==3) goto 100

           100 Continue
               Vcor_getHmat0=Ekin+Epot

        ! [ 1-mode exc. ]-------------------------------------------------------- 
            Case(1) 
               mi=mode(1)
               Ekin=Vcor_getKinmat(mi,mCnf(mi),nCnf(mi))

             ! 1-mode terms
               m=0
               Do n=1,nS1
                  mi=mS1(n)
                  if(mi/=mode(1)) cycle

                  Call PEGrid_getnGrid1(nGi,n)
                  Allocate(Vi(nGi),Di(nGi))

                  Call PEGrid_getV1(n,Vi)
                  Call Vcor_getXwfn(mi,nGi,mCnf(mi),nCnf(mi),Di)
                  tmp=0.D+00
                  Do qi=1,nGi
                     tmp=tmp + Di(qi)*Vi(qi)
                  End do
                  Epot=Epot+tmp

                  Deallocate(Vi,Di)
                  m=1
                  exit
               End do

               if(m==0) then

                  Do n=1,nQ1
                     mi=mQ1(n)
                     if(mi/=mode(1)) cycle

                     Epot=Epot + qav(1,mi)*coeff1(1,n) &
                               + qav(2,mi)*coeff1(2,n) &
                               + qav(3,mi)*coeff1(3,n) &
                               + qav(4,mi)*coeff1(4,n)
                     exit
                  End do

               endif

               if(MR==1) goto 200

             ! 2-mode terms
               m=0; m0=Nfree-1
               Do n=1,nS2
                  mi=mS2(1,n)
                  mj=mS2(2,n)
                  if(mi==mode(1) .or. mj==mode(1)) then
                     Call PEGrid_getnGrid2(nGi,nGj,n)
                     Allocate(Vij(nGj,nGi),Di(nGi),Dj(nGj))
                     Call PEGrid_getV2(n,Vij)

                     if(mi==mode(1)) then
                        Call Vcor_getXwfn(mi,nGi,mCnf(mi),nCnf(mi),Di)
                        Call Vcor_getDwfn(mj,nGj,mCnf(mj),Dj)
                     else
                        Call Vcor_getDwfn(mi,nGi,mCnf(mi),Di)
                        Call Vcor_getXwfn(mj,nGj,mCnf(mj),nCnf(mj),Dj)
                     endif

                     tmp=0.D+00
                     Do qi=1,nGi
                        tmp2=0.D+00
                        Do qj=1,nGj
                           tmp2=tmp2 + Dj(qj)*Vij(qj,qi)
                        End do
                        tmp=tmp + Di(qi)*tmp2
                     End do
                     Epot=Epot+tmp

                     Deallocate(Vij,Di,Dj)
                     m=m+1
                     if(m==m0) exit
                  endif
               End do

               if(m<m0) then
                  Do n=1,nQ2
                     mi=mQ2(1,n)
                     mj=mQ2(2,n)
                     if(mi==mode(1) .or. mj==mode(1)) then
                        Epot=Epot + qav(1,mi)*qav(1,mj)*coeff2(1,n) &
                                  + qav(2,mi)*qav(1,mj)*coeff2(2,n) &
                                  + qav(1,mi)*qav(2,mj)*coeff2(3,n) &
                                  + qav(2,mi)*qav(2,mj)*coeff2(4,n) &
                                  + qav(3,mi)*qav(1,mj)*coeff2(5,n) &
                                  + qav(1,mi)*qav(3,mj)*coeff2(6,n) 
                        m=m+1
                        if(m==m0) exit
                     endif
                  End do
               endif

               if(MR==2) goto 200

             ! 3-mode terms
               m=0; m0=(Nfree-1)*(Nfree-2)/2
               Do n=1,nS3
                  mi=mS3(1,n)
                  mj=mS3(2,n)
                  mk=mS3(3,n)
                  if(mi==mode(1) .or. mj==mode(1) .or. mk==mode(1)) then
                     Call PEGrid_getnGrid3(nGi,nGj,nGk,n)
                     Allocate(Vijk(nGk,nGj,nGi))
                     Allocate(Di(nGi),Dj(nGj),Dk(nGk))
                     Call PEGrid_getV3(n,Vijk)
                     if(mi==mode(1)) then
                        Call Vcor_getXwfn(mi,nGi,mCnf(mi),nCnf(mi),Di)
                        Call Vcor_getDwfn(mj,nGj,mCnf(mj),Dj)
                        Call Vcor_getDwfn(mk,nGk,mCnf(mk),Dk)
                     elseif(mj==mode(1)) then
                        Call Vcor_getDwfn(mi,nGi,mCnf(mi),Di)
                        Call Vcor_getXwfn(mj,nGj,mCnf(mj),nCnf(mj),Dj)
                        Call Vcor_getDwfn(mk,nGk,mCnf(mk),Dk)
                     elseif(mk==mode(1)) then
                        Call Vcor_getDwfn(mi,nGi,mCnf(mi),Di)
                        Call Vcor_getDwfn(mj,nGj,mCnf(mj),Dj)
                        Call Vcor_getXwfn(mk,nGk,mCnf(mk),nCnf(mk),Dk)
                     endif

                     tmp=0.D+00
                     Do qi=1,nGi
                        tmp2=0.D+00
                        Do qj=1,nGj
                           tmp3=0.D+00
                           Do qk=1,nGk
                              tmp3=tmp3 + Dk(qk)*Vijk(qk,qj,qi)
                           End do
                           tmp2=tmp2 + Dj(qj)*tmp3
                        End do
                        tmp=tmp + Di(qi)*tmp2
                     End do
                     Epot=Epot+tmp

                     Deallocate(Vijk,Di,Dj,Dk)
                     m=m+1
                     if(m==m0) exit

                  endif
               End do

               if(m<m0) then
                  Do n=1,nQ3
                     mi=mQ3(1,n)
                     mj=mQ3(2,n)
                     mk=mQ3(3,n)
                     if(mi==mode(1) .or. mj==mode(1) .or. mk==mode(1)) then
                        Epot=Epot + qav(1,mi)*qav(1,mj)*qav(1,mk)*coeff3(1,n) &
                                  + qav(2,mi)*qav(1,mj)*qav(1,mk)*coeff3(2,n) &
                                  + qav(1,mi)*qav(2,mj)*qav(1,mk)*coeff3(3,n) &
                                  + qav(1,mi)*qav(1,mj)*qav(2,mk)*coeff3(4,n) 
                        m=m+1
                        if(m==m0) exit
                     endif
                  End do
               endif

               if(MR==3) goto 200

           200 Continue
               Vcor_getHmat0=Ekin+Epot

        ! [ 2-mode exc. ]-------------------------------------------------------- 
            Case(2) 

             ! 2-mode terms
               m=0
               Do n=1,nS2
                  mi=mS2(1,n)
                  mj=mS2(2,n)
                  if(mi==mode(2) .and. mj==mode(1)) then
                     Call PEGrid_getnGrid2(nGi,nGj,n)
                     Allocate(Vij(nGj,nGi),Di(nGi),Dj(nGj))
                     Call PEGrid_getV2(n,Vij)
                     Call Vcor_getXwfn(mi,nGi,mCnf(mi),nCnf(mi),Di)
                     Call Vcor_getXwfn(mj,nGj,mCnf(mj),nCnf(mj),Dj)

                     tmp=0.D+00
                     Do qi=1,nGi
                        tmp2=0.D+00
                        Do qj=1,nGj
                           tmp2=tmp2 + Dj(qj)*Vij(qj,qi)
                        End do
                        tmp=tmp + Di(qi)*tmp2
                     End do
                     Epot=Epot+tmp

                     Deallocate(Vij,Di,Dj)
                     m=1
                     exit
                  endif
               End do

               if(m==0) then
                  Do n=1,nQ2
                     mi=mQ2(1,n)
                     mj=mQ2(2,n)
                     if(mi==mode(2) .and. mj==mode(1)) then
                        Epot=Epot + qav(1,mi)*qav(1,mj)*coeff2(1,n) &
                                  + qav(2,mi)*qav(1,mj)*coeff2(2,n) &
                                  + qav(1,mi)*qav(2,mj)*coeff2(3,n) &
                                  + qav(2,mi)*qav(2,mj)*coeff2(4,n) &
                                  + qav(3,mi)*qav(1,mj)*coeff2(5,n) &
                                  + qav(1,mi)*qav(3,mj)*coeff2(6,n) 
                        m=1
                        exit
                     endif
                  End do
               endif

               if(MR==2) goto 300

             ! 3-mode terms
               m=0; m0=(Nfree-2)
               Do n=1,nS3
                  mi=mS3(1,n)
                  mj=mS3(2,n)
                  mk=mS3(3,n)
                  if((mi==mode(2) .and. mj==mode(1)) .or. &
                     (mi==mode(2) .and. mk==mode(1)) .or. &
                     (mj==mode(2) .and. mk==mode(1))) then

                     Call PEGrid_getnGrid3(nGi,nGj,nGk,n)
                     Allocate(Vijk(nGk,nGj,nGi))
                     Allocate(Di(nGi),Dj(nGj),Dk(nGk))
                     Call PEGrid_getV3(n,Vijk)

                     if(mi==mode(2) .and. mj==mode(1)) then
                        Call Vcor_getXwfn(mi,nGi,mCnf(mi),nCnf(mi),Di)
                        Call Vcor_getXwfn(mj,nGj,mCnf(mj),nCnf(mj),Dj)
                        Call Vcor_getDwfn(mk,nGk,mCnf(mk),Dk)
                     elseif(mi==mode(2) .and. mk==mode(1)) then
                        Call Vcor_getXwfn(mi,nGi,mCnf(mi),nCnf(mi),Di)
                        Call Vcor_getDwfn(mj,nGj,mCnf(mj),Dj)
                        Call Vcor_getXwfn(mk,nGk,mCnf(mk),nCnf(mk),Dk)
                     else
                        Call Vcor_getDwfn(mi,nGi,mCnf(mi),Di)
                        Call Vcor_getXwfn(mj,nGj,mCnf(mj),nCnf(mj),Dj)
                        Call Vcor_getXwfn(mk,nGk,mCnf(mk),nCnf(mk),Dk)
                     endif

                     tmp=0.D+00
                     Do qi=1,nGi
                        tmp2=0.D+00
                        Do qj=1,nGj
                           tmp3=0.D+00
                           Do qk=1,nGk
                              tmp3=tmp3 + Dk(qk)*Vijk(qk,qj,qi)
                           End do
                           tmp2=tmp2 + Dj(qj)*tmp3
                        End do
                        tmp=tmp + Di(qi)*tmp2
                     End do
                     Epot=Epot+tmp

                     Deallocate(Vijk,Di,Dj,Dk)
                     m=m+1
                     if(m==m0) exit

                  endif
               End do

               if(m<m0) then
                  Do n=1,nQ3
                     mi=mQ3(1,n)
                     mj=mQ3(2,n)
                     mk=mQ3(3,n)
                     if((mi==mode(2) .and. mj==mode(1)) .or. &
                        (mi==mode(2) .and. mk==mode(1)) .or. &
                        (mj==mode(2) .and. mk==mode(1))) then

                        Epot=Epot + qav(1,mi)*qav(1,mj)*qav(1,mk)*coeff3(1,n) &
                                  + qav(2,mi)*qav(1,mj)*qav(1,mk)*coeff3(2,n) &
                                  + qav(1,mi)*qav(2,mj)*qav(1,mk)*coeff3(3,n) &
                                  + qav(1,mi)*qav(1,mj)*qav(2,mk)*coeff3(4,n) 

                        m=m+1
                        if(m==m0) exit
                     endif
                  End do
               endif

               if(MR==3) goto 300

           300 Continue
               Vcor_getHmat0=Epot

        !-[ 3-mode exc. ]-------------------------------------------------------- 
            Case(3) 

             ! 3-mode terms
               m=0
               Do n=1,nS3
                  mi=mS3(1,n)
                  mj=mS3(2,n)
                  mk=mS3(3,n)
                  if(mi==mode(3) .and. mj==mode(2) .and. mk==mode(1)) then

                     Call PEGrid_getnGrid3(nGi,nGj,nGk,n)
                     Allocate(Vijk(nGk,nGj,nGi))
                     Allocate(Di(nGi),Dj(nGj),Dk(nGk))
                     Call PEGrid_getV3(n,Vijk)

                     Call Vcor_getXwfn(mi,nGi,mCnf(mi),nCnf(mi),Di)
                     Call Vcor_getXwfn(mj,nGj,mCnf(mj),nCnf(mj),Dj)
                     Call Vcor_getXwfn(mk,nGk,mCnf(mk),nCnf(mk),Dk)

                     tmp=0.D+00
                     Do qi=1,nGi
                        tmp2=0.D+00
                        Do qj=1,nGj
                           tmp3=0.D+00
                           Do qk=1,nGk
                              tmp3=tmp3 + Dk(qk)*Vijk(qk,qj,qi)
                           End do
                           tmp2=tmp2 + Dj(qj)*tmp3
                        End do
                        tmp=tmp + Di(qi)*tmp2
                     End do
                     Epot=Epot+tmp

                     Deallocate(Vijk,Di,Dj,Dk)
                     m=1
                     exit

                  endif
               End do

               if(m==0) then
                  Do n=1,nQ3
                     mi=mQ3(1,n)
                     mj=mQ3(2,n)
                     mk=mQ3(3,n)
                     if(mi==mode(3) .and. mj==mode(2) .and. mk==mode(1)) then

                        Epot=Epot + qav(1,mi)*qav(1,mj)*qav(1,mk)*coeff3(1,n) &
                                  + qav(2,mi)*qav(1,mj)*qav(1,mk)*coeff3(2,n) &
                                  + qav(1,mi)*qav(2,mj)*qav(1,mk)*coeff3(3,n) &
                                  + qav(1,mi)*qav(1,mj)*qav(2,mk)*coeff3(4,n) 

                        m=1
                        exit
                     endif
                  End do
               endif

               if(MR==3) goto 400


           400 Continue
               Vcor_getHmat0=Epot

         End select

         if(nQ1/=0 .or. nQ2/=0 .or. nQ3/=0) then
            Call spr_memdealloc(dble(Nfree)*4.D+00*8.D+00)
            Deallocate(qav)
         endif

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Vcor_rVSCF(ierr)

      USE Vcor_mod

      Implicit None

         Integer :: i,j,k,l,nGi,ierr
         Real(8), allocatable :: Cwfn(:,:),Ene(:)

         Character :: wftyp*4

         ierr=0
         Read(rwf,end=100) wftyp
         Read(rwf) ref
         Do i=1,Nfree
            Read(rwf) nGi
            if(nGi==0) cycle
            Allocate(Cwfn(nGi,nGi),Ene(nGi))
            Read(rwf) Cwfn
            Read(rwf) Ene
            Call Modal_setCwfn(i,Cwfn)
            Call Modal_setEne(i,Ene)
            Deallocate(Cwfn,Ene)
         End do

         zp=.true.
         Do i=1,Nfree
            if(ref(i)/=0) then
               zp=.false. 
               exit
            endif
         End do

         tar=ref

         return

     100 Continue
         ierr=-1
         return

      !--------------------------------------------------------------------------

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
