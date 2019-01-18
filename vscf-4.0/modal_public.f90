!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/04/03
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_Construct(ierr)

      USE Modal_mod

      Implicit None

         Integer :: ierr,spr_memalloc
         Integer :: rCHO,i,j,k,mHO,mCHO
         Real(8), allocatable :: CHOi(:,:)

         if(.not. CHODVR) return

       ! Setup ptCHO
         ierr=spr_memalloc(1,dble(Nfree)*4.D+00)
         if(ierr<0) return
         Allocate(ptCHO(Nfree))

       ! Read CHO
         memsz=0
         Call file_indicator(50,rCHO)
         Open(rCHO,file='cho-r.wfn',status='OLD',form='UNFORMATTED')
         Do i=1,Nfree
            Read(rCHO) mHO,mCHO
            if(mCHO==0) cycle
            Allocate(CHOi(mHO,mCHO))
            Read(rCHO) CHOi
            Deallocate(CHOi)

            nHO(i)=mHO
            if(nCHO(i)>mCHO) nCHO(i)=mCHO
            ptCHO(i)=memsz
            memsz=memsz + nHO(i)*nCHO(i)

         End do

       ! Setup CHO
         ierr=spr_memalloc(1,dble(memsz)*8.D+00)
         if(ierr<0) return
         Allocate(CHO(memsz))

         Rewind(rCHO)

         Do i=1,Nfree
            Read(rCHO) mHO,mCHO
            if(mCHO==0) cycle
            Allocate(CHOi(mHO,mCHO))
            Read(rCHO) CHOi
            Call Modal_setCHO(i,CHOi(:,1:nCHO(i)))
            Deallocate(CHOi)

         End do
         Close(rCHO)

         !Write(6,'(3i4)') nCHO
         !Write(6,'(3i4)') nHO
         !Write(6,'(3i4)') ptCHO
         !Do i=1,Nfree
         !   Allocate(CHOi(nHO(i),nCHO(i)))
         !   Call Modal_getCHO(i,CHOi)
         !   Do j=1,nCHO(i)
         !      write(6,'(11f8.4)') CHOi(:,j)
         !   End do
         !   write(6,*)
         !   Deallocate(CHOi)
         !End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Modal_open(ierr)

      USE Modal_mod

      Implicit None

         Integer :: ierr

         Integer :: i,j,k,l
         Integer :: isz,spr_memalloc

       ! Allocate typeTable1, Ene
         i=idx1(Nfree)+nCHO(Nfree)
         ierr=spr_memalloc(1,dble(i)*12.D+00)
         if(ierr<0) return
         allocate(type_Table1(i),Ene(i))

       ! Initialize
         Ntype=0
         memsz=0
         type_Table1=-1

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Modal_add(mode,nGrid)

      USE Modal_mod

      Implicit None

         Integer :: mode,nGrid,pt,nG2

         pt=idx1(mode)+nGrid
         if(type_Table1(pt) >= 0 ) return

         type_Table1(pt)=memsz

         nG2=nGrid*nGrid
         memsz=memsz + nG2

         if(nGrid /= nCHO(mode)) then
            memsz=memsz + nGrid
         else
            memsz=memsz + nGrid*4
         endif

         Ntype=Ntype+1

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Modal_Close(ierr)

      USE Modal_mod

      Implicit None

         Integer :: ierr
         Integer :: spr_memalloc
         Integer :: i,j,k, pt
         Real(8), dimension(:,:), allocatable :: Cwfn

         !dbg Do i=1,Nfree-1
         !dbg    write(6,'(11i5)') type_Table1(idx1(i)+1:idx1(i+1))
         !dbg End do
         !dbg write(6,'(11i5)') type_Table1(idx1(Nfree)+1:)
         !dbg write(6,'(i5)') memsz

         ierr=spr_memalloc(1,dble(memsz)*8.D+00)
         if(ierr<0) return
         allocate(block(memsz))

         i=idx2(Nfree)+nCHO(Nfree)*nCHO(Nfree)
         ierr=spr_memalloc(1,dble(i)*8.D+00)
         if(ierr<0) return
         allocate(Cwfn0(i))

         i=3*Ntype
         ierr=spr_memalloc(1,dble(i)*4.D+00)
         if(ierr<0) return
         allocate(type_Table2(3,Ntype))

         k=1
         Do i=1,Nfree
            pt=idx1(i)
            Do j=1,nCHO(i)
               if(type_Table1(pt+j)<0) cycle

               type_Table2(1,k)=i
               type_Table2(2,k)=j
               type_Table2(3,k)=type_Table1(pt+j)
               k=k+1

            End do

         End do
         !dbg write(6,*) Ntype
         !dbg write(6,'(4i5)') (i,type_Table2(:,i),i=1,Ntype)

         Do i=1,Nfree
            Allocate(Cwfn(nCHO(i),nCHO(i)))
            Cwfn=0.D+00
            Do j=1,nCHO(i)
               Cwfn(j,j)=1.D+00
            End do
            Call Modal_setCwfn(i,Cwfn)
            Deallocate(Cwfn)
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_Destruct()

      USE Modal_mod

      Implicit None

         Real(8) :: rsz

         Call Modal_DumpCHO()

         rsz=dble(size(type_Table1))*4.D+00
         Call spr_memdealloc(rsz)
         Deallocate(type_Table1)

         rsz=dble(size(type_Table2))*4.D+00
         Call spr_memdealloc(rsz)
         Deallocate(type_Table2)

         rsz=dble(size(block))*8.D+00
         Call spr_memdealloc(rsz)
         Deallocate(block)

         rsz=dble(size(Ene))*8.D+00
         Call spr_memdealloc(rsz)
         Deallocate(Ene)

         rsz=dble(size(Cwfn0))*8.D+00
         Call spr_memdealloc(rsz)
         Deallocate(Cwfn0)

         if(CHODVR) then
            rsz=dble(size(CHO))*8.D+00
            Call spr_memdealloc(rsz)
            Deallocate(CHO)
         endif

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_genDVR_ALL()

      USE Modal_mod

      Implicit None

         Integer :: i,j,k,pt

         Do k=1,Ntype
            i=type_Table2(1,k)
            j=type_Table2(2,k)
            Call genDVR(i,j)
         End do

         Contains

         Subroutine genDVR(mode,nGrid)
   
         Implicit None

            Integer :: mode,nGrid
            Real(8) :: xdvr(nGrid*nGrid),qq(nGrid),q2(nGrid),q3(nGrid),q4(nGrid)

            Integer :: i,j,k,l


            Call Modal_genDVR(mode,nGrid,xdvr,qq)

            i=type_Table1(idx1(mode)+nGrid)
            j=i+nGrid*nGrid
            block(i+1:j)=xdvr

            if(nGrid/= nCHO(mode)) then
               block(j+1:j+nGrid)=qq
            else
               Do i=1,nGrid
                  q2(i)=qq(i)*qq(i)
                  q3(i)=q2(i)*qq(i)
                  q4(i)=q2(i)*q2(i)
               End do
               block(j+1:j+nGrid)=qq
               j=j+nGrid
               block(j+1:j+nGrid)=q2
               j=j+nGrid
               block(j+1:j+nGrid)=q3
               j=j+nGrid
               block(j+1:j+nGrid)=q4
               !dbg write(6,'(i4)') j+nGrid
            endif

         End subroutine

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Modal_genDVR(mode,nGrid,xdvr,qq)

      USE Modal_mod

      Implicit None

         Integer :: mode,nGrid
         Real(8) :: xdvr(nGrid,nGrid),qq(nGrid)

         Integer :: i,j,k,l

         Real(8), dimension(nGrid,nGrid) :: xx
         Real(8), dimension(:,:), allocatable :: xho,CHOm

        !-------------------------------------
        ! >>  Quadrature points / Angs(amu)1/2
        !-------------------------------------

         if(.not. CHODVR) then
            Call HO_xmat(nGrid,omegaf(mode),xx)

         else
            allocate(xho(nHO(mode),nHO(mode)))
            Call HO_xmat(nHO(mode),omegaf(mode),xho)

            allocate(CHOm(nHO(mode),nCHO(mode)))
            Call Modal_getCHO(mode,CHOm)

            ! xx = CHO*xho*CHO
            Do i=1,nGrid
               xx(i,i)=0.D+00
               Do k=1,nHO(mode)
               Do l=1,nHO(mode)
                  xx(i,i)=xx(i,i) + CHOm(l,i)*xho(l,k)*CHOm(k,i)
               End do
               End do
            End do

            Do i=1,nGrid
            Do j=1,i-1
               xx(j,i)=0.D+00
               Do k=1,nHO(mode)
               Do l=1,nHO(mode)
                  xx(j,i)=xx(j,i) + CHOm(l,j)*xho(l,k)*CHOm(k,i)
               End do
               End do
               xx(i,j)=xx(j,i)
            End do
            End do
            Deallocate(xho,CHOm)

         endif
         !Do j=0,l-1
         !   write(6,'(20f8.4)') xx(:,j)
         !End do

         Call diag2(nGrid,nGrid,xx,xdvr,qq)
         !Call huckeler(nGrid,nGrid,xx,qq,xdvr)

         !dbg write(6,*)
         !dbg write(6,*) mode,nGrid
         !dbg Do i=1,nGrid
         !dbg    write(6,'(12f10.4)') qq(i),xdvr(:,i)
         !dbg End do

        !-------------------------------------


      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_setCwfn(mode,Cwfn)

      USE Modal_mod

      Implicit None

         Integer :: mode, nG, nG2, pt
         Real(8) :: Cwfn(nCHO(mode)*nCHO(mode))

         nG=nCHO(mode)
         nG2=nG*nG
         pt=idx2(mode)
         Cwfn0(pt+1:pt+nG2)=Cwfn

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_getCwfn(mode,Cwfn)

      USE Modal_mod

      Implicit None

         Integer :: mode,nG,nG2,pt
         Real(8) :: Cwfn(nCHO(mode)*nCHO(mode))

         nG=nCHO(mode)
         nG2=nG*nG
         pt=idx2(mode)
         Cwfn=Cwfn0(pt+1:pt+nG2)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_setEne(mode,Ein)

      USE Modal_mod

      Implicit None

         Integer :: mode, pt
         Real(8) :: Ein(nCHO(mode))

         pt=idx1(mode)
         Ene(pt+1:pt+nCHO(mode))=Ein

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_getEne(mode,Ein)

      USE Modal_mod

      Implicit None

         Integer :: mode, pt
         Real(8) :: Ein(nCHO(mode))

         pt=idx1(mode)
         Ein=Ene(pt+1:pt+nCHO(mode))

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_getEne_i(mode,qi,Ein)

      USE Modal_mod

      Implicit None

         Integer :: mode, qi
         Real(8) :: Ein

         Ein=Ene(idx1(mode)+qi+1)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_getxdvr(mode,nGrid,xdvr)

      USE Modal_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: xdvr(nGrid*nGrid)

         nGrid2=nGrid*nGrid
         pt=type_Table1(idx1(mode)+nGrid)
         xdvr=block(pt+1:pt+nGrid2)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_getQ(mode,nGrid,qq)

      USE Modal_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: qq(nGrid)

         nGrid2=nGrid*nGrid
         pt=type_Table1(idx1(mode)+nGrid)+nGrid2
         qq=block(pt+1:pt+nGrid)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Modal_getQ4(mode,qq)

      USE Modal_mod

      Implicit None

         Integer :: mode,nGrid,nGrid2,pt
         Real(8) :: qq(4*nCHO(mode))

         nGrid=nCHO(mode)
         nGrid2=nGrid*nGrid
         pt=type_Table1(idx1(mode)+nGrid)+nGrid2
         qq=block(pt+1:pt+4*nGrid)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Modal_setCHO(m,CHOm)

      USE Modal_mod

      Implicit None

         Integer :: m
         Real(8) :: CHOm(nHO(m)*nCHO(m))

         CHO(ptCHO(m)+1:ptCHO(m)+nHO(m)*nCHO(m))=CHOm

      End subroutine
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Modal_getCHO(m,CHOm)

      USE Modal_mod

      Implicit None

         Integer :: m
         Real(8) :: CHOm(nHO(m)*nCHO(m))

         CHOm=CHO(ptCHO(m)+1:ptCHO(m)+nHO(m)*nCHO(m))

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Modal_gettypeTable1(table)

      USE Modal_mod

      Implicit None

         Integer :: table(idx1(Nfree)+nCHO(Nfree))

         table=type_Table1

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Modal_DumpCHO()

      USE Modal_mod

      Implicit None

         Integer :: i,j,k,l
         Integer :: wCHO

         Real(8), allocatable :: Cwfn(:,:),CHOn(:,:),CHOm(:,:)

         Call file_indicator(40,wCHO)
         Open(wCHO,file='cho-w.wfn',status='unknown',form='UNFORMATTED')

         if(CHODVR) then
            Do i=1,Nfree
               write(wCHO) nHO(i),nCHO(i)

               if(nCHO(i)==0) cycle
               Allocate(Cwfn(nCHO(i),nCHO(i)))
               Allocate(CHOn(nHO(i),nCHO(i)),CHOm(nHO(i),nCHO(i)))
               Call Modal_getCwfn(i,Cwfn)
               Call Modal_getCHO(i,CHOn)

               Do j=1,nCHO(i)
               Do k=1,nHO(i)
                  CHOm(k,j)=0.D+00
                  Do l=1,nCHO(i)
                     CHOm(k,j)=CHOm(k,j) + CHOn(k,l)*Cwfn(l,j)
                  End do
               End do
               End do

               write(wCHO) CHOm

               Deallocate(Cwfn,CHOn,CHOm)

            End do
         else
            Do i=1,Nfree
               write(wCHO) nCHO(i),nCHO(i)

               if(nCHO(i)==0) cycle
               Allocate(Cwfn(nCHO(i),nCHO(i)))
               Call Modal_getCwfn(i,Cwfn)

               write(wCHO) Cwfn

               Deallocate(Cwfn)

            End do

         endif

         Close(wCHO)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
