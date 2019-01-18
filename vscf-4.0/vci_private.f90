!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/07/19
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Module Vci_mod

        USE Vcor_mod

        ! -- (FILE description) -------------------------------------------------
        !     dump :: if true, dump VCI wfn

        Logical :: dump

        ! -- (CI basis parameters) ----------------------------------------------
        !     Nst                 :: The number of CI states to solve
        !     CI space selection
        !      o maxExc(Nfree)    :: Max quantum number of excitation
        !      o maxSum           :: Max sum of quantum number
        !      o nCUP             :: Max mode coupling
        !      o nCI              :: Max number of VSCF configurations

        Integer :: Nst, maxSum, nCUP, nCI
        Integer, dimension(:),   allocatable:: maxExc
        Double Precision, parameter :: pi=3.14159265358979323846

        ! -- (VCI wfn) ----------------------------------------------------------
        !     CIwfn(nCI*Nstate)    :: CI coefficient matrix
        !     CIene(nCI)           :: CI energies

        Real(8), dimension(:), allocatable :: CIwfn, CIene

      CONTAINS


      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine setCIbasis(ierr)

      Implicit None

         Integer :: ierr,in,io,spr_memalloc

         Integer :: i,j,k,ist
         Integer :: jtmp,ktmp,jold,kold
         Real(8) :: tmp,Ei,Ej,memsz
         Real(8), parameter :: H2wvn=2.194746D+05

         Integer :: lb(Nfree)

         Integer, dimension(:), allocatable :: im,Lsum,Lcup,nP
         Integer, dimension(:,:), allocatable :: Label
         Real(8), dimension(:), allocatable :: Ene,DelE
         Real(8), dimension(:,:), allocatable :: E1wfn

         Integer, dimension(:), allocatable :: mode,vv


         ierr=spr_memalloc(1,dble(Nfree)*dble(maxCHO)*8.D+00)
         if(ierr<0) return
         allocate(E1wfn(maxCHO,Nfree))

         if(lvscf) then

           ! ONE-TO-ALL VCI/VSCF

           Call Vcor_rVSCF(i)
           Do i=1,Nfree
              Call Modal_getEne(i,E1wfn(:,i))
           End do
           E1wfn=E1wfn*H2wvn
           !dbg Do i=1,Nfree
           !dbg    write(6,*)
           !dbg    write(6,'(f12.3)') (E1wfn(j,i),j=1,nCHO(i))
           !dbg End do

           Do i=2,Nfree
              tmp=E1wfn(1,i)
              Do j=1,i-1
                 if(abs(tmp-E1wfn(1,j))<1.D-02) then
                    Do k=1,nCHO(j)
                       E1wfn(k,j)=E1wfn(k,j)+1.D+00*dble(k)
                    End do
                 endif
              End do
           End do


         else

           ! ONE-TO-ALL VCI/HO

           ! Set the reference
           ref=0; tar=0

           Do i=2,Nfree
              tmp=omegaf(i)
              Do j=1,i-1
                 if(abs(tmp-omegaf(j))<1.D-02) then
                    omegaf(j)=omegaf(j)+1.D+00
                 endif
              End do
           End do

           Do i=1,Nfree
              tmp=omegaf(i)
              Do j=1,nCHO(i)
                 E1wfn(j,i)=(dble(j)-0.5D+00)*tmp
              End do
              Do j=nCHO(i)+1,maxCHO
                 E1wfn(j,i)=0.D+00
              End do
           End do

         endif
         !dbg Do i=1,Nfree
         !dbg Do j=1,nCHO(i)
         !dbg    write(6,'(2i4,f10.2)') i,j-1,E1wfn(j,i)-E1wfn(1,i)
         !dbg End do
         !dbg End do

         memsz=dble(nCI)*16.D+00+dble(Nfree)*12D+00
         ierr=spr_memalloc(1,memsz)
         if(ierr < 0) return
         Allocate(Ene(nCI),im(Nfree),DelE(Nfree),lsum(nCI),lcup(nCI))

         ierr=spr_memalloc(1,dble(nCI)*dble(Nfree)*4.D+00)
         if(ierr<0) return
         Allocate(Label(Nfree,nCI))

         Allocate(nP(nCUP))

         label(:,1)=0
         lsum(1)=0
         lcup(1)=0
         Ene(1)=0.D+00
         Do i=1,Nfree
            Ene(1)=Ene(1)+E1wfn(1,i)
         End do
         nP=0

         jold=-1;kold=-1
         ist=2
         Do while(ist<=nCI) 
            Ei=1.D+15
            Do j=1,ist-1
               if(lsum(j)>=maxSum) cycle
               im=label(:,j)+1

               if(lcup(j) < nCUP) then
                  Do k=1,Nfree
                     if(im(k)<maxExc(k)) then
                        DelE(k)=E1wfn(im(k)+1,k)-E1wfn(im(k),k)
                     else
                        DelE(k)=0.D+00
                     endif
                  End do
               else
                  Do k=1,Nfree
                     if(im(k)<maxExc(k) .and. Label(k,j)/=ref(k)) then
                        DelE(k)=E1wfn(im(k)+1,k)-E1wfn(im(k),k)
                     else
                        DelE(k)=0.D+00
                     endif
                  End do
               endif

               tmp=Ene(j)+MaxVal(DelE)
               if(tmp<Ene(ist-1)) cycle

               Do k=1,Nfree
                  Ej=Ene(j)+DelE(k)
                  if(Ej<Ei.and.(Ej-Ene(ist-1))>1.D-04) then 
                     Ei=Ej
                     jtmp=j
                     ktmp=k
                  endif
               End do
            End do
            im=label(:,jtmp)+1
            label(:,ist)=label(:,jtmp)
            label(ktmp,ist)=label(ktmp,jtmp)+1

            if(jtmp==jold .and. ktmp==kold) then
                nCI=ist-1
                exit
            endif
            jold=jtmp
            kold=ktmp

            lsum(ist)=0
            lcup(ist)=0
            Do j=1,Nfree
               lsum(ist)=lsum(ist)+label(j,ist)
               if(label(j,ist)/=ref(j)) lcup(ist)=lcup(ist)+1
            End do
            Ene(ist)=Ene(jtmp)+E1wfn(im(ktmp)+1,ktmp)-E1wfn(im(ktmp),ktmp)
            !dbg write(6,'(3i5)') im
            !dbg write(6,'(2i5,''|'',3i5,f12.4)') jtmp,ktmp,label(:,i),Ene(i)

            nP(lcup(ist))=nP(lcup(ist))+1
            ist=ist+1
         End do

         if(lvscf .and. .not. zp) then
            Do i=1,nCI
               Do j=1,Nfree
                  k=0
                  if(Label(j,i)/=ref(j)) then
                     k=1
                     exit
                  endif
                  if(k==0) then
                     lb=Label(:,i)
                     Label(:,i)=Label(:,1)
                     Label(:,1)=lb
                  endif
               End do
            End do
         endif

         !dbg  write(6,*) nCI
         !dbg  Do i=1,nCI
         !dbg     write(6,'(f12.1,7i3)') Ene(i),Label(:,i)
         !dbg  End do
         !dbg  write(6,*) nP
         !dbg  write(6,'(6i2)') ref

         Allocate(mode(nCUP),vv(nCUP))
         Call Vcor_setup_qCnf(nCUP,nCI-1)
         Do i=1,Nci

            k=1
            lcup(i)=0
            Do j=Nfree,1,-1
               if(Label(j,i)/=ref(j)) then
                  mode(k)=j
                  vv(k)=Label(j,i)-ref(j)
                  k=k+1
                  lcup(i)=lcup(i)+1
               endif
            End do
            if(lcup(i)==0) cycle

            Call Vcor_add_qCnf(lcup(i),mode(1:lcup(i)),vv(1:lcup(i)))
            !dbg write(6,'(10i3)') lcup(i),mode(1:lcup(i)),vv(1:lcup(i))
         End do
         Deallocate(mode,vv)

         Call spr_memdealloc(memsz)
         Deallocate(Ene,im,DelE,lsum,lcup)
         Call spr_memdealloc(dble(maxCHO)*dble(Nfree)*4.D+00)
         Deallocate(E1wfn)
         Call spr_memdealloc(dble(nCI)*dble(Nfree)*4.D+00)
         Deallocate(Label,nP)

      End subroutine

      !--1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      End module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
