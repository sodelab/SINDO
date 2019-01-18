!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Module mkhs_mod

        Integer :: Nat,Nfree,MR,Nf,base,m1,m2
        Integer, dimension(:), allocatable :: Qmode
        Real(8), dimension(:), allocatable :: qrf,delq
        Character :: title*80

      End Module

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      PROGRAM main

      Implicit None

      Integer :: i
         Call v2ai_Const()
         Call Read_hs()
         Call mkqff_main
         Call Finalize

      End 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Read_hs

      USE mkhs_mod

      Implicit None

        Integer, parameter :: nmax=100
        Integer :: i,j,k,Get_Nfree,Get_Nat
        Integer :: Mode(nmax),fMode(nmax)
        Real(8) :: delx
        Real(8) :: Qref(nmax)
        Real(8), dimension(:), allocatable :: Omega,Omg

        Namelist /mkqff/Qref,delx,Mode,fMode,m1,m2,MR,base,title

!----------------------------------------------------------------------

        Write(6,100)
  100   Format(/,'--------------------(    MKHS MODULE   )--------------------',/)

! --    Read parameters

        Nat=Get_Nat()
        Nfree=Get_Nfree()
        Nf=Get_Nfree()
        MR=3
        delx=0.5D+00
        Qref=0.D+00
        base=0
        title=''
        Mode=-1
        fMode=-1
        m1=-1; m2=-1

        Rewind(5)
        Read(5,mkqff,end=10000)
      10000 Continue

! --    Set the modes of the system
        if(Mode(1)/=-1) then
           Nfree=0
           Do while(Mode(Nfree+1)/=-1)
              Nfree=Nfree+1
           End do
           allocate(Qmode(Nfree))
           Do i=1,Nfree
              Qmode(i)=mode(i)
           End do

        elseif(fMode(1)/=-1) then
           i=0
           Do while(fMode(i+1)/=-1)
              i=i+1
           End do
           Call isort(i,fMode(1:i))
           write(6,*) fMode(1:i)

           Nfree=Nf-i
           allocate(Qmode(Nfree))

           j=1;k=1
           Do i=1,Nf
              if(fMode(j)/=i) then
                 Qmode(k)=i
                 k=k+1
              else
                 j=j+1
              endif
           End do

        else
           allocate(Qmode(Nfree))
           Do i=1,Nfree
              Qmode(i)=i
           End do

        endif

        Write(6,105) Nfree
  105   Format(' -> Number of active mode:',i10)

! --    qq
        allocate(qrf(Nfree))
        qrf=Qref(1:Nfree)
        Write(6,*) '-> The reference point is set at:'
        Write(6,110) qrf
  110   Format('    ',3f12.6)

!
! --    delx
        Write(6,130) delx
        Allocate(Omega(Nfree),delq(Nfree))
        if(Nfree==Nf) then
           Call Get_Omega(Omega)
        else
           Allocate(Omg(Nf))
           Call Get_Omega(Omg)
           Do i=1,Nfree
              Omega(i)=Omg(Qmode(i))
           End do
        endif
        Call GetStpsz(Nfree,Omega,delx,delq)
        Deallocate(Omega)
        Do i=1,Nfree
           Write(6,135) Qmode(i),delq(i)
        End do
  130   Format(' -> Step size is set to: ',f10.4)
  135   Format('  > Mode= ',i3,':  ',f8.3,' / Angs(amu)^1/2')
 

! --    base
        if(base==0) then 
           Write(6,140) 
        else
           Write(6,141) 
        endif
  140   Format(' -> Numerical derivatives based on [  ENERGY  ]')
  141   Format(' -> Numerical derivatives based on [  HESSIAN  ]')


      Contains


         Subroutine isort(N,C)

         Implicit None

            Integer :: N,C(N)

            Integer :: i,j(1),k,tmp

               Do i=1,N
                  j=MinLoc(C(i:N))
                  k=j(1)+i-1

                  tmp=C(i)
                  C(i)=C(k)
                  C(k)=tmp

               End do

         End Subroutine

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mkqff_main()

      USE mkhs_mod

      Implicit None

      Integer :: Isym(3)
      Logical :: dpl,Get_dpl

        Isym=0
        dpl=Get_dpl()
        if(dpl) then 
           Write(6,*) 'warning warning warning warning'
           Write(6,*) 'Sorry dipole surface is not ready'
           Write(6,*) 'warning warning warning warning'
           dpl=.false.
        endif


        if(m1==-1) then
           Open(100,file='tmp1',status='unknown')
           Open(155,file='PES',status='unknown')!MK
           if(base==0) then
              Call eBase_ref
              Call eBase_1MR
              if(MR>1) Call eBase_2MR
              if(MR>2) Call eBase_3MR
              Close(100)
              Call read_eBase_data(MR,Nfree,delq,qrf,Qmode,title)
              Call eBase_hs(Isym,dpl)
           else
              Call hBase_ref
              Call hBase_1MR
              Close(100)
              Call read_hBase_data(Nfree,delq,qrf,title)
              Call hBase_hs
           endif

        else
           Open(100,file='tmp1',status='new')
           if(base==0) then
             ! Energy-only algorithm

           else
              if(m1==0) then 
                 Call hBase_ref
                 m1=1
              endif
              Call hBase_1MR
              Close(100)
           endif
        endif

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Finalize

      USE mkhs_mod

      Implicit None

        Deallocate(qrf,delq,Qmode)

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine GetStpsz(N,Omega,delx,delq)

      Implicit None

        Integer :: N,M
        Real(8) :: delx,Omega(N),delq(N)

        Integer :: i
        Real(8) :: PI,const  ! change planck's constant
        Real(8), parameter :: h=6.62606896D-34,    &  ! Planck's constant in J.s !Change it
                              c=2.99792458D+10, &  ! velocity of light in cm.s-1
                              Na=6.02214D+23       ! Avogadoro's constant

           PI=Acos(-1.0D+00)
           Do i=1,N
              const=2.0D+00*PI*SQRT(c*Omega(i)/h/Na*1.0D-23)
              delq(i)=delx/const
           End do
        return

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
! Generate ref.com
!
      Subroutine eBase_ref

      USE mkhs_mod

      Implicit None

        Integer :: i
        Real(8) :: Ene,di(3)

          Read(100,*,end=100)
          return

 !     100 Continue
      100 Backspace(100)
          Call Run0(qrf,Ene,di)
          Write(100,'(i4,f22.12,3f14.8)') 0,Ene,di
          Call myFlsh(100)

      End Subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Run0(qin,Ene,di)

      USE mkhs_mod

      Implicit None

        Integer :: i
        Real(8) :: qin(Nfree),q1(Nf),Ene,di(3)
        Logical :: Get_dpl

          q1=0.D+00
          Do i=1,Nfree
             q1(Qmode(i))=qin(i)
          End do
          if(Get_dpl()) then
             Call Run_ed(q1,Ene,di)
          else
             di=0.D+00
             Call Run_e(q1,Ene)
          endif


      End Subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
! Generate qi-a.com
!
!          a=1   +3di,       =4   -1di,
!            2   +2di,       =5   -2di,
!            3   +1di,       =6   -3di.
!
      Subroutine eBase_1MR

      USE mkhs_mod

        Implicit None

        Integer :: i,j,k,i1,j1
        Integer :: Msym,s(3)
        Real(8) :: dxi,x(3,Nat),qq(Nfree)
        Real(8) :: Ene(6),di(3,6)

          !i1=1; j1=0
          Do i=1,Nfree
             i1=i
             Read(100,*,end=110,err=110)
             Do j=1,6
                Read(100,*)
             End do
          End do
          return

      110 Continue

          qq=qrf

          Do i=i1,Nfree

             Call SymOp(Qmode(i),Msym,s)
             dxi=delq(i)

             if(j1/=1) Write(100,'(''MODE='',i4)') Qmode(i)
             j1=0

             if(Msym==0) then
                Do j=1,6
                   Call MkGrid(i,j,Nfree,dxi,qrf,qq)
                   Call Run0(qq,Ene(j),di(:,j))
                End do

             else
                Do j=1,3
                   Call MkGrid(i,j,Nfree,dxi,qrf,qq)
                   Call Run0(qq,Ene(j),di(:,j))
                   Ene(7-j)=Ene(j)
                   Do k=1,3
                      di(k,7-j)=di(k,j)*dble(s(k))
                   End do
                End do

             endif

             Write(100,'(i4,f22.12,3f14.8)') (j,Ene(j),di(:,j),j=1,6)
             Call myFlsh(100)

             qq(i)=qrf(i)

          End do 

          Contains 

          Subroutine MkGrid(i,j,Nfree,dxi,qref,qq)

             Implicit None

             Integer :: i,j,Nfree
             Real(8) :: dxi,qref(Nfree),qq(Nfree)

                Select case(j)

                   case(1)
                      qq(i)= qref(i)+dxi*3.D+00
                   case(2)
                      qq(i)= qref(i)+dxi*2.D+00
                   case(3)
                      qq(i)= qref(i)+dxi
                   case(4)
                      qq(i)= qref(i)-dxi
                   case(5)
                      qq(i)= qref(i)-dxi*2.D+00
                   case(6)
                      qq(i)= qref(i)-dxi*3.D+00

                End select

             return

          End Subroutine

      End Subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine eBase_2MR

      USE mkhs_mod

        Implicit None

        Integer :: i,j,k,k1,k2,l,ii,jj,kk
        Integer :: Mi,si(3),Mj,sj(3),t(3)
        Real(8) :: dxi,dxj,x(3,Nat),qq(Nfree)
        Real(8) :: Ene(12),di(3,12)

        Integer :: i1(12),i2(12),i3(12),i4(12)
        data i1 / 1,2,3,7,8,9, 4,5,6,10,11,12 /
        data i2 / 1,2,3,4,5,6, 9,8,7,12,11,10 /
        data i3 / 1,2,3,4,5,6, 12,11,10,9,8,7 /
        data i4 / 1,2,3, 4,5,6, 9,8,7, 12,11,10 /

          !ii=2; jj=1; kk=1
          Do i=2,Nfree
          ii=i
          Do j=1,i-1
             jj=j
             Read(100,*,end=110,err=110)
             Do k=1,12
                Read(100,*)
             End do
          End do
          End do
          return

      110 Continue

          qq=qrf

          Do i=ii,Nfree

             Call SymOp(Qmode(i),Mi,si)
             dxi=delq(i)

             Do j=jj,i-1

                Call SymOp(Qmode(j),Mj,sj)
                dxj=delq(j)

                !if(kk/=1) Write(100,'(''MODE='',2i4)') Qmode(i),Qmode(j)
                !kk=0
                Write(100,'(''MODE='',2i4)') Qmode(i),Qmode(j)

                if(Mi==0) then 
                   if(Mj==0) then
                      Do k=1,12
                         Call MkGrid(i,j,k,Nfree,dxi,dxj,qrf,qq)
                         Call Run0(qq,Ene(k),di(:,k))
                      End do

                   else
                      Do k=1,6
                         k1=i1(k)
                         Call MkGrid(i,j,k1,Nfree,dxi,dxj,qrf,qq)
                         Call Run0(qq,Ene(k1),di(:,k1))
                         k2=i1(6+k)
                         Call cpy(k1,k2,Ene,di,sj)
                      End do

                   endif

                else
                   if(Mj==0) then
                      Do k=1,6
                         k1=i2(k)
                         Call MkGrid(i,j,k1,Nfree,dxi,dxj,qrf,qq)
                         Call Run0(qq,Ene(k1),di(:,k1))
                         k2=i2(6+k)
                         Call cpy(k1,k2,Ene,di,si)
                      End do

                   else
                      if(Mi==Mj) then 
                         Do k=1,6
                            k1=i3(k)
                            Call MkGrid(i,j,k1,Nfree,dxi,dxj,qrf,qq)
                            Call Run0(qq,Ene(k1),di(:,k1))
                            k2=i3(6+k)
                            Call cpy(k1,k2,Ene,di,si)
                         End do

                      else
                         Do k=1,3
                            k1=i4(k)
                            Call MkGrid(i,j,k1,Nfree,dxi,dxj,qrf,qq)
                            Call Run0(qq,Ene(k1),di(:,k1))
                            k2=i4(3+k)
                            Call cpy(k1,k2,Ene,di,sj)
                            k2=i4(6+k)
                            Call cpy(k1,k2,Ene,di,si)
                            k2=i4(9+k)
                            Do l=1,3
                               t(l)=si(l)*sj(l)
                            End do
                            Call cpy(k1,k2,Ene,di,t)
                         End do
                      endif
                   endif

                end if

                Write(100,'(i4,f22.12,3f14.8)') (k,Ene(k),di(:,k),k=1,12)
                Call myFlsh(100)

                qq(j)=qrf(j)
             End do 

             qq(i)=qrf(i)
             jj=1
          End do 


          Contains

          !-   1  (+3di,+1dj), 4  (+3di,-1dj), 7  (-1di,+3dj), 10  (-1di,-3dj), 
          !-   2  (+1di,+1dj), 5  (+1di,-1dj), 8  (-1di,+1dj), 11  (-1di,-1dj), 
          !-   3  (+1di,+3dj), 6  (+1di,-3dj), 9  (-3di,+1dj), 12  (-3di,-1dj), 
          ! 
          Subroutine MkGrid(i,j,k,Nfree,dxi,dxj,qref,qq)

             Implicit None

             Integer :: i,j,k,Nfree
             Real(8) :: dxi,dxj,qref(Nfree),qq(Nfree)

                Select case(k)

                   case(1)
                      qq(i)= qref(i)+dxi*3.D+00
                      qq(j)= qref(j)+dxj 
                   case(2)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)+dxj 
                   case(3)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)+dxj*3.D+00 
                   case(4)
                      qq(i)= qref(i)+dxi*3.D+00
                      qq(j)= qref(j)-dxj 
                   case(5)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)-dxj 
                   case(6)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)-dxj*3.D+00 
                   case(7)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)+dxj*3.D+00 
                   case(8)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)+dxj 
                   case(9)
                      qq(i)= qref(i)-dxi*3.D+00
                      qq(j)= qref(j)+dxj 
                   case(10)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)-dxj*3.D+00 
                   case(11)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)-dxj 
                   case(12)
                      qq(i)= qref(i)-dxi*3.D+00
                      qq(j)= qref(j)-dxj 

                End select

             return

          End Subroutine

      End Subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine eBase_3MR

      USE mkhs_mod

        Implicit None

        Integer :: i,j,k,l,l1,l2,m,ii,jj,kk,ll
        Integer :: Mi,Mj,Mk,si(3),sj(3),sk(3),t(3)
        Real(8) :: dxi,dxj,dxk,x(3,Nat),qq(Nfree)
        Real(8) :: Ene(8),di(3,8)

        Integer :: n1(8),n2(8),n3(8),n4(8),n5(8),n6(8),n7(8),n8(8),n9(8), &
                   n10(8),n11(8),n12(8) 

        Integer :: getChr
        Character :: Sym*3,getSym*3

        data n1 / 1,2,5,6, 3,4,7,8 /
        data n2 / 1,2,5,6, 7,8,3,4 /
        data n3 / 1,2, 3,4, 5,6, 7,8 /
        data n4 / 1,3,5,7, 2,4,6,8 /
        data n5 / 1,3,5,7, 6,8,2,4 /
        data n6 / 1,3, 2,4, 5,7, 6,8 /
        data n7 / 1,2,5,6, 4,3,8,7 /
        data n8 / 1,2,5,6, 8,7,4,3 /
        data n9 / 1,2, 4,3, 5,6, 8,7 /
        data n10 / 1,5, 2,6, 3,7, 4,8 /
        data n11 / 1,5, 6,2, 3,7, 8,4 /
        data n12 / 1,5, 2,6, 7,3, 8,4 /


          !ii=3; jj=2; kk=1; ll=0
          Do i=3,Nfree
          ii=i
          Do j=2,i-1
          jj=j
          Do k=1,j-1
          kk=k
             Read(100,*,end=110,err=110)
             Do l=1,8
                Read(100,*)
             End do
          End do
          End do
          End do
          return

      110 Continue

          Sym=getSym()

          qq=qrf
          Do i=ii,Nfree

             Call SymOp(Qmode(i),Mi,si)
             dxi=delq(i)

             Do j=jj,i-1

                Call SymOp(Qmode(j),Mj,sj)
                dxj=delq(j)

                Do k=kk,j-1

                   Call SymOp(Qmode(k),Mk,sk)
                   dxk=delq(k)

                   !if(ll/=1) Write(100,'(''MODE='',3i4)') Qmode(i),Qmode(j),Qmode(k)
                   !ll=0
                   Write(100,'(''MODE='',3i4)') Qmode(i),Qmode(j),Qmode(k)

                   !---
                   if(Mi==0) then
                      if(Mj==0) then
                         if(Mk==0) then
                            Do l=1,8
                               Call MkGrid(i,j,k,l,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l),di(:,l))
                            End do

                         else
                            Do l=1,4
                               Call MkGrid(i,j,k,l,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l),di(:,l))
                               l1=l+4
                               Call cpy(l,l1,Ene,di,sk)
                            End do

                         endif

                      else
                         if(Mk==0) then
                            Do l=1,4
                               l1=n1(l)
                               Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l1),di(:,l1))
                               l2=n1(4+l)
                               Call cpy(l1,l2,Ene,di,sj)
                            End do
                         elseif(Mj==Mk) then
                            Do l=1,4
                               l1=n2(l)
                               Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l1),di(:,l1))
                               l2=n2(4+l)
                               Call cpy(l1,l2,Ene,di,sj)
                            End do
                         else
                            Do l=1,2
                               l1=n3(l)
                               Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l1),di(:,l1))
                               l2=n3(2+l)
                               Call cpy(l1,l2,Ene,di,sj)
                               l2=n3(4+l)
                               Call cpy(l1,l2,Ene,di,sk)
                               l2=n3(6+l)
                               Do m=1,3
                                  t(m)=sj(m)*sk(m)
                               End do
                               Call cpy(l1,l2,Ene,di,t)
                            End do
                         endif
                      endif
                   else
                      if(Mj==0) then
                         if(Mk==0) then
                            Do l=1,4
                               l1=n4(l)
                               Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l1),di(:,l1))
                               l2=n4(4+l)
                               Call cpy(l1,l2,Ene,di,si)
                            End do
                         elseif(Mk==Mi) then
                            Do l=1,4
                               l1=n5(l)
                               Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l1),di(:,l1))
                               l2=n5(4+l)
                               Call cpy(l1,l2,Ene,di,si)
                            End do
                         else
                            Do l=1,2
                               l1=n6(l)
                               Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(l1),di(:,l1))
                               l2=n6(2+l)
                               Call cpy(l1,l2,Ene,di,si)
                               l2=n6(4+l)
                               Call cpy(l1,l2,Ene,di,sk)
                               l2=n6(6+l)
                               Do m=1,3
                                  t(m)=si(m)*sk(m)
                               End do
                               Call cpy(l1,l2,Ene,di,t)
                            End do
                         endif
                      else 
                         if(Mj==Mi) then
                            if(Mk==0) then
                               Do l=1,4
                                  l1=n7(l)
                                  Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                                  Call Run0(qq,Ene(l1),di(:,l1))
                                  l2=n7(4+l)
                                  Call cpy(l1,l2,Ene,di,si)
                               End do
                            elseif(Mi==Mk) then
                               Do l=1,4
                                  l1=n8(l)
                                  Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                                  Call Run0(qq,Ene(l1),di(:,l1))
                                  l2=n8(4+l)
                                  Call cpy(l1,l2,Ene,di,si)
                               End do
                            else
                               Do l=1,2
                                  l1=n9(l)
                                  Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                                  Call Run0(qq,Ene(l1),di(:,l1))
                                  l2=n9(2+l)
                                  Call cpy(l1,l2,Ene,di,si)
                                  l2=n9(4+l)
                                  Call cpy(l1,l2,Ene,di,sk)
                                  l2=n9(6+l)
                                  Do m=1,3
                                     t(m)=si(m)*sk(m)
                                  End do
                                  Call cpy(l1,l2,Ene,di,t)
                               End do
                            endif
                         else
                            if(Mk==0) then
                               Do l=1,2
                                  l1=n10(l)
                                  Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                                  Call Run0(qq,Ene(l1),di(:,l1))
                                  l2=n10(2+l)
                                  Call cpy(l1,l2,Ene,di,si)
                                  l2=n10(4+l)
                                  Call cpy(l1,l2,Ene,di,sj)
                                  l2=n10(6+l)
                                  Do m=1,3
                                     t(m)=si(m)*sj(m)
                                  End do
                                  Call cpy(l1,l2,Ene,di,t)
                               End do
                            elseif(Mi==Mk) then
                               Do l=1,2
                                  l1=n11(l)
                                  Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                                  Call Run0(qq,Ene(l1),di(:,l1))
                                  l2=n11(2+l)
                                  Call cpy(l1,l2,Ene,di,si)
                                  l2=n11(4+l)
                                  Call cpy(l1,l2,Ene,di,sj)
                                  l2=n11(6+l)
                                  Do m=1,3
                                     t(m)=si(m)*sj(m)
                                  End do
                                  Call cpy(l1,l2,Ene,di,t)
                               End do
                            elseif(Mj==Mk) then
                               Do l=1,2
                                  l1=n12(l)
                                  Call MkGrid(i,j,k,l1,Nfree,dxi,dxj,dxk,qrf,qq)
                                  Call Run0(qq,Ene(l1),di(:,l1))
                                  l2=n12(2+l)
                                  Call cpy(l1,l2,Ene,di,si)
                                  l2=n12(4+l)
                                  Call cpy(l1,l2,Ene,di,sj)
                                  l2=n12(6+l)
                                  Do m=1,3
                                     t(m)=si(m)*sj(m)
                                  End do
                                  Call cpy(l1,l2,Ene,di,t)
                               End do
                            else
                               !<<<<
                               if(Sym=='D2H' .and. getChr(Mi,Mj,Mk)==1) then
                               Call MkGrid(i,j,k,1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(1),di(:,1))
                               Call cpy(1,2,Ene,di,t)

                               !<<<<
                               else
                               Call MkGrid(i,j,k,1,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(1),di(:,1))
                               Call MkGrid(i,j,k,2,Nfree,dxi,dxj,dxk,qrf,qq)
                               Call Run0(qq,Ene(2),di(:,2))

                               !<<<<
                               endif

                               Do m=1,3
                                  t(m)=si(m)*sj(m)
                               End do
                               Call cpy(1,4,Ene,di,t)
                               Call cpy(2,3,Ene,di,t)

                               Do m=1,3
                                  t(m)=si(m)*sk(m)
                               End do
                               Call cpy(1,6,Ene,di,t)
                               Call cpy(2,5,Ene,di,t)

                               Do m=1,3
                                  t(m)=sj(m)*sk(m)
                               End do
                               Call cpy(1,7,Ene,di,t)
                               Call cpy(2,8,Ene,di,t)

                            endif
                         endif
                      endif
                   endif

                   Write(100,'(i4,f22.12,3f14.8)') (l,Ene(l),di(:,l),l=1,8)
                   Call myFlsh(100)

                   qq(k)=qrf(k)
                End do 

                qq(j)=qrf(j)
                kk=1
             End do 

             qq(i)=qrf(i)
             jj=2
          End do 


          Contains

          !-  1  (+1di,+1dj,+1dk), 3  (+1di,-1dj,+1dk),
          !-  2  (-1di,+1dj,+1dk), 4  (-1di,-1dj,+1dk),
          !
          !-  5  (+1di,+1dj,-1dk), 7  (+1di,-1dj,-1dk), 
          !-  6  (-1di,+1dj,-1dk), 8  (-1di,-1dj,-1dk). 
          Subroutine MkGrid(i,j,k,l,Nfree,dxi,dxj,dxk,qref,qq)

             Implicit None

             Integer :: i,j,k,l,Nfree
             Real(8) :: dxi,dxj,dxk,qref(Nfree),qq(Nfree)

                Select case(l)

                   case(1)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)+dxj 
                      qq(k)= qref(k)+dxk 
                   case(2)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)+dxj 
                      qq(k)= qref(k)+dxk 
                   case(3)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)-dxj 
                      qq(k)= qref(k)+dxk 
                   case(4)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)-dxj 
                      qq(k)= qref(k)+dxk 
                   case(5)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)+dxj 
                      qq(k)= qref(k)-dxk 
                   case(6)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)+dxj 
                      qq(k)= qref(k)-dxk 
                   case(7)
                      qq(i)= qref(i)+dxi
                      qq(j)= qref(j)-dxj 
                      qq(k)= qref(k)-dxk 
                   case(8)
                      qq(i)= qref(i)-dxi
                      qq(j)= qref(j)-dxj 
                      qq(k)= qref(k)-dxk 

                End select

             return

          End Subroutine

      End Subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine cpy(i,j,Ene,di,s)

      Implicit None

        Integer :: i,j,s(3),m
        Real(8) :: Ene(*),di(3,*)

          Ene(j)=Ene(i)
          Do m=1,3
             di(m,j)=di(m,i)*s(m)
          End do

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
