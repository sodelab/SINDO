!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
! Ene, Grad, and Hess for the reference point
!
      Subroutine hBase_ref

      USE mkhs_mod

      Implicit None

        Integer :: i,idx
        Real(8) :: Ene,grd(Nfree),hess(Nfree,Nfree)

          
          Read(100,*,end=100)
          idx=2+Nfree/3+Nfree*Nfree/3
          Do i=1,idx
             Read(100,*)
          End do
          return
      100 Continue

          idx=-1
          Do i=1,Nfree
             if(abs(qrf(i)) > 1.D-09) idx=0
          End do
          Call RunHess(qrf,Ene,grd,hess,idx)
          Write(100,'(i4,f22.12)') 0,Ene
          Write(100,*) 'Grad / Hartree Angs-1 amu-1/2'
          Write(100,'(3f18.12)') grd
          Write(100,*) 'Hess / Hartree Angs-2 amu-1'
          Write(100,'(3f18.12)') hess
          Call myFlsh(100)

      End Subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
! Generate qi-a.com
!
!          a=1   +1di,       =2   -1di.
!
      Subroutine hBase_1MR

      USE mkhs_mod

        Implicit None

        Integer :: i,j,k,i1,i2,j1
        Integer :: Msym,Nir,Get_Num_of_Irrp,s(3)
        Integer, dimension(:), allocatable :: Lmode
        Real(8) :: dxi,x(3,Nat),qq(Nf)
        Real(8) :: Ene(3),grd(Nf,3),hess(Nf,Nf,3)

          if(m1 == -1) then
             Write(6,'('' -> Loop over modes'')')

             !i1=1; j1=0
             j1=(3+Nf/3+Nf*Nf/3)*2
             Do i=1,Nfree
                i1=i
                Read(100,*,end=110,err=110)
                Do j=1,j1
                   Read(100,*)
                End do
                Write(6,'(''  > Done for Mode '',i4)') Qmode(i) 
             End do
             Write(6,*)
             return

      110    Continue
             i2=Nfree

          else
             Write(6,'('' -> Loop over modes from ('',i3,'') to ('',i3,'')'')') &
             m1,m2
             i1=m1; i2=m2

          endif

          qq=0.D+00
          Do i=1,Nfree
             qq(Qmode(i))=qrf(i)
          End do

          Do i=i1,i2

             Call SymOp(i,Msym,s)
             dxi=delq(i)

             if(j1/=1) Write(100,'(''MODE='',i4)') Qmode(i)
             j1=0

             if(Msym==0) then
                Do j=1,2
                   Call MkGrid(i,j,dxi,qrf,qq)
                   Call RunHess(qq,Ene(j),grd(:,j),hess(:,:,j),0)
                End do

             else
                Call MkGrid(i,1,dxi,qrf,qq)
                Call RunHess(qq,Ene(1),grd(:,1),hess(:,:,1),0)
                !dbg Write(6,'(2i3)') i,Nir
                !dbg Write(6,'(5i3)') Lmode
                Call hcpy(i,Nf,Ene,grd,hess)

             endif

             Do j=1,2
                Write(100,'(i4,f22.12)') j,Ene(j)
                Write(100,*) 'Grad / Hartree Angs-1 amu-1/2'
                Write(100,'(3f18.12)') grd(:,j)
                Write(100,*) 'Hess / Hartree Angs-2 amu-1'
                Write(100,'(3f18.12)') hess(:,:,j)
                Call myFlsh(100)
             End do

             qq(Qmode(i))=qrf(i)
             Write(6,'(''  > Done for Mode= '',i4)') Qmode(i) 

          End do 
          Write(6,*)

          Contains 

          Subroutine MkGrid(i,j,dxi,qref,qq)

             Implicit None

             Integer :: i,j
             Real(8) :: dxi,qref(Nfree),qq(Nf)

                Select case(j)

                   case(1)
                      qq(i)= qref(i)+dxi
                   case(2)
                      qq(i)= qref(i)-dxi

                End select

             return

          End Subroutine

      End Subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine hcpy(mode,Nf,E,g,h)

      Implicit None

        Integer :: mode,Nf
        Real(8) :: E(2),g(Nf,2),h(Nf,Nf,2)

        Integer :: i,j
        Integer :: Ms0,Msi,Msij,getIrrep1,getIrrep2

          E(2)=E(1)
          g(:,2)=g(:,1)
          h(:,:,2)=h(:,:,1)

          Ms0=getIrrep1(mode)
          !dbg write(6,'(2i2)') mode,Ms0
          Do i=1,Nf
             Msi=getIrrep1(mode)
             if(Msi==Ms0) g(i,2)=-g(i,2)
          End do
          Do i=1,Nf
          Do j=1,i-1
             Msij=getIrrep2(i,j)
             if(Msij==Ms0) then
                h(i,j,2)=-h(i,j,2)
                h(j,i,2)=h(i,j,2)
             endif
          End do
          End do

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
