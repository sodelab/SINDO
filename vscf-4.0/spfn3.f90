!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/08
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Module Spfn3_mod 

        USE Vib_mod, ONLY: mS3,get_fname3

        !  ni,nj   :: Number of grid points
        Integer :: ni,nj,nk

        !  title  :: Title of the PEF
        Character :: title*40
 
        !  qi,qj,qk  :: Grid points
        !  V3  :: Potential energy values at the grid points
        Real(8), dimension(:), allocatable :: qi,qj,qk
        Real(8), dimension(:,:,:), allocatable :: V3,V3s

      End module
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Spfn3_Gen1(Iopt,nn)

      USE Spfn3_mod

      Implicit None

         Integer :: Iopt,nn,mi,mj,mk

         mi=mS3(1,nn)
         mj=mS3(2,nn)
         mk=mS3(3,nn)
         Call Spfn3_Gen2(Iopt,mi,mj,mk)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Spfn3_Gen2(Iopt,mi,mj,mk)

      USE Spfn3_mod

      Implicit None

         Integer :: Iopt,mi,mj,mk

         Integer :: i,j,k,ifl,spr_memalloc
         Character :: fp*80

         if(mi<mj) then
            i=mi; mi=mj; mj=i
         endif
         if(mj<mk) then
            j=mj; mj=mk; mk=j
         endif
         if(mi<mj) then
            i=mi; mi=mj; mj=i
         endif

         Call get_fname3(mi,mj,mk,fp)
         Call file_indicator(12,ifl)
         Open(unit=ifl,file=fp,status='OLD')
         Read(ifl,'(a)') title

         Read(ifl,*)
         Read(ifl,*) nk,nj,ni
         i=spr_memalloc(-1,dble(ni*nj*nk+ni+nj+nk)*8.D+00)
         Allocate(qi(ni),qj(nj),qk(nk),V3(nk,nj,ni))
         Read(ifl,*)
         Do i=1,ni
         Do j=1,nj
         Do k=1,nk
            Read(ifl,*) qk(k),qj(j),qi(i),V3(k,j,i)
         End do
         End do
         End do
         Close(ifl)

         if(Iopt==1) then
            i=spr_memalloc(-1,dble(ni*nj*nk)*8.D+00)
            Allocate(V3s(nk,nj,ni))
            Call DSplie3(qk,qj,qi,V3,nk,nj,ni,V3s)
         endif

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Spfn3_Dispose()

      USE Spfn3_mod

      Implicit None

      Real(8) :: memsz

         memsz=dble(ni*nj*nk+ni+nj+nk)*8.D+00
         Call spr_memdealloc(memsz)
         Deallocate(qi,qj,qk,V3)

         if(allocated(V3s)) then
            memsz=dble(ni*nj*nk)*8.D+00
            Call spr_memdealloc(memsz)
            Deallocate(V3s)
         endif

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function Spfn3(qqi,qqj,qqk)

      USE Spfn3_mod

      Implicit None

      Real(8) :: Spfn3,qqi,qqj,qqk

         Call DSplin3(qk,qj,qi,V3,V3s,nk,nj,ni,qqk,qqj,qqi,Spfn3)

      End Function
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn3_getnGrid(nni,nnj,nnk)

      USE Spfn3_mod

      Implicit None

      Integer :: nni,nnj,nnk

         nni=ni; nnj=nj; nnk=nk

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn3_getGrid(qqi,qqj,qqk)

      USE Spfn3_mod

      Implicit None

      Real(8) :: qqi(ni),qqj(nj),qqk(nk)

         qqi=qi
         qqj=qj
         qqk=qk

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn3_getTitle(ch)

      USE Spfn3_mod

      Implicit None

      Character(*) :: ch

         ch=Title

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
