!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/03
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Module Spfn2_mod 

        USE Vib_mod, ONLY: mS2,get_fname2

        !  ni,nj   :: Number of grid points
        Integer :: ni,nj

        !  title  :: Title of the PEF
        Character :: title*40
 
        !  qi,qj  :: Grid points
        !  V2  :: Potential energy values at the grid points
        Real(8), dimension(:), allocatable :: qi,qj
        Real(8), dimension(:,:), allocatable :: V2,V2s

      End module
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Spfn2_Gen1(Iopt,nn)

      USE Spfn2_mod

      Implicit None

         Integer :: Iopt,nn,mi,mj

         mi=mS2(1,nn)
         mj=mS2(2,nn)
         Call Spfn2_Gen2(Iopt,mi,mj)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Spfn2_Gen2(Iopt,mi,mj)

      USE Spfn2_mod

      Implicit None

         Integer :: Iopt,mi,mj

         Integer :: i,j,ifl,spr_memalloc
         Character :: fp*80

         if(mi<mj) then
            i=mi; mi=mj; mj=i
         endif

         Call get_fname2(mi,mj,fp)
         Call file_indicator(12,ifl)
         Open(unit=ifl,file=fp,status='OLD')
         Read(ifl,'(a)') title

         Read(ifl,*)
         Read(ifl,*) nj,ni
         i=spr_memalloc(-1,dble(ni*nj+ni+nj)*8.D+00)
         Allocate(qi(ni),qj(nj),V2(nj,ni))
         Read(ifl,*)
         Do i=1,ni
         Do j=1,nj
            Read(ifl,*) qj(j),qi(i),V2(j,i)
         End do
         End do
         Close(ifl)

         if(Iopt==1) then
            i=spr_memalloc(-1,dble(ni*nj)*8.D+00)
            Allocate(V2s(nj,ni))
            Call DSplie2(qj,qi,V2,nj,ni,V2s)
         endif

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Spfn2_Dispose()

      USE Spfn2_mod

      Implicit None

      Real(8) :: memsz

         memsz=dble(ni*nj+ni+nj)*8.D+00
         Call spr_memdealloc(memsz)
         Deallocate(qi,qj,V2)

         if(allocated(V2s)) then
            memsz=dble(ni*nj)*8.D+00
            Call spr_memdealloc(memsz)
            Deallocate(V2s)
         endif

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function Spfn2(qqi,qqj)

      USE Spfn2_mod

      Implicit None

      Real(8) :: Spfn2,qqi,qqj

         Call DSplin2(qj,qi,V2,V2s,nj,ni,qqj,qqi,Spfn2)

      End Function
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn2_getnGrid(nni,nnj)

      USE Spfn2_mod

      Implicit None

      Integer :: nni,nnj

         nni=ni; nnj=nj

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn2_getGrid(qqi,qqj)

      USE Spfn2_mod

      Implicit None

      Real(8) :: qqi(ni),qqj(nj)

         qqi=qi
         qqj=qj

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn2_getTitle(ch)

      USE Spfn2_mod

      Implicit None

      Character(*) :: ch

         ch=Title

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
