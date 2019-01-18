!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/03
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Module Spfn1_mod 

        USE Vib_mod, ONLY: mS1,get_fname1

        !  ni   :: Number of grid points
        Integer :: ni

        !  title  :: Title of the PEF
        Character :: title*40
 
        !  qi  :: Grid points
        !  V1  :: Potential energy values at the grid points
        Real(8), dimension(:), allocatable :: qi,V1,V1s

      End module
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Spfn1_Gen1(Iopt,nn)

      USE Spfn1_mod

      Implicit None

         Integer :: Iopt,nn,mode

         mode=mS1(nn)
         Call Spfn1_Gen2(Iopt,mode)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Spfn1_Gen2(Iopt,mode)

      USE Spfn1_mod

      Implicit None

         Integer :: Iopt,mode

         Integer :: i,ifl,spr_memalloc
         Character :: fp*80

         Call get_fname1(mode,fp)
         Call file_indicator(12,ifl)
         Open(unit=ifl,file=fp,status='OLD')
         Read(ifl,'(a)') title

         Read(ifl,*)
         Read(ifl,*) ni
         i=spr_memalloc(-1,dble(ni*2)*8.D+00)
         Allocate(qi(ni),V1(ni))
         Read(ifl,*)
         Do i=1,ni
            Read(ifl,*) qi(i),V1(i)
         End do
         Close(ifl)

         if(Iopt==1) then
            i=spr_memalloc(-1,dble(ni)*8.D+00)
            Allocate(V1s(ni))
            Call DSpline(qi,V1,ni,1.D+32,1.D+32,V1s)
         endif

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Spfn1_Dispose()

      USE Spfn1_mod

      Implicit None

      Real(8) :: memsz

         memsz=dble(ni*2)*8.D+00
         Call spr_memdealloc(memsz)
         Deallocate(qi,V1)

         if(allocated(V1s)) then
            memsz=dble(ni)*8.D+00
            Call spr_memdealloc(memsz)
            Deallocate(V1s)
         endif

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Function Spfn1(qq)

      USE Spfn1_mod

      Implicit None

      Real(8) :: Spfn1,qq

         Call DSplint(qi,V1,V1s,ni,qq,Spfn1) 

      End Function
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn1_getnGrid(nni)

      USE Spfn1_mod

      Implicit None

      Integer :: nni

         nni=ni

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn1_getGrid(qq)

      USE Spfn1_mod

      Implicit None

      Real(8) :: qq(ni)

         qq=qi

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Spfn1_getTitle(ch)

      USE Spfn1_mod

      Implicit None

      Character(*) :: ch

         ch=Title

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
