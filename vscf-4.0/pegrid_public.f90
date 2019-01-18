!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/15
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_Construct(ierr)

      USE PEGrid_mod

      Implicit None

         Integer :: ierr

         Integer :: spr_memalloc
         Real(8) :: msz

       ! 1MR
         msz=dble(nS1)*8.D+0 + 4.D+000
         ierr=spr_memalloc(1,msz)
         if(ierr<0) return
         Allocate(nGrid1(nS1),ptV1(nS1+1))

       ! 2MR
         if(nS2 /=0) then
            msz=dble(nS2)*12.D+00 + 4.D+00
            ierr=spr_memalloc(1,msz)
            if(ierr<0) return
            Allocate(nGrid2(2,nS2),ptV2(nS2+1))
         endif

       ! 3MR
         if(nS3 /=0) then
            msz=dble(nS3)*12.D+00 + 4.D+00
            ierr=spr_memalloc(1,msz)
            if(ierr<0) return
            Allocate(nGrid3(3,nS3),ptV3(nS3+1))
         endif

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine PEGrid_Destruct()

      USE PEGrid_mod

      Implicit None

         Real(8) :: rsz1,rsz2

         rsz1=dble(size(nGrid1))+dble(size(nGrid2))+dble(size(nGrid3))
         rsz1=rsz1+dble(size(ptV1))+dble(size(ptV2))+dble(size(ptV3))
         rsz2=dble(size(V1))+dble(size(V2))+dble(size(V3))

         Call spr_memdealloc(rsz1*4.D+00+rsz2*8.D+00)
         if(nS1 /=0) then
            Deallocate(nGrid1,ptV1,V1)
         endif
         if(nS2 /=0) then
            Deallocate(nGrid2,ptV2,V2)
         endif
         if(nS3 /=0) then
            Deallocate(nGrid3,ptV3,V3)
         endif

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine PEGrid_onGrid(ierr)

      USE PEGrid_mod

      Implicit None

         Integer :: ierr, spr_memalloc
         Integer :: i,j,k,l

         if(nS1/=0) then
            ptV1(1)=0
            Do i=2,nS1+1
               ptV1(i)=ptV1(i-1) + nGrid1(i-1)
            End do
            i=ptV1(nS1+1)
            ierr=spr_memalloc(1,dble(i)*8.D+00)
            if(ierr<0) return
            Allocate(V1(i))

            Do i=1,nS1
               Call setV1(i,nGrid1(i),V1(ptV1(i)+1:ptV1(i+1)))
            End do

         endif

         if(nS2/=0) then
            ptV2(1)=0
            Do i=2,nS2+1
               ptV2(i)=ptV2(i-1) + nGrid2(1,i-1)*nGrid2(2,i-1)
            End do
            i=ptV2(nS2+1)
            ierr=spr_memalloc(1,dble(i)*8.D+00)
            if(ierr<0) return
            Allocate(V2(i))

            Do i=1,nS2
               Call setV2(i,nGrid2(1,i),nGrid2(2,i),V2(ptV2(i)+1:ptV2(i+1)))
            End do

         endif

         if(nS3/=0) then
            ptV3(1)=0
            Do i=2,nS3+1
               ptV3(i)=ptV3(i-1) + nGrid3(1,i-1)*nGrid3(2,i-1)*nGrid3(3,i-1)
            End do
            i=ptV3(nS3+1)
            ierr=spr_memalloc(1,dble(i)*8.D+00)
            if(ierr<0) return
            Allocate(V3(i))

            Do i=1,nS3
               Call setV3(i,nGrid3(1,i),nGrid3(2,i),nGrid3(3,i), &
                          V3(ptV3(i)+1:ptV3(i+1)))
            End do

         endif

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine PEGrid_setnGrid1(nGrid,n)

      USE PEGrid_mod

      Implicit None

         Integer :: n,nGrid

         nGrid1(n)=nGrid

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_setnGrid2(nGi,nGj,n)

      USE PEGrid_mod

      Implicit None

         Integer :: nGi,nGj,n

         nGrid2(1,n)=nGi
         nGrid2(2,n)=nGj

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_setnGrid3(nGi,nGj,nGk,n)

      USE PEGrid_mod

      Implicit None

         Integer :: nGi,nGj,nGk,n

         nGrid3(1,n)=nGi
         nGrid3(2,n)=nGj
         nGrid3(3,n)=nGk

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine PEGrid_getnGrid1(nGrid,n)

      USE PEGrid_mod

      Implicit None

         Integer :: n,nGrid

         nGrid=nGrid1(n)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_getnGrid2(nGi,nGj,n)

      USE PEGrid_mod

      Implicit None

         Integer :: nGi,nGj,n

         nGi=nGrid2(1,n)
         nGj=nGrid2(2,n)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_getnGrid3(nGi,nGj,nGk,n)

      USE PEGrid_mod

      Implicit None

         Integer :: nGi,nGj,nGk,n

         nGi=nGrid3(1,n)
         nGj=nGrid3(2,n)
         nGk=nGrid3(3,n)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_getV1(n,VV1)

      USE PEGrid_mod

      Implicit None

         Integer :: n,i,j
         Real(8) :: VV1(*)

         j=ptV1(n)
         Do i=1,nGrid1(n)
            j=j+1
            VV1(i)=V1(j)
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_getV2(n,VV2)

      USE PEGrid_mod

      Implicit None

         Integer :: n,i,j
         Real(8) :: VV2(nGrid2(1,n)*nGrid2(2,n))

         VV2=V2(ptV2(n)+1:ptV2(n)+nGrid2(1,n)*nGrid2(2,n))
         !j=ptV2(n)
         !Do i=1,nGrid2(1,n)*nGrid2(2,n)
         !   j=j+1
         !   VV2(i)=V2(j)
         !End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine PEGrid_getV3(n,VV3)

      USE PEGrid_mod

      Implicit None

         Integer :: n,i,j
         Real(8) :: VV3(nGrid3(1,n)*nGrid3(2,n)*nGrid3(3,n))

         VV3=V3(ptV3(n)+1:ptV3(n)+nGrid3(1,n)*nGrid3(2,n)*nGrid3(3,n))
         !j=ptV3(n)
         !Do i=1,nGrid3(1,n)*nGrid3(2,n)*nGrid3(3,n)
         !   j=j+1
         !   VV3(i)=V3(j)
         !End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function PEGrid_sizeofV1()

      USE PEGrid_mod

      Implicit None

         Integer :: i,j,k,l,n
         Real(8) :: ri,PEGrid_sizeofV1

         ri=0.D+00
         Do n=1,nS1
            ri=ri + dble(nGrid1(n))
         End do

         PEGrid_sizeofV1=ri*8.D+00

      End function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function PEGrid_sizeofV2()

      USE PEGrid_mod

      Implicit None

         Integer :: i,j,k,l,n
         Real(8) :: ri,PEGrid_sizeofV2

         ri=0.D+00
         Do n=1,nS2
            ri=ri + dble(nGrid2(1,n))*dble(nGrid2(2,n))
         End do

         PEGrid_sizeofV2=ri*8.D+00

      End function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function PEGrid_sizeofV3()

      USE PEGrid_mod

      Implicit None

         Integer :: i,j,k,l,n
         Real(8) :: ri,PEGrid_sizeofV3

         ri=0.D+00
         Do n=1,nS3
            ri=ri + dble(nGrid3(1,n))*dble(nGrid3(2,n))*dble(nGrid3(3,n))
         End do

         PEGrid_sizeofV3=ri*8.D+00

      End function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
