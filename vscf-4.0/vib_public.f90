!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/07/19
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine Vib_Construct(ierr)

      USE Vib_mod

      Implicit None

         Integer, parameter :: nmax=1000

         Integer :: i,j,k,l,n
         Integer :: ierr,in,iout
         Integer :: spr_GetNfree,spr_memalloc

         Integer :: vmax_base
         Integer, dimension(nmax) :: vmax
         Logical :: vscf,vci,vcc,vpt,tdh,tdv,vav

         Character :: title*40

         Integer :: ni,nj,nk

         Namelist /vib/MR,vmax,vmax_base,CHODVR, &
                       vscf,vci,vcc,vpt,tdh,tdv,vav,summary

         ierr=0
         Call spr_Getio(in,iout)
         Nfree=spr_GetNfree()
         if(Nfree>nmax) then
            ierr=-1
            Write(iout,'(''  ERROR: MAXIMUM NUMBER OF MODE IS '',i4)') nmax
            Write(iout,'(''  ERROR: TERMINATED IN VIB_CONSTRUCT '')')
            return
         endif

         ierr=spr_memalloc(1,dble(Nfree*4)+dble(Nfree*8))
         if(ierr<0) return
         allocate(nCHO(Nfree),omegaf(Nfree))

         ! Frequencies for harmonic oscillator basis functions.
         Call spr_GetFreq(omegaf)

         ! >> Read input

         ! --- DEFAULT ---
         ! - Run option 
         summary=.true. !MK
         vscf=.false. 
         vci=.false.
         vcc=.false.
         vpt=.false.
         tdh=.false.
         tdv=.false.
         vav=.false.

         ! - Mode Coupling 
         MR = 3

         ! - Maximum quanta of basis functions.
         vmax=-2
         vmax_base=10
         CHODVR=.false.

         Rewind(in)
         Read(in,vib,end=10)
      10 Continue

         ! MR=1,2,3
         if(MR>3) then
            Write(iout,*) 'ERROR: '
            Write(iout,*) 'ERROR: SORRY! MR>3 IS NOT READY.'
            Write(iout,*) 'ERROR: MR =',MR
            Write(iout,*) 'ERROR: '
            ierr=-1
            return
         endif

         !  Run option
         lvscf=vscf
         lvci=vci
         lvcc=vcc
         lvpt=vpt
         ltdh=tdh
         ltdv=tdv
         lvav=vav
         
         write(iout,100)
         write(iout,110) vscf,vci,vcc,vpt,tdh,tdv,vav
     100 Format(//,'(  ENTER VIB MODULE  )',//)
     110 Format(/,3x,'>> RUN OPTIONS',/ &
         9x,'VSCF =',l6,/, &
         9x,' VCI =',l6,/, &
         9x,' VCC =',l6,/, &
         9x,' VPT =',l6,/, &
         9x,' TDH =',l6,/, &
         9x,' TDV =',l6,/, &
         9x,' VAV =',l6,/)

         !  Num. of basis function
         nCHO=vmax_base
         Do i=1,Nfree
            if(vmax(i)==0) then
               nCHO(i)=-1
               cycle
            endif
            if(vmax(i)/=-2) nCHO(i)=vmax(i)
         End do
         nCHO=nCHO+1
         maxCHO=maxVal(nCHO)

         write(iout,120)
         i=mod(Nfree,6)
         j=(Nfree-i)/6
         Do k=1,j
            write(iout,125) (6*(k-1)+l,l=1,6)
            write(iout,126) nCHO(6*(k-1)+1:6*k)-1
            write(iout,128) omegaf(6*(k-1)+1:6*k)
         End do
         if(i/=0) then
            write(iout,125) (6*j+k,k=1,i)
            write(iout,126) nCHO(6*j+1:Nfree)-1
            write(iout,128) omegaf(6*j+1:Nfree)
         endif

         if(CHODVR) then
         !  CHODVR 
            ierr=spr_memalloc(1,dble(Nfree*4))
            if(ierr<0) return
            allocate(nHO(Nfree))

         endif
         write(iout,'(/)')

     120 Format(/,3x,'>> BASIS FUNCTIONS')
     125 Format(7x,'MODE : ',6i9)
     126 Format(7x,'MAXV : ',6i9)
     128 Format(7x,'FREQ : ',6f9.2,/)


         Call setupPEF(ierr)
         if(ierr<0) return

         write(iout,130) MR
     130 Format(/,3x,'>> POTENTIAL',/, &
                /,6x,'[  OPTIONS  ]',//, &
              9x,'MR     = ',i4)

         if(Len_PotDir /= 0) then
            write(iout,140) PotDir
         else
            write(iout,141) 
         endif
     140 Format(9x,'POTDIR = ',a80,/)
     141 Format(9x,'POTDIR =  PWD',/)

         Do i=1,Nfree
            if(nCHO(i)>0) cycle
            write(iout,145) i
         End do
     145 Format(9x,'DISABLED: MODE=',i4)

         write(iout,'(/,9x,''1MR-PEF'',/)')
         if(nS1/=0) then
            write(iout,150) 
            Do i=1,nS1
               Call Spfn1_Gen1(0,i)
               Call Spfn1_GetnGrid(j)
               Call Spfn1_GetTitle(title)
               Call Spfn1_Dispose()
               write(iout,151) mS1(i),j,title
            End do
            write(iout,*)

         endif

         if(nQ1/=0) then
            write(iout,160) tQ1 
            Do i=1,nQ1
               write(iout,161) mQ1(i)
            End do
            write(iout,*)

         endif
         if(MR==1) goto 1000

         write(iout,'(/,9x,''2MR-PEF'',/)')
         if(nS2/=0) then
            write(iout,150) 
            Do i=1,nS2
               Call Spfn2_Gen1(0,i)
               Call Spfn2_GetnGrid(j,k)
               Call Spfn2_GetTitle(title)
               Call Spfn2_Dispose()
               write(iout,152) mS2(:,i),j,k,title
            End do
            write(iout,*)
         endif

         if(nQ2/=0) then
            write(iout,160) tQ2 
            Do i=1,nQ2
               write(iout,162) mQ2(:,i)
            End do
            write(iout,*)

         endif
         if(MR==2) goto 1000

         write(iout,'(/,9x,''3MR-PEF'',/)')
         if(nS3/=0) then
            write(iout,150) 
            Do i=1,nS3
               Call Spfn3_Gen1(0,i)
               Call Spfn3_GetnGrid(j,k,l)
               Call Spfn3_GetTitle(title)
               Call Spfn3_Dispose()
               write(iout,153) mS3(:,i),j,k,l,title
            End do
            write(iout,*)
         endif

         if(nQ3/=0) then
            write(iout,160) tQ3
            Do i=1,nQ3
               write(iout,163) mQ3(:,i)
            End do
            write(iout,*)

         endif
         if(MR==3) goto 1000

     150 Format(9x,'o SPLINE FNC.',/)
     151 Format(11x,'MODE=',i4,', GRID=',i4,3x,a40)
     152 Format(11x,'MODE=',2i4,', GRID=',2i4,3x,a40)
     153 Format(11x,'MODE=',3i4,', GRID=',3i4,3x,a40)
     160 Format(9x,'o QFF: ',a,/)
     161 Format(11x,'MODE=',i4)
     162 Format(11x,'MODE=',2i4)
     163 Format(11x,'MODE=',3i4)

    1000 Continue

         write(iout,500)
         Call spr_meminfo
         Call timer(1)

     500 Format(/,'(  SETUP OF VIB COMPLETED  )',//)
    ! if (summary) write(666,666)vmax_base,MR ! MK
     666 Format('vmax_base=',i2,' MR=',i1)
      End subroutine

!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Vib_finalz()

      USE Vib_mod

      Implicit None 

         Integer :: in,io
         Real(8) :: memsz

         Call Modal_Destruct()
         Call PEGrid_Destruct()

         memsz=dble(size(omegaf)*8)+dble(size(nCHO)*4)
         Call spr_memdealloc(memsz)
         Deallocate(omegaf,nCHO)

         if(nS1/=0) then
            memsz=dble(size(mS1)*4)
            Call spr_memdealloc(memsz)
            Deallocate(mS1)
         endif

         if(nQ1/=0) then
            memsz=dble(size(mQ1)*4)
            Call spr_memdealloc(memsz)
            Deallocate(mQ1)
         endif

         Call spr_Getio(in,io)
         Write(io,100)
     100 Format(/,'(  FINALIZE VIB MODULE  )',/)

      End subroutine 
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine Vib_setupGrid(ierr)

      USE Vib_mod

      Implicit None

         Integer :: ierr,in,iout
         Integer :: i,j,k,l,n,ni,nj,nk

         Real(8) :: ri,PEGrid_sizeofV1,PEGrid_sizeofV2,PEGrid_sizeofV3

         Call spr_Getio(in,iout)

         write(iout,200) CHODVR 
     200 Format(//,'(  SETUP DVR BASIS  )',//, &
                /,6x,'[  OPTIONS  ]',//, &
              9x,'CHODVR     = ',l4)

         Call genGrid(ierr)
         if(ierr<0) return

         if(CHODVR) then
            Do i=1,Nfree
               write(iout,205) i,nHO(i)
            End do
         endif
     205 Format(9x,'> MODE  = ',i4,',  NHO = ',i4)

         write(iout,230)
         Do i=1,Nfree
            if(nCHO(i)/=0) then
               write(iout,210) i,nCHO(i)
               Call printDVR(i,nCHO(i))
            endif
         End do

         if(nS2>0) then
            write(iout,240) 
            Do n=1,nS2
               Call PEGrid_getnGrid2(ni,nj,n)
               write(iout,211) mS2(:,n),ni,nj
               !Call printDVR(mS2(1,n),ni)
               !Call printDVR(mS2(2,n),nj)
            End do
         endif

         if(nS3>0) then
            write(iout,250) 
            Do n=1,nS3
               Call PEGrid_getnGrid3(ni,nj,nk,n)
               write(iout,212) mS3(:,n),ni,nj,nk
               !Call printDVR(mS3(1,n),ni)
               !Call printDVR(mS3(2,n),nj)
               !Call printDVR(mS3(3,n),nk)
            End do
         endif

     210 Format(9x,'> MODE  = ',i4,',  NGRID = ',i4,/)
     211 Format(9x,'> MODE  = ',2i4,',  NGRID = ',2i4,/)
     212 Format(9x,'> MODE  = ',3i4,',  NGRID = ',3i4,/)
     230 Format(/,9x,'1MR TERMS',/)
     240 Format(/,9x,'2MR TERMS',/)
     250 Format(/,9x,'3MR TERMS',/)

         write(iout,300)
         ri=PEGrid_sizeofV1()*1.D-03
         write(iout,310) ri
         if(MR>1) then
            ri=PEGrid_sizeofV2()*1.D-06
            write(iout,320) ri
         endif
         if(MR>2) then
            ri=PEGrid_sizeofV3()*1.D-06
            write(iout,330) ri
         endif
         write(iout,*)

     300 Format(/,9x,'POTENTIAL ENERGY VALUES GENERATED INCORE',/)
     310 Format(9x,'> 1MR :  ',f8.2,' KB')
     320 Format(9x,'> 2MR :  ',f8.2,' MB')
     330 Format(9x,'> 3MR :  ',f8.2,' MB')

         write(iout,500)
         Call spr_meminfo
         Call timer(1)

     500 Format(/,'(  SETUP OF DVR COMPLETED  )',//)

         CONTAINS

         Subroutine printDVR(i,ni)

             Integer :: i,j,k,l,n,ni
             Real(8), dimension(:), allocatable   :: qq
             Real(8), dimension(:,:), allocatable :: xdvr

             Allocate(qq(ni),xdvr(ni,ni))
             Call Modal_getxdvr(i,ni,xdvr)
             Call Modal_getQ(i,ni,qq)
             j=mod(ni,9)
             if(j/=0) then
                j=(ni-j)/9
             else
                j=ni/9-1
             endif
             Do k=1,j
                write(iout,100) (l,l=(k-1)*10,k*10-1)
                Do l=1,ni
                   write(iout,200) qq(l),(xdvr(n,l),n=(k-1)*10+1,k*10)
                End do
                write(iout,*)
             End do
             write(iout,100) (l,l=j*10,ni-1)
             Do k=1,ni
                write(iout,200) qq(k),(xdvr(l,k),l=j*10+1,ni)
             End do
             write(iout,*)
             Deallocate(qq,xdvr)

         100 Format(11x,'   Q',10i8)
         200 Format(9x,11f8.4)
	 
         End subroutine


      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Function Vib_Getlvscf()

      USE Vib_mod

      Implicit None

           Logical Vib_Getlvscf

           Vib_Getlvscf=lvscf

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vib_Getlvci()

      USE Vib_mod

      Implicit None

           Logical Vib_Getlvci

           Vib_Getlvci=lvci

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vib_Getlvcc()

      USE Vib_mod

      Implicit None

           Logical Vib_Getlvcc

           Vib_Getlvcc=lvcc

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vib_Getlvpt()

      USE Vib_mod

      Implicit None

           Logical Vib_Getlvpt

           Vib_Getlvpt=lvpt

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vib_Getltdh()

      USE Vib_mod

      Implicit None

           Logical Vib_Getltdh

           Vib_Getltdh=ltdh

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vib_Getltdv()

      USE Vib_mod

      Implicit None

           Logical Vib_Getltdv

           Vib_Getltdv=ltdv

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vib_Getlvav()

      USE Vib_mod

      Implicit None

           Logical Vib_Getlvav

           Vib_Getlvav=lvav

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
