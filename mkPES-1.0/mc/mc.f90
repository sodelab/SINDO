!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Module mc_mod


        !     int               :: Max number of Hermitte-Gauss grid
        !     quad(int,Nfree)     :: Quadrature points / Angs(amu)1/2

        Integer :: Nfree
        Integer, parameter :: int=2
        Real(8), dimension(:,:), allocatable :: quad

        !     Potential energy
        Real(8) :: V0
        Real(8), allocatable :: V1(:,:),V2(:,:),V3(:,:)

        !     Irst,Irst2 :: File indicators
        Integer :: Irst,Irst2

      End Module
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      PROGRAM main

      Implicit None

          Call v2AI_Const()
          Call mc_Inp()
          Call mc_nMR_PES()
          Call mc_main()
          Call mc_finalz
          Call v2AI_Dest()

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine mc_main()

      USE mc_mod

      Implicit None


        Integer :: i,j,k,l,ij,ijk,nf2,nf3,i2(2),i3(3)
        Real(8) :: x,H2wvn,res2,res3

        Integer, allocatable :: map2(:,:),map3(:,:), &
                                is2(:),is3(:),a(:,:)
        Real(8), allocatable :: d2(:),s2(:),x2(:),Om(:), &
                                d3(:),s3(:),x3(:)

          write(6,100) 
      100 Format(//,' -> Strength of Mode Coupling terms: ',/)

          Allocate(Om(Nfree))
          Call Get_omega(Om)
          Om=Om/H2wvn()

          nf2=Nfree*(Nfree-1)/2
          nf3=Nfree*(Nfree-1)*(Nfree-2)/6

          Allocate(map2(2,nf2), is2(nf2), d2(nf2), s2(nf2), x2(nf2), &
                   map3(3,nf3), is3(nf3), d3(nf3), s3(nf3), x3(nf3))

          ij=1
          Do i=2,Nfree
          Do j=1,i-1
             map2(1,ij)=i 
             map2(2,ij)=j 

             x=0.D+00
             Do k=1,4
                !x=x + abs(V2(k,ij))
                x=x + V2(k,ij)
             End do
             !x=x/4D+00
             x=abs(x)/4D+00
             s2(ij)=x
             d2(ij)=res2(Om(i),Om(j))
             x2(ij)=s2(ij)*d2(ij)

             ij=ij+1
          End do
          End do

          Call sort(nf2,is2,x2)
          Write(6,'(''  > 2MR'')')
          Do i=1,nf2
             Write(6,'(2i4,4x,f12.6,f12.4,f12.6)') &
                          !map2(:,i),s2(i),d2(i),x2(i)
                          map2(:,is2(i)),s2(is2(i)),d2(is2(i)),x2(i)
          End do

          ijk=1
          Do i=3,Nfree
          Do j=2,i-1
          Do k=1,j-1
             map3(1,ijk)=i
             map3(2,ijk)=j
             map3(3,ijk)=k
             x=0.D+00
             Do l=1,8
                !x=x + abs(V3(l,ijk))
                x=x + V3(l,ijk)
             End do
             !x=x/8D+00
             x=abs(x)/8D+00
             s3(ijk)=x
             d3(ijk)=res3(Om(i),Om(j),Om(k))
             x3(ijk)=s3(ijk)*d3(ijk)

             ijk=ijk+1
          End do
          End do
          End do

          Call sort(nf3,is3,x3)
          Write(6,*)
          Write(6,'(''  > 3MR'')')
          Do i=1,nf3
             Write(6,'(3i4,f12.6,f12.4,f12.6)') &
                      !map3(:,i),s3(i),d3(i),x3(i)
                      map3(:,is3(i)),s3(is3(i)),d3(is3(i)),x3(i)
          End do

          Write(6,*)
          Write(6,*)

          !Call sort(nf2,is2,x2)
          !Call sort(nf3,is3,x3)
          i=1; j=1

        ! very Strong
          Write(6,'(''  > vS (>1)'')')
          k=1; l=1
          Do while(x2(i)>1.0D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,2i4,$)') map2(:,is2(i))
                l=2
             else
                Write(6,'(2x,2i4)') map2(:,is2(i))
                l=3
             endif
             i=i+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

          k=1; l=1
          Do while(x3(j)>1.0D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,3i4,$)') map3(:,is3(j))
                l=2
             else
                Write(6,'(2x,3i4)') map3(:,is3(j))
                l=3
             endif
             j=j+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

        ! Strong
          Write(6,'(''  > S (0.25 - 1)'')')
          k=1; l=1
          Do while(x2(i)>0.25D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,2i4,$)') map2(:,is2(i))
                l=2
             else
                Write(6,'(2x,2i4)') map2(:,is2(i))
                l=3
             endif
             i=i+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

          k=1; l=1
          Do while(x3(j)>0.25D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,3i4,$)') map3(:,is3(j))
                l=2
             else
                Write(6,'(2x,3i4)') map3(:,is3(j))
                l=3
             endif
             j=j+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

        ! Medium
          Write(6,'(''  > M (0.1 - 0.25)'')')
          k=1; l=1
          Do while(x2(i)>0.1D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,2i4,$)') map2(:,is2(i))
                l=2
             else
                Write(6,'(2x,2i4)') map2(:,is2(i))
                l=3
             endif
             i=i+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

          k=1; l=1
          Do while(x3(j)>0.1D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,3i4,$)') map3(:,is3(j))
                l=2
             else
                Write(6,'(2x,3i4)') map3(:,is3(j))
                l=3
             endif
             j=j+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

        ! Weak
          Write(6,'(''  > W (0.05 - 0.1)'')')
          k=1; l=1
          Do while(x2(i)>0.05D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,2i4,$)') map2(:,is2(i))
                l=2
             else
                Write(6,'(2x,2i4)') map2(:,is2(i))
                l=3
             endif
             i=i+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

          k=1; l=1
          Do while(x3(j)>0.05D+00) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,3i4,$)') map3(:,is3(j))
                l=2
             else
                Write(6,'(2x,3i4)') map3(:,is3(j))
                l=3
             endif
             j=j+1
             k=k+1
          End do
          if(l==2) Write(6,'(/)')
          if(l==3) Write(6,*)

        ! very Weak
          Write(6,'(''  > vW (<0.05)'')')
          k=1; l=1
          Do while(i<=nf2) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,2i4,$)') map2(:,is2(i))
                l=-1
             else
                Write(6,'(2x,2i4)') map2(:,is2(i))
                !l=1
             endif
             i=i+1
             k=k+1
          End do
          if(l<0) Write(6,'(/)')

          k=1
          Do while(j<=nf3) 
             if(mod(k,5)/=0) then
                Write(6,'(2x,3i4,$)') map3(:,is3(j))
             else
                Write(6,'(2x,3i4)') map3(:,is3(j))
             endif
             j=j+1
             k=k+1
          End do
          Write(6,*)

          Deallocate(map2,map3,x2,x3,is2,is3)

      Contains

      Subroutine sort(N,L,C)

      Implicit None

         Integer :: N,L(N)
         Real(8) :: C(N)

         Integer :: i,j(1),k,itmp
         Real(8) :: tmp

            Do i=1,N
               L(i)=i
            End do

            !dbg write(6,*)
            !dbg write(6,'(i3,f9.4)') (i,C(i),i=1,N)
            Do i=1,N
               j=MaxLoc(C(i:N))
               k=j(1)+i-1

               tmp=C(i)
               C(i)=C(k)
               C(k)=tmp

               itmp=L(i)
               L(i)=L(k)
               L(k)=itmp
            End do
            !dbg write(6,*)
            !dbg write(6,'(i3,f9.4)') (L(i),C(i),i=1,N)

      End subroutine
      !-----------------------------------------------------------------

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
      
      Function res2(wi,wj)

      Implicit None

         Real(8) :: res2,wi,wj,wii,wjj,del,H2wvn

            del=10.D+00/H2wvn()

            if(wi>wj) then 
               wii=wi; wjj=wj
            else
               wii=wj; wjj=wi
            endif
            res2=13.D+00/24.D+00/(abs(2.D+00*wj-wi)+del) + &
                  1.D+00/6.D+00 /(abs(wi-wj)+del)        + &
                  1.D+00/6.D+00 /(abs(3.D+00*wj-wi)+del) + &
                  1.D+00/24.D+00/(abs(4.D+00*wj-wi)+del)

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function res3(wi,wj,wk)

      Implicit None

         Real(8) :: res3,wi,wj,wk,wii,wjj,wkk,del,H2wvn

         del=10.D+00/H2wvn()

         if(wi>wj) then
            if(wk>wi) then
               wii=wk; wjj=wi; wkk=wj
            elseif(wk>wj) then
               wii=wi; wjj=wk; wkk=wj
            else
               wii=wi; wjj=wj; wkk=wk
            endif
         else
            if(wk>wj) then
               wii=wk; wjj=wj; wkk=wi
            elseif(wk>wi) then
               wii=wj; wjj=wk; wkk=wi
            else
               wii=wj; wjj=wi; wkk=wk
            endif
         endif
         res3=4.D+00/3.D+00/(abs(wii-wjj-wkk)+del) + &
              0.5D+00/(abs(wii-wjj-2.D+00*wkk)+del) + &
              0.5D+00/(abs(wii-2.D+00*wjj-wkk)+del) + &
              1.D+00/6.D+00/(abs(wii-wjj-3.D+00*wkk)+del) + &
              1.D+00/6.D+00/(abs(wii-3.D+00*wjj-wkk)+del) + &
              1.D+00/4.D+00/(abs(wii-2.D+00*wjj-2.D+00*wkk)+del)

      End

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mc_Inp()

      USE mc_mod

      Implicit None


        Integer :: i,j,k,l
        Integer :: Get_Nfree
        Real(8), dimension(:), allocatable :: Om
        Real(8) :: yy,H2wvn,Elmass,B2A
        Logical :: op

!----------------------------------------------------------------------

        Write(6,100)
  100   Format(/,'--------------------(    MC MODULE    )--------------------',/)

! --    Read parameters
        Nfree=Get_Nfree()
        yy=2.D+00

        Allocate(V1(2,Nfree),V2(4,Nfree*(Nfree-1)/2), &
                 V3(8,Nfree*(Nfree-1)*(Nfree-2)/6))

!--     Restart Files
        Irst=210
        Do while(.true.)
           Inquire(Irst,opened=op)
           if(op) then
             Irst=Irst+1
           else
             exit
           endif
        End do
        Open(Irst,file='restart-r',status='unknown')

        Irst2=220
        Do while(.true.)
           Inquire(Irst2,opened=op)
           if(op) then
             Irst2=Irst2+1
           else
             exit
           endif
        End do
        Open(Irst2,file='restart-w',status='unknown')

        write(6,150)
    150 Format(/,' -> Restart(read)  FILE :  [ restart-r ]',/, &
                 ' -> Restart(write) FILE :  [ restart-w ]',/)

!--     Quadrature Grid Points
        Allocate(quad(int,Nfree))
        Allocate(Om(Nfree))
        Call Get_omega(Om)
        Om=Om/H2wvn()
        Do i=1,Nfree
           quad(2,i)=sqrt(2.D+00/Om(i)/Elmass())*B2A()
           !quad(2,i)=1/sqrt(Om(i)*Elmass())*B2A()*1.225D+00
           quad(1,i)=-quad(2,i)
        End do
        write(6,200)
        Do i=1,Nfree
           write(6,210) i
           write(6,220) quad(:,i)
        End do
        Write(6,*)

    200 Format(2x,'> Grids')
    210 Format(6x,'Mode :',i3)
    220 Format(6x,6f9.4)

        End
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!

      Subroutine mc_finalz

      USE mc_mod

      Implicit None

      Integer :: in,io

          Call spr_Getio(in,io)
          Write(io,100)
      100 Format(/,'(  FINALIZE MC MODULE  )',/)

          Deallocate(quad,V1,V2,V3)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mc_nMR_PES()

      USE mc_mod

      Implicit None

      ! ------------------------------------------------------------------------

      Integer :: io,in

      ! ------------------------------------------------------------------------

          Call spr_GetIO(in,io)

          !----------------------
          ! >>  Q=0
          !----------------------
          write(io,100) 
      100 Format(' -> Energies on grids ',/, &
                2x,'> Q=0')
          Call mc_e_0(io)

          !----------------------
          ! >>  1MR PES 
          !----------------------
          Call mc_e_1MR(io)

          !----------------------
          ! >>  2MR PES 
          !----------------------
          Call mc_e_2MR(io)

          !----------------------
          ! >>  3MR PES 
          !----------------------
          Call mc_e_3MR(io)

          
          ! --------------------------------------------------------------------

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mc_e_0(io)

      USE mc_mod

      Implicit None

      Integer :: io
      Real(8), dimension(:), allocatable   :: qq

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree))
          qq=0.D+00
          Call mc_airun_e(0,qq,V0)
          Write(io,2000) 0,V0
          Deallocate(qq)

     2000 Format(i8,f18.10)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine mc_e_1MR(io)

      USE mc_mod

      Implicit None

      Integer :: io

      Integer :: i,im,idx,Mi,si(3)
      Real(8), dimension(:), allocatable :: qq,qu

      ! ------------------------------------------------------------------------

          Allocate(qq(Nfree))
          qq=0.D+00
          idx=0

          Do i=1,Nfree

             Write(io,1000) i

             Allocate(qu(int))
             qu=quad(:,i)

             Call SymOp(i,Mi,si)

             if(Mi==0) then
                qq(i)=qu(1)
                idx=idx+1
                Call mc_airun_e(idx,qq,V1(1,i))
                Write(io,2000) idx,V1(1,i),qq(i)

                qq(i)=qu(2)
                idx=idx+1
                Call mc_airun_e(idx,qq,V1(2,i))
                Write(io,2000) idx,V1(2,i),qq(i)

             else
                qq(i)=qu(1)
                idx=idx+1
                Call mc_airun_e(idx,qq,V1(1,i))
                Write(io,2000) idx,V1(1,i),qq(i)

                V1(2,i)=V1(1,i)
                Write(io,2000) 0,V1(2,i),qu(2)

             endif
             qq(i)=0.D+00

             V1(:,i)=V1(:,i)-V0

             Deallocate(qu)

          End do

          Deallocate(qq)

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,f8.3)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine mc_e_2MR(io)

      USE mc_mod

      Implicit None

      Integer :: io

      Integer :: i,j,k,idx,Mi,Mj,si(3),sj(3)
      Real(8), dimension(:), allocatable :: qui,quj

      ! ------------------------------------------------------------------------

          idx=0
          k=0
          Do i=2,Nfree
          Do j=1,i-1

             k=k+1
             Write(io,1000) i,j

             Call SymOp(i,Mi,si)
             Call SymOp(j,Mj,sj)
             Allocate(qui(int),quj(int))
             qui=quad(:,i)
             quj=quad(:,j)

             if(Mi==0) then

                if(Mj==0) then

                   Call PEF2D(0,0,Nfree,i,j,idx,qui,quj,V2(:,k))

                else

                   Call PEF2D(0,1,Nfree,i,j,idx,qui,quj,V2(:,k))

                   ! qj>0 region : generate by symmetry
                   V2(2,k)=V2(1,k)
                   Write(io,2000) 0,V2(2,k),qui(1),quj(2)
                   V2(4,k)=V2(3,k)
                   Write(io,2000) 0,V2(2,k),qui(2),quj(2)

                endif

             else
                if(Mj==0) then

                   Call PEF2D(1,0,Nfree,i,j,idx,qui,quj,V2(:,k))

                   ! qi>0 region : generate by symmetry
                   V2(3,k)=V2(1,k)
                   Write(io,2000) 0,V2(3,k),qui(2),quj(1)
                   V2(4,k)=V2(2,k)
                   Write(io,2000) 0,V2(4,k),qui(2),quj(2)

                elseif(Mi==Mj) then
                   Call PEF2D(1,0,Nfree,i,j,idx,qui,quj,V2(:,k))

                   ! qi>0 region : generate by symmetry
                   V2(3,k)=V2(2,k)
                   Write(io,2000) 0,V2(3,k),qui(2),quj(1)
                   V2(4,k)=V2(1,k)
                   Write(io,2000) 0,V2(4,k),qui(2),quj(2)

                else
                   Call PEF2D(1,1,Nfree,i,j,idx,qui,quj,V2(:,k))

                   V2(2,k)=V2(1,k)
                   Write(io,2000) 0,V2(2,k),qui(1),quj(2)
                   V2(3,k)=V2(1,k)
                   Write(io,2000) 0,V2(3,k),qui(2),quj(1)
                   V2(4,k)=V2(1,k)
                   Write(io,2000) 0,V2(4,k),qui(2),quj(2)

                endif
             endif

             V2(:,k)=V2(:,k)-V0
             V2(1,k)=V2(1,k)-V1(1,i)-V1(1,j)
             V2(2,k)=V2(2,k)-V1(1,i)-V1(2,j)
             V2(3,k)=V2(3,k)-V1(2,i)-V1(1,j)
             V2(4,k)=V2(4,k)-V1(2,i)-V1(2,j)

             Deallocate(qui,quj)

          End do
          End do

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,2f8.3)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     ir>0  ... qi<0
!     ir<0  ... qi>0
!     ir=0  ... disabled

!     jr>0  ... qj<0
!     jr<0  ... qj>0
!     jr=0  ... disabled 

      Subroutine PEF2D(ir,jr,Nfree,ii,jj,idx,qqi,qqj,V2)

      Implicit None

      Integer :: ir,jr,Nfree,ii,jj,idx
      Real(8) :: qqi(2),qqj(2),V2(4)

      Real(8) :: qq(Nfree)

         qq=0.D+00

         ! qi<0
         if(ir>=0) then

         ! qi<0, qj<0 region
         if(jr>=0) then
            qq(ii)=qqi(1)
            qq(jj)=qqj(1)
            idx=idx+1
            Call mc_airun_e(idx,qq,V2(1))
            Write(6,2000) idx,V2(1),qq(ii),qq(jj)
         endif

         ! qi<0, qj>0 region
         if(jr<=0) then
            qq(ii)=qqi(1)
            qq(jj)=qqj(2)
            idx=idx+1
            Call mc_airun_e(idx,qq,V2(2))
            Write(6,2000) idx,V2(2),qq(ii),qq(jj)
         endif
         endif

         if(ir<=0) then

         ! qi>0, qj<0 region
         if(jr>=0) then
            qq(ii)=qqi(2)
            qq(jj)=qqj(1)
            idx=idx+1
            Call mc_airun_e(idx,qq,V2(3))
            Write(6,2000) idx,V2(3),qq(ii),qq(jj)
         endif

         ! qi>0, qj>0 region
         if(jr<=0) then
            qq(ii)=qqi(2)
            qq(jj)=qqj(2)
            idx=idx+1
            Call mc_airun_e(idx,qq,V2(4))
            Write(6,2000) idx,V2(4),qq(ii),qq(jj)
         endif

         endif

    2000 Format(i8,f18.10,100f8.3)

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

      Subroutine mc_e_3MR(io)

      USE mc_mod

      Implicit None

      Integer :: io

      Integer :: i,j,k,l,ij,jk,ik,ijk,idx,jdx, &
                 Mi,Mj,Mk,si(3),sj(3),sk(3),getChr
      Real(8), dimension(:), allocatable :: qui,quj,quk
      Character :: Sym*3,getSym*3

      ! ------------------------------------------------------------------------

          Sym=getSym()

          ! --------------------------------------------------------------------
          ! >>>>>>>

          idx=0
          l=0
          Do i=3,Nfree
          Do j=2,i-1
          Do k=1,j-1

             l=l+1
             Write(io,1000) i,j,k

             Call SymOp(i,Mi,si)
             Call SymOp(j,Mj,sj)
             Call SymOp(k,Mk,sk)
             Allocate(qui(2),quj(2),quk(2))
             qui=quad(:,i)
             quj=quad(:,j)
             quk=quad(:,k)

             if(Mi==0) then
                if(Mj==0) then
                   if(Mk==0) then
                      ! OK
                      Call PEF3D(0,0,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                   else
                      ! OK
                      ! qk<0 region
                      Call PEF3D(0,0,1,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                      ! qk>0 region : generate by symmetry
                      V3(2,l)=V3(1,l)
                      Write(io,2000) 0,V3(2,l),qui(1),quj(1),quk(2)
                      V3(4,l)=V3(3,l)
                      Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                      V3(6,l)=V3(5,l)
                      Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                      V3(8,l)=V3(7,l)
                      Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)
                
                   endif

                else
                   if(Mk==0) then
                      ! OK
                      ! qj<0 region
                      Call PEF3D(0,1,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                      ! qj>0 region : generate by symmetry
                      V3(3,l)=V3(1,l)
                      Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                      V3(4,l)=V3(2,l)
                      Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                      V3(7,l)=V3(5,l)
                      Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                      V3(8,l)=V3(6,l)
                      Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                   elseif(Mk==Mj) then
                      ! OK
                      ! qj<0 region
                      Call PEF3D(0,1,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                      ! qj>0 region : generate by symmetry
                      V3(3,l)=V3(2,l)
                      Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                      V3(4,l)=V3(1,l)
                      Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                      V3(7,l)=V3(6,l)
                      Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                      V3(8,l)=V3(5,l)
                      Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                   else !(Mk/=Mj)
                      ! OK
                      ! qj<0,qk<0 region
                      Call PEF3D(0,1,1,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))
 
                      ! qj<0,qk>0 region : generate by symmetry
                      V3(2,l)=V3(1,l)
                      Write(io,2000) 0,V3(2,l),qui(1),quj(1),quk(2)
                      V3(3,l)=V3(1,l)
                      Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                      V3(4,l)=V3(1,l)
                      Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                      V3(6,l)=V3(5,l)
                      Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                      V3(7,l)=V3(5,l)
                      Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                      V3(8,l)=V3(5,l)
                      Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                   endif

                endif

             else

                if(Mj==0) then
                   if(Mk==0 .or. Mk==Mi) then
                      ! qi<0 region
                      Call PEF3D(1,0,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                      if(Mk==0) then
                         ! OK
                         ! qi>0 region
                         V3(5,l)=V3(1,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(6,l)=V3(2,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(7,l)=V3(3,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(8,l)=V3(4,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                      else
                         ! Mk==Mi
                         ! OK
                         ! qi>0 region
                         V3(5,l)=V3(2,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(6,l)=V3(1,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(7,l)=V3(4,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(8,l)=V3(3,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                      endif

                   else
                      ! OK
                      ! Mk /=0 .and. Mi /= Mk
                      ! qi<0,qk<0 region
                      Call PEF3D(1,0,1,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                      ! qi<0,qk>0 region
                      V3(2,l)=V3(1,l)
                      Write(io,2000) 0,V3(2,l),qui(1),quj(1),quk(2)
                      V3(5,l)=V3(1,l)
                      Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                      V3(6,l)=V3(1,l)
                      Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                      V3(4,l)=V3(3,l)
                      Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                      V3(7,l)=V3(3,l)
                      Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                      V3(8,l)=V3(3,l)
                      Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                   endif
                else
                   if(Mk==0) then
                      if(Mi==Mj) then
                         ! OK
                         ! Mi=Mj /= 0, Mk=0
                         ! qi<0 region
                         Call PEF3D(1,0,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                         ! qi>0 region
                         V3(5,l)=V3(3,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(6,l)=V3(4,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(7,l)=V3(1,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(8,l)=V3(2,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                      else
                         ! OK
                         ! Mi/=Mj/= 0, Mk=0
                         ! qi<0,qj<0 region
                         Call PEF3D(1,1,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                         V3(3,l)=V3(1,l)
                         Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                         V3(5,l)=V3(1,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(7,l)=V3(1,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(4,l)=V3(2,l)
                         Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                         V3(6,l)=V3(2,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(8,l)=V3(2,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                      endif
                   else
                      if(Mi==Mj .and. Mj/=Mk) then
                         ! OK
                         ! qi<0,qk<0 region
                         Call PEF3D(1,0,1,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                         V3(2,l)=V3(1,l)
                         Write(io,2000) 0,V3(2,l),qui(1),quj(1),quk(2)
                         V3(7,l)=V3(1,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(8,l)=V3(1,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)
                         V3(4,l)=V3(3,l)
                         Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                         V3(5,l)=V3(3,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(6,l)=V3(3,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)

                      elseif(Mi==Mk .and. Mk/=Mj) then
                         ! OK
                         ! qj<0,qk<0 region
                         Call PEF3D(0,1,1,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))
                         V3(3,l)=V3(1,l)
                         Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                         V3(6,l)=V3(1,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(8,l)=V3(1,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)
                         V3(2,l)=V3(5,l)
                         Write(io,2000) 0,V3(2,l),qui(1),quj(1),quk(2)
                         V3(7,l)=V3(5,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(4,l)=V3(5,l)
                         Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)

                      elseif(Mi/=Mj .and. Mj==Mk) then
                         ! OK
                         ! qi<0,qj<0 region
                         Call PEF3D(1,1,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))
                         V3(3,l)=V3(2,l)
                         Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                         V3(6,l)=V3(2,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(7,l)=V3(2,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(4,l)=V3(1,l)
                         Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                         V3(5,l)=V3(1,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(8,l)=V3(1,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                      elseif(Mi==Mj .and. Mj==Mk) then
                         ! OK
                         ! qk<0 region
                         Call PEF3D(0,0,1,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                         V3(2,l)=V3(7,l)
                         Write(io,2000) 0,V3(2,l),qui(1),quj(1),quk(2)
                         V3(4,l)=V3(5,l)
                         Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                         V3(6,l)=V3(3,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(8,l)=V3(1,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                      else !(Mi/=Mj/=Mk)
                         !<<<<
                         if(Sym=='D2H' .and. getChr(Mi,Mj,Mk)==1) then

                         ! OK
                         ! qi<0,qj<0,qk<0 region
                         Call PEF3D(1,1,1,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                         V3(2,l)=V3(1,l)
                         Write(io,2000) 0,V3(2,l),qui(1),quj(1),quk(2)
                         V3(3,l)=V3(1,l)
                         Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                         V3(4,l)=V3(2,l)
                         Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                         V3(5,l)=V3(1,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(6,l)=V3(2,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(7,l)=V3(1,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)
                         V3(8,l)=V3(2,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)

                         !<<<<
                         else

                         ! OK
                         ! qi<0,qj<0 region
                         Call PEF3D(1,1,0,Nfree,i,j,k,idx,qui,quj,quk,V3(:,l))

                         V3(3,l)=V3(2,l)
                         Write(io,2000) 0,V3(3,l),qui(1),quj(2),quk(1)
                         V3(5,l)=V3(2,l)
                         Write(io,2000) 0,V3(5,l),qui(2),quj(1),quk(1)
                         V3(8,l)=V3(2,l)
                         Write(io,2000) 0,V3(8,l),qui(2),quj(2),quk(2)
                         V3(4,l)=V3(1,l)
                         Write(io,2000) 0,V3(4,l),qui(1),quj(2),quk(2)
                         V3(6,l)=V3(1,l)
                         Write(io,2000) 0,V3(6,l),qui(2),quj(1),quk(2)
                         V3(7,l)=V3(1,l)
                         Write(io,2000) 0,V3(7,l),qui(2),quj(2),quk(1)

                         endif

                      endif
                   endif
                endif

             endif

             ij= i*(i-1)/2 -(i-1)+j
             jk= j*(j-1)/2 -(j-1)+k
             ik= i*(i-1)/2 -(i-1)+k

             V3(:,l)=V3(:,l)-V0
             V3(1,l)=V3(1,l)-V2(1,ij)-V2(1,jk)-V2(1,ik) &
                            -V1(1,i)-V1(1,j)-V1(1,k)
             V3(2,l)=V3(2,l)-V2(1,ij)-V2(2,jk)-V2(2,ik) &
                            -V1(1,i)-V1(1,j)-V1(2,k)
             V3(3,l)=V3(3,l)-V2(2,ij)-V2(3,jk)-V2(1,ik) &
                            -V1(1,i)-V1(2,j)-V1(1,k)
             V3(4,l)=V3(4,l)-V2(2,ij)-V2(4,jk)-V2(2,ik) &
                            -V1(1,i)-V1(2,j)-V1(2,k)
             V3(5,l)=V3(5,l)-V2(3,ij)-V2(1,jk)-V2(3,ik) &
                            -V1(2,i)-V1(1,j)-V1(1,k)
             V3(6,l)=V3(6,l)-V2(3,ij)-V2(2,jk)-V2(4,ik) &
                            -V1(2,i)-V1(1,j)-V1(2,k)
             V3(7,l)=V3(7,l)-V2(4,ij)-V2(3,jk)-V2(3,ik) &
                            -V1(2,i)-V1(2,j)-V1(1,k)
             V3(8,l)=V3(8,l)-V2(4,ij)-V2(4,jk)-V2(4,ik) &
                            -V1(2,i)-V1(2,j)-V1(2,k)

             Deallocate(qui,quj,quk)

          End do
          End do
          End do

     1000 Format(2x,'> MODE=',10i4)
     2000 Format(i8,f18.10,3f8.3)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     ir,jr,kr>0  ... q<0
!     ir,jr,kr<0  ... q>0
!     ir,jr,kr=0  ... disabled

      Subroutine PEF3D(ir,jr,kr,Nfree,ii,jj,kk,idx,qqi,qqj,qqk,V3)

      Implicit None

      Integer :: ir,jr,kr,Nfree,ii,jj,kk,idx,ni,nj,nk
      Real(8) :: qqi(2),qqj(2),qqk(2),V3(8),qq(Nfree)

         qq=0.D+00

         ! qi<0
         if(ir>=0) then
            qq(ii)=qqi(1)

         ! qi<0, qj<0 region
         if(jr>=0) then
            qq(jj)=qqj(1)

            ! qi<0, qj<0, qk<0 region
            if(kr>=0) then
               qq(kk)=qqk(1)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(1))
               Write(6,2000) idx,V3(1),qq(ii),qq(jj),qq(kk)
            endif

            ! qi<0, qj<0, qk>0 region
            if(kr<=0) then
               qq(kk)=qqk(2)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(2))
               Write(6,2000) idx,V3(2),qq(ii),qq(jj),qq(kk)
            endif

         endif

         ! qi<0, qj>0 region
         if(jr<=0) then
            qq(jj)=qqj(2)

            ! qi<0, qj>0, qk<0 region
            if(kr>=0) then
               qq(kk)=qqk(1)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(3))
               Write(6,2000) idx,V3(3),qq(ii),qq(jj),qq(kk)
            endif

            ! qi<0, qj>0, qk>0 region
            if(kr<=0) then
               qq(kk)=qqk(2)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(4))
               Write(6,2000) idx,V3(4),qq(ii),qq(jj),qq(kk)
            endif

         endif

         endif

         ! qi>0
         if(ir<=0) then
            qq(ii)=qqi(2)

         ! qi>0, qj<0 region
         if(jr>=0) then
            qq(jj)=qqj(1)

            ! qi>0, qj<0, qk<0 region
            if(kr>=0) then
               qq(kk)=qqk(1)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(5))
               Write(6,2000) idx,V3(5),qq(ii),qq(jj),qq(kk)

            endif

            ! qi>0, qj<0, qk>0 region
            if(kr<=0) then
               qq(kk)=qqk(2)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(6))
               Write(6,2000) idx,V3(6),qq(ii),qq(jj),qq(kk)
            endif

         endif

         ! qi>0, qj>0 region
         if(jr<=0) then
            qq(jj)=qqj(2)

            ! qi>0, qj>0, qk<0 region
            if(kr>=0) then
               qq(kk)=qqk(1)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(7))
               Write(6,2000) idx,V3(7),qq(ii),qq(jj),qq(kk)
            endif

            ! qi>0, qj>0, qk>0 region
            if(kr<=0) then
               qq(kk)=qqk(2)
               idx=idx+1
               Call mc_airun_e(idx,qq,V3(8))
               Write(6,2000) idx,V3(8),qq(ii),qq(jj),qq(kk)
            endif

         endif

         endif

      2000 Format(i8,f18.10,3f8.3)

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine mc_airun_e(indx0,qq,E)

      USE mc_mod

      Implicit None

      Integer :: indx0,indx1
      Real(8) :: qq(Nfree),E
      Integer :: i

         Read(Irst,*,end=100,err=100) indx1,E
         if(indx0 .ne. indx1) then 
            Write(6,*) 'Invalid line is detected in the restart file:'
            !Write(6,'(i8,f20.10)') indx1,E
            Call write_e(6,indx1,E)
            Write(6,'('' It must have an index '',i8)') indx0
            Stop
         endif
         Call write_e(Irst2,indx0,E)
         Call myFlsh(Irst2)
         return

  100    Continue
         Call Run_e(qq,E)
         Call write_e(Irst2,indx0,E)
         Call myFlsh(Irst2)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      Subroutine write_e(ifl,indx,E)

      Implicit None

          Integer :: ifl,indx
          Real(8) :: E

          Write(ifl,500) indx,E
      500 Format(i8,f20.10)

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
