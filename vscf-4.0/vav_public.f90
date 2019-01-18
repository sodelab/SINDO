!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/12/09
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vav_Construct()

      USE Vav_mod

      Implicit None

         Integer :: in,iout,i
         Logical :: op

         Namelist /vav/r_v,NP,MR2

          Call spr_Getio(in,iout)

          ! >> Read input

          ! --- default ---
          r_v=.true. 
          NP=0
          MR2=1

          Rewind(in)
          Read(in,vav,end=10)
       10 Continue

          write(iout,100)
      100 Format(//,'(  ENTER VAV MODULE  )',//)

          if(r_v) then
             write(iout,110)
             Call Vav_setup_nma()
          endif
      110 Format(/,3x,'>> VIBRATIONALLY AVERAGED STRUCTURE',/)

          if(NP/=0) then
             if(MR2>1) then
                write(iout,*)' WARNING: VAV MODULE IS CURRENTLY AVAILABLE ONLY'
                write(iout,*)' WARNING: FOR MR2=1.'
                MR2=1
             endif

             write(iout,120) MR2,NP
             write(iout,121) 
             Call setupPRPT(iout)
             write(iout,122)
             Do i=1,NP
                write(iout,123) i,P0(i)
             End do

          endif
      120 Format(/,3x,'>> VIBRATIONAL AVERAGE OVER PROPERTY SURFACES',//, &
                   3x,'    MR =',i4,/&
                   3x,'    NP =',i4)
      121 Format(/,3x,'    [ PROPERTY SURFACES ]',/)
      122 Format(/,3x,'    [ PROPERTY_VALUE@EQ ]')
      123 Format(  3x,'     o P[',i2.2,']= ',f15.4)


          write(iout,200)
          Call spr_meminfo
          Call timer(1)

  200 Format(//,'(  SETUP OF VAV COMPLETED  )',//)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vav_finalz()

      USE Vav_mod

      Implicit None

      Integer :: in,io

          Call spr_Getio(in,io)
          Write(io,100)
      100 Format(/,'(  FINALIZE VAV MODULE  )',/)

          if(NP>0) then
             if(allocated(P0)) then
                Call spr_memdealloc(dble(size(P0)*8))
                Deallocate(P0)
             endif
             if(allocated(P1)) then
                Call spr_memdealloc(dble(size(P1)*8))
                Deallocate(P1)
             endif
          endif
          if(r_v) Call Vav_finlz_nma

          Call spr_meminfo
          Call timer(1)

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vav_Getr_v()

      USE Vav_mod, ONLY : r_v

        Implicit None

           Logical Vav_Getr_v

           Vav_Getr_v=r_v

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Function Vav_GetNP()

      USE Vav_mod, ONLY : NP

        Implicit None

           Integer :: Vav_GetNP

           Vav_GetNP=NP

      End Function

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vav_setup_nma()

      Implicit None

         Integer :: in,io,ierr,spr_memalloc
         Integer :: Na,Nf,spr_GetNat,spr_GetNfree
         Real(8) :: rNa,rNf
         Real(8), allocatable :: x0(:),CL(:,:),Ms(:)

         Na=spr_GetNat()
         Nf=spr_GetNfree()

         rNa=dble(Na); rNf=dble(Nf)
         ierr=spr_memalloc(-1,(rNa*3.D+00+rNa*rNf+rNa)*8.D+00)
         Allocate(x0(Na*3),CL(Na*3,Nf),Ms(Na))

         Call spr_Getio(in,io)
         Call spr_Getxin(x0)
         Call spr_GetMass(Ms)
         Call spr_GetL(CL)
         Call nma_Construct(Na,Nf,Ms,x0,CL)
         Call nma_Print(io)

         Deallocate(x0,CL,Ms)

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine Vav_finlz_nma()

        Call nma_destruct

      End Subroutine 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
      Subroutine Vav_q2x(Na,Nf,x0,QQ)

        Implicit None

        Integer :: Na,Nf
        Real(8) :: x0(3,Na),QQ(Nf),CL(Na*3,Nf),Ms(Na)

        Call nma_q2x(x0,QQ)
        !Write(6,'(3f12.6)') x0

      End Subroutine Vav_q2x
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
