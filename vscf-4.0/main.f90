!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/02/04
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80


      PROGRAM main

      Implicit None

        Integer :: i

          Call Launch
          Call timer(1)

          Call Version !MK just startup till now

          i=0
          Call spr_Const(i)!/mol/Nat,Nfree,mass,x,omega,L and /sys/Maxmem
          if(i/=0) Goto 10

          Call Sindo_main !/vib/MR,vmax,vmax_base,CHODVR,vscf,vci,vcc,vpt,tdh,tdv,vav,summary

          Call spr_Dest
   10     Continue

          Call timer(2)

      End 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!     JOB 05:
!
      Subroutine Sindo_main

      Implicit None

         Integer :: ierr
         Logical :: lvscf,Vib_Getlvscf, &
                    lvci ,Vib_Getlvci,  &
                    lvpt ,Vib_Getlvpt,  &
                    lvav ,Vib_Getlvav


         Call Vib_Construct(ierr) !/vib/MR,vmax,vmax_base,CHODVR,vscf,vci,vcc,vpt,tdh,tdv,vav,summary
         if(ierr/=0) Goto 10
 
         lvscf= Vib_Getlvscf()
         lvci= Vib_Getlvci()
         lvpt= Vib_Getlvpt()
         lvav= Vib_Getlvav()

         Call Vib_setupGrid(ierr)
         if(ierr<0) return

         if(lvav) then
           Call vav_Construct()
         endif

         if(lvscf) then
           Call Vscf_Construct(ierr)
           if(ierr<0) return
           Call Vscf_main()
           Call Vscf_finalz()
         endif

         if(lvpt) then
           Call Vpt_Construct(ierr)
           if(ierr<0) return
           Call Vpt_main()
           Call Vpt_finalz()
         endif

         if(lvci) then
           Call Vci_Construct(ierr)
           if(ierr<0) return
           Call Vci_main()
           Call Vci_finalz()
         endif

         if(lvav) then
           Call vav_finalz
         endif

         Call Vib_finalz
   10    Continue

         return

      End subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
