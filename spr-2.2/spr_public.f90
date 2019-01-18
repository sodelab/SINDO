!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!   Last modified  2007/05/09
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!-----   Constructer reads the input file.
!
  Subroutine spr_Const(ierr)

    USE SPR_PRIVATE

    Implicit None
    Integer :: ierr,i
    Integer, parameter :: n=100,m=3*n-6
    Integer :: Maxmem
    Double Precision :: mass(n), x(3*n)
    Double Precision :: omega(m)
    Double Precision, dimension(3*n*m) :: L

    Namelist /mol/Nat,Nfree,mass,x,omega,L
    Namelist /sys/Maxmem


      !  >> Read &sys
      !  >> default memory: 1500MB.!MK
       Maxmem=1500
       Read(Inp,sys,end=10)
    10 Continue
       if(Maxmem>0) then
          rMaxmem=dble(Maxmem)*10**6
       else
          rMaxmem=-1
       endif
       rmem=0.D+00
       smem=0.D+00

      !  >> Read &mol
      !  >> initialize
       Nat=0; Nfree=0
       mass=0.D+00
       x=0.D+00
       omega=0.D+00
       L=0.D+00

       Rewind(Inp)
       Read(Inp,mol,end=20)
    20 Continue

       if(Nat>n) then 
          Write(Iout,*) 'ERROR IN SPR_CONST'
          Write(Iout,*) 'ERROR:: LIMITATION EXCEEDED'
          Write(Iout,*) 'ERROR:: NAT =',Nat
          Write(Iout,*) 'ERROR:: MAX NAT =',n
          ierr=-1
          return
       endif
       if(Nfree>m) then 
          Write(Iout,*) 'ERROR IN SPR_CONST'
          Write(Iout,*) 'ERROR:: LIMITATION EXCEEDED'
          Write(Iout,*) 'ERROR:: NFREE =',Nfree
          Write(Iout,*) 'ERROR:: MAX NFREE =',m
          ierr=-1
          return
       endif
       if(Nat==0.and.Nfree==0) then 
          Write(Iout,*) 'ERROR IN SPR_CONST'
          Write(Iout,*) 'ERROR::  NO INPUT GIVEN FOR EITHER NAT OR NFREE.'
          Write(Iout,*) 'ERROR::  >   NAT=',Nat
          Write(Iout,*) 'ERROR::  > NFREE=',Nfree
          ierr=-1
          return
       endif
       Nat3=Nat*3
       if(Nfree==0) then 
          Nfree=Nat3-6
       endif

       if(mass(1)/=0.D+00) then
          if(Nat==0) then 
             Write(Iout,*) 'ERROR IN SPR_CONST'
             Write(Iout,*) 'ERROR::  "MASS" REQUIRES NAT FOR INPUT.'
             Write(Iout,*) 'ERROR::  >   NAT=',Nat
             ierr=-1
          else
             Call spr_SetMass(mass)
          endif
       endif

       if(any(x/=0.D+00)) then
          if(Nat==0) then 
             Write(Iout,*) 'ERROR IN SPR_CONST'
             Write(Iout,*) 'ERROR::  "X" REQUIRES NAT FOR INPUT.'
             Write(Iout,*) 'ERROR::  >   NAT=',Nat
             ierr=-1
          else
             Call spr_Setxin(x)
          endif
       endif

       if(any(omega/=0.D+00)) then
          if(Nfree==0) then 
             Write(Iout,*) 'ERROR IN SPR_CONST'
             Write(Iout,*) 'ERROR::  "OMEGA" REQUIRES NFREE FOR INPUT.'
             Write(Iout,*) 'ERROR::  > NFREE=',Nfree
             ierr=-1
          else
             Call spr_SetFreq(omega)
          endif
       endif

       if(any(L/=0.D+00)) then
          if(Nat==0.or.Nfree==0) then 
             Write(Iout,*) 'ERROR IN SPR_CONST'
             Write(Iout,*) 'ERROR::  "L" REQUIRES NAT AND NFREE FOR INPUT.'
             Write(Iout,*) 'ERROR::  >   NAT=',Nat
             Write(Iout,*) 'ERROR::  > NFREE=',Nfree
             ierr=-1
          else
             Call spr_SetL(L)
          endif
       endif

       return

  End Subroutine spr_Const
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_Dest

    USE SPR_PRIVATE

       if(allocated(Zmass)) then
          Deallocate(Zmass)
          rmem=rmem-dble(Nat*8)
       endif
       if(allocated(xin)) then 
          Deallocate(xin)
          rmem=rmem-dble(Nat3*8)
       endif
       if(allocated(Freq)) then 
          Deallocate(Freq)
          rmem=rmem-dble(Nfree*8)
       endif
       if(allocated(CL)) then 
          Deallocate(CL)
          rmem=rmem-dble(Nat3*Nfree*8)
       endif

       Call spr_mem_final

  End Subroutine spr_Dest
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10

!-----   These routines set the parameters individually.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_SetNat(N)

    USE SPR_PRIVATE 

    Implicit None
    Integer  N

      Nat=N

    return

  End Subroutine spr_SetNat
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_SetNfree(N)

    USE SPR_PRIVATE

    Implicit None
    Integer  N

      Nfree=N

    return

  End Subroutine spr_SetNfree
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_SetMaxmem(N)

    USE SPR_PRIVATE

    Implicit None
    Integer  N

      rMaxmem=dble(N)*10**6

    return

  End Subroutine spr_SetMaxmem
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_SetMass(M)

    USE SPR_PRIVATE

    Implicit None
    Integer :: i
    Double Precision, dimension(*) :: M

      if(.not. allocated(Zmass)) then 
         allocate(Zmass(Nat))
         rmem=rmem+dble(Nat*8)
      endif
      Do i=1,Nat
         Zmass(i)=M(i)
      End do

    return

  End Subroutine spr_SetMass
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_Setxin(x)

    USE SPR_PRIVATE

    Implicit None
    Integer :: i,j,k
    Double Precision, dimension(*) :: x

      if(.not. allocated(xin)) then 
          allocate(xin(Nat3))
          rmem=rmem+dble(Nat3*8)
      endif
      k=1
      Do i=1,Nat
         Do j=1,3
            xin(k)=x(k)
            k=k+1
         End do
      End do

    return

  End Subroutine spr_Setxin
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_SetFreq(F)

    USE SPR_PRIVATE

    Implicit None
    Integer :: i
    Double Precision, dimension(*) :: F

      if(.not. allocated(Freq)) then 
         allocate(Freq(Nfree))
         rmem=rmem+dble(Nfree*8)
      endif
      Do i=1,Nfree
         Freq(i)=F(i)
      End do

    return

  End Subroutine spr_SetFreq
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_SetL(L)

    USE SPR_PRIVATE

    Implicit None
    Integer :: i,j,k
    Double Precision, dimension(*) :: L

      if(.not. allocated(CL)) then 
         allocate(CL(Nat3,Nfree))
         rmem=rmem+dble(Nat3*Nfree*8)
      endif
      k=1
      Do i=1,Nfree
         Do j=1,Nat3
            CL(j,i)=L(k)
            k=k+1
         End do
      End do

    return

  End Subroutine spr_SetL
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10

!-----   These routines provides the parameters.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Function spr_GetNat()

    USE SPR_PRIVATE 

    Implicit None
    Integer  spr_GetNat

      spr_GetNat=Nat

    return

  End Function spr_GetNat
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Function spr_GetNfree()

    USE SPR_PRIVATE 

    Implicit None
    Integer  spr_GetNfree

      spr_GetNfree=Nfree

    return

  End Function spr_GetNfree
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_GetMass(M)

    USE SPR_PRIVATE

    Implicit None
    Double Precision, dimension(Nat) :: M

      if(.not. allocated(Zmass)) then 
         Write(Iout,*) 'ERROR IN SPR_GETMASS:: MASS OF EACH ATOMS IS NOT GIVEN'
         Stop
      endif
      M = Zmass

    return

  End Subroutine spr_GetMass
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_Getxin(x)

    USE SPR_PRIVATE

    Implicit None
    Double Precision, dimension(Nat3) :: x

      if(.not. allocated(xin)) then 
         Write(Iout,*) 'ERROR IN SPR_GETXIN:: CARTESIAN COORDINATES ARE NOT GIVEN'
         Stop
      endif
      x=xin

    return

  End Subroutine spr_Getxin
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_GetFreq(F)

    USE SPR_PRIVATE

    Implicit None
    Double Precision, dimension(Nfree) :: F

      if(.not. allocated(Freq)) then 
         Write(Iout,*) 'ERROR IN SPR_GETFREQ:: FREQUENCIES ARE NOT GIVEN'
         Stop
      endif
      F=Freq

    return

  End Subroutine spr_GetFreq
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_GetL(L)

    USE SPR_PRIVATE

    Implicit None
    Double Precision, dimension(Nat3,Nfree) :: L

      if(.not. allocated(CL)) then 
         Write(Iout,*) 'ERROR IN SPR_GETL:: CL VECTORS ARE NOT GIVEN'
         Stop
      endif
      L = CL

    return

  End Subroutine spr_GetL
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_GetIO(In,Io)

    USE SPR_PRIVATE 

      In=Inp
      Io=Iout

    return

  End Subroutine spr_GetIO
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_SetIO(In,Io)

    USE SPR_PRIVATE 

       Inp=In
       Iout=Io

    return

  End Subroutine spr_SetIO
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! -- Print memory information
!
  Subroutine spr_meminfo

    USE SPR_PRIVATE

       Write(Iout,100) rmem/1.0E+06,rmaxmem/1.0D+06,rmem/rmaxmem*1.0E+02

    100 Format(7x,' -- MEMORY USAGE INFO [ ',f8.2,' MB / ',f8.2,' MB (',f5.1,' % ) ]  -- ',/)

    return

  End Subroutine spr_meminfo
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
  Subroutine spr_meminfo2(rm,rmm)

    USE SPR_PRIVATE

    Real(8) :: rm,rmm

       rm=rmem
       rmm=rMaxmem

    return

  End Subroutine spr_meminfo2
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!
  Subroutine spr_mem_final

    USE SPR_PRIVATE

       Write(Iout,100) smem/1.0E+06,rmaxmem/1.0D+06,smem/rmaxmem*1.0E+02

    100 Format(7x,' -- TOTAL MEMORY USAGE [ ',f8.2,' MB / ',f8.2,' MB (',f5.1,' % ) ]  -- ',/)

    return

  End Subroutine spr_mem_final
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! iopt:: =0, Silent mode
!        =1, Output (but don't stop)
!        =-1, Output and stop when error.
!  ri :: Request of the memory allocation in MB.
! spr_memalloc::  = 0, Success
!                 =-1, Fail (not enough space)
!
  Function spr_memalloc(iopt,ri)

    USE SPR_PRIVATE

    Integer :: iopt,spr_memalloc
    Real(8) :: ri

       if(rmem+ri < rMaxmem .or. rMaxmem<0) then
          rmem=rmem+ri
          spr_memalloc=0
          if(rmem>smem) smem=rmem

       else
          spr_memalloc=-1
          if(iopt/=0) then 
             Write(Iout,*) 'ERROR IN SPR_MEMALLOC'
             Write(Iout,*) 
             Call spr_meminfo
             Write(Iout,*) 'ERROR::  NOT ENOUGH MEMORY SPACE TO ALLOCATE'
             Write(Iout,'('' ERROR::'',f12.1)') ri*1.D-06
             Write(Iout,*) 'ERROR::  MEGA BYTE.'
             Write(Iout,*) 
             if(iopt==-1) Stop
          endif
           
       endif

    return

  End Function
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! ri :: Request of the memory deallocation in byte.
!
  Subroutine spr_memdealloc(ri)

    USE SPR_PRIVATE

    !Integer :: i
    Real(8) :: ri

      rmem=rmem-ri

    return

  End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
