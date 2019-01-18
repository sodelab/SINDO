!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/07/19
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
  SUBROUTINE qff_Construct(ier)
!
  USE qff_mod
  USE Vib_mod, ONLY : PotDir,Len_PotDir

  Implicit None
!
    Integer :: mm,ier,in,spr_GetNfree,spr_memalloc
!
!
      Nfree=spr_GetNfree()
      Call spr_Getio(in,Iout)

      !quiet Write(Iout,100) Nfree
!
!=======================(  Set up reference points  ) =======================
!
! Set up Matrix
!
      ier=spr_memalloc(1,dble(Nfree*8))
         if(ier<0) goto 10
      Allocate (EC(Nfree))

      ! -- 1MR --
      mm=Nfree
      ier=spr_memalloc(1,dble(mm*4*8))
         if(ier<0) goto 10
      Allocate (Gi(mm),Hii(mm),Tiii(mm),Uiiii(mm))

      ! -- 2MR --
      mm=Nfree*(Nfree-1)/2
      ier=spr_memalloc(1,dble(mm*2*8))
         if(ier<0) goto 10
      Allocate (Hij(mm),Uiijj(mm))

      mm=Nfree*(Nfree-1)
      ier=spr_memalloc(1,dble(mm*2*8))
         if(ier<0) goto 10
      Allocate (Tiij(mm),Uiiij(mm))

      ! -- 3MR --
      mm=Nfree*(Nfree-1)*(Nfree-2)/6
      ier=spr_memalloc(1,dble(mm*2*8))
         if(ier<0) goto 10
      Allocate (Tijk(mm))

      mm=mm*3
      ier=spr_memalloc(1,dble(mm*2*8))
         if(ier<0) goto 10
      Allocate (Uiijk(mm))

!
! Read data
!
      Open(unit=7,file=PotDir(:Len_PotDir)//'001.hs',status='OLD')
      Call qff_Readhs
      Close(7)
      E0=0.D+00
!
!===============================================================================
!
      !quiet Write(Iout,300)
      !quiet Call spr_meminfo
      !quiet Call timer(1)
      return

      10 Write(Iout,*) 'ALLOCATION ERROR OF QFF DATA. IER=',ier; return
!
!
  100 Format(//,'(  ENTER QFF MODULE  )',//, &
             3x,'o NUMBER OF MODES = ',i8)
  300 Format(/,'(  SETUP OF QFF MODULE COMPLETED  )',//)
!
  END SUBROUTINE  qff_Construct
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
  SUBROUTINE qff_Destruct
!
  USE qff_mod
!
  Implicit None
!
  Real(8) :: memsz
!
      memsz=dble(size(EC))+dble(size(Gi))+dble(size(Hii)) &
           +dble(size(Tiii))+dble(size(Uiiii)) &
           +dble(size(Hij))+dble(size(Uiijj))  &
           +dble(size(Tiij))+dble(size(Uiiij)) &
           +dble(size(Tijk))+dble(size(Uiijk))
      memsz=memsz*8.D+00
      Call spr_memdealloc(memsz)
      Deallocate (EC)
      Deallocate (Gi,Hii,Tiii,Uiiii)
      Deallocate (Hij,Uiijj,Tiij,Uiiij)
      Deallocate (Tijk,Uiijk)

      !quiet Write(Iout,100) 
      !quiet Call spr_meminfo
      !quiet Call timer(1)
!
      return
!
!
  100 Format(/,'(  EXIT QFF MODULE... EXECUTED NORMALLY  ) ',/)
!
  END SUBROUTINE qff_Destruct
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
  SUBROUTINE qff_Readhs
!
  USE qff_mod
!
  Implicit None
!
!--------------------------------------------------------------------------------
!
    Integer :: i,j,k,l,m,n,n1,n2,n3
    Real(8) :: xx
    Character :: ch*80
!
!--------------------------------------------------------------------------------
!
!   >> >> Potential energy and Geometry << <<
!
      Read(7,*)
      Read(7,*) E0
      Read(7,*)

      j=mod(Nfree,3)
      if(j==0) then
         n=Nfree/3
         Do i=1,n
            Read(7,*) (EC(3*(i-1)+j),j=1,3)
         End do
      else
         n=(Nfree-j)/3
         Do i=1,n
            Read(7,*) (EC(3*(i-1)+j),j=1,3)
         End do
         Read(7,*) (EC(j),j=3*n+1,Nfree)
      endif

!
!   >> >> 1 MR << <<
!
      Read(7,'(5x,a)') tl1

      Read(7,*)
      Do i=1,Nfree
         Read(7,*) n,xx
         Gi(i)=xx
      End do
      Read(7,*)
      Do i=1,Nfree
         Read(7,*) n,xx
         Hii(i)=xx
      End do
      Hii=Hii*0.5D+00

      Read(7,*)
      Do i=1,Nfree
         Read(7,*) n,xx
         Tiii(i)=xx
      End do
      Tiii=Tiii/6.D+00

      Read(7,*)
      Do i=1,Nfree
         Read(7,*) n,xx
         Uiiii(i)=xx
      End do
      Uiiii=Uiiii/24.D+00

!
!   >> >> 2 MR << <<
!
      Read(7,'(5x,a)',end=110) tl2

      Read(7,*)
      k=1
      Do i=2,Nfree
      Do j=1,i-1
         Read(7,*) n1,n2,xx
         Hij(k)=xx
         k=k+1
      End do
      End do

      Read(7,*)
      k=1
      Do i=2,Nfree
      Do j=1,i-1
         Read(7,*) n1,n2,xx
         Uiijj(k)=xx
         k=k+1
      End do
      End do
      Uiijj=Uiijj*0.25D+00

      Read(7,*)
      k=1
      Do i=2,Nfree
      Do j=1,i-1
         Read(7,*) n1,n2,xx
         Tiij(k)=xx
         Read(7,*) n1,n2,xx
         Tiij(k+1)=xx
         k=k+2
      End do
      End do
      Tiij=Tiij*0.5D+00

      Read(7,*)
      k=1
      Do i=2,Nfree
      Do j=1,i-1
         Read(7,*) n1,n2,xx
         Uiiij(k)=xx
         Read(7,*) n1,n2,xx
         Uiiij(k+1)=xx
         k=k+2
      End do
      End do
      Uiiij=Uiiij/6.D+00

!   >> >> 3 MR << <<
!
      Read(7,'(5x,a)',end=100) tl3

      Read(7,*)
      l=1
      Do i=3,Nfree
      Do j=2,i-1
      Do k=1,j-1
         Read(7,*) n1,n2,n3,xx
         Tijk(l)=xx
         l=l+1
      End do
      End do
      End do

      Read(7,*)
      l=1
      Do i=3,Nfree
      Do j=2,i-1
      Do k=1,j-1
         Read(7,*) n1,n2,n3,xx
         Uiijk(l)=xx
         Read(7,*) n1,n2,n3,xx
         Uiijk(l+1)=xx
         Read(7,*) n1,n2,n3,xx
         Uiijk(l+2)=xx
         l=l+3
      End do
      End do
      End do
      Uiijk=Uiijk*0.5D+00

      MR=3
      return

!--------------------------------------------------------------------------------
!
  110 Continue
      Write(*,*) 'No 2MR QFF found'
      Hij = 0.0D+00
      Uiijj = 0.0D+00
      Tiij = 0.0D+00
      Uiiij = 0.0D+00
      MR=1
      return

  100 Continue
      Write(*,*) 'No 3MR QFF found'
      Tijk = 0.0D+00
      Uiijk = 0.0D+00
      MR=2
      return
!
  END SUBROUTINE qff_Readhs
!


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!MK Added by MK
!
! SUBROUTINE qff_Readhrm
!
! USE qff_mod
!
! Implicit None
!
!--------------------------------------------------------------------------------
!
!   Integer :: i,j,k,l,m,n,n1,n2,n3
!   Real(8) :: xx
!   Character :: ch*80
!
!--------------------------------------------------------------------------------
!
!!  >> >> Potential energy and Geometry << <<
!
!     Read(7,*)
!     Read(7,*) E0
!     Read(7,*)

!     j=mod(Nfree,3)
!     if(j==0) then
!        n=Nfree/3
!        Do i=1,n
!           Read(7,*) (EC(3*(i-1)+j),j=1,3)
!        End do
!     else
!        n=(Nfree-j)/3
!        Do i=1,n
!           Read(7,*) (EC(3*(i-1)+j),j=1,3)
!        End do
!        Read(7,*) (EC(j),j=3*n+1,Nfree)
!     endif

!
!!  >> >> 1 MR << <<
!
!     Read(7,'(5x,a)') tl1

!     Read(7,*)
!     Do i=1,Nfree
!        Read(7,*) 
!        Gi(i)=0.D+00
!     End do
!     Read(7,*)
!     Do i=1,Nfree
!        Read(7,*) n,xx
!        Hii(i)=xx
!     End do
!     Hii=Hii*0.5D+00

!     Read(7,*)
!     Do i=1,Nfree
!        Read(7,*) 
!        Tiii(i)=0.D+00
!     End do
!     Tiii=Tiii/6.D+00

!     Read(7,*)
!     Do i=1,Nfree
!        Read(7,*) 
!        Uiiii(i)=0.D+00
!     End do
!     Uiiii=Uiiii/24.D+00

!     Write(*,*) 'Harmonic approximation is used.(Only 1MR quadratic force constants are nonzero.)'
!     Hij = 0.0D+00
!     Uiijj = 0.0D+00
!     Tiij = 0.0D+00
!     Uiiij = 0.0D+00
!     MR=1
!     return

!
! END SUBROUTINE qff_Readhrm



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
  SUBROUTINE qff_PES(QQ,V)
!
  USE qff_mod
!
  Implicit None
!
!--------------------------------------------------------------------------------
!
! Geometry in normal coordinate (Angs(amu)^1/2)
!
    Double Precision, dimension(Nfree), intent(in):: QQ
    Double Precision, intent(out):: V

!--------------------------------------------------------------------------------
!
!
    Integer :: i,j,k,l,n,a1,a2

    Double Precision :: Di,Dj,Dk,Dl,G1,H1,T1,U1,H2,T2,U21,U22,T3,U3,U4
!
!--------------------------------------------------------------------------------
!
!   >>  initialize  <<
!
!
        G1=0.D+00; H1=0.D+00; T1=0.D+00; U1=0.D+00
        H2=0.D+00; T2=0.D+00; U21=0.D+00; U22=0.D+00
        T3=0.D+00; U3=0.D+00
        U4=0.D+00

        Do i=1,Nfree
           Di=QQ(i) - EC(i)

           G1 = G1 + Di * Gi(i)
           H1 = H1 + Di*Di * Hii(i)
           T1 = T1 + Di*Di*Di * Tiii(i)
           U1 = U1 + Di*Di*Di*Di * Uiiii(i)

        End do

        a1=1
        a2=1
        Do i=2,Nfree
           Di=QQ(i) - EC(i)
           Do j=1,i-1
              Dj = QQ(j) - EC(j)

              H2 = H2 + Di*Dj * Hij(a1)
              U21= U21+ Di*Di*Dj*Dj * Uiijj(a1)
              a1=a1+1
              T2 = T2 + Di*Di*Dj * Tiij(a2)
              T2 = T2 + Dj*Dj*Di * Tiij(a2+1)
              U22= U22+ Di*Di*Di*Dj * Uiiij(a2)
              U22= U22+ Dj*Dj*Dj*Di * Uiiij(a2+1)
              a2=a2+2

           End do
        End do

        if(MR==3) then
           a1=1
           a2=1
           Do i=3,Nfree
              Di=QQ(i) - EC(i)
              Do j=2,i-1
                 Dj=QQ(j) - EC(j)
                 Do k=1,j-1
                    Dk=QQ(k) - EC(k)

                    T3 = T3 + Di*Dj*Dk * Tijk(a1)
                    a1=a1+1
                    U3 = U3 + Di*Di*Dj*Dk * Uiijk(a2)
                    U3 = U3 + Dj*Dj*Dk*Di * Uiijk(a2+1)
                    U3 = U3 + Dk*Dk*Di*Dj * Uiijk(a2+2)
                    a2=a2+3
   
                 End do
              End do
           End do
        endif

        V = E0 + & 
            G1 + H1 + T1 + U1  &
          + H2 + T2 + U21+ U22 &
          + T3 + U3

!
!--------------------------------------------------------------------------------
!
     return
!
!
  END SUBROUTINE qff_PES
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
  SUBROUTINE qff_PES1(Ni,Qi,V)
!
  USE qff_mod
!
  Implicit None
!
!--------------------------------------------------------------------------------
!
! Geometry in normal coordinate (Angs(amu)^1/2)
!
    Integer, intent(in):: Ni
    Double Precision, intent(in):: Qi
    Double Precision, intent(out):: V

!--------------------------------------------------------------------------------
!
!
    Integer :: i
    Double Precision :: Di,G1,H1,T1,U1
!
!--------------------------------------------------------------------------------
!
!   >>  initialize  <<
!
!
        G1=0.D+00; H1=0.D+00; T1=0.D+00; U1=0.D+00

        i=Ni
        Di=Qi - EC(i)
        G1 = G1 + Di * Gi(i)
        H1 = H1 + Di*Di * Hii(i)
        T1 = T1 + Di*Di*Di * Tiii(i)
        U1 = U1 + Di*Di*Di*Di * Uiiii(i)

        V = E0 + & 
            G1 + H1 + T1 + U1

!
!--------------------------------------------------------------------------------
!
     return
!
!
  END SUBROUTINE qff_PES1
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
  SUBROUTINE qff_PES2(Ni,Nj,Qi,Qj,V)
!
  USE qff_mod
!
  Implicit None
!
!--------------------------------------------------------------------------------
!
!
! Index of normal modes (Ni>Nj)
    Integer, intent(in):: Ni,Nj
! Geometry in normal coordinate (Angs(amu)^1/2)
    Double Precision, intent(in):: Qi,Qj
! 2MR-PE
    Double Precision, intent(out):: V

!--------------------------------------------------------------------------------
!
!
    Integer :: i,j,k,l,n,a1,a2

    Double Precision :: H2,T2,U21,U22
!
!--------------------------------------------------------------------------------
!
!   >>  initialize  <<
!
!
        H2=0.D+00; 
        U21=0.D+00; U22=0.D+00
        T2 =0.D+00

        a1=(Ni-2)*(Ni-1)/2 + Nj
        a2=a1*2 -1

        H2 = H2 + Qi*Qj * Hij(a1)
        U21= U21+ Qi*Qi*Qj*Qj * Uiijj(a1)
        T2 = T2 + Qi*Qi*Qj * Tiij(a2)
        T2 = T2 + Qj*Qj*Qi * Tiij(a2+1)
        U22= U22+ Qi*Qi*Qi*Qj * Uiiij(a2)
        U22= U22+ Qj*Qj*Qj*Qi * Uiiij(a2+1)
   
        V = H2 + T2 + U21 + U22

!
!--------------------------------------------------------------------------------
!
     return
!
!
  END SUBROUTINE qff_PES2
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
  SUBROUTINE qff_PES3(Ni,Nj,Nk,Qi,Qj,Qk,V)
!
  USE qff_mod
!
  Implicit None
!
!--------------------------------------------------------------------------------
!
!
! Index of normal modes (Ni>Nj>Nk)
    Integer, intent(in):: Ni,Nj,Nk
! Geometry in normal coordinate (Angs(amu)^1/2)
    Double Precision, intent(in):: Qi,Qj,Qk
! 3MR-PE
    Double Precision, intent(out):: V

!--------------------------------------------------------------------------------
!
!
    Integer :: i,j,k,l,n,a1,a2

    Double Precision :: T3,U3,U4
!
!--------------------------------------------------------------------------------
!
!   >>  initialize  <<
!
!
        T3=0.D+00; U3=0.D+00

        a1=(Ni-3)*(Ni-2)*(Ni-1)/6 + (Nj-2)*(Nj-1)/2 + Nk
        a2=(Ni-3)*(Ni-2)*(Ni-1)/2 + (Nj-2)*(Nj-1)/2*3 + Nk*3 - 2

        T3 = T3 + Qi*Qj*Qk * Tijk(a1)
        U3 = U3 + Qi*Qi*Qj*Qk * Uiijk(a2)
        U3 = U3 + Qj*Qj*Qk*Qi * Uiijk(a2+1)
        U3 = U3 + Qk*Qk*Qi*Qj * Uiijk(a2+2)
   
        V =T3 + U3

!
!--------------------------------------------------------------------------------
!
     return
!
!
  END SUBROUTINE qff_PES3
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetGi(i)
!
  USE qff_mod

  Implicit None

    Integer :: i

    Double Precision :: qff_GetGi

       qff_GetGi=Gi(i)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetHii(i)
!
  USE qff_mod

  Implicit None

    Integer :: i

    Double Precision :: qff_GetHii

       qff_GetHii=Hii(i)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetTiii(i)
!
  USE qff_mod

  Implicit None

    Integer :: i

    Double Precision :: qff_GetTiii

       qff_GetTiii=Tiii(i)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetUiiii(i)
!
  USE qff_mod

  Implicit None

    Integer :: i

    Double Precision :: qff_GetUiiii

       qff_GetUiiii=Uiiii(i)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetHij(k)
!
  USE qff_mod

  Implicit None

    ! k = (i-1)*(i-2)/2 + j  (i>j)
    Integer :: k

    Double Precision :: qff_GetHij

       qff_GetHij=Hij(k)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetUiijj(k)
!
  USE qff_mod

  Implicit None

    ! k = (i-1)*(i-2)/2 + j  (i>j)
    Integer :: k

    Double Precision :: qff_GetUiijj

       qff_GetUiijj=Uiijj(k)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetTiij(k,iopt)
!
  USE qff_mod

  Implicit None

    ! k = (i-1)*(i-2)/2 + j  (i>j)
    Integer :: k
    ! iopt =0 : tiij, =1: tjji
    Integer :: iopt

    Double Precision :: qff_GetTiij

       qff_GetTiij=Tiij(2*k-1+iopt)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetUiiij(k,iopt)
!
  USE qff_mod

  Implicit None

    ! k = (i-1)*(i-2)/2 + j  (i>j)
    Integer :: k
    ! iopt =0 : uiiij, =1: ujjji
    Integer :: iopt

    Double Precision :: qff_GetUiiij

       qff_GetUiiij=Uiiij(2*k-1+iopt)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetTijk(l)
!
  USE qff_mod

  Implicit None

    ! l = (i-1)*(i-2)*(i-3)/6 + (j-1)*(j-2)/2 + k  (i>j>k)
    Integer :: l

    Double Precision :: qff_GetTijk

       qff_GetTijk=Tijk(l)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  FUNCTION qff_GetUiijk(l,iopt)
!
  USE qff_mod

  Implicit None

    ! l = (i-1)*(i-2)*(i-3)/6 + (j-1)*(j-2)/2 + k  (i>j>k)
    Integer :: l
    ! iopt =0 : uiijk, =1: uijjk, =2: uijkk
    Integer :: iopt

    Double Precision :: qff_GetUiijk

       qff_GetUiijk=Uiijk(3*l-2+iopt)

  END FUNCTION
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE qff_GetTitle(ch1,ch2,ch3)
!
  USE qff_mod

  Implicit None

    Character :: ch1*(*),ch2*(*),ch3*(*)

       ch1=tl1; ch2=tl2; ch3=tl3

  END SUBROUTINE
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
