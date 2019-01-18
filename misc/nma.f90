!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!   Code description by K.Yagi:  yagi@qcl.t.u-tokyo.ac.jp
!
!   Last modified : 2002/02/28
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!=======================================================================================================
!  The information of equilibrium structure.
!
!        x0   : Equilibrium structure (Angs)
!        CL   : Normal displacement vector at EQ
!
!=======================================================================================================
!
!
Module nma_private
!
!     ----------------------------------(Parameters)----------------------------------
!    |    Nat  : The number of atoms                                                  |
!    |    Nat3 : The number of atoms * 3                                              |
!    |    Nfree: The number of internal degrees of freedom, i.e. Nat3-6(or Nat3-5)    |
!    |    Zmass: Mass of each atom                                                    |
!     --------------------------------------------------------------------------------
!
  Integer :: Nat, Nat3, Nfree
!
! MASS
  Double Precision, dimension(:), Allocatable :: Zmass
!
! GEOMETRY
  Double Precision, dimension(:,:), Allocatable:: x0
!
! NORMAL DISPLACEMENT VECTOR
  Double Precision, dimension(:,:), Allocatable:: CL
!
  Logical :: setup
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
End Module
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!
  SUBROUTINE nma_Construct(N,Nf,Z,x,L)

  USE nma_private

  Implicit None
!
    Integer, intent(in):: N,Nf
    Double Precision, dimension(N) :: Z
    Double Precision, dimension(3,N) :: x
    Double Precision, dimension(N*3,Nf) :: L
!
      Nat=N
      Nat3=N*3
      Nfree=Nf
!
      Allocate(Zmass(Nat),x0(3,Nat),CL(Nat3,Nfree))
 
      Zmass=Z
      x0=x
      CL=L
      setup=.true. 
!
  End SUBROUTINE nma_Construct
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!
  SUBROUTINE nma_Destruct

  USE nma_private

  Implicit None
!

       Deallocate(Zmass,x0,CL)
       setup=.false.
       return

!
  End SUBROUTINE nma_Destruct
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!
  SUBROUTINE nma_Print(io)

  USE nma_private

  Implicit None
 
    Integer :: io,i

      if(.NOT.setup) Call nma_Error(1)

      Write(io,'(6x,''o ATOMIC MASS'')')
      Write(io,'(6x,5f12.4)') Zmass
      Write(io,*)
      Write(io,'(6x,''o REFERENCE GEOMETRY'')')
      Write(io,'(6x,3f12.6)') x0
      Write(io,*)
      Write(io,'(6x,''o NORMAL DISPLACEMENT VECTOR'')')
      Do i=1,Nat
         Write(io,'(6x,3f12.6)') CL(:,i)
         Write(io,*)
      End do

  End SUBROUTINE nma_Print
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!
  SUBROUTINE nma_ChkCL(Ierr)

  USE nma_private

  Implicit None

    Integer, intent(in):: Ierr

    Integer :: i,j,k
    Double Precision:: x,tmp,prcn, err

      if(.NOT.setup) Call nma_Error(1)
!
!
! ===== Check : precision criteria is 'prcn' = 1.d-05 =====
!
      prcn=1.D-05
      err =0.D+00
      Write(IERR,*)
      Write(IERR,*) '====( NMA_CHKCL::  CHECKING THE PRECISION OF CL )===='

      Do i=1,Nfree
         Do j=1,i

            tmp=0.D+00
            Do k=1,Nat3
               tmp = tmp + CL(k,j)*CL(k,i)
            End do

            if(i==j) then
               x=abs(tmp-1.D+00)
               if(x>prcn) then
                  Write(Ierr,*) '== WARNING == NORMAL DISPLACEMENT VECTOR IS NOT NORMALIZED'
                  Write(Ierr,'(3x,''|L('',i2,'')| = '',f8.6)') i,tmp 
               endif
            else
               x=abs(tmp)
               if(x>prcn) then
                  Write(Ierr,*) '== WARNING == NORMAL DISPLACEMENT VECTOR IS NOT ORTHOGONAL'
                  Write(Ierr,'(3x,''L('',i2,'')*L('',i2,'') = '',f8.6)') i,j,tmp 
               endif
            endif
            if(x>err) err = x

         End do
      End do
!
!
      Write(Ierr,*)
      if(err<1.D-04) then
         Write(Ierr,100) err
      else
         Write(Ierr,200) 
      endif
      Write(Ierr,*) '====( NMA_CHKCL::  END )===='
      Write(IERR,*)
!
      return

  100 Format(' MAX ERRROR IS ',e10.3,' :    [PASSED]',/)
  200 Format(' MAX ERRROR IS LARGER THAN 1.D-04:    [ERROR]  ',/)
!
!
  END SUBROUTINE
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! x :  Geometry in Cartesian coordinate (length)
! q :  Geometry in Normal coordinate (length amu1/2)
!
  SUBROUTINE nma_x2q(x,q)

  USE nma_private

  Implicit None
!
    Integer :: i,j,k,kk
    Double Precision  :: AMS,TMP

    Double Precision, intent(in):: x(3,Nat)
    Double Precision, intent(out):: q(Nfree)
!
!
    if(.NOT.setup) Call nma_Error(1)
    q = 0.D+00

    kk=0
    Do i=1,Nat
       AMS = SQRT(Zmass(i))
       Do j=1,3
          kk=kk+1
          TMP = (x(j,i)-x0(j,i))*AMS
          Do k=1,Nfree
             q(k) = q(k)+TMP*CL(kk,k)
          End do
       End do
    End do

    return

  END SUBROUTINE nma_x2q
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! x :  Geometry in Cartesian coordinate (length)
! q :  Geometry in Normal coordinate (length amu1/2)
!
  SUBROUTINE nma_q2x(x,q)

  USE nma_private

  Implicit None
!
    Integer :: i,j,k,kk
    Double Precision  :: AMS

    Double Precision, intent(out):: x(3,Nat)
    Double Precision, intent(in):: q(Nfree)
!
    if(.NOT.setup) Call nma_Error(1)
    x=x0
!
    kk=0
    Do i=1,Nat
       AMS = SQRT(Zmass(i))
       Do j=1,3
          kk=kk+1
          Do k=1,Nfree
             x(j,i)=x(j,i) + q(k)*CL(kk,k)/AMS
          End do
       End do
    End do
    return

  END SUBROUTINE nma_q2x
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!  v: velocity in Cartesian ( length/time )
!  p: velocity/momentum in normal Coordinate ( length amu1/2/time )
!
  SUBROUTINE nma_v2p(v,p)

  USE nma_private

  Implicit None
!
    Integer :: i,j,k,kk
    Double Precision  :: AMS,TMP

    Double Precision, intent(in):: v(3,Nat)
    Double Precision, intent(out):: p(Nfree)
!
!
    if(.NOT.setup) Call nma_Error(1)
    p = 0.D+00

    kk=0
    Do i=1,Nat
       AMS = SQRT(Zmass(i))
       Do j=1,3
          kk=kk+1
          TMP = v(j,i)*AMS
          Do k=1,Nfree
             p(k) = p(k)+TMP*CL(kk,k)
          End do
       End do
    End do
!
    return
!
  END SUBROUTINE nma_v2p
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!  v: velocity in Cartesian ( length/time )
!  p: velocity/momentum in normal Coordinate ( length amu1/2/time )
!
  SUBROUTINE nma_p2v(v,p)

  USE nma_private

  Implicit None
!
    Integer :: i,j,k,kk
    Double Precision  :: AMS

    Double Precision, intent(out):: v(3,Nat)
    Double Precision, intent(in):: p(Nfree)
!
!
    if(.NOT.setup) Call nma_Error(1)
    v = 0.D+00
!
    kk=0
    Do i=1,Nat
       AMS = SQRT(Zmass(i))
       Do j=1,3
          kk=kk+1
          Do k=1,Nfree
             v(j,i)=v(j,i) + p(k)*CL(kk,k)/AMS
          End do
       End do
    End do
!
    return
!
  END SUBROUTINE nma_p2v
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! g : Gradient in Cartesian coordinate ( energy/length )
! qg: Gradient in Normal coordinate ( energy/length amu1/2 )
!
  SUBROUTINE nma_g2qg(g,qg)

  USE nma_private

  Implicit None
!
    Integer :: i,j,k,kk
    Double Precision  :: AMS,TMP

    Double Precision, intent(in):: g(3,Nat)
    Double Precision, intent(out):: qg(Nfree)
!
    if(.NOT.setup) Call nma_Error(1)
!
! initialize
    qg=0.D+00
!
    kk=0
    Do i=1,Nat
       AMS = SQRT(Zmass(i))
       Do j=1,3
          kk=kk+1
          TMP = g(j,i)/AMS
          Do k=1,Nfree
             qg(k)= qg(k)+TMP*CL(kk,k)
          End do
       End do
    End do
 
    return
!
  END SUBROUTINE nma_g2qg
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! g : Gradient in Cartesian coordinate ( energy/length )
! qg: Gradient in Normal coordinate ( energy/length amu1/2 )
!
  SUBROUTINE nma_qg2g(g,qg)

  USE nma_private

  Implicit None
!
    Integer :: i,j,k,kk
    Double Precision  :: AMS

    Double Precision, intent(out):: g(3,Nat)
    Double Precision, intent(in):: qg(Nfree)
!
    if(.NOT.setup) Call nma_Error(1)
!
! initialize
    g=0.D+00
!
    kk=1
    Do i=1,Nat
       AMS=SQRT(ZMASS(i))
       Do j=1,3
          Do k=1,Nfree
             g(j,i)=g(j,i)+qg(k)*CL(kk,k)
          End Do
          g(j,i) = g(j,i)*AMS
          kk=kk+1
       End do
    End do
 
    return
! 
  END SUBROUTINE nma_qg2g
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
! h : Hessian in Cartesian coordinate ( energy/length^2 )
! qh: Hessian in Normal coordinate ( energy/length^2(amu) )
!
  SUBROUTINE nma_h2qh(h,qh)

  USE nma_private

  Implicit None
!
    Integer :: i,ii,j,jj,k,kk
    Double Precision  :: AMS

    Double Precision, intent(in):: h(Nat3,Nat3)
    Double Precision, intent(out):: qh(Nfree,Nfree)
    Double Precision, dimension(Nat3,Nfree) :: CL2
!
    if(.NOT.setup) Call nma_Error(1)
!
! initialize
!
    qh = 0.D+00
!
    kk=1
    Do i=1,Nat
       AMS = SQRT(Zmass(i))
       Do j=1,3
          Do k=1,Nfree
            CL2(kk,k) = CL(kk,k)/AMS
          End do
          kk=kk+1
       End do
    End do
!
    Do ii=1,Nfree
       Do jj=1,ii
          Do i=1,Nat3
             Do j=1,Nat3
                qh(jj,ii) = qh(jj,ii)+CL2(i,ii)*h(i,j)*CL2(j,jj)
  
             End do
          End do
          qh(ii,jj) = qh(jj,ii)
       End do
    End do
!
  return
!
!
  END SUBROUTINE nma_h2qh
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
!
!  q        : A (unit) vector of 3N-6 dimension in Normal Coordinate
!  rm    : Reduced Mass
!
!       rm = (|q|/|x|)^2
!           x = M^(-1/2)*CL^(-1)*q
!
  SUBROUTINE nma_getrm(q,rm)

  USE nma_private

  Implicit None
!
    Integer :: i,j,k,l
    Double Precision  :: qq,xx,tmp,AMS

    Double Precision, dimension(Nfree), intent(in):: q
    Double Precision, intent(out):: rm
 
    Double Precision, dimension(Nat3):: x
!
    if(.NOT.setup) Call nma_Error(1)
!
    x=0.D+00
    qq=0.D+00
    qq=Dot_Product(q,q)
!
!   Check 
!
    tmp=ABS(qq-1.D+00)
    if(tmp>=1.D-08) Write(6,100)
!
    l=0
    Do i=1,Nat
       AMS = SQRT(Zmass(i))
       Do j=1,3
          l=l+1
          Do k=1,Nfree
             x(l)=x(l) + q(k)*CL(l,k)/AMS
          End do
       End do
    End do
!
    xx=0.D+00
    xx=Dot_Product(x,x)
!
    rm=qq/xx
!
    return
!
    100 Format(10x,'=== WARNING (GET_RED_MASS) == GIVEN Q IS NOT A UNIT VECTOR')
!
!
  END SUBROUTINE nma_getrm
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10

  SUBROUTINE nma_Error(ierr)

    Integer :: ierr

       Write(6,*) '   >(ERROR)> ERROR IN NMA_MODULE.'

       Select case(ierr)
          case(1)
             Write(6,*) '   >(ERROR)> NMA_MODULE IS NOT SETUP.'
             Write(6,*) '   >(ERROR)> CONSTRUCT THE MODULE FIRST!'
       End Select
       Stop

  END SUBROUTINE nma_Error
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----10
