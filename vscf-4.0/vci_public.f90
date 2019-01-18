!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!   Last modified  2007/12/09
!   Code description by K.Yagi
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
SUBROUTINE Vci_Construct(ierr)

  USE Vci_mod

  IMPLICIT NONE

  INTEGER :: ierr
  INTEGER :: i,j,k,l,in,iout,spr_memalloc
  INTEGER :: maxEx(100),maxExAll

  INTEGER, PARAMETER :: nmax=10000,m=60
  INTEGER, DIMENSION(m*nmax):: state
  INTEGER, DIMENSION(m,nmax):: state_label

  NAMELIST /vci/Nstate,nCI,maxEx,maxExAll,maxSum,nCUP, &
       nLvl,nGen,pth0,pth1,pth2,dump,ss,state,state_label

  ! ----------------------------------------------------------------

  ierr=0
  CALL spr_Getio(in,iout)

  ! ------------------------------------------------------------
  ! >> Read input 

  ! --- default ---

  ! - VCI space selection 
  !
  ! > One-to-All CI
  Nstate=-1
  nCI=0
  maxEx=-1
  maxExAll=-1
  maxSum=-1
  nCUP=4

  ! > State-specific CI
  nLvl=0
  nGen=3
  pth0=500
  pth1=1.D-01
  pth2=0.9D+00
  ss=.FALSE.

  REWIND(in)
  READ(in,vci,END=10)
10 CONTINUE

  WRITE(iout,100)
100 FORMAT(//,'(  ENTER VCI MODULE  )',//, &
       3x,'>> VCI OPTIONS',/)

  IF(nLvl==0) THEN

     WRITE(iout,105)

     WRITE(iout,110) Nstate
     Nst=Nstate

     WRITE(iout,120)

     IF(nCUP < 1 .OR. nCUP > Nfree) nCUP=Nfree
     WRITE(iout,130) nCUP

     IF(maxSum>0) THEN 
        WRITE(iout,131) maxSum
        !write(666,131) maxSum!MK
     ELSE
        WRITE(iout,132) 
     ENDIF

     ierr=spr_memalloc(1,DBLE(Nfree)*4.D+00)
     IF(ierr<0) RETURN
     ALLOCATE(maxExc(Nfree))
     DO i=1,Nfree
        IF(maxEx(i)/=-1) THEN 
           maxExc(i)=maxEx(i)+1
        ELSEIF(maxExAll/=-1) THEN
           maxExc(i)=maxExAll+1
        ELSE
           maxExc(i)=nCHO(i)
        ENDIF
        IF(maxExc(i)>nCHO(i)) maxExc(i)=nCHO(i)
        IF(maxSum>0 .AND. maxExc(i)>maxSum+1) maxExc(i)=maxSum+1
     END DO
     IF(Nfree<4) THEN 
        WRITE(iout,121) (maxExc(i)-1,i=1,Nfree)
     ELSE
        WRITE(iout,121) (maxExc(i)-1,i=1,3)
        WRITE(iout,122) (maxExc(i)-1,i=4,Nfree)
     ENDIF

     IF(maxSum>0) THEN
        k=MAXVAL(maxExc)-1
        j=1; l=1
        DO i=1,nCUP
           l=l*(Nfree+1-i)/i*(maxSum+1-i)/i
           j=j+l
        END DO

     ELSE
        maxSum=0
        DO i=1,Nfree
           IF(maxExc(i)>0) maxSum=maxSum + maxExc(i)-1
        END DO
        !maxSum=maxSum-Nfree

        k=MAXVAL(maxExc)-1
        j=1; l=1
        DO i=1,nCUP
           l=l*(Nfree+1-i)/i*k
           j=j+l
        END DO
     ENDIF
     IF(nCI>j .OR. nCI==0) nCI=j

     WRITE(iout,133) nCI
     WRITE(iout,*)

105  FORMAT(5x,'[  ONE-TO-ALL CI  ]',/)
110  FORMAT(7x,'o NUM_OF_STATES  : ',i10)
120  FORMAT(7x,'o CI SPACE')
121  FORMAT(9x,'- MAX Q        :  ',3i3)
122  FORMAT(9x,'                  ',3i3)
130  FORMAT(9x,'- MODE COUPLING:    ',i7)
131  FORMAT(9x,'- MAX SUM_OF_Q :    ',i7)
132  FORMAT(9x,'- MAX SUM_OF_Q :  UNLIMITED')
133  FORMAT(9x,'- MAX N_CONFIG.:    ',i7)

  ELSE

     IF(ss) THEN 
        WRITE(iout,145)
     ELSE
        WRITE(iout,146)
        CALL Vcor_setTarget(state,state_label(1:Nfree,1:Nstate))

        IF(Nfree<=30) THEN
           DO j=1,Nstate
              WRITE(iout,115) j,(label0(k+(j-1)*Nfree),k=1,Nfree)
           END DO
        ELSE
           DO j=1,Nstate
              WRITE(iout,115) j,(label0(k+(j-1)*Nfree),k=1,30)
              WRITE(iout,116) (label0(k+(j-1)*Nfree),k=31,Nfree)
           END DO
        ENDIF
        WRITE(iout,*)

     ENDIF
115  FORMAT(7x,'TARGET STATE ',i3.3,':',3x,6(5i1,1x))
116  FORMAT(27x,6(5i1,1x))

     IF(.NOT. lvscf) THEN
        WRITE(iout,*)
        WRITE(iout,'(8x,''ERROR:  STATE-SPECIFIC CI REQUIRES VSCF REFERENCES '')')
        WRITE(iout,'(8x,''ERROR:  RE-RUN WITH VSCF=.TRUE. '')')
        WRITE(iout,*)
        STOP
     ENDIF

     WRITE(iout,150) nLvl
     WRITE(iout,154) nGen
     WRITE(iout,155) pth0
     WRITE(iout,156) pth1
     WRITE(iout,157) pth2
     WRITE(iout,*)

145  FORMAT(5x,'[  STATE-SPECIFIC CI  ]',/)
146  FORMAT(5x,'[  STATE-SPECIFIC CI WITH ZERO-POINT VSCF REFERENCE  ]',/)
150  FORMAT(7x,'o LEVEL    : ',i10)
154  FORMAT(7x,'o N_GEN    : ',i10)
155  FORMAT(7x,'o PTH0     : ',f10.1,' [CM-1]')
156  FORMAT(7x,'o PTH1     : ',f10.4)
157  FORMAT(7x,'o PTH2     : ',f10.4)

  ENDIF

  ! ------------------------------------------------------------

  IF(lvscf) THEN
     WRITE(iout,170)
     CALL file_indicator(10,rwf)
     OPEN(rwf,file='vscf-w.wfn',form='UNFORMATTED')
  ENDIF

170 FORMAT(7x,'o READ VSCF WFN  : vscf-w.wfn')

  ! ------------------------------------------------------------

  IF(dump) WRITE(iout,180)
180 FORMAT(7x,'o DUMP VCI WFN   : vci-w.wfn')

  ! ------------------------------------------------------------

  WRITE(iout,200)
  CALL spr_meminfo
  CALL timer(1)
200 FORMAT(//,'(  SETUP OF VCI COMPLETED  )',//)

END SUBROUTINE Vci_Construct

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

SUBROUTINE Vci_finalz()

  USE Vci_mod

  IMPLICIT NONE

  IF(lvscf) CLOSE(rwf)

  IF(ALLOCATED(maxExc)) THEN
     CALL spr_memdealloc(DBLE(Nfree)*4.D+00)
     DEALLOCATE(maxExc)
  ENDIF

END SUBROUTINE Vci_finalz

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

SUBROUTINE Vci_main()

  USE Vci_mod

  IMPLICIT NONE

  INTEGER :: i,j,k,l
  INTEGER :: in,io,ierr,spr_memalloc,Ist
  INTEGER :: lb(Nfree)

  CALL spr_Getio(in,io)
  WRITE(io,100)

  CALL Vcor_Construct(ierr)
  IF(ierr<0) THEN
     WRITE(io,*) 'ERROR IN VCI_MAIN WHILE SETTING UP VCOR_MODULE'
     STOP
  ENDIF

  IF(nLvl==0) THEN
     CALL setCIbasis(ierr)
     IF(ierr<0) THEN 
        WRITE(io,*) 'ERROR IN VCI_MAIN WHILE SETTING UP CI BASIS' 
        STOP
     ENDIF
     CALL Vcor_setKinmat()
     CALL Vcor_setWfn()
     CALL Vcor_setQmat()

     IF(Nst>nCI) Nst=nCI
     WRITE(io,110) Nst,nCI
     !write(666,110) Nst,nCI !MK
     IF(Nfree<31) THEN
        WRITE(io,115) ref
     ELSE
        WRITE(io,115) (ref(i),i=1,30)
        WRITE(io,116) (ref(i),i=31,Nfree)
     ENDIF
     WRITE(io,*)

     ierr=spr_memalloc(-1,DBLE(nCI)*DBLE(Nst+1)*8.D+00)
     ALLOCATE(CIwfn(nCI*Nst),CIene(nCI))
     CALL Vci_boost(io)
     CALL Vci_print(io)
     IF(dump) CALL Vci_Dump()
     CALL spr_memdealloc(DBLE(nCI)*DBLE(Nst+1)*8.D+00)
     DEALLOCATE(CIwfn,CIene)


  ELSE
     WRITE(io,120) 

     Ist=1
     DO WHILE (.TRUE.) 

        IF(ss .OR. Ist==1) THEN 
           CALL Vcor_rVSCF(ierr)
           IF(ierr<0) EXIT
           IF(Ist==1 .AND. .NOT. zp) THEN
              WRITE(io,*) 'ERROR:  THE FIRST VSCF REFERENCE MUST BE THE ZERO POINT STATE'
              WRITE(io,*) 'ERROR:  RESET state IN &VSCF'
              STOP
           ENDIF
        ELSE
           IF(Ist==(Nstate+1)) EXIT
           tar=label0(Nfree*(Ist-1)+1:Nfree*Ist)
           zp=.FALSE.
        ENDIF

        CALL Vcor_setKinmat()
        CALL Vcor_setWfn()
        CALL Vcor_setQmat()

        IF(Nfree<31) THEN
           WRITE(io,125) Ist,tar
        ELSE
           WRITE(io,125) Ist,(tar(i),i=1,30)
           WRITE(io,126) (tar(i),i=31,Nfree)
        ENDIF
        WRITE(io,*)

        CALL Vcor_set_pCnf()

        IF(npCnf>0) THEN
           WRITE(io,128) npCnf
           IF(Nfree<31) THEN
              DO i=1,npCnf
                 CALL Vcor_getLabel(i+1,lb)
                 WRITE(io,130) i,lb
              END DO
           ELSE
              DO i=1,npCnf
                 CALL Vcor_getLabel(i+1,lb)
                 WRITE(io,130) j,(lb(j),j=1,30)
                 WRITE(io,131) (lb(j),j=31,Nfree)
              END DO
           ENDIF
           WRITE(io,*)
        ENDIF

        CALL Vcor_set_qCnf()
        WRITE(io,135)
        WRITE(io,136) nqCnf
        DO i=1,nqCUP
           WRITE(io,137) i,nqCnfi(i)
        END DO

        nCI=1+npCnf+nqCnf
        Nst=1+npCnf+nqCnf

        ierr=spr_memalloc(-1,DBLE(nCI)*DBLE(Nst+1)*8.D+00)
        ALLOCATE(CIwfn(nCI*Nst),CIene(nCI))
        CALL Vci_boost(io)
        CALL Vci_print2(io)
        CALL spr_memdealloc(DBLE(nCI)*DBLE(Nst+1)*8.D+00)
        DEALLOCATE(CIwfn,CIene)

        Ist=Ist+1
        CALL Vcor_ClearCnf

        CALL timer(1)
        WRITE(io,*)

     END DO

  ENDIF

  CALL Vcor_finalz()


  WRITE(io,200)
  CALL spr_meminfo
  CALL timer(1)

100 FORMAT(/'(  ENTER VCI MAIN MODULE  )',//)
110 FORMAT(/,3x,'>> ONE-TO-ALL CI-BASIS SELECTION',//, &
       9x,'o NUM_OF_STATES  : ',i6,/, &
       9x,'o NUM_OF_CONFIG. : ',i6,/)
115 FORMAT(9x,'o REFERENCE VSCF :',3x,6(5i3,1x))
116 FORMAT(17x,                ':',3x,6(5i3,1x))
120 FORMAT(/,3x,'[  LOOP OVER VCI STATES  ]',//)
125 FORMAT(5x,'>> STATE ',i3.3,':',3x,6(5i3,1x))
126 FORMAT(21x,6(5i3,1x))
128 FORMAT(9x,'o RESONANT STATES   (P-SPACE) : ',i6,/)
130 FORMAT(12x,'(',i2.2,'):  ',6(5i3,1x))
131 FORMAT(19x,6(5i3,1x))
135 FORMAT(9x,'o CORRELATED STATES (Q-SPACE)',/)
136 FORMAT(12x,'=== TOTAL ===  : ',i10)
137 FORMAT(19x,i2,'-MODE : ',i10)
200 FORMAT(/'(  EXIT VCI MAIN MODULE  )',//)

END SUBROUTINE Vci_main

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

SUBROUTINE Vci_boost(io)

  USE Vci_mod

  IMPLICIT NONE

  INTEGER :: i,j,k,l,io,spr_memalloc
  INTEGER :: lb(Nfree),lb2(Nfree)
  INTEGER :: cp,Vcor_getCUP
  REAL(8) :: Ene,Vcor_getHmat,H2wvn
  REAL(8), ALLOCATABLE :: Hmat(:)

  PRINT*,"Nci = ",nCI
  i=spr_memalloc(-1,DBLE(nCI*(nCI+1)/2)*8.D+00)
  ALLOCATE(Hmat(nCI*(nCI+1)/2))

  ! Generate Hamiltonian matrix
  IF(io>0) WRITE(io,100)
  k=1
  DO i=1,nCI
     cp=Vcor_getCUP(1,i)
     IF(cp>MR) THEN
        Hmat(k)=0.D+00
     ELSEIF(cp/=1) THEN
        Hmat(k)=Vcor_getHmat(cp,1,i)
     ELSE
        IF(.NOT. ss .OR. .NOT. lvscf) THEN 
           Hmat(k)=Vcor_getHmat(cp,1,i)
        ELSE
           Hmat(k)=0.D+00
        ENDIF
     ENDIF
     k=k+1

     DO j=2,i
        cp=Vcor_getCUP(j,i)
        IF(cp>MR) THEN
           Hmat(k)=0.D+00
        ELSE
           Hmat(k)=Vcor_getHmat(cp,j,i)
        ENDIF
        k=k+1
     END DO
  END DO
  IF(io>0) WRITE(io,200)

  IF(io>0) WRITE(io,110)
  !
  ! vci matrix output //shiozaki
  !>>>
  !        open(40,file='hmat',form='unformatted')
  !        write(40) Hmat
  !        stop
  !<<<
  !
  CALL diag(nCI,Nst,Hmat,CIwfn,CIene)
  IF(io>0) WRITE(io,200)

  CALL spr_memdealloc(DBLE(nCI*(nCI+1)/2)*8.D+00)
  DEALLOCATE(Hmat)

100 FORMAT(/,9x,'o FORMING HAMILTONIAN MATRIX')
110 FORMAT(/,9x,'o DIAGONALIZING HAMILTONIAN MATRIX')
200 FORMAT(9x,'  ... DONE')

  !dbg Ene=Vcor_getHmat(0,1,1)
  !dbg write(6,'(f12.4)') Ene*H2wvn()
  !dbg Do i=2,nCI
  !dbg    cp=Vcor_getCUP(1,i)
  !dbg    if(cp<=3) then
  !dbg       Call Vcor_getLabel(1,lb)
  !dbg       Call Vcor_getLabel(i,lb2)
  !dbg       Ene=Vcor_getHmat(cp,1,i)
  !dbg       write(6,'(i2,2x,4i3,2x,4i3,f12.4)') cp,lb,lb2,Ene*H2wvn()
  !dbg    endif
  !dbg End do

  !dbg Do i=1,nCI
  !dbg   Call Vcor_getLabel(i,lb)
  !dbg   Do j=1,i
  !dbg      Call Vcor_getLabel(j,lb2)
  !dbg      cp=Vcor_getCUP(i,j)
  !dbg      write(6,'(2x,3i5)') i,j,cp
  !dbg      write(6,'(3x,10i2)') lb
  !dbg      write(6,'(3x,10i2)') lb2
  !dbg   End do
  !dbg End do

END SUBROUTINE Vci_boost
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!    ist :: an index for the current state

SUBROUTINE Vci_ave(ist)

  USE Vci_mod

  IMPLICIT NONE

  INTEGER :: spr_GetNat,Nat,ie,spr_memalloc
  INTEGER :: i,j,k,l,m,n,i1,i2,n1,n2,j1,j2

  INTEGER :: ist,jst
  INTEGER :: Label1(Nfree),Label2(Nfree)
  Double Precision :: wgt,D1(maxCHO),A0
  !       Real(8), allocatable :: 
  Double Precision , PARAMETER :: thresh=1.D-06
  DOUBLE PRECISION, ALLOCATABLE :: QQ(:),quad(:,:),x0(:,:)
  LOGICAL :: ABC=.true. ,ABC2=.false.  !MK
  DOUBLE PRECISION, ALLOCATABLE :: R0(:,:),Masses(:)!,xfromquad(:,:,:),x0exp(:,:)!MK
  DOUBLE PRECISION :: rot_ABC(3),eValues(3),eVectors(3,3),Inertia(6),work(9)!MK
  DOUBLE PRECISION :: COMx, COMy, COMz
  INTEGER :: info
  Double Precision, PARAMETER :: u_const=1.660538783162726D-027
  Double Precision, PARAMETER :: h_const=6.62606896D-34 
  Double Precision, PARAMETER :: c_const=2.99792458D+10 !in cm/s  
  Nat=spr_GetNat()
  ie=spr_memalloc(1,DBLE(Nfree+3*Nat+maxCHO*Nfree)*16.D+00)
  ALLOCATE(QQ(Nfree),quad(maxCHO,Nfree),x0(3,Nat),R0(Nat,Nat),Masses(Nat))!,xfromquad(maxCHO,3,Nat),x0exp(3,Nat))

  DO i=1,Nfree
     CALL Modal_getQ(i,nCHO(i),quad(:,i))
  END DO

 ! if (ABC2) then ! expectation value of x !MK
 !   Do i=1,maxCHO
 !       Call vav_q2x(Nat,Nfree,xfromquad(i,:,:),quad(i,:))
 !   enddo
 !  endif! end expectation

  jst=(ist-1)*nCI

  QQ=0.D+00
 ! x0exp=0.d0
  DO i=1,nCI
     wgt=CIwfn(i+jst)*CIwfn(i+jst)
     IF(ABS(wgt)<thresh) CYCLE
     CALL Vcor_getLabel(i,Label1)
     DO j=1,Nfree
        n=Label1(j)
        CALL Vcor_getDwfn(j,nCHO(j),n,D1)
        DO k=1,nCHO(j)
           QQ(j)=QQ(j) + wgt*quad(k,j)*D1(k)
!           if (ABC2) then ! expectation value of x !MK
!              Do j1=1,Nat
!                Do j2=1,3
!                 x0exp(j2,j1)=x0exp(j2,j1)+ wgt*xfromquad(k,j2,j1)*D1(k)
!                enddo
!              enddo
!           endif ! end expectation
        END DO
     END DO
  END DO

  DO i1=1,nCI
     CALL Vcor_getLabel(i1,Label1)
     DO i2=1,i1-1
        CALL Vcor_getLabel(i2,Label2)
        wgt=2.D+00*CIwfn(i1+jst)*CIwfn(i2+jst)
        IF(ABS(wgt)<thresh) CYCLE
             
        DO j=1,Nfree
           DO k=1,j-1
              IF(Label1(k)/=Label2(k)) GOTO 1000
           END DO
           DO k=j+1,Nfree
              IF(Label1(k)/=Label2(k)) GOTO 1000
           END DO

           n1=Label1(j)
           n2=Label2(j)
           CALL Vcor_getXwfn(j,nCHO(j),n1,n2,D1)
           DO k=1,nCHO(j)
              QQ(j)=QQ(j) + wgt*quad(k,j)*D1(k)
 !             if (ABC2) then ! expectation value of x !MK
 !             Do j1=1,Nat
 !               Do j2=1,3
 !                x0exp(j2,j1)=x0exp(j2,j1)+ wgt*xfromquad(k,j2,j1)*D1(k)
 !               enddo
 !             enddo
 !             endif ! end expectation
           END DO

1000       CONTINUE
        END DO
     END DO
  END DO

  CALL vav_q2x(Nat,Nfree,x0,QQ)
  CALL spr_GetIO(i,j)
  WRITE(j,100) 
  WRITE(j,200) QQ
  WRITE(j,110) 
  WRITE(j,200) x0
  !MK
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

  WRITE(j,120)
  R0=0.D0;   
  DO j=1,Nat
     DO k=j+1,Nat
        R0(j,k)=SQRT((x0(1,j)-x0(1,k))**2.d0+(x0(2,j)-x0(2,k))**2.d0+(x0(3,j)-x0(3,k))**2.d0)
        Write(6,250)j,k,R0(j,k)  
        R0(k,j)=R0(j,k)
     END DO
  END DO
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  WRITE(6,130)!MK
  A0=0.0;   
  If (Nat .ge. 3) then
    DO i=1,Nat
     DO j=1,Nat
        IF(i==j)CYCLE
        DO k=1,Nat
           IF(j==k .OR. i==k .OR. k>j)CYCLE 
           A0=ACOS((R0(i,j)**2.d0+R0(i,k)**2.d0-R0(k,j)**2.d0)/(2.d0*R0(i,j)*R0(i,k)))*180.0D0/pi
           WRITE(6,300)j,i,k,A0
        END DO
     END DO
    END DO
   Endif
  !MK
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  !08/11/08 MK rotational constants
  IF(ABC) then
    CALL spr_GetMass(Masses)
    COMx = sum( Masses*x0(1,:) ) / sum(Masses)
    COMy = sum( Masses*x0(2,:) ) / sum(Masses)
    COMz = sum( Masses*x0(3,:) ) / sum(Masses)
   Inertia=0.d0
   do k=1,Nat
    Inertia(1)= Inertia(1) + Masses(k) * ( (x0(2,k)-COMy)**2.d0 + (x0(3,k)-COMz)**2.d0 )!Ixx
    Inertia(2)= Inertia(2) - Masses(k) * ( (x0(1,k)-COMx)*(x0(2,k)-COMy) )!Ixy
    Inertia(3)= Inertia(3) + Masses(k) * ( (x0(1,k)-COMx)**2.d0 + (x0(3,k)-COMz)**2.d0 )!Iyy
    Inertia(4)= Inertia(4) - Masses(k) * ( (x0(1,k)-COMx) * (x0(3,k)-COMz) )!Ixz
    Inertia(5)= Inertia(5) - Masses(k) * ( (x0(2,k)-COMy) * (x0(3,k)-COMz) )!Iyz
    Inertia(6)= Inertia(6) + Masses(k) * ( (x0(2,k)-COMy)**2.d0 + (x0(1,k)-COMx)**2.d0 )!Izz
   enddo
   call dspev('N','U',3,Inertia,eValues,eVectors,3,Work,Info)
   evalues=evalues*u_const*1.d-20
   rot_abc=h_const/(8.d0*pi*pi*c_const*eValues)
  WRITE(6,140)!MK
  WRITE(6,200) rot_ABC!MK

  endif
  !end rotational constants


  !        

  CALL spr_memdealloc(DBLE(SIZE(QQ)+SIZE(quad)+SIZE(x0))*8)
  DEALLOCATE(QQ,quad,x0,R0)

100 FORMAT(/,10x,'o VIBRATIONALLY AVERAGED STRUCTURE',//,10x,'  - Q0 -')
110 FORMAT(/,10x,'  - X0 -')
120 FORMAT(/,10x,'  - R0 -')
130 FORMAT(/,10x,'  - A0 -')
140 FORMAT(/,10x,'  - A B C -')
200 FORMAT(12x,3f12.6)
250 FORMAT(12x,2i2,f9.4)
300 FORMAT(12x,3i2,f9.4)

END SUBROUTINE Vci_ave

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!    ist :: an index for the current state

SUBROUTINE Vci_prpt(ist)

  USE Vci_mod
  USE Vav_mod

  IMPLICIT NONE

  INTEGER :: spr_memalloc
  INTEGER :: i,j,k,l,m,n,i1,i2

  INTEGER :: ist,jst,Lbl1(Nfree),Lbl2(Nfree)
  REAL(8) :: wgt
  REAL(8), DIMENSION(NP) :: P,tmp
  REAL(8), DIMENSION(maxCHO) :: D1

  jst=(ist-1)*nCI

  P=0.D+00
  DO i=1,nCI
     wgt=CIwfn(i+jst)*CIwfn(i+jst)
     CALL Vcor_getLabel(i,Lbl1)

     tmp=0.D+00
     DO j=1,Nfree
        CALL Vcor_getDwfn(j,nCHO(j),Lbl1(j),D1)

        DO k=1,NP
           DO l=1,nCHO(j)
              !P(k)=P(k) + wgt*P1(l,k,j)*D1(l)
              tmp(k)=tmp(k) + P1(l,k,j)*D1(l)
           END DO

        END DO
     END DO

     P=P+tmp*wgt

  END DO

  DO i1=1,Nci
     DO i2=1,i1-1
        wgt=2.D+00*CIwfn(i1+jst)*CIwfn(i2+jst)
        CALL Vcor_getLabel(i1,Lbl1)
        CALL Vcor_getLabel(i2,Lbl2)

        tmp=0.D+00
        DO j=1,Nfree
           DO k=1,j-1
              IF(Lbl1(k)/=Lbl2(k)) GOTO 1000
           END DO
           DO k=j+1,Nfree
              IF(Lbl1(k)/=Lbl2(k)) GOTO 1000
           END DO

           CALL Vcor_getXwfn(j,nCHO(j),Lbl1(j),Lbl2(j),D1)
           DO k=1,NP
              DO l=1,nCHO(j)
                 !P(k)=P(k) + wgt*P1(l,k,j)*D1(l)
                 tmp(k)=tmp(k) + P1(l,k,j)*D1(l)
              END DO

           END DO
1000       CONTINUE

        END DO

        P=P+tmp*wgt

     END DO
  END DO

  IF(MR2==1) GOTO 999

  !=========================
  !  MR=2,3,.. is not ready
  !=========================

999 CONTINUE
  CALL spr_GetIO(i,j)
  WRITE(j,100) 
  DO i=1,NP
     WRITE(j,200) i,P0(i),P(i),P0(i)+P(i)
  END DO

100 FORMAT(/,10x,'o VIBRATIONAL AVERAGE OVER PROPERTY SURFACES',/, &
       10x,'     [i]      P0             dP             <P>')
200 FORMAT(15x,i2,1x,3f15.6)

END SUBROUTINE Vci_prpt

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

SUBROUTINE Vci_print(io)

  USE Vci_mod

  IMPLICIT NONE

  INTEGER :: i,j,k 
  INTEGER :: io,spr_memalloc

  REAL(8) :: E0,E,realE0
  INTEGER, DIMENSION(:), ALLOCATABLE :: Lbl,Label
  REAL(8), DIMENSION(:), ALLOCATABLE :: Ci,Wi
  REAL(8), PARAMETER :: H2wvn=2.194746D+05, thsh=1.0D-03 !1.0D-03

  LOGICAL :: r_v,vav_Getr_v
  INTEGER :: NP,vav_GetNP

  r_v=vav_Getr_v()
  NP=vav_GetNP()

  i=spr_memalloc(1,DBLE(nCI*2)*8.D+00+DBLE(nCI)*4.D+00+DBLE(Nfree)*4.D+00)
  ALLOCATE(Ci(nCI),Wi(nCI),Lbl(nCI),Label(Nfree))

  ! Print Ground state
  Ci=CIwfn(1:nCI)
  DO i=1,nCI
     Wi(i)=Ci(i)*Ci(i)
  END DO
  CALL sort(nCI,Lbl,Wi)

  E0=CIene(1)*H2wvn
  CALL Vcor_getLabel(Lbl(1),Label)
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  !MK
  IF (SUM(label)==0)THEN
     WRITE(666,666)MR,nci,maxsum,nCUP,nGen,nLvl,maxexc          
     WRITE(666,'(f10.3)',advance='no')E0
666  FORMAT(/,"VCI(cm^-1)  for MR=",i1," nci=",i5," maxsum=",i2," nCUP=",i1," nGEN=", i2," nLvl=",i1," maxexc=",10i3)
     RealE0=E0
  ENDIF
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  !MK              
  IF(Nfree<31) THEN
     WRITE(io,100) 0,Label
     WRITE(io,110) E0
     WRITE(io,120)
     i=1
     DO WHILE(Wi(i) > thsh) 
        CALL Vcor_getLabel(Lbl(i),Label)
        WRITE(io,125) Ci(Lbl(i)),Wi(i),Label
        i=i+1
     END DO

  ELSE
     WRITE(io,100) 0,Label(1:30)
     WRITE(io,101) 0,Label(31:)
     WRITE(io,110) E0

     WRITE(io,120)
     i=1
     DO WHILE(Wi(i) > thsh) 
        CALL Vcor_getLabel(Lbl(i),Label)
        WRITE(io,125) Ci(Lbl(i)),Wi(i),Label(1:30)
        WRITE(io,126) Label(31:)
        i=i+1
     END DO
  ENDIF

  IF(r_v) CALL Vci_ave(1)
  IF(NP>0) CALL Vci_prpt(1)

  ! Print Excited state
  DO i=2,Nst
     Ci=CIwfn(nCI*(i-1)+1:nCI*i)
     DO j=1,nCI
        Wi(j)=Ci(j)*Ci(j)
     END DO
     CALL sort(nCI,Lbl,Wi)

     E=CIene(i)*H2wvn
     CALL Vcor_getLabel(Lbl(1),Label)
     !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
     !MK
     IF (SUM(label)==0)THEN
        WRITE(666,666)MR,nci,maxsum,nCUP,nGen,nLvl,maxexc          
        WRITE(666,*)"(Ghost ground state warning!)"             
        WRITE(666,'(f10.3)',advance='no')E
        realE0=E
     ENDIF
     IF (SUM(label)==1)THEN
        WRITE(666,'(f10.3)',advance='no')E-realE0

        IF(E0/=realE0)WRITE(777,'(f10.3)',advance='no')E-E0
     ENDIF
     !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
     !MK
     IF(Nfree<31) THEN
        WRITE(io,100) i-1,Label
        WRITE(io,115) E,E-E0
        WRITE(io,120)
        j=1
        DO WHILE(Wi(j) > thsh) 
           CALL Vcor_getLabel(Lbl(j),Label)
           WRITE(io,125) Ci(Lbl(j)),Wi(j),Label
           j=j+1
        END DO
     ELSE
        WRITE(io,100) i-1,Label(1:30)
        WRITE(io,101) i-1,Label(31:)
        WRITE(io,115) E,E-E0
        WRITE(io,120)
        j=1
        DO WHILE(Wi(j) > thsh) 
           CALL Vcor_getLabel(Lbl(j),Label)
           WRITE(io,125) Ci(Lbl(j)),Wi(j),Label(1:30)
           WRITE(io,126) Label(31:)
           j=j+1
        END DO
     ENDIF
     IF(r_v) CALL Vci_ave(i)
     IF(NP>0) CALL Vci_prpt(i)

  END DO

  CALL spr_memdealloc(DBLE(SIZE(Ci)+SIZE(Wi))*8.D+00+DBLE(SIZE(Lbl))*4.D+00)
  CALL spr_memdealloc(DBLE(SIZE(Label))*4.D+00)
  DEALLOCATE(Ci,Wi,Lbl,Label)

100 FORMAT(/,9x,'> STATE ',i3.3,': ',6(5i3,1x))!MK i1>i2
101 FORMAT(17x,': ',6(5i3,1x))!MK i1>i2
110 FORMAT(/,10x,'   E(VCI)   =',f12.2)
115 FORMAT(/,10x,'   E(VCI)   =',f12.2,/ &
       10x,'   E(VCI)-E0=',f12.2)
120 FORMAT(/,10x,'   COEFF.  WEIGHT      CONFIG.') 
125 FORMAT(13x,f6.3,2x,f6.3,6x,6(5i3,1x))!MK i1>i2
126 FORMAT(33x,6(5i3,1x))!MK i1>i2

END SUBROUTINE Vci_print

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!

SUBROUTINE Vci_print2(io)

  USE Vci_mod

  IMPLICIT NONE

  INTEGER :: i,j,k,kk
  INTEGER :: io,spr_memalloc

  REAL(8) :: E0,E,cc
  INTEGER, DIMENSION(:), ALLOCATABLE :: Lbl,Label
  REAL(8), DIMENSION(:), ALLOCATABLE :: Ci,Wi
  REAL(8), PARAMETER :: H2wvn=2.194746D+05, thsh=1.0D-03

  LOGICAL :: r_v,vav_Getr_v
  INTEGER :: NP,vav_GetNP

  SAVE :: E0

  r_v=vav_Getr_v()
  NP=vav_GetNP()

  i=spr_memalloc(-1,DBLE(nCI*2)*8.D+00+DBLE(nCI)*4.D+00+DBLE(Nfree)*4.D+00)
  ALLOCATE(Ci(nCI),Wi(nCI),Lbl(nCI),Label(Nfree))

  !dbg Do i=1,nCI
  !dbg    write(6,'(i3,10f8.3)') i,CIwfn((i-1)*nCI+1:(i-1)*nCI+10)
  !dbg End do

  ! Print Reference state
  j=1; kk=1; cc=0.D+00
  DO i=1,nCI
     IF(cc<ABS(CIwfn(j))) THEN 
        kk=i
        cc=ABS(CIwfn(j))
     ENDIF
     j=j+nCI
  END DO
  !dbg write(6,*) kk
  Ci=CIwfn((kk-1)*nCI+1:kk*nCI)
  DO i=1,nCI
     Wi(i)=Ci(i)*Ci(i)
  END DO
  CALL sort(nCI,Lbl,Wi)

  CALL Vcor_getLabel(Lbl(1),Label)
  IF(Nfree<31) THEN
     WRITE(io,100) Label
  ELSE
     WRITE(io,100) Label(1:30)
     WRITE(io,102) Label(31:)
  ENDIF
  IF(.NOT. zp) THEN
     E=CIene(kk)*H2wvn
     WRITE(io,115) E,E-E0
     !mk
     IF (SUM(label)==1) WRITE(666,'(f10.3)',advance='no')E-E0
     !MK
     
  ELSE
     E0=CIene(kk)*H2wvn
     WRITE(io,110) E0
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  !MK
  IF (SUM(label)==0)THEN
     WRITE(666,666)MR,nci,maxsum,nCUP,nGen,nLvl,maxexc          
     WRITE(666,'(f10.3)',advance='no')E0
666  FORMAT(/,"SS.VCI(cm^-1)  for MR=",i1," nci=",i5," maxsum=",i2," nCUP=",i1," nGEN=", i2," nLvl=",i1," maxexc=",10i3)
  ENDIF
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  !MK              
  ENDIF
  WRITE(io,120)
  i=1
  DO WHILE(Wi(i) > thsh) 
     CALL Vcor_getLabel(Lbl(i),Label)
     IF(Nfree<31) THEN 
        WRITE(io,125) Ci(Lbl(i)),Wi(i),Label
     ELSE
        WRITE(io,125) Ci(Lbl(i)),Wi(i),Label(1:30)
        WRITE(io,126) Label(31:)
     ENDIF
     i=i+1
  END DO
  IF(r_v) CALL Vci_ave(kk)
  IF(NP>0) CALL Vci_prpt(kk)

  IF(npCnf==0) RETURN

  ! Print P-states
  WRITE(io,101)
  DO i=1,npCnf
     k=1+i; kk=1; cc=0.D+00
     DO j=1,nCI
        IF(cc<ABS(CIwfn(k))) THEN 
           kk=j
           cc=ABS(CIwfn(k))
        ENDIF
        k=k+nCI
     END DO
     !dbg write(6,*) kk
     Ci=CIwfn(nCI*(kk-1)+1:nCI*kk)
     DO j=1,nCI
        Wi(j)=Ci(j)*Ci(j)
     END DO
     CALL sort(nCI,Lbl,Wi)

     CALL Vcor_getLabel(Lbl(1),Label)
     E=CIene(kk)*H2wvn
     !mk
     IF (SUM(label)==1) WRITE(666,'(f10.3)',advance='no')E-E0
     !MK
     IF(Nfree<31) THEN
        WRITE(io,103) i,Label
        WRITE(io,115) E,E-E0
     ELSE
        WRITE(io,103) i,Label(1:30)
        WRITE(io,102) Label(31:)
        WRITE(io,115) E,E-E0
     ENDIF

     WRITE(io,120)
     j=1
     DO WHILE(Wi(j) > thsh)
        CALL Vcor_getLabel(Lbl(j),Label)
        IF(Nfree<31) THEN
           WRITE(io,125) Ci(Lbl(j)),Wi(j),Label
        ELSE
           WRITE(io,125) Ci(Lbl(j)),Wi(j),Label(1:30)
           WRITE(io,126) Label(31:)
        ENDIF
        j=j+1
     END DO
     IF(r_v) CALL Vci_ave(i)
     IF(NP>0) CALL Vci_prpt(i)

  END DO

  CALL spr_memdealloc(DBLE(SIZE(Ci)+SIZE(Wi))*8.D+00+DBLE(SIZE(Lbl))*4.D+00)
  CALL spr_memdealloc(DBLE(SIZE(Label))*4.D+00)
  DEALLOCATE(Ci,Wi,Lbl,Label)

100 FORMAT(/,9x,'o REFERENCE STATE : ',6(5i3,1x))!MK i1>i2
101 FORMAT(/,9x,'*****  RESONANT  STATE  *****')
102 FORMAT(27x,': ',6(5i3,1x))
103 FORMAT(/,9x,'ooo STATE ',i2.2,8x,' : ',6(5i3,1x))!MK i1>i2
110 FORMAT(/,10x,'   E(VCI)   =',f12.2)
115 FORMAT(/,10x,'   E(VCI)   =',f12.2,/ &
       10x,'   E(VCI)-E0=',f12.2)
120 FORMAT(/,10x,'   COEFF.  WEIGHT      CONFIG.') 
125 FORMAT(13x,f6.3,2x,f6.3,6x,6(5i3,1x))!MK i1>i2
126 FORMAT(33x,6(5i3,1x))!MK i1>i2

END SUBROUTINE Vci_print2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

SUBROUTINE Vci_Dump()

  USE Vci_mod

  IMPLICIT NONE

  INTEGER :: wwf,i,j
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Label
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: C1wfn,E1wfn
  CHARACTER :: wftyp*4

  ! ------------------------------------------------------------
  ! Dump VCI wavefunction
  CALL file_indicator(30,wwf)
  OPEN(wwf,file='vci-w.wfn',status='unknown',form='UNFORMATTED')
  WRITE(wwf) Nfree
  WRITE(wwf) nCHO         ! nCHO(Nfree)
  WRITE(wwf) omegaf       ! omegaf(Nfree)

  IF(lvscf) THEN
     ALLOCATE(C1wfn(maxCHO*maxCHO,Nfree),E1wfn(maxCHO,Nfree))
     C1wfn=0; E1wfn=0
     DO i=1,Nfree
        CALL Modal_getCwfn(i,C1wfn(:,i))
        CALL Modal_getEne(i,E1wfn(:,i))
     END DO
     WRITE(wwf) C1wfn     ! C1wfn(int*int,Nfree)
     WRITE(wwf) E1wfn     ! E1wfn(int,Nfree)
     DEALLOCATE(C1wfn,E1wfn)
  ENDIF

  ALLOCATE(Label(Nfree,Nci))
  DO i=1,Nci
     CALL Vcor_getLabel(i,Label(:,i))
  END DO

  wftyp='VCI'
  WRITE(wwf) wftyp
  WRITE(wwf) Nstate,Nci
  WRITE(wwf) Label        ! Label(Nfree,Nci)
  WRITE(wwf) CIwfn        ! CIwfn(Nci*Nstate)
  WRITE(wwf) CIene        ! CIene(Nci)
  CLOSE(wwf)

  !Open(wwf,file='vci-w.wfn',status='unknown',form='FORMATTED')
  !Write(wwf,'(i4)') Nfree
  !Write(wwf,'(3i4)') nho
  !Write(wwf,'(3f8.2)') omegaf
  !Do i=1,Nfree
  !    Write(6,'(11f8.4)') C1wfn(:,i)
  !    Write(6,*)
  !End do
  !Write(wwf,'(2i4)') Nstate,Nci
  !Write(wwf,'(3i4)') Label
  !Do i=1,Nstate
  !   Write(wwf,'(10f12.6)') CIwfn(Nci*(i-1)+1:Nci*i)
  !   Write(wwf,*)
  !End do

  DEALLOCATE(Label)

END SUBROUTINE Vci_Dump

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
