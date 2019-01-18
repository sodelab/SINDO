!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
  SUBROUTINE Launch
!

   Integer, Parameter :: Iout=6
   Integer:: IsHost, Hostnm

   Character(len=24) Date
   Character(len=32) Name
!

     Ishost=Hostnm(Name)
     Call Fdate(Date)

     Write(Iout,100) Date
     Write(Iout,200) Name
!
   return
!
   100 Format(/,5x,'JOB STARTED AT: ',a24)
   200 Format(5x,'RUNNING ON ',a32)
!
  END SUBROUTINE
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE timer(int)
!
  Integer, Parameter :: Iout=6
  Real(4) Delta,dtime,Elasped,etime,tarry(2)
!
     if(int==1) then

       Delta=dtime(tarry)

       Write(Iout,1000) tarry(1)
     elseif(int==2) then

       Elasped=etime(tarry)

       Write(Iout,1001) tarry(1)
     End if
!
   1000 Format(/,3x,'(CLOCK) ---->  LAST STEP ',F8.2,' SECS  <----',/)
   1001 Format(/,3x,'(CLOCK) ---->  TOTAL ELAPSED TIME', F14.2,' SECS  <----',/)
!
  END SUBROUTINE
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
  SUBROUTINE myFlsh(Lunit)
!
  Integer :: Lunit

      Call Flush(Lunit)

  END SUBROUTINE
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
