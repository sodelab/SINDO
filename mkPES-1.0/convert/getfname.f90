!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine get_fname1(imod,fname)

      Implicit None

         Integer :: i,imod
         Character(*) :: fname

            if(imod<10) then
               write(fname,'(''q'',i1,''.pot'')') imod
            elseif(imod<100) then
               write(fname,'(''q'',i2,''.pot'')') imod
            else
               write(fname,'(''q'',i3,''.pot'')') imod
            endif

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      ! imod>jmod
      Subroutine get_fname2(imod,jmod,fname)

      Implicit None

         Integer :: imod,jmod,i,j
         Character(*) :: fname

         i=0
         Call get_header(imod,i,fname(i+1:))
         Call get_header(jmod,i,fname(i+1:))
         fname(i+1:)='.pot'

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      ! i>j>k
      Subroutine get_fname3(imod,jmod,kmod,fname)

      Implicit None

         Integer :: i,j,k,imod,jmod,kmod
         Character(*) :: fname

         i=0
         Call get_header(imod,i,fname(i+1:))
         Call get_header(jmod,i,fname(i+1:))
         Call get_header(kmod,i,fname(i+1:))
         fname(i+1:)='.pot'

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine get_header(imod,ipt,fp)

      Implicit None

         Integer :: imod,ipt
         Character(*) :: fp
         Character :: cm1*1,cm2*2,cm3*3

         if(imod<10) then
            write(cm1,'(i1)') imod
            fp='q'//cm1
            ipt=ipt+2
         elseif(imod<100) then
            write(cm2,'(i2)') imod
            fp='q'//cm2
            ipt=ipt+3
         else
            write(cm3,'(i3)') imod
            fp='q'//cm3
            ipt=ipt+4
         endif

      End Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
