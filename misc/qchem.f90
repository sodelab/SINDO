!
!         Read(readfile,500) index,E
!      if (E>minLimit && E<maxLimit) then
!            E=get_qc_energy(filenumber)
!      endif
!      Write(ifl,500) indx,E
!      500 Format(i8,f20.10)
!      function get_qc_energy(filenumber)
!      double precision :: get_qc_energy
!      integer :: filenumber
!      character(25) :: qc_filename,key_string='Convergence criterion met'
!
!      qc_filename=(filenumber)
!      do while (position==0)
!      read(21,'(A)')line
!      position=index(line,key_string)
!      read(line())
!      end do
