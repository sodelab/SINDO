!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

  Subroutine VERSION()

      Integer, parameter :: Iout=6

      Write(Iout,*) '-------------------------------------------------------'
      Write(Iout,*) '                                                     '
      Write(Iout,*) '   ***   WELCOME TO SINDO PROGRAM                    '
      Write(Iout,*) '    ***                (  VER. 2.0     2007/12/10 )  '
      Write(Iout,*) '     ***                                             '
      Write(Iout,*) '      ***                  KIYOSHI YAGI              '
      Write(Iout,*) '       ***                 YAGI@QCL.T.U-TOKYO.AC.JP  '
      Write(Iout,*) '                                                     '
      Write(Iout,*) '-------------------------------------------------------'
  END

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine diag(n, m, H, C, E)

      Implicit None

         Integer :: n,m,m1,in,io
         Real(8), dimension(n*(n+1)/2) :: H
         Real(8), dimension(n,m) :: C
         Real(8), dimension(n) :: E

         Integer :: i,j,k,spr_memalloc

         Character :: jobz,range,uplo
         Real(8)   :: vl,vu,abstol
         Integer   :: il,iu,ldz,info
         Integer, dimension(:), allocatable :: ifail,iwork
         Real(8), dimension(:), allocatable :: work


            !write(6,*) n,m
            !Do i=1,n
            !   write(6,'(11f12.6)') (H(j),j=i*(i-1)/2+1,i*(i+1)/2)
            !End do

            i=spr_memalloc(-1,dble(10*n)*8.D+00) 
            i=spr_memalloc(-1,dble(n+10*n)*4.D+00) 
            Allocate(work(10*n),ifail(n),iwork(10*n))

            jobz='V'
            uplo='U' 
            vl=0.D+00
            vu=0.D+00
            il=0
            iu=0
            if(n==m) then
               range='A'
            else
               range='I'; il=1; iu=m
            endif

            abstol=0.D+00
            ldz=n

            m1=0
            ifail=0; info=0
            Call dspevx(jobz,range,uplo,n,H,vl,vu,il,iu,abstol,m1,E, &
                        C,ldz,work,iwork,ifail,info)

            !write(6,*) info,ifail
            Call spr_memdealloc(dble(size(work))*8.D+00+ &
                            dble(size(ifail)+size(iwork))*4.D+00)
            Deallocate(work,ifail,iwork)

            !write(6,'(11f12.6)') E

            if(info==0) return

            Call spr_Getio(in,io)
            if(info<0) then 
               write(io,'(''ERROR IN '',i3,''TH PARAMETER'')') info
            else
               write(io,'(3x,i3,''EIGENVECTORS FAILED TO CONVERGE'')') info
            endif

      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
      Subroutine diag2(n, m, H, C, E)

      Implicit None

         Integer :: n,m,in,io
         Real(8), dimension(n,n) :: H
         Real(8), dimension(n,m) :: C
         Real(8), dimension(n) :: E

         Integer :: i,j,k,spr_memalloc
         Real(8), dimension(n*(n+1)/2) :: H0

            H0=0.D+00; k=1
            Do i=1,n
               Do j=1,i
                  H0(k)=H(j,i)
                  k=k+1
               End do
            End do

            !write(6,'(11f12.6)') H
            !write(6,*)
            !Do i=1,n
            !   write(6,'(11f12.6)') (H0(j),j=i*(i-1)/2+1,i*(i+1)/2)
            !End do
            Call diag(n,m,H0,C,E)

            return

      End Subroutine
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
!
!     Huckeler ver. 0.2
      Subroutine Huckeler(ncol,mcol,a,d,v)
      implicit double precision(a-h,o-z)
!      parameter (mcol=3)
      dimension a(mcol,mcol),d(mcol),v(mcol,mcol), &
                tv(mcol)
!     Initialization
      nrot=0
      do i=1,ncol
         d(i)=0.
         do j=1,ncol
            v(i,j)=0.
         enddo
      enddo
!     Calling jacobi
      call jacobi(a,ncol,mcol,d,v,nrot)
!     Printing out results
      do i=1,ncol
         if (v(1,i).lt.0.) then
            do j=1,ncol
               v(j,i)=-v(j,i)
            enddo
         endif
      enddo
      do k=1,ncol-1
         do i=k+1,ncol
            if (d(i).gt.d(k)) then
               td=d(k)
               d(k)=d(i)
               d(i)=td
               do j=1,ncol
                  tv(j)=v(j,k)
                  v(j,k)=v(j,i)
                  v(j,i)=tv(j)
               enddo
            endif
         enddo
      enddo
!      do i=1,ncol
!         write(6,99) d(i)
!         write(6,98) (v(j,i),j=1,ncol)
!         write(6,*)
!      enddo
!99   format(f8.4)
!98   format(10f8.4)
      end

      SUBROUTINE jacobi(adum,n,np,d,v,nrot)
      implicit double precision(a-h,o-z)
      PARAMETER (NMAX=50000)
      dimension a(np,np),adum(np,np),d(np),v(np,np)
      dimension b(NMAX),z(NMAX)
!
!      do i=1,n
!         do j=1,n
!            write(6,*) 'a(',i,',',j,')=',a(i,j)
!         enddo
!      enddo
!
      do j=1,n
         do i=1,n
            a(i,j)=adum(i,j)
         end do
      end do
!
      do ip=1,n
         do iq=1,n
            v(ip,iq)=0.0d0
         end do
         v(ip,ip)=1.0d0
      end do
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0.0d0
      end do
      nrot=0
      do i=1,50
         sm=0.0d0
         do ip=1,n-1
            do iq=ip+1,n
               sm=sm+abs(a(ip,iq))
            end do
         end do
        if (sm.eq.0.0d0) return
        if (i.lt.4) then
           tresh=0.2d0*sm/n**2
        else
           tresh=0.0d0
        end if
        do ip=1,n-1
           do iq=ip+1,n
              g=100.0d0*abs(a(ip,iq))
              if ((i.gt.4).and.(abs(d(ip))+ &
                 g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq)))) then
                 a(ip,iq)=0.0d0
              else if (abs(a(ip,iq)).gt.tresh) then
                 h=d(iq)-d(ip)
                 if (abs(h)+g.eq.abs(h)) then
                    t=a(ip,iq)/h
                 else
                    theta=0.5d0*h/a(ip,iq)
                    t=1./(abs(theta)+sqrt(1.+theta**2))
                    if (theta.lt.0.) t=-t
                 end if
                 c=1.0d0/sqrt(1+t**2)
                 s=t*c
                 tau=s/(1.0d0+c)
                 h=t*a(ip,iq)
                 z(ip)=z(ip)-h
                 z(iq)=z(iq)+h
                 d(ip)=d(ip)-h
                 d(iq)=d(iq)+h
                 a(ip,iq)=0.0d0
                 do j=1,ip-1
                    g=a(j,ip)
                    h=a(j,iq)
                    a(j,ip)=g-s*(h+g*tau)
                    a(j,iq)=h+s*(g-h*tau)
                 end do
                 do j=ip+1,iq-1
                    g=a(ip,j)
                    h=a(j,iq)
                    a(ip,j)=g-s*(h+g*tau)
                    a(j,iq)=h+s*(g-h*tau)
                 end do
                 do j=iq+1,n
                    g=a(ip,j)
                    h=a(iq,j)
                    a(ip,j)=g-s*(h+g*tau)
                    a(iq,j)=h+s*(g-h*tau)
                 end do
                 do j=1,n
                    g=v(j,ip)
                    h=v(j,iq)
                    v(j,ip)=g-s*(h+g*tau)
                    v(j,iq)=h+s*(g-h*tau)
                 end do
                 nrot=nrot+1
              end if
           end do
        end do
        do ip=1,n
           b(ip)=b(ip)+z(ip)
           d(ip)=b(ip)
           z(ip)=0.0d0
        end do
      end do
      Stop 'too many iterations in jacobi'
      return
      END

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine lw2up(buf)
! **************************************************
! *  convert lower-case letters to uppercase ones  *
! *                  (a-z)               (A-Z)     *
! *    coded by Naoki MAEDA on 2000/10/26          *
! **************************************************
      character buf*(*)
      n2=len(buf)
      do 200 i=1,n2
      if(ichar(buf(i:i)).ge.97.and.ichar(buf(i:i)).le.122) &
            buf(i:i)=char(ichar(buf(i:i))-32)
  200 continue
      return
      end Subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      subroutine up2lw(buf)
! **************************************************
! *  convert uppercase letters to lower-case ones  *
! *                 (A-Z)               (a-z)      *
! *    coded by Naoki MAEDA on 2000/10/26          *
! **************************************************
      character buf*(*)
      n2=len(buf)
      do 200 i=1,n2
      if(ichar(buf(i:i)).ge.65.and.ichar(buf(i:i)).le.90) &
            buf(i:i)=char(ichar(buf(i:i))+32)
  200 continue
      return
      end

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine file_indicator(in,io)

      Implicit None

         Integer :: in,io
         Logical :: op

         io=in
         Do while(.true.)
            Inquire(io,opened=op)
            if(op) then 
              io=io+1
            else
              exit
            endif
         End do

      End subroutine

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80

      Subroutine sort(N,L,C)

      Implicit None

         Integer :: N,L(N)
         Real(8) :: C(N)

         Integer :: i,j(1),k,itmp  !MK ???j(1)???
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

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
     Function B2A()

     Double Precision :: B2A

!       Bohr -> Angs
        B2A = 0.52917724924d+00

     End Function B2A
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
     Function Elmass()

     Double Precision :: Elmass

!       emu -> amu
        Elmass = 1822.88853D+00

     End Function Elmass
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
!
     Function H2wvn()

     Double Precision :: H2wvn

!       hartree -> cm-1
        H2wvn = 2.194746E+05

     End Function H2wvn
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----80
