!
!      Program test
!
!      Integer, parameter :: Nat=3
!      Real(8) :: x(3,Nat),ms(Nat),g(3,Nat),h(Nat*3,Nat*3),w(Nat*3-6), &
!                 L(Nat*3,Nat*3-6)
!      Character :: lbl(Nat)*2
!
!         Call rdftb2(Nat,x,ms,lbl,w,L)
!
!         Write(6,'(3f12.6)') x
!         Write(6,*)
!         Do i=1,Nat*3-6
!            Write(6,'(3f12.6)') L(:,i)
!            Write(6,*)
!         End do
!
!         Call nma_Construct(3,3,ms,x,L)
!         Call nma_ChkCL(6)
!
!      End

      Subroutine rdftb(Nat,e,g,h)

      Implicit None

      Integer :: Nat,Nat3
      Real(8) :: x(3,Nat),e,g(3,Nat),h(Nat*3,Nat*3)

      Integer :: i,j,k,l
      Character :: lbl(Nat)*2
      Real(8) :: ms(Nat)

         Nat3=Nat*3

         Open(10,file='sccoutput',status='old')

         Read(10,*)
         Read(10,*) e

         Read(10,*)
         Read(10,*)
         Do i=1,Nat
            Read(10,*) j,lbl(i),ms(i),x(:,i)
         End do

         Read(10,*)
         Read(10,*)
         Do i=1,Nat
            Read(10,*) j,lbl(i),g(:,i)
         End do

         Read(10,*)
         Do i=1,Nat3
            Read(10,*) h(:,i)
         End do

         Close(10)

      End

      Subroutine rdftb2(Nat,x,ms,lbl,w,L)

      Implicit None
      Integer :: Nat,Nat3,Nfree

      Real(8) :: x(3,Nat),ms(Nat),e,g(3,Nat),h(Nat*3,Nat*3),w(Nat*3-6), &
                 L(Nat*3,Nat*3-6)
      Character :: lbl(Nat)*2

      Integer :: i,j,k,ll
      Real(8) :: L2(Nat*3,Nat*3),tmp

         Nat3=Nat*3
         Nfree=Nat3-6

         Open(10,file='sccoutput',status='old')

         Read(10,*)
         Read(10,*) e

         Read(10,*)
         Read(10,*)
         Do i=1,Nat
            Read(10,*) j,lbl(i),ms(i),x(:,i)
         End do

         Read(10,*)
         Read(10,*)
         Do i=1,Nat
            Read(10,*) j,lbl(i),g(:,i)
         End do

         Read(10,*)
         Do i=1,Nat3
            Read(10,*) h(:,i)
         End do

         Read(10,*)
         Do i=1,6
            Read(10,*) 
         End do
         Do i=1,Nfree
            Read(10,*) j,w(i)
         End do

         Read(10,*)
         Do i=1,Nat3
            Read(10,*) L2(:,i)
         End do
         Close(10)

         Do i=1,Nfree
            L(:,i)=L2(6+i,:)
            ll=1
            Do j=1,Nat
               tmp=SQRT(ms(j))
               Do k=1,3
                  L(ll,i)=L(ll,i)*tmp
                  ll=ll+1
               End do
            End do
         End do

         Do i=1,Nfree
            tmp=0
            Do j=1,Nat3
               tmp = tmp + L(j,i)**2
            End do
            tmp=SQRT(tmp)
            Do j=1,Nat3
               L(j,i)=L(j,i)/tmp
            End do
         End do

      End
