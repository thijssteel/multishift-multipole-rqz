      module z_rqz_aed
         use common
         use z_rqz_small
         implicit none

      contains

      subroutine aed_rqz_complex(wantS,wantQ,wantZ,n,ilo,ihi,ws,we,A,
     $   ldA,B,ldB,Q,ldQ,Z,ldZ,ns,nds,nde,alpha,beta,Qc,ldQc,Zc,ldZc,
     $   work,lwork)
      implicit none
*     Arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,ldQc,ldZc,lwork
      integer,intent(inout) :: ws,we

      double complex,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,
     $   *),alpha(*),beta(*)
      integer,intent(out) :: ns,nds,nde
      double complex :: Qc(ldQc,*),Zc(ldZc,*),work(*)

*     Local Scalars
      integer :: kwetop,kwebot,kwstop,kwsbot,istopm,istartm,k,k2,
     $   ztgexc_info,ifst,ilst,lworkreq,i,ishift,istartb,istopb,npoles
      double complex :: s_top1,s_top2,s_bot1,s_bot2,temp
      double precision :: smlnum,ulp,safmin,safmax

*     External Functions
      double precision,external :: dlamch

*     External Subroutines
      external :: zgemm,zlacpy,ztgexc

      lworkreq = max(ws,we+1)*n
      if (lwork .eq.-1) then
*     workspace query, quick return
         work(1) = lworkreq
         return
      else if (lwork .lt. lworkreq) then
         write (*,*) "workspace provided to aed is too small",lworkreq,
     $      lwork
         return
      end if

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      CALL dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

*     Top AED part

*     Set up deflation windows
      if( (we .ge. ihi-ilo+1) .or. (we .ge. ihi-ilo+1) ) then
         ws = ihi-ilo+1
         we = 0
         kwetop = ihi
         kwebot = ihi
         kwstop = ilo
         kwsbot = ilo+ws-1
      else
         kwetop = ihi-we+1
         ws = min(ws,kwetop-ilo+1)
         kwstop = ilo
         kwsbot = ilo+ws-1
      end if
      if (kwetop .eq. ilo) then
         s_top1 = czero
         s_top2 = czero
      else
         s_top1 = A(kwetop,kwetop-1)
         s_top2 = B(kwetop,kwetop-1)
      end if
      if (kwsbot .eq. ihi) then
         s_bot1 = czero
         s_bot2 = czero
      else
         s_bot1 = A(kwsbot+1,kwsbot)
         s_bot2 = B(kwsbot+1,kwsbot)
      end if

      if (ihi .eq. kwstop) then
*        1 by 1 deflation window, just try a regular deflation
         alpha(kwstop) = A(kwstop,kwstop)
         beta(kwstop) = B(kwstop,kwstop)
         nds = 0
      if (abs(s_bot1) .le. max(smlnum,ulp*abs(A(kwstop,
     $   kwstop))) .and.abs(s_bot2) .le. max(smlnum,ulp*abs(B(kwstop,
     $   kwstop)))) then
            nds = 1
            if (kwstop .gt. ilo) then
               A(kwstop,kwstop-1) = czero
               B(kwstop,kwstop-1) = czero
            end if
         end if
      end if

      call zlaset('Full',ws,ws,czero,cone,Qc,ldQc)
      call zlaset('Full',ws,ws,czero,cone,Zc,ldZc)

*     Transform window to real schur form
      call z_rqz(.true.,.true.,.true.,ws,1,ws,A(kwstop,kwstop),ldA,
     $   B(kwstop,kwstop),ldB,Qc,ldQc,Zc,ldZc,alpha(kwstop),
     $   beta(kwstop))


*     Deflation detection loop
      if (kwsbot .eq. ihi) then
         kwstop = kwsbot+1
      else
         k = 1
         k2 = ws
         do while (k .le. ws)

*           Try to deflate single eigenvalue
            if ((abs(s_bot1*Zc(ws,kwstop-ilo+1))).le.ulp*(abs(A(kwstop,
     $         kwstop))+abs(A(kwstop+1,kwstop+1))).and.(abs(s_bot2*
     $         Zc(ws,kwstop-ilo+1))).le.ulp*(abs(B(kwstop,
     $         kwstop))+abs(B(kwstop+1,kwstop+1)))) then
*              Deflatable
               kwstop = kwstop+1
            else
*              Not deflatable, move out of the way
               ifst = kwstop-ilo+1
               ilst = k2
               call ztgexc(.true.,.true.,ws,A(ilo,ilo),ldA,B(ilo,ilo),
     $            ldB,Qc,ldQc,Zc,ldZc,ifst,ilst,ztgexc_info)
               if (ztgexc_info .ne. 0) then
                  write (*,*) "swap warning",ztgexc_info
               end if
               k2 = k2-1
            end if

            k = k+1
         end do
      end if

*     Store eigenvalues
      nds = kwstop-ilo
      do k = kwstop,kwsbot
         alpha(k) = A(k,k)
         beta(k) = B(k,k)
      end do

      if (kwsbot .ne. ihi) then
*        Reorder the submatrix to allow for single swap with the spike afterwards
         do k = kwstop,kwsbot
            temp = A(kwsbot,k)
            A(kwstop+2:kwsbot,k) = A(kwstop+1:kwsbot-1,k)
            A(kwstop+1,k) = temp

            temp = B(kwsbot,k)
            B(kwstop+2:kwsbot,k) = B(kwstop+1:kwsbot-1,k)
            B(kwstop+1,k) = temp

         end do
         do k = kwsbot-1,kwstop+1,-1
            do k2 = 1,ws
               temp = Qc(k2,k-ilo+1)
               Qc(k2,k-ilo+1) = Qc(k2,k+1-ilo+1)
               Qc(k2,k+1-ilo+1) = temp
            end do
         end do

      end if

*     Apply Qc and Zc to rest of the matrix
      if (wantS) then
         istartm = 1
         istopm = n
      else
         istartm = ilo
         istopm = ihi
      end if

      if (istopm-kwsbot > 0) then
         call zgemm('C','N',ws,istopm-kwsbot,ws,cone,Qc,ldQc,A(ilo,
     $      kwsbot+1),ldA,czero,work,ws)
         call zlacpy('ALL',ws,istopm-kwsbot,work,ws,A(ilo,kwsbot+1),ldA)
         call zgemm('C','N',ws,istopm-kwsbot,ws,cone,Qc,ldQc,B(ilo,
     $      kwsbot+1),ldB,czero,work,ws)
         call zlacpy('ALL',ws,istopm-kwsbot,work,ws,B(ilo,kwsbot+1),ldB)
      end if
      if (wantQ) then
         call zgemm('N','N',n,ws,ws,cone,Q(1,ilo),ldQ,Qc,ldQc,czero,
     $      work,n)
         call zlacpy('ALL',n,ws,work,n,Q(1,ilo),ldQ)
      end if

      if (ilo-1-istartm+1 > 0) then
         call zgemm('N','N',ilo-istartm,ws,ws,cone,A(istartm,ilo),ldA,
     $      Zc,ldZc,czero,work,ilo-istartm)
        call zlacpy('ALL',ilo-istartm,ws,work,ilo-istartm,A(istartm,
     $     ilo),ldA)
         call zgemm('N','N',ilo-istartm,ws,ws,cone,B(istartm,ilo),ldB,
     $      Zc,ldZc,czero,work,ilo-istartm)
        call zlacpy('ALL',ilo-istartm,ws,work,ilo-istartm,B(istartm,
     $     ilo),ldB)
      end if
      if (wantZ) then
         call zgemm('N','N',n,ws,ws,cone,Z(1,ilo),ldZ,Zc,ldZc,czero,
     $      work,n)
         call zlacpy('ALL',n,ws,work,n,Z(1,ilo),ldZ)
      end if

      if (ws .ne. ihi-ilo+1) then
*        Calculate the spike
         A(kwsbot+1,ilo+nds:kwsbot) = A(kwsbot+1,kwsbot)*Zc(ws,1+nds:ws)
         B(kwsbot+1,ilo+nds:kwsbot) = B(kwsbot+1,kwsbot)*Zc(ws,1+nds:ws)

*        Swap the spike up to get a Hessenberg-Hessenberg pencil
         do k = istartm,istopm
            temp = A(kwsbot+1,k)
            A(kwsbot+1,k) = A(kwstop+1,k)
            A(kwstop+1,k) = temp
         end do
         do k = istartm,istopm
            temp = B(kwsbot+1,k)
            B(kwsbot+1,k) = B(kwstop+1,k)
            B(kwstop+1,k) = temp
         end do
         do k = 1,n
            temp = Q(k,kwsbot+1)
            Q(k,kwsbot+1) = Q(k,kwstop+1)
            Q(k,kwstop+1) = temp
         end do
      end if

*     Bottom AED part

      if(we .eq. 0) return

      if (ihi .eq. kwetop) then
*     1 by 1 deflation window, just try a regular deflation
         alpha(kwetop) = A(kwetop,kwetop)
         beta(kwetop) = B(kwetop,kwetop)
         ns = 1
         nde = 0
      if (abs(s_top1) .le. max(smlnum,ulp*abs(A(kwetop,
     $   kwetop))) .and.abs(s_top2) .le. max(smlnum,ulp*abs(B(kwetop,
     $   kwetop)))) then
            ns = 0
            nde = 1
            if (kwetop .gt. ilo) then
               A(kwetop,kwetop-1) = czero
               B(kwetop,kwetop-1) = czero
            end if
         end if
      end if

      call zlaset('Full',we,we,czero,cone,Qc,ldQc)
      call zlaset('Full',we+1,we+1,czero,cone,Zc,ldZc)

*     Transform window to real schur form
      call z_rqz(.true.,.true.,.true.,we,1,we,A(kwetop,kwetop),ldA,
     $   B(kwetop,kwetop),ldB,Qc,ldQc,Zc(2,2),ldZc,alpha(kwetop),
     $   beta(kwetop))

      kwebot = ihi

*     Deflation detection loop
      if (kwetop .eq. ilo) then
         kwebot = kwetop-1
      else
         kwebot = ihi
         k = 1
         k2 = 1
         do while (k .le. we)

*           Try to deflate single eigenvalue
            if ((abs(s_top1*Qc(1,kwebot-kwetop+1))) .le.ulp*
     $         (abs(A(kwebot,kwebot))+abs(A(kwebot-1,
     $         kwebot-1))).and.(abs(s_top2*Qc(1,kwebot-kwetop+
     $         1))) .le.ulp*(abs(B(kwebot,kwebot))+abs(B(kwebot-1,
     $         kwebot-1)))) then
*              Deflatable
               kwebot = kwebot-1
            else
*              Not deflatable, move out of the way
               ifst = kwebot-kwetop+1
               ilst = k2
               call ztgexc(.true.,.true.,we,A(kwetop,kwetop),ldA,
     $            B(kwetop,kwetop),ldB,Qc,ldQc,Zc(2,2),ldZc,ifst,ilst,
     $            ztgexc_info)
               if (ztgexc_info .ne. 0) then
                  write (*,*) "swap warning",ztgexc_info
            write (*,*) "after swap: ",A(kwebot,kwebot)/B(kwebot,kwebot)
               end if
               k2 = k2+1
            end if

            k = k+1
         end do
      end if

*     Store eigenvalues
      nde = ihi-kwebot
      ns = we-nde
      k = kwetop
      do k = kwetop,ihi
         alpha(k) = A(k,k)
         beta(k) = B(k,k)
      end do

      if (kwetop .ne. ilo) then
*        Calculate the spike
         A(kwetop:kwebot,kwetop-1) = A(kwetop,kwetop-1)*dconjg(Qc(1,
     $      1:we-nde))
         B(kwetop:kwebot,kwetop-1) = B(kwetop,kwetop-1)*dconjg(Qc(1,
     $      1:we-nde))

*        Reorder the subblock to get back Hessenberg-Hesserberg pencil
         do k = kwetop-1,kwebot-2
            do k2 = kwetop,kwebot
               temp = A(k2,k)
               A(k2,k) = A(k2,k+1)
               A(k2,k+1) = temp
            end do
            do k2 = kwetop,kwebot
               temp = B(k2,k)
               B(k2,k) = B(k2,k+1)
               B(k2,k+1) = temp
            end do
            do k2 = 1,we+1
               temp = Zc(k2,k-(kwetop-1)+1)
               Zc(k2,k-(kwetop-1)+1) = Zc(k2,k+1-(kwetop-1)+1)
               Zc(k2,k+1-(kwetop-1)+1) = temp
            end do
         end do

*        Introduce undeflated eigenvalues from top AED as poles in the
*        top of the bottom AED window. Other poles in the window
*        are set to infinity.

*        istartb points to the first row we will be updating
         istartb = kwetop
*        istopb points to the last column we will be updating
         istopb = ihi
         npoles = min(we,kwsbot-kwstop+1)

         do i = 1,npoles
*           Set the pole
            call z_setpole_end(istartb,.true.,kwebot,alpha(kwstop+i-1),
     $         beta(kwstop+i-1),A,ldA,B,ldB,we+1,kwetop-1,Zc,ldZc)
*           Chase the pole up
            do ishift = kwebot-2,kwetop-2+i,-1
               call z_swap11(istartb,istopb,.true.,.true.,ishift,1,A,
     $            ldA,B,ldB,we,kwetop,Qc,ldQc,we+1,kwetop-1,Zc,ldZc)
            end do
         end do
         do i = npoles+1,we
*           Set the pole to infinity
            call z_setpole_end(istartb,.true.,kwebot,cone,czero,A,ldA,B,
     $         ldB,we+1,kwetop-1,Zc,ldZc)
*           Chase the pole up
            do ishift = kwebot-2,kwetop-2+i,-1
               call z_swap11(istartb,istopb,.true.,.true.,ishift,1,A,
     $            ldA,B,ldB,we,kwetop,Qc,ldQc,we+1,kwetop-1,Zc,ldZc)
            end do
         end do

      end if

*     Apply Qc and Zc to rest of the matrix
      if (wantS) then
         istartm = 1
         istopm = n
      else
         istartm = ilo
         istopm = ihi
      end if

      if (istopm-ihi > 0) then
         call zgemm('C','N',we,istopm-ihi,we,cone,Qc,ldQc,A(kwetop,
     $      ihi+1),ldA,czero,work,we)
         call zlacpy('ALL',we,istopm-ihi,work,we,A(kwetop,ihi+1),ldA)
         call zgemm('C','N',we,istopm-ihi,we,cone,Qc,ldQc,B(kwetop,
     $      ihi+1),ldB,czero,work,we)
         call zlacpy('ALL',we,istopm-ihi,work,we,B(kwetop,ihi+1),ldB)
      end if
      if (wantQ) then
         call zgemm('N','N',n,we,we,cone,Q(1,kwetop),ldQ,Qc,ldQc,czero,
     $      work,n)
         call zlacpy('ALL',n,we,work,n,Q(1,kwetop),ldQ)
      end if

      if (kwetop-1-istartm+1 > 0) then
         call zgemm('N','N',kwetop-istartm,we+1,we+1,cone,A(istartm,
     $      kwetop-1),ldA,Zc,ldZc,czero,work,kwetop-istartm)
        call zlacpy('ALL',kwetop-istartm,we+1,work,kwetop-istartm,
     $     A(istartm,kwetop-1),ldA)
         call zgemm('N','N',kwetop-istartm,we+1,we+1,cone,B(istartm,
     $      kwetop-1),ldB,Zc,ldZc,czero,work,kwetop-istartm)
        call zlacpy('ALL',kwetop-istartm,we+1,work,kwetop-istartm,
     $     B(istartm,kwetop-1),ldB)
      end if
      if (wantZ) then
         call zgemm('N','N',n,we+1,we+1,cone,Z(1,kwetop-1),ldZ,Zc,ldZc,
     $      czero,work,n)
         call zlacpy('ALL',n,we+1,work,n,Z(1,kwetop-1),ldZ)
      end if


      end subroutine aed_rqz_complex

      end module z_rqz_aed
