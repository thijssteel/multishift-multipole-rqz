      module d_rqz_aed
         use common
         use d_rqz_small
         use d_rqz_small2
         implicit none

      contains

      subroutine aed_rqz_real(wantS,wantQ,wantZ,n,ilo,ihi,ws,we,A,ldA,B,
     $   ldB,Q,ldQ,Z,ldZ,ns,nds,nde,alphar,alphai,beta,Qc,ldQc,Zc,ldZc,
     $   work,lwork)
      implicit none
*     Arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,ldQc,ldZc,lwork
      integer,intent(inout) :: ws,we

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),alphar(*),alphai(*),beta(*)
      integer,intent(out) :: ns,nds,nde
      double precision :: Qc(ldQc,*),Zc(ldZc,*),work(*)

*     Local Scalars
      integer :: kwetop,kwebot,kwstop,kwsbot,istopm,istartm,k,k2,
     $   dtgexc_info,ifst,ilst,lworkreq,i,n_shifts,swapinfo,n1,n2,
     $   qz_small_info
      double precision :: s_top1,s_top2,s_top3,s_bot1,s_bot2,s_bot3,
     $   temp,c,s
      double precision :: smlnum,ulp,safmin,safmax
      logical :: bulge,defl

*     External Functions
      double precision,external :: dlamch

*     External Subroutines
      external :: dgemm,dlacpy,ztgexc

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
      call dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

      ns = 0
      nds = 0
      nde = 0

*     Set up deflation windows
      if( (we .ge. ihi-ilo+1) .or. (we .ge. ihi-ilo+1) ) then
         ws = ihi-ilo+1
         we = 0
         kwetop = ihi
         kwebot = ihi
         kwstop = ilo
         kwsbot = ilo+ws-1
      else
         if( A(ihi-we+2,ihi-we) .ne. zero ) then
            we = we-1
         end if
         kwetop = ihi-we+1
         ws = min(ws,kwetop-ilo+1)
         kwebot = ihi
         if( A(ilo+ws,ilo+ws-2) .ne. zero) then
            ws = ws-1
         end if
         kwstop = ilo
         kwsbot = ilo+ws-1
      end if
      if (kwsbot .eq. ihi) then
         s_top1 = zero
         s_top2 = zero
         s_top3 = zero
      else
         s_top1 = A(kwsbot+1,kwsbot)
         s_top2 = B(kwsbot+1,kwsbot)
         if(kwsbot+2 .le. ihi) then
            s_top3 = A(kwsbot+2,kwsbot)
         else
            s_top3 = zero
         end if
      end if
      if (kwetop .eq. ilo) then
         s_bot1 = zero
         s_bot2 = zero
         s_bot3 = zero
      else
         s_bot1 = A(kwetop,kwetop-1)
         s_bot2 = B(kwetop,kwetop-1)
         if(kwetop-2 .ge. ilo) then
            s_bot3 = A(kwetop,kwetop-2)
         else
            s_bot3 = zero
         end if
      end if

*     Top AED part

      if (ihi .eq. kwstop) then
*        1 by 1 deflation window, just try a regular deflation
         nds = 0
      if ((abs(s_top1)+abs(s_top3)) .le. max(smlnum,ulp*abs(A(kwstop,
     $   kwstop))) .and.abs(s_top2) .le. max(smlnum,ulp*abs(B(kwstop,
     $   kwstop)))) then
            nds = 1
            if (kwstop .gt. ilo) then
               A(kwstop,kwstop-1) = zero
               B(kwstop,kwstop-1) = zero
            end if
         end if
      end if


*     Store window in case of convergence failure
      call dlacpy('ALL',ws,ws,A(kwstop,kwstop),ldA,work,ws)
      call dlacpy('ALL',ws,ws,B(kwstop,kwstop),ldB,work(ws**2+1),ws)

*     Transform window to real schur form
      call dlaset('Full',ws,ws,zero,one,Qc,ldQc)
      call dlaset('Full',ws,ws,zero,one,Zc,ldZc)
      call d_rqz2(.true.,.true.,.true.,ws,1,ws,A(kwstop,kwstop),ldA,
     $   B(kwstop,kwstop),ldB,Qc,ldQc,Zc,ldZc,alphar(kwstop),
     $   alphai(kwstop),beta(kwstop),work(2*ws**2+1),lwork-2*ws**2,
     $   qz_small_info,n_shifts)

      if(qz_small_info .ne. 0) then
*        Convergence failure, restore the window and skip to bottom aed
         write(*,*) "top aed window failed to converge"
         nds = 0
         call dlacpy('ALL',ws,ws,work,ws,A(kwstop,kwstop),ldA)
         call dlacpy('ALL',ws,ws,work(ws**2+1),ws,B(kwstop,kwstop),ldB)
      else

*     Deflation detection loop
      if (kwsbot .eq. ihi) then
         kwstop = kwsbot+1
      else
         kwebot = ihi
         k = 1
         k2 = ws
         do while (k .le. ws)
            defl = .false.
            bulge = .false.
            if (kwsbot-kwstop+1 .ge. 2) then
               bulge = A(kwstop+1,kwstop) .ne. zero
            end if
            if (bulge) then
               n1 = 2
*              Try to deflate complex conjugate eigenvalue pair
               if ((abs(s_top1*Zc(ws,kwstop-ilo+2))+abs(s_top3*Zc(ws,
     $            kwstop-ilo+2))).le.ulp*(abs(A(kwstop+1,
     $            kwstop+1))+abs(A(kwstop+2,kwstop+2))).and.(abs(s_top2*
     $            Zc(ws,kwstop-ilo+2))).le.ulp*(abs(B(kwstop+1,
     $            kwstop+1))+abs(B(kwstop+2,kwstop+2)))) then
*                 Deflatable
                  kwstop = kwstop+2
                  defl = .true.
               end if
            else
*              Try to deflate real eigenvalue
               n1 = 1
               if ((abs(s_top1*Zc(ws,kwstop-ilo+1))+abs(s_top3*Zc(ws,
     $            kwstop-ilo+1))).le.ulp*(abs(A(kwstop,
     $            kwstop))+abs(A(kwstop+1,kwstop+1))).and.(abs(s_top2*
     $            Zc(ws,kwstop-ilo+1))).le.ulp*(abs(B(kwstop,
     $            kwstop))+abs(B(kwstop+1,kwstop+1)))) then
*                 Deflatable
                  kwstop = kwstop+1
                  defl = .true.
               end if
            end if

            if(.not. defl) then
               i = kwstop-ilo+1
               do while(i .lt. ws-n1)
                  n2 = 1
                  if(A(ilo+i-1+n1+1,ilo+i-1+n1) .ne. zero) then
                     n2 = 2
                  end if
                  call d_swap(1,ws,.true.,.true.,i,n1,n2,0,A(ilo,ilo),
     $               ldA,B(ilo,ilo),ldB,ws,1,Qc,ldQc,ws,1,Zc,ldZc,work,
     $               lwork,dtgexc_info)
                  i = i+n2
               end do


!             ifst = kwstop-ilo+1
!             ilst = ws-n1
!             if(ilst .lt. ws .and. ilst .gt. 1) then
!                call dtgexc(.true.,.true.,ws,A(ilo,ilo),ldA,B(ilo,
!   $               ilo),ldB,Qc,ldQc,Zc,ldZc,ifst,ilst,work,lwork,
!   $               dtgexc_info)
!             end if
            end if
            k = k+n1

         end do
      end if

*     Store eigenvalues
      nds = kwstop-ilo
      k = ilo
      do while (k .le. kwsbot)
         bulge = .false.
         if (k .lt. kwsbot) then
            if (A(k+1,k) .ne. zero) then
               bulge = .true.
            end if
         end if
         if (bulge) then
*           2x2 eigenvalue block
            call eig22(A(k:k+1,k:k+1),B(k:k+1,k:k+1),alphar(k),
     $         alphar(k+1),alphai(k),beta(k))
            alphai(k+1) =-alphai(k)
            beta(k+1) = beta(k)

            ! alphar(k) = one
            ! alphai(k) = zero
            ! beta(k) = zero
            ! alphar(k+1) = one
            ! alphai(k+1) = zero
            ! beta(k+1) = zero

            k = k+2

         else
*           1x1 eigenvalue block
            alphar(k) = A(k,k)
            alphai(k) = zero
            beta(k) = B(k,k)

            ! alphar(k) = one
            ! alphai(k) = zero
            ! beta(k) = zero
            k = k+1
         end if
      end do

*     Handle the spike
      if (ws .ne. ihi-ilo+1) then
         do i=1+nds,ws-1
            call dlartg(Zc(ws,i+1),Zc(ws,i),c,s,temp)
            call drot(ws,Zc(1,i+1),1,Zc(1,i),1,c,s)
            call drot(ws,A(kwsbot-ws+1,kwsbot-ws+1+i),1,A(kwsbot-ws+1,
     $         kwsbot-ws+1+i-1),1,c,s)
            call drot(ws,B(kwsbot-ws+1,kwsbot-ws+1+i),1,B(kwsbot-ws+1,
     $         kwsbot-ws+1+i-1),1,c,s)
         end do

         A(kwsbot+1,kwsbot) = Zc(ws,ws)*A(kwsbot+1,kwsbot)
         if(kwsbot+2 .le. ihi) then
            A(kwsbot+2,kwsbot) = Zc(ws,ws)*A(kwsbot+2,kwsbot)
         end if
         B(kwsbot+1,kwsbot) = Zc(ws,ws)*B(kwsbot+1,kwsbot)

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
         call dgemm('C','N',ws,istopm-kwsbot,ws,one,Qc,ldQc,A(ilo,
     $      kwsbot+1),ldA,zero,work,ws)
         call dlacpy('ALL',ws,istopm-kwsbot,work,ws,A(ilo,kwsbot+1),ldA)
         call dgemm('C','N',ws,istopm-kwsbot,ws,one,Qc,ldQc,B(ilo,
     $      kwsbot+1),ldB,zero,work,ws)
         call dlacpy('ALL',ws,istopm-kwsbot,work,ws,B(ilo,kwsbot+1),ldB)
      end if
      if (wantQ) then
         call dgemm('N','N',n,ws,ws,one,Q(1,ilo),ldQ,Qc,ldQc,zero,work,
     $      n)
         call dlacpy('ALL',n,ws,work,n,Q(1,ilo),ldQ)
      end if

      if (ilo-1-istartm+1 > 0) then
         call dgemm('N','N',ilo-istartm,ws,ws,one,A(istartm,ilo),ldA,Zc,
     $      ldZc,zero,work,ilo-istartm)
        call dlacpy('ALL',ilo-istartm,ws,work,ilo-istartm,A(istartm,
     $     ilo),ldA)
         call dgemm('N','N',ilo-istartm,ws,ws,one,B(istartm,ilo),ldB,Zc,
     $      ldZc,zero,work,ilo-istartm)
        call dlacpy('ALL',ilo-istartm,ws,work,ilo-istartm,B(istartm,
     $     ilo),ldB)
      end if
      if (wantZ) then
         call dgemm('N','N',n,ws,ws,one,Z(1,ilo),ldZ,Zc,ldZc,zero,work,
     $      n)
         call dlacpy('ALL',n,ws,work,n,Z(1,ilo),ldZ)
      end if
         
      end if

*     Bottom AED part

      if(we .eq. 0) return

*     Store window in case of convergence failure
  100 call dlacpy('ALL',we,we,A(kwetop,kwetop),ldA,work,we)
      call dlacpy('ALL',we,we,B(kwetop,kwetop),ldB,work(we**2+1),we)

*     Transform window to real schur form
      call dlaset('Full',we,we,zero,one,Qc,ldQc)
      call dlaset('Full',we,we,zero,one,Zc,ldZc)
      call d_rqz2(.true.,.true.,.true.,we,1,we,A(kwetop,kwetop),ldA,
     $   B(kwetop,kwetop),ldB,Qc,ldQc,Zc,ldZc,alphar(kwetop),
     $   alphai(kwetop),beta(kwetop),work(2*we**2+1),lwork-2*we**2,
     $   qz_small_info,n_shifts)

      if(qz_small_info .ne. 0) then
*        Convergence failure, restore the window and exit
         write(*,*) "bottom aed window failed to converge"
         nde = 0
         ns = we-qz_small_info
         call dlacpy('ALL',we,we,work,we,A(kwetop,kwetop),ldA)
         call dlacpy('ALL',we,we,work(we**2+1),we,B(kwetop,kwetop),ldB)
         return
      end if 

*     Deflation detection loop
      if (kwetop .eq. ilo .or. s_bot1 .eq. zero) then
         kwebot = kwetop-1
      else
         k = 1
         k2 = 1
         do while (k .le. we)
            bulge = .false.
            defl = .false.
            if (kwebot-kwetop+1 .ge. 2) then
               bulge = A(kwebot,kwebot-1) .ne. zero
            end if
            if (bulge) then
               n2 = 2
*              Try to deflate complex conjugate eigenvalue pair
               if ((abs(s_bot1*Qc(1,kwebot-kwetop))+abs(s_bot3*Qc(1,
     $            kwebot-kwetop))) .le. ulp*(abs(A(kwebot-1,
     $            kwebot-1))+abs(A(kwebot-2,
     $            kwebot-2))) .and. (abs(s_bot2*Qc(1,
     $            kwebot-kwetop))) .le. ulp*(abs(B(kwebot-1,
     $            kwebot-1))+abs(B(kwebot-2,kwebot-2))) ) then
*                 Deflatable
                  kwebot = kwebot-2
                  defl = .true.
               end if
            else
               n2 = 1
*              Try to deflate real eigenvalue
               if ((abs(s_bot1*Qc(1,kwebot-kwetop+1))+abs(s_bot3*Qc(1,
     $            kwebot-kwetop+1))) .le. ulp*(abs(A(kwebot,
     $            kwebot))+abs(A(kwebot-1,kwebot-1))) .and. (abs(s_bot2*
     $            Qc(1,kwebot-kwetop+1))) .le. ulp*(abs(B(kwebot,
     $            kwebot))+abs(B(kwebot-1,kwebot-1))) ) then
*                 Deflatable
                  kwebot = kwebot-1
                  defl = .true.
               end if
            end if

            if(.not. defl) then
               i = kwebot-kwetop+2-n2
               do while(i .gt. k2)
                  n1 = 1
                  if(i .ge. 3) then
                     if(A(kwetop+i-2,kwetop+i-3) .ne. zero) then
                        n1 = 2
                     end if
                  end if
                  call d_swap(1,we,.true.,.true.,i-n1,n1,n2,0,A(kwetop,
     $               kwetop),ldA,B(kwetop,kwetop),ldB,we,1,Qc,ldQc,we,1,
     $               Zc,ldZc,work,lwork,dtgexc_info)
                  i = i-n1
               end do
!             ifst = kwebot-kwetop+1
!             ilst = k2
!             call dtgexc(.true.,.true.,we,A(kwetop,kwetop),ldA,
!   $            B(kwetop,kwetop),ldB,Qc,ldQc,Zc,ldZc,ifst,ilst,work,
!   $            lwork,dtgexc_info)
               k2 = k2+n2
            end if


            k = k+n2
         end do
      end if

*     Store eigenvalues
   50 nde = ihi-kwebot
      ns = we-nde
      k = kwetop
      do while (k .le. ihi)
         bulge = .false.
         if (k .lt. ihi) then
            if (A(k+1,k) .ne. zero) then
               bulge = .true.
            end if
         end if
         if (bulge) then
*           2x2 eigenvalue block
            call eig22(A(k:k+1,k:k+1),B(k:k+1,k:k+1),alphar(k),
     $         alphar(k+1),alphai(k),beta(k))
            alphai(k+1) =-alphai(k)
            beta(k+1) = beta(k)
            k = k+2
         else
*           1x1 eigenvalue block
            alphar(k) = A(k,k)
            alphai(k) = zero
            beta(k) = B(k,k)
            k = k+1
         end if
      end do

*     Handle the spike
   60 if (kwetop .ne. ilo) then
         do i=we-1-nde,1,-1
            call dlartg(Qc(1,i),Qc(1,i+1),c,s,temp)
            call drot(we,Qc(1,i),1,Qc(1,i+1),1,c,s)
            call drot(we,A(kwetop+i-1,kwetop),ldA,A(kwetop+i,kwetop),
     $         ldA,c,s)
            call drot(we,B(kwetop+i-1,kwetop),ldB,B(kwetop+i,kwetop),
     $         ldB,c,s)
         end do
         
         A(kwetop,kwetop-1) = Qc(1,1)*A(kwetop,kwetop-1)
         if(kwetop-2 .ge. ilo) then
            A(kwetop,kwetop-2) = Qc(1,1)*A(kwetop,kwetop-2)
         end if
         B(kwetop,kwetop-1) = Qc(1,1)*B(kwetop,kwetop-1)

      end if

*     Set the poles
   70 k=1
      do while(k .le. ws-nds)
         if ( alphai(kwstop+k-1) .eq. zero .or. kwebot-
     $      2 .lt. kwetop ) then
            bulge = .false.
            if( kwebot-2 .ge. kwetop ) then
               bulge = A(kwebot,kwebot-2) .ne. zero
            end if
            if(bulge) then
               call d_remove_double_shift(.true.,.true.,we,we-nde,
     $            A(kwetop,kwetop),ldA,B(kwetop,kwetop),ldB,Qc,ldQc,Zc,
     $            ldZc,1,we)
            end if
            call d_setpole_end(1,.true.,we-nde,alphar(kwstop+k-1),
     $         beta(kwstop+k-1),A(kwetop,kwetop),ldA,B(kwetop,kwetop),
     $         ldB,we,1,Zc,ldZc)
            n2 = 1
            k2 = we-nde-n2
            do while(k2 .gt. k)
               n1 = 1
               if(k2 .gt. 2) then
                  if(A(kwetop+k2-1,kwetop+k2-3) .ne. zero) then
                     n1 = 2
                  end if
               end if
               call d_swap(1,we,.true.,.true.,k2-n1,n1,n2,1,A(kwetop,
     $            kwetop),ldA,B(kwetop,kwetop),ldB,we,1,Qc,ldQc,we,1,Zc,
     $            ldZc,work,lwork,swapinfo)
               k2 = k2-n1
            end do
            k = k+1
         else
            bulge = .false.
            if( kwebot-3 .ge. kwetop ) then
               bulge = A(kwebot-1,kwebot-3) .ne. zero
            end if
            if(bulge) then
*              There is a double pole block in second to last position,
*              move it down to make room.
               call d_swap(1,we,.true.,.true.,kwebot-3-kwetop+1,2,1,1,
     $            A(kwetop,kwetop),ldA,B(kwetop,kwetop),ldB,we,1,Qc,
     $            ldQc,we,1,Zc,ldZc,work,lwork,swapinfo)
            end if
            call d_remove_double_shift(.true.,.true.,we,we-nde,A(kwetop,
     $         kwetop),ldA,B(kwetop,kwetop),ldB,Qc,ldQc,Zc,ldZc,1,we)
            call d_setdoublepole_end(1,.true.,we-nde,alphar(kwstop+k-1),
     $         alphar(kwstop+k),alphai(kwstop+k),beta(kwstop+k-1),
     $         A(kwetop,kwetop),ldA,B(kwetop,kwetop),ldB,we,1,Zc,ldZc)
            n2 = 2
            k2 = we-nde-n2
            do while(k2 .gt. k)
               n1 = 1
               if(k2 .gt. 2) then
                  if(A(kwetop+k2-1,kwetop+k2-3) .ne. zero) then
                     n1 = 2
                  end if
               end if
               call d_swap(1,we,.true.,.true.,k2-n1,n1,n2,1,A(kwetop,
     $            kwetop),ldA,B(kwetop,kwetop),ldB,we,1,Qc,ldQc,we,1,Zc,
     $            ldZc,work,lwork,swapinfo)
               k2 = k2-n1
            end do
            k = k+2
         end if
      end do

*     Apply Qc and Zc to rest of the matrix
   80 if (wantS) then
         istartm = 1
         istopm = n
      else
         istartm = ilo
         istopm = ihi
      end if

      if (istopm-ihi > 0) then
         call dgemm('T','N',we,istopm-ihi,we,one,Qc,ldQc,A(kwetop,
     $      ihi+1),ldA,zero,work,we)
      call dlacpy('ALL',we,istopm-ihi,work,we,A(kwetop,ihi+1),ldA)
         call dgemm('T','N',we,istopm-ihi,we,one,Qc,ldQc,B(kwetop,
     $      ihi+1),ldB,zero,work,we)
      call dlacpy('ALL',we,istopm-ihi,work,we,B(kwetop,ihi+1),ldB)
      end if
      if (wantQ) then
         call dgemm('N','N',n,we,we,one,Q(1,kwetop),ldQ,Qc,ldQc,zero,
     $      work,n)
         call dlacpy('ALL',n,we,work,n,Q(1,kwetop),ldQ)
      end if

      if (kwetop-1-istartm+1 > 0) then
         call dgemm('N','N',kwetop-istartm,we,we,one,A(istartm,kwetop),
     $      ldA,Zc,ldZc,zero,work,kwetop-istartm)
        call dlacpy('ALL',kwetop-istartm,we,work,kwetop-istartm,
     $     A(istartm,kwetop),ldA)
         call dgemm('N','N',kwetop-istartm,we,we,one,B(istartm,kwetop),
     $      ldB,Zc,ldZc,zero,work,kwetop-istartm)
        call dlacpy('ALL',kwetop-istartm,we,work,kwetop-istartm,
     $     B(istartm,kwetop),ldB)
      end if
      if (wantZ) then
         call dgemm('N','N',n,we,we,one,Z(1,kwetop),ldZ,Zc,ldZc,zero,
     $      work,n)
         call dlacpy('ALL',n,we,work,n,Z(1,kwetop),ldZ)
      end if
     

      end subroutine aed_rqz_real

      end module d_rqz_aed
