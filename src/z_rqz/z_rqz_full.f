      module z_rqz_full
         use common
         use z_rqz_small
         use z_rqz_aed
         use z_rqz_sweep
         implicit none

      contains

      subroutine z_rqz_f(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,
     $   Z,ldZ,alpha,beta,work,lwork,n_aed,n_sweep,n_shifts)

*     Arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,lwork

      double complex,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,
     $   *),alpha(*),beta(*),work(*)

      integer,intent(out) :: n_aed,n_sweep,n_shifts

*     Local scalars
      double precision :: smlnum,ulp,safmin,safmax
      double complex :: eshift
      integer :: istart,istop,iiter,maxit,istart2,k,ld,nshifts,nblock,
     $   we,ws,nmin,nibble,n_undeflated,n_deflated_start,n_deflated_end,
     $   ns,sweep_info,shiftpos,lworkreq,nit_total,nshifts_total

*     External Functions
      double precision,external :: dlamch

*     Get the parameters to find out required workspace
      call z_rqz_getparameters(n,nshifts,nblock,ws,we,nmin,nibble)
      k = max(we+1,ws,nmin)
      lworkreq = max(n*nblock+2*nblock**2,n*k+2*k**2)
      if (lwork .eq.-1) then
         work(1) = dble(lworkreq)
         return
      else if (lwork .lt. lworkreq) then
         write (*,*) "workspace provided to multishift qz is too small"
         return
      end if

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      call dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

      istart = ilo
      istop = ihi
      maxit = 3*(ihi-ilo+1)
*     maxit = 10
      ld = 0

      n_aed = 0
      n_sweep = 0
      n_shifts = 0

      do iiter = 1,maxit
         ! write (*,*) iiter,istop-istart
         if (istart+1 .ge. istop) then
            exit
         end if

*        Check deflations at the end
         if (abs(A(istop,istop-1)) .le. max(smlnum,ulp*(abs(A(istop,
     $      istop))+abs(A(istop-1,istop-1)))) .and. abs(B(istop,
     $      istop-1)) .le. max(smlnum,ulp*(abs(B(istop,
     $      istop))+abs(B(istop-1,istop-1))))) then
            A(istop,istop-1) = zero
            B(istop,istop-1) = zero
            istop = istop-1
            ld = 0
            eshift = zero
         end if

         if (istart .ge. istop) then
            exit
         end if

*        Check interior deflations
         istart2 = istart
         do k = istop,istart+1,-1
            if (abs(A(k,k-1)) .le. max(smlnum,ulp*(abs(A(k,k))+abs(A(k-
     $         1,k-1)))) .and. abs(B(k,k-1)) .le. max(smlnum,
     $         ulp*(abs(B(k,k))+abs(B(k-1,k-1))))) then
               A(k,k-1) = zero
               B(k,k-1) = zero
               istart2 = k
               exit
            end if
         end do

*        istart2 now points to the top of the bottom right
*        unreduced Hessenberg block
         if (istart2 .ge. istop) then
            istop = istart2-1
            cycle
         end if

*        Get the parameters, we do this for each iteration because
*        The blocks might get smaller leading to a better choice
*        of parameter being available
         call z_rqz_getparameters(n,nshifts,nblock,ws,we,nmin,nibble)

         if (istop-istart2+1 .lt. nmin) then
*           Setting ws to the size of the subblock will make AED deflate
*           all the eigenvalues. This is slightly more efficient than just
*           using qz_small because the off diagonal part gets updated via BLAS.
            if (istop-istart+1 .lt. nmin) then
               ws = istop-istart+1
               istart2 = istart
            else
               ws = istop-istart2+1
            end if
         end if

         k = max(ws,we+1)
*        Time for aed
         call aed_rqz_complex(wantS,wantQ,wantZ,n,istart2,istop,ws,we,A,
     $      ldA,B,ldB,Q,ldQ,Z,ldZ,n_undeflated,n_deflated_start,
     $      n_deflated_end,alpha,beta,work,k,work(k**2+1),k,work(2*k**2+
     $      1),lwork-2*k**2)

         n_aed = n_aed+1

         if (n_deflated_start > 0) then
            istart2 = istart2+n_deflated_start
         end if

         if (n_deflated_end > 0) then
            istop = istop-n_deflated_end
            ld = 0
            eshift = zero
         end if

         if (100*n_deflated_end > nibble*(n_deflated_end+
     $      n_undeflated).or.istop-istart2+1 .lt. nmin) then
*           write (*,*) "qz step skipped"
*           AED has uncovered many eigenvalues. Skip a QZ sweep and run
*           AED again.
            cycle
         end if

         ld = ld+1

         ns = min(nshifts,istop-istart2,n_undeflated)
         shiftpos = istop-n_deflated_end-n_undeflated+1

         if (mod(ld,6) .eq. 0) then
            call z_eig22(A(istop-1,istop-1),ldA,B(istop-1,istop-1),ldB,
     $         alpha(shiftpos),alpha(shiftpos+1),beta(shiftpos))
            beta(shiftpos+1) = beta(shiftpos)
            ns = 2
         end if

         call rqz_sweep_complex(wantS,wantQ,wantZ,n,istart2,istop,ns,
     $      nblock,alpha(shiftpos),beta(shiftpos),A,ldA,B,ldB,Q,ldQ,Z,
     $      ldZ,work,nblock,work(nblock**2+1),nblock,work(2*nblock**2+
     $      1),lwork-2*nblock**2,sweep_info)

         n_sweep = n_sweep+1

         n_shifts = n_shifts+ns

      end do

*     TODO: store eigenvalues and standardize blocks

      end subroutine z_rqz_f

*     Subroutine: getparameters
*    
*    ------------------------------------------------------------------------------
*     DESCRIPTION:
*    >  Returns some tuning parameters for a given pencil size
*    
*       This subroutine currently only accounts for the size of the unreduced
*       subblock, higher performance might be attained by also accounting for the
*       size of the full matrix.
*    
*    ------------------------------------------------------------------------------
*     ARGUMENTS
*    >  n          integer [IN]
*                     The size of the subblock
*    >  nshifts    integer [OUT]
*                     The desired number of simultanious shifts.
*                     This should always be a multiple of two.
*    >  nblock     integer [OUT]
*                     The desired block size to use during the qz sweep.
*                     Good values for this are a multiple of the blas kernel size.
*    >  nw         integer [OUT]
*                     The desired deflation window size
*    >  nmin       integer [OUT]
*                     Threshold value to switch between double and multishift code.
*    >  nibble     integer [OUT]
*                     Threshold value to skip a QZ sweep, expressed as a percentage.
*    
*    ------------------------------------------------------------------------------
      subroutine z_rqz_getparameters(n,nshifts,nblock,ws,we,nmin,nibble)
      integer,intent(in) :: n
      integer,intent(out) :: nshifts,nblock,ws,we,nmin,nibble

      real,parameter :: c = 0.1d0

      integer :: k

      nmin = 75
      nibble = 8

      if (n .lt. 30) then
         nshifts = 2
         ws = 2
         we = 4
      else if (n .lt. 150) then
         nshifts = 4
         ws = 4
         we = 8
      else if (n .lt. 590) then
         nshifts = 32
         ws = 32
         we = 48
      else if (n .lt. 3000) then
         nshifts = 40
         ws = 40
         we = 96
      else if (n .lt. 6000) then
         nshifts = 64
         ws = 64
         we = 96
      else
         nshifts = 64
         ws = 64
         we = 96
      end if


      k = int(nshifts/sqrt( 1+2*nshifts/(c*n) ))
      k = ((k-1)/4)*4+4
      nblock = nshifts+k

      end subroutine z_rqz_getparameters

      end module z_rqz_full
