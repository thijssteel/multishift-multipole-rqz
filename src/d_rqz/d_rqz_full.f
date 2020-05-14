      module d_rqz_full
         use common
         use d_rqz_small
         use d_rqz_aed
         use d_rqz_sweep
         implicit none

      contains

      subroutine d_rqz_f(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,
     $   Z,ldZ,alphar,alphai,beta,work,lwork,n_aed,n_sweep,n_shifts)

*     Arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,lwork

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),alphar(*),alphai(*),beta(*),work(*)

      integer,intent(out) :: n_aed,n_sweep,n_shifts

*     Local scalars
      double precision :: smlnum,ulp,safmin,safmax
      double precision :: eshift
      integer :: istart,istop,iiter,maxit,istart2,k,ld,nshifts,nblock,
     $   we,ws,nmin,nibble,n_undeflated,n_deflated_start,n_deflated_end,
     $   ns,sweep_info,shiftpos,lworkreq,nit_total,nshifts_total,
     $   istartm,istopm
      logical :: defl

*     External Functions
      double precision,external :: dlamch

*     Get the parameters to find out required workspace
      call d_rqz_getparameters(n,nshifts,nblock,ws,we,nmin,nibble)
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
         if (istart+1 .ge. istop) then
            exit
         end if

         if (wantS) then
            istartm = 1
            istopm = n
         else
            istartm = istart
            istopm = istop
         end if

         call d_rqz_deflations(wantQ,wantZ,n,istartm,istopm,istart,
     $      istart2,istop,ulp,smlnum,defl,A,ldA,B,ldB,Q,ldQ,Z,ldZ)

         if (istart2+1 .ge. istop) then
            cycle
         end if

*        Get the parameters, we do this for each iteration because
*        The blocks might get smaller leading to a better choice
*        of parameter being available
         call d_rqz_getparameters(n,nshifts,nblock,ws,we,nmin,nibble)

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
         call aed_rqz_real(wantS,wantQ,wantZ,n,istart2,istop,ws,we,A,
     $      ldA,B,ldB,Q,ldQ,Z,ldZ,n_undeflated,n_deflated_start,
     $      n_deflated_end,alphar,alphai,beta,work,k,work(k**2+1),k,
     $      work(2*k**2+1),lwork-2*k**2)
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
   
*           write (*,*) "qz step skipped",n_deflated_end,n_deflated_end+
*     $        n_deflated_end
*           AED has uncovered many eigenvalues. Skip a QZ sweep and run
*           AED again.
            cycle
         end if

         ld = ld+1

         ns = min(nshifts,istop-istart2-1)
         ns = min(ns,n_undeflated)
         shiftpos = istop-n_deflated_end-n_undeflated+1

         if (mod(ld,6) .eq. 0) then
            call eig22(A(istop-1:istop,istop-1:istop),B(istop-1:istop,
     $         istop-1:istop),alphar(shiftpos),alphar(shiftpos+1),
     $         alphai(shiftpos),beta(shiftpos))
            alphai(shiftpos+1) =-alphai(shiftpos)
            beta(shiftpos+1) = beta(shiftpos)
            ns = 2
         end if
         call rqz_sweep_real(wantS,wantQ,wantZ,n,istart2,istop,ns,
     $      nblock,alphar(shiftpos),alphai(shiftpos),beta(shiftpos),A,
     $      ldA,B,ldB,Q,ldQ,Z,ldZ,work,nblock,work(nblock**2+1),nblock,
     $      work(2*nblock**2+1),lwork-2*nblock**2,sweep_info)

         n_sweep = n_sweep+1

         n_shifts = n_shifts+ns

      end do

*     TODO: store eigenvalues and standardize blocks

      end subroutine d_rqz_f

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
      subroutine d_rqz_getparameters(n,nshifts,nblock,ws,we,nmin,nibble)
      integer,intent(in) :: n
      integer,intent(out) :: nshifts,nblock,ws,we,nmin,nibble

      real,parameter :: c = 0.1

      integer :: k

      nmin = 80
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

      end subroutine d_rqz_getparameters

      end module d_rqz_full
