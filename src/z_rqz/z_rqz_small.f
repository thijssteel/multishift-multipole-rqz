      module z_rqz_small
         use common
         use z_swaps
         implicit none

      contains

      subroutine z_rqz(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,Z,
     $   ldZ,alpha,beta)

*     Arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ

      double complex,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,
     $   *),alpha(*),beta(*)

*     Local scalars
      double precision :: smlnum,ulp,safmin,safmax
      double complex :: s1,s2,ss,eshift
      integer :: istartm,istopm,istart,istop,iiter,maxit,istart2,k,ld

*     External Functions
      double precision,external :: dlamch

*     External Subroutines
      external            :: dlartg

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      CALL dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

      istart = ilo
      istop = ihi
      maxit = 30*(ihi-ilo+1)
      ld = 0

      do iiter = 1,maxit
         if (istart .ge. istop) then
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

         ld = ld+1
*        Calculate shift
         if (mod(ld,6) .eq. 0) then
*           TODO: new shift strategy
            call z_eig22(A(istop-1,istop-1),ldA,B(istop-1,istop-1),ldB,
     $         s1,s2,ss)
         else
            call z_eig22(A(istop-1,istop-1),ldA,B(istop-1,istop-1),ldB,
     $         s1,s2,ss)
         end if

*        Get range to apply rotations to
         if (wantS) then
            istartm = 1
            istopm = n
         else
            istartm = istart2
            istopm = istop
         end if

*        Introduce shift
         call z_setpole_start(istopm,wantQ,istart2,s1,ss,A,ldA,B,ldB,n,
     $      1,Q,ldQ)

*        The actual RQZ sweep
         do k = istart2,istop-2

*           Apply transformation from the right
            call z_swap11(istartm,istopm,wantQ,wantZ,k,1,A,ldA,B,ldB,n,
     $         1,Q,ldQ,n,1,Z,ldZ)

         end do

         call z_eig22(A(istart2,istart2),ldA,B(istart2,istart2),ldB,s1,
     $      s2,ss)

*        Set the pole
         call z_setpole_end(istartm,wantZ,istop,cone,czero,A,ldA,B,ldB,
     $      n,1,Z,ldZ)

      end do

      end subroutine z_rqz

      end module z_rqz_small
