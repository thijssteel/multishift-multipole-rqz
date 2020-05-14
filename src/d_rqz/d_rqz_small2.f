      module d_rqz_small2
         use common
         use d_swaps2
         use d_swaps
         use d_manipulate_poles
         implicit none

      contains

*     Although a bit slower than rqz, removing the poles and using LAPACK
*     solves some problems relating to the normalization of the resulting schur form
      subroutine d_rqz2(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,Z,
     $   ldZ,alphar,alphai,beta,work,lwork,info,n_shifts)

*     Arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,lwork
      integer,intent(inout) :: n_shifts,info

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),alphar(*),alphai(*),beta(*),work(*)

*     Local scalars
      double precision :: smlnum,ulp,safmin,safmax,sr1,sr2,si,sb,temp,
     $   v(3),c1,s1,c2,s2,eshift
      integer :: istartm,istopm,istart,istop,iiter,maxit,istart2,k,k2,
     $   ld,n1,n2,swapinfo
      logical :: defl,bulge

*     External Functions
      double precision,external :: dlamch

*     External Subroutines
      external            :: dlartg

*     Remove the poles in the pencil

      k = ihi
      do while(k.gt.ilo)
         bulge = .false.
         if(k-2 .ge. ilo) then
            bulge = A(k,k-2) .ne. zero
         end if
         if (bulge) then
            do k2 = k-2,ihi-3
               call d_swap_21_inf(.true.,.true.,k2,A,ldA,B,ldB,n,Q,ldQ,
     $            n,Z,ldZ,1,n)
            end do
            call d_remove_double_shift(.true.,.true.,n,ihi,A,ldA,B,ldB,
     $         Q,ldQ,Z,ldZ,1,n)
            k = k-2
         else
            do k2 = k-1,ihi-2
               call d_swap_11_inf(.true.,.true.,n,k2,A,ldA,B,ldB,Q,ldQ,
     $            Z,ldZ,1,n)
            end do
            call d_remove_single_shift(.true.,n,ihi,A,ldA,B,ldB,Z,ldZ,
     $         1)
            k = k-1
            end if
      end do

*     Use LAPACK

      call dhgeqz('S','V','V',n,ilo,ihi,A,ldA,B,ldB,alphar,alphai,beta,
     $   Q,ldQ,Z,ldZ,work,lwork,info)

      n_shifts = 0

      end subroutine d_rqz2

      end module d_rqz_small2
