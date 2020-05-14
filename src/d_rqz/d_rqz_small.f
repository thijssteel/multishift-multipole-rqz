      module d_rqz_small
         use common
         use d_swaps2
         use d_swaps
         use d_manipulate_poles
         implicit none

      contains

      subroutine d_rqz(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,Z,
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
      integer :: istartm,istopm,istart,istop,iiter,maxit,istart2,k,ld,
     $   n1,n2,swapinfo
      logical :: defl

*     External Functions
      double precision,external :: dlamch

*     External Subroutines
      external            :: dlartg

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      call dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

      istart = ilo
      istop = ihi
      maxit = 30*(ihi-ilo+1)
      ld = 0
      n_shifts = 0
      
      do iiter = 1,maxit+1
         if(iiter .ge. maxit) then
            write(*,*) "warning,rqz did not converge"
            info = istop
            return
         end if
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

         if (defl) then
            ld = 0
            eshift = zero
         end if
         ld = ld+1

         if (istart+1 .ge. istop) then
            exit
         end if
         if (istart2+1 .ge. istop) then
            cycle
         end if

*        Calculate shift
         if (mod(ld,6) .eq. 0) then
            call random_number(sr1)
            sr2 = sr1
            call random_number(si)
            call random_number(sb)
         else
            call eig22(A(istop-1:istop,istop-1:istop),B(istop-1:istop,
     $         istop-1:istop),sr1,sr2,si,sb)
         end if

*        Get range to apply rotations to
         if (wantS) then
            istartm = 1
            istopm = n
         else
            istartm = istart2
            istopm = istop
         end if

         if(si .ne. zero .or. use_double_pole .ne. 1) then
*           Swap poles
            if( A(istart2+3,istart2+1) .ne. zero ) then
               call d_swap(istartm,istopm,wantQ,wantZ,istart2,1,2,1,A,
     $            ldA,B,ldB,n,1,Q,ldQ,n,1,Z,ldZ,work,lwork,swapinfo)
            end if

*           Remove the poles in the pencil
            call dlartg(B(istart2,istart2),B(istart2+1,istart2),c1,s1,
     $         temp)
            B(istart2,istart2) = temp
            B(istart2+1,istart2) = zero
            call drot(istopm-(istart2+1)+1,B(istart2,istart2+1),ldB,
     $         B(istart2+1,istart2+1),ldB,c1,s1)
            call drot(istopm-istart2+1,A(istart2,istart2),ldA,
     $         A(istart2+1,istart2),ldA,c1,s1)
            if (wantQ) then
               call drot(n,Q(1,istart2),1,Q(1,istart2+1),1,c1,s1)
            end if

            call dlartg(B(istart2+1,istart2+1),B(istart2+2,istart2+1),
     $         c2,s2,temp)
            B(istart2+1,istart2+1) = temp
            B(istart2+2,istart2+1) = zero
            call drot(istopm-(istart2+2)+1,B(istart2+1,istart2+2),ldB,
     $         B(istart2+2,istart2+2),ldB,c2,s2)
            call drot(istopm-istart2+1,A(istart2+1,istart2),ldA,
     $         A(istart2+2,istart2),ldA,c2,s2)
            if (wantQ) then
               call drot(n,Q(1,istart2+1),1,Q(1,istart2+2),1,c2,s2)
            end if

            call dlartg(A(istart2+1,istart2),A(istart2+2,istart2),c1,s1,
     $         temp)
            A(istart2+1,istart2) = temp
            A(istart2+2,istart2) = zero
            call drot(istopm-(istart2+1)+1,A(istart2+1,istart2+1),ldA,
     $         A(istart2+2,istart2+1),ldA,c1,s1)
            call drot(istopm-istart2+1,B(istart2+1,istart2),ldB,
     $         B(istart2+2,istart2),ldB,c1,s1)
            if (wantQ) then
               call drot(n,Q(1,istart2+1),1,Q(1,istart2+2),1,c1,s1)
            end if

            call dlartg(B(istart2,istart2),B(istart2+1,istart2),c1,s1,
     $         temp)
            B(istart2,istart2) = temp
            B(istart2+1,istart2) = zero
            call drot(istopm-(istart2+1)+1,B(istart2,istart2+1),ldB,
     $         B(istart2+1,istart2+1),ldB,c1,s1)
            call drot(istopm-istart2+1,A(istart2,istart2),ldA,
     $         A(istart2+1,istart2),ldA,c1,s1)
            if (wantQ) then
               call drot(n,Q(1,istart2),1,Q(1,istart2+1),1,c1,s1)
            end if

            call d_swap(istartm,istopm,wantQ,wantZ,istart2,1,1,1,A,ldA,
     $         B,ldB,n,1,Q,ldQ,n,1,Z,ldZ,work,lwork,swapinfo)


            call d_setpole_start(istartm,istopm,wantQ,wantZ,istart2,one,
     $         zero,A,ldA,B,ldB,n,1,Q,ldQ,n,1,Z,ldZ)

            B(istart2+1,istart2) = zero
            B(istart2+2,istart2+1) = zero

*           Introduce double shift
            n1 = 2
            call shiftcolumn(A(istart2,istart2),ldA,B(istart2,istart2),
     $         ldB,sr1,sr2,si,sb,v)

            call dlartg(v(2),v(3),c1,s1,temp)
            v(2) = temp
            call dlartg(v(1),v(2),c2,s2,temp)

*           Apply rotations from the left
            call drot(istopm-istart2+1,A(istart2+1,istart2),ldA,
     $         A(istart2+2,istart2),ldA,c1,s1)
            call drot(istopm-istart2+1,A(istart2,istart2),ldA,
     $         A(istart2+1,istart2),ldA,c2,s2)
            call drot(istopm-istart2+1,B(istart2+1,istart2),ldB,
     $         B(istart2+2,istart2),ldB,c1,s1)
            call drot(istopm-istart2+1,B(istart2,istart2),ldB,
     $         B(istart2+1,istart2),ldB,c2,s2)
            if (wantQ) then
               call drot(n,Q(1,istart2+1),1,Q(1,istart2+2),1,c1,s1)
               call drot(n,Q(1,istart2),1,Q(1,istart2+1),1,c2,s2)
            end if

*           Reduce to normal form
            call dlartg(B(istart2+1,istart2),B(istart2+2,istart2),c1,s1,
     $         temp)
            B(istart2+1,istart2) = temp
            B(istart2+2,istart2) = zero
            call drot(istopm-(istart2+1)+1,B(istart2+1,istart2+1),ldB,
     $         B(istart2+2,istart2+1),ldB,c1,s1)
            call drot(istopm-istart2+1,A(istart2+1,istart2),ldA,
     $         A(istart2+2,istart2),ldA,c1,s1)
            if(wantQ) then
               call drot(n,Q(1,istart2+1),1,Q(1,istart2+2),1,c1,s1)
            end if
         else
*           Introduce single shift

            call d_setpole_start(istartm,istopm,wantQ,wantZ,istart2,sr1,
     $         sb,A,ldA,B,ldB,n,1,Q,ldQ,n,1,Z,ldZ)

            n1 = 1
         end if

         k = istart2
*        The actual QZ sweep
         do while(k .lt. istop-n1)

            if( k+n1+2 .le. istop ) then
               if( A(k+n1+2,k+n1) .ne. zero ) then
                  n2 = 2
               else
                  n2 = 1
               end if
            else
               n2 = 1
            end if

            call d_swap(istartm,istopm,wantQ,wantZ,k,n1,n2,1,A,ldA,B,
     $         ldB,n,1,Q,ldQ,n,1,Z,ldZ,work,lwork,swapinfo)

            k = k+n2

         end do

*        Remove the shift
         if(n1 .eq. 1) then
            call dlartg(B(istop,istop),B(istop,istop-1),c1,s1,temp)
            B(istop,istop) = temp
            B(istop,istop-1) = zero
            call drot(istop-istartm,B(istartm,istop),1,B(istartm,
     $         istop-1),1,c1,s1)
            call drot(istop-istartm+1,A(istartm,istop),1,A(istartm,
     $         istop-1),1,c1,s1)
            if (wantZ) then
               call drot(n,Z(1,istop),1,Z(1,istop-1),1,c1,s1)
            end if
         else
            call d_remove_double_shift(wantQ,wantZ,n,istop,A,ldA,B,ldB,
     $         Q,ldQ,Z,ldZ,istartm,istopm)
         end if

*        Set the pole
         call eig22(A(istart2:istart2+1,istart2:istart2+1),
     $      B(istart2:istart2+1,istart2:istart2+1),sr1,sr2,si,sb)
         if( si .eq. zero .and. use_double_pole .ne. 2 ) then
            
*           Introduce single pole
            v(1) = sb*A(istop,istop-1)-sr1*B(istop,istop-1)
            v(2) = sb*A(istop,istop)-sr1*B(istop,istop)

            call dlartg(v(2),v(1),c1,s1,temp)
            call drot(istop-istartm+1,A(istartm,istop),1,A(istartm,
     $         istop-1),1,c1,s1)
            call drot(istop-istartm+1,B(istartm,istop),1,B(istartm,
     $         istop-1),1,c1,s1)
            if(wantZ) then
               call drot(n,Z(1,istop),1,Z(1,istop-1),1,c1,s1)
            end if

         else

*           Introduce double pole if there is room
            if( istop .gt. 3 ) then
               if( A(istop-1,istop-3) .eq. zero ) then

                  call d_remove_double_shift(wantQ,wantZ,n,istop,A,ldA,
     $               B,ldB,Q,ldQ,Z,ldZ,istartm,istopm)
                  
                  if(use_double_pole .ne. 0) then

                     call polecolumn(A(istop-2,istop-2),ldA,B(istop-2,
     $                  istop-2),ldB,sr1,sr2,si,sb,v)

                     call dlartg(v(2),v(1),c1,s1,temp)
                     v(2) = temp
                     call dlartg(v(3),v(2),c2,s2,temp)

                     call drot(istop-istartm+1,A(istartm,istop-1),1,
     $                  A(istartm,istop-2),1,c1,s1)
                     call drot(istop-istartm+1,A(istartm,istop),1,
     $                  A(istartm,istop-1),1,c2,s2)
                     call drot(istop-istartm+1,B(istartm,istop-1),1,
     $                  B(istartm,istop-2),1,c1,s1)
                     call drot(istop-istartm+1,B(istartm,istop),1,
     $                  B(istartm,istop-1),1,c2,s2)
                     if(wantZ) then
                        call drot(n,Z(1,istop-1),1,Z(1,istop-2),1,c1,s1)
                        call drot(n,Z(1,istop),1,Z(1,istop-1),1,c2,s2)
                     end if
                  end if

               end if
            end if

         end if

         n_shifts = n_shifts+n1

      end do

      call d_normalize_eigenvalues(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,
     $   ldB,Q,ldQ,Z,ldZ,alphar,alphai,beta,ulp,smlnum)

      info = 0

      end subroutine d_rqz

      end module d_rqz_small
