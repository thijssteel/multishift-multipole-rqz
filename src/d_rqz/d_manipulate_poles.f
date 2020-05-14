      module d_manipulate_poles
         use common
         implicit none


      contains

      subroutine d_setpole_end(istartm,wantZ,ihi,alpha,beta,A,ldA,B,ldB,
     $   nz,zStart,Z,ldZ)
*        Arguments
         integer,intent(in) :: istartm,ihi,ldA,ldB,ldZ,zStart,nz
         logical,intent(in) :: wantZ
         double precision,intent(in) :: alpha,beta
         double precision,intent(inout) :: A(ldA,*),B(ldB,*),Z(ldZ,*)

*        Local variables
         double precision :: temp,temp2,temp3,c,s

         temp2 = beta*A(ihi,ihi-1)-alpha*B(ihi,ihi-1)
         temp3 = beta*A(ihi,ihi)-alpha*B(ihi,ihi)

         call dlartg(temp3,temp2,c,s,temp)
         call drot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c,
     $      s)
         call drot(ihi-istartm+1,B(istartm,ihi),1,B(istartm,ihi-1),1,c,
     $      s)
         if (wantZ) then
            call drot(nz,Z(1,ihi-zStart+1),1,Z(1,ihi-1-zStart+1),1,c,s)
         end if

      end subroutine

      subroutine d_setdoublepole_end(istartm,wantZ,ihi,sr1,sr2,si,sb,A,
     $   ldA,B,ldB,nz,zStart,Z,ldZ)
*        Arguments
         integer,intent(in) :: istartm,ihi,ldA,ldB,ldZ,zStart,nz
         logical,intent(in) :: wantZ
         double precision,intent(in) :: sr1,sr2,si,sb
         double precision,intent(inout) :: A(ldA,*),B(ldB,*),Z(ldZ,*)

*        Local variables
         double precision :: temp,temp2,temp3,c1,s1,c2,s2,v(3)


         call polecolumn(A(ihi-2,ihi-2),ldA,B(ihi-2,ihi-2),ldB,sr1,sr2,
     $      si,sb,v)

         call dlartg(v(2),v(1),c1,s1,temp)
         v(2) = temp
         call dlartg(v(3),v(2),c2,s2,temp)

         call drot(ihi-istartm+1,A(istartm,ihi-1),1,A(istartm,ihi-2),1,
     $      c1,s1)
         call drot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c2,
     $      s2)
         call drot(ihi-istartm+1,B(istartm,ihi-1),1,B(istartm,ihi-2),1,
     $      c1,s1)
         call drot(ihi-istartm+1,B(istartm,ihi),1,B(istartm,ihi-1),1,c2,
     $      s2)
         if(wantZ) then
            call drot(nz,Z(1,ihi-1-zStart+1),1,Z(1,ihi-2-zStart+1),1,c1,
     $         s1)
            call drot(nz,Z(1,ihi-zStart+1),1,Z(1,ihi-1-zStart+1),1,c2,
     $         s2)
         end if

      end subroutine

      subroutine d_setpole_start(istartm,istopm,wantQ,wantZ,ilo,alpha,
     $   beta,A,ldA,B,ldB,nq,qStart,Q,ldQ,nz,zStart,Z,ldZ)
*        Arguments
         integer,intent(in) :: istartm,istopm,ilo,ldA,ldB,ldQ,qStart,nq,
     $      zStart,nz,ldZ
         logical,intent(in) :: wantQ,wantZ
         double precision,intent(in) :: alpha,beta
         double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $      Z(ldZ,*)

*        Local variables
         double precision :: c1,s1,c2,s2,temp,temp2,temp3

         if (A(ilo+2,ilo) .ne. zero) then
*           There is a double pole block in the first position => clear it first
            call dlartg(B(ilo+1,ilo),B(ilo+2,ilo),c1,s1,temp)
            B(ilo+2,ilo) = zero
            B(ilo+1,ilo) = temp
            call drot(istopm-ilo+1,A(ilo+1,ilo),ldA,A(ilo+2,ilo),ldA,c1,
     $         s1)
            call drot(istopm-(ilo+1)+1,B(ilo+1,ilo+1),ldB,B(ilo+2,
     $         ilo+1),ldB,c1,s1)
            if (wantQ) then
               call drot(nq,Q(1,ilo+1-qStart+1),1,Q(1,ilo+2-qStart+1),1,
     $            c1,s1)
            end if

            call dlartg(B(ilo,ilo),B(ilo+1,ilo),c1,s1,temp)
            B(ilo+1,ilo) = zero
            B(ilo,ilo) = temp
            call drot(istopm-ilo+1,A(ilo,ilo),ldA,A(ilo+1,ilo),ldA,c1,
     $         s1)
            call drot(istopm-(ilo+1)+1,B(ilo,ilo+1),ldB,B(ilo+1,ilo+1),
     $         ldB,c1,s1)
            if (wantQ) then
               call drot(nq,Q(1,ilo-qStart+1),1,Q(1,ilo+1-qStart+1),1,
     $            c1,s1)
            end if

            call dlartg(A(ilo+1,ilo),A(ilo+2,ilo),c1,s1,temp)
            A(ilo+2,ilo) = zero
            A(ilo+1,ilo) = temp
            call drot(istopm-(ilo+1)+1,A(ilo+1,ilo+1),ldA,A(ilo+2,
     $         ilo+1),ldA,c1,s1)
            call drot(istopm-ilo+1,B(ilo+1,ilo),ldB,B(ilo+2,ilo),ldB,c1,
     $         s1)
            if (wantQ) then
               call drot(nq,Q(1,ilo+1-qStart+1),1,Q(1,ilo+2-qStart+1),1,
     $            c1,s1)
            end if
         end if

         temp2 = beta*A(ilo,ilo)-alpha*B(ilo,ilo)
         temp3 = beta*A(ilo+1,ilo)-alpha*B(ilo+1,ilo)

         call dlartg(temp2,temp3,c1,s1,temp)
         call drot(istopm-ilo+1,A(ilo,ilo),ldA,A(ilo+1,ilo),ldA,c1,s1)
         call drot(istopm-ilo+1,B(ilo,ilo),ldB,B(ilo+1,ilo),ldB,c1,s1)
         if (wantQ) then
            call drot(nq,Q(1,ilo-qStart+1),1,Q(1,ilo+1-qStart+1),1,c1,
     $         s1)
         end if

      end subroutine

      subroutine d_remove_double_pole_start(wantQ,wantZ,ilo,A,ldA,B,ldB,
     $   nq,qStart,Q,ldQ,nz,zStart,Z,ldZ,istartm,istopm)
*        Arguments
         logical,intent(in) :: wantQ,wantZ
         integer,intent(in) :: ilo,ldA,ldB,ldQ,ldZ,istartm,istopm,nq,
     $      qStart,nz,zStart
         double precision :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,*)

*        Local variables
         double precision :: c1,s1,c2,s2,temp

*        Make B upper triangular
         call dlartg(B(ilo+1,ilo),B(ilo+2,ilo),c1,s1,temp)
         B(ilo+1,ilo) = temp
         B(ilo+2,ilo) = zero
         call dlartg(B(ilo,ilo),B(ilo+1,ilo),c2,s2,temp)
         B(ilo,ilo) = temp
         B(ilo+1,ilo) = zero

         call drot(istopm-(ilo+1)+1,B(ilo+1,ilo+1),ldB,B(ilo+2,ilo+1),
     $      ldB,c1,s1)
         call drot(istopm-(ilo+1)+1,B(ilo,ilo+1),ldB,B(ilo+1,ilo+1),ldB,
     $      c2,s2)
         call drot(istopm-ilo+1,A(ilo+1,ilo),ldA,A(ilo+2,ilo),ldA,c1,s1)
         call drot(istopm-ilo+1,A(ilo,ilo),ldA,A(ilo+1,ilo),ldA,c2,s2)
         if(wantQ) then
            call drot(nq,Q(1,ilo+1-qStart+1),1,Q(1,ilo+2-qStart+1),1,c1,
     $         s1)
            call drot(nq,Q(1,ilo-qStart+1),1,Q(1,ilo+1-qStart+1),1,c2,
     $         s2)
         end if

         call dlartg(B(ilo+1,ilo+1),B(ilo+2,ilo+1),c1,s1,temp)
         B(ilo+1,ilo+1) = temp
         B(ilo+2,ilo+1) = zero
         call drot(istopm-(ilo+2)+1,B(ilo+1,ilo+2),ldB,B(ilo+2,ilo+2),
     $      ldB,c1,s1)
         call drot(istopm-ilo+1,A(ilo+1,ilo),ldA,A(ilo+2,ilo),ldA,c1,s1)
         if (wantQ) then
            call drot(nq,Q(1,ilo+1-qStart+1),1,Q(1,ilo+2-qStart+1),1,c1,
     $         s1)
         end if

         call dlartg(A(ilo+2,ilo+1),A(ilo+2,ilo),c1,s1,temp)
         A(ilo+2,ilo+1) = temp
         A(ilo+2,ilo) = zero
         call drot(ilo+1-istartm+1,A(istartm,ilo+1),1,A(istartm,ilo),1,
     $      c1,s1)
         call drot(ilo+2-istartm+1,B(istartm,ilo+1),1,B(istartm,ilo),1,
     $      c1,s1)
         if (wantZ) then
            call drot(nz,Z(1,ilo+1-qStart+1),1,Z(1,ilo-qStart+1),1,c1,
     $         s1)
         end if
         
         call dlartg(B(ilo,ilo),B(ilo+1,ilo),c2,s2,temp)
         B(ilo,ilo) = temp
         B(ilo+1,ilo) = zero
         call drot(istopm-(ilo+1)+1,B(ilo,ilo+1),ldB,B(ilo+1,ilo+1),ldB,
     $      c2,s2)
         call drot(istopm-ilo+1,A(ilo,ilo),ldA,A(ilo+1,ilo),ldA,c2,s2)
         if(wantQ) then
            call drot(nq,Q(1,ilo-qStart+1),1,Q(1,ilo+1-qStart+1),1,c2,
     $         s2)
         end if

      end subroutine

      subroutine d_rqz_deflations(wantQ,wantZ,n,istartm,istopm,istart,
     $   istart2,istop,ulp,smlnum,defl,A,ldA,B,ldB,Q,ldQ,Z,ldZ)
         double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $      Z(ldZ,*)
         integer,intent(in) :: ldA,ldB,ldQ,ldZ
         double precision :: ulp,smlnum,c,s,temp,At(3,3),Bt(3,3)
         integer,intent(in) :: istartm,istopm,n
         integer,intent(inout) :: istart,istop
         integer,intent(out) :: istart2
         logical,intent(in) :: wantQ,wantZ
         logical,intent(out) :: defl

         integer :: k
         logical :: bulge
         
         defl = .false.

         if (istop-2 .lt. istart) then
            defl = .true.
            istop = istart
            return
         end if

*        Check deflations at the end
         bulge = .false.
         if (istop-3 .ge. istart) then
            if ( A(istop-1,istop-3) .ne. zero ) then
               bulge = .true.
            end if
         end if
         if ( bulge ) then

            k = istop-1
            if (abs(A(istop,istop-1)) .le. max(smlnum,ulp*(abs(A(istop,
     $         istop))+abs(A(istop-1,istop-1))) ) .and.abs(B(istop,
     $         istop-1)) .le. max(smlnum,ulp*(abs(B(istop,
     $         istop))+abs(B(istop-1,istop-1))) )) then
               A(istop,istop-1) = zero
               B(istop,istop-1) = zero
               istop = istop-1
               defl = .true.
            end if

         else

            if( A(istop,istop-2) .ne. zero ) then
               k = istop-2
            else
               k = istop-1
            end if

            if ((abs(A(istop,istop-2))+abs(A(istop-1,
     $         istop-2))) .le. max(smlnum,ulp*(abs(A(istop-1,
     $         istop-1))+abs(A(istop-2,istop-2)))) .and. (abs(B(istop,
     $         istop-2))+abs(B(istop-1,istop-2))) .le. max(smlnum,
     $         ulp*(abs(B(istop-1,istop-1))+abs(B(istop-2,
     $         istop-2))))) then
               A(istop,istop-2) = zero
               A(istop-1,istop-2) = zero
               B(istop-1,istop-2) = zero
               B(istop,istop-2) = zero
               istop = istop-2
               defl = .true.
            else if ((abs(A(istop,istop-2))+abs(A(istop,
     $         istop-1))) .le. max(smlnum,ulp*(abs(A(istop,
     $         istop))+abs(A(istop-1,istop-1))) ) .and. (abs(B(istop,
     $         istop-2))+abs(B(istop,istop-1))) .le. max(smlnum,
     $         ulp*(abs(B(istop,istop))+abs(B(istop-1,
     $         istop-1))) )) then
               A(istop,istop-2) = zero
               A(istop,istop-1) = zero
               B(istop,istop-2) = zero
               B(istop,istop-1) = zero
               istop = istop-1
               defl = .true.
            end if

         end if

         if (istart+1 .ge. istop) then
            return
         end if

*        If no direct deflation has been found, check for a deflating rotation
         if(.not. defl) then
            bulge = .false.
            if (istop-3 .ge. istart) then
               if ( A(istop-1,istop-3) .ne. zero ) then
                  bulge = .true.
               end if
            end if
            if ( bulge ) then

               At(1:2,1:2) = A(istop-1:istop,istop-1:istop)
               Bt(1:2,1:2) = B(istop-1:istop,istop-1:istop)

               if( abs(A(istop,istop-1)) .ge. abs(B(istop,
     $            istop-1)) ) then
                  call dlartg(At(2,2),At(2,1),c,s,temp)
                  At(2,2) = temp
                  At(2,1) = zero
                  call drot(1,At(1,2),1,At(1,1),1,c,s)
                  call drot(2,Bt(1,2),1,Bt(1,1),1,c,s)
               else
                  call dlartg(Bt(2,2),Bt(2,1),c,s,temp)
                  Bt(2,2) = temp
                  Bt(2,1) = zero
                  call drot(1,Bt(1,2),1,Bt(1,1),1,c,s)
                  call drot(2,At(1,2),1,At(1,1),1,c,s)
               end if

               if( abs(At(2,1)) .lt. ulp*( abs(At(1,1))+abs(At(2,
     $            2)) ) .and. abs(Bt(2,1)) .lt. ulp*( abs(Bt(1,
     $            1))+abs(Bt(2,2)) ) ) then
                  call drot(istop-istartm+1,A(istartm,istop),1,
     $               A(istartm,istop-1),1,c,s)
                  call drot(istop-istartm+1,B(istartm,istop),1,
     $               B(istartm,istop-1),1,c,s)
                  if(wantZ) then
                     call drot(n,Z(1,istop),1,Z(1,istop-1),1,c,s)
                  end if
                  A(istop,istop-1) = zero
                  B(istop,istop-1) = zero
                  istop = istop-1
                  defl = .true.
               end if

            end if
         end if

         if (istart+1 .ge. istop) then
            return
         end if

*        Check interior deflations
         istart2 = istart
         do while (k .ge. istart+1)
            bulge = .false.
            if ( k-2 .ge. istart ) then
               if ( A(k,k-2) .ne. zero ) then
                  bulge = .true.
               end if
            end if

            if (bulge) then
               if ((abs(A(k,k-2))+abs(A(k,k-1))) .le. ulp*(abs(A(k,
     $            k))+abs(A(k-1,k-1))) .and. (abs(B(k,k-2))+abs(B(k,
     $            k-1))) .le.ulp*(abs(B(k,k))+abs(B(k-1,k-1)))) then
                  A(k,k-2) = zero
                  B(k,k-2) = zero
                  A(k,k-1) = zero
                  B(k,k-1) = zero
                  istart2 = k
                  ! defl = .true.
                  exit
               end if
               if ((abs(A(k,k-2))+abs(A(k-1,k-2))) .le. ulp*(abs(A(k-1,
     $            k-1))+abs(A(k-2,k-2))) .and. (abs(B(k,k-2))+abs(B(k-1,
     $            k-2))) .le. ulp*(abs(B(k-1,k-1))+abs(B(k-2,
     $            k-2)))) then
                  A(k,k-2) = zero
                  B(k,k-2) = zero
                  A(k-1,k-2) = zero
                  B(k-1,k-2) = zero
                  istart2 = k-1
                  ! defl = .true.
                  exit
               end if
               k = k-2
            else
               if (abs(A(k,k-1)) .le. ulp*(abs(A(k,k))+abs(A(k-1,
     $            k-1))) .and. abs(B(k,k-1)) .le.ulp*(abs(B(k,
     $            k))+abs(B(k-1,k-1)))) then
                  A(k,k-1) = zero
                  B(k,k-1) = zero
                  istart2 = k
                  ! defl = .true.
                  exit
               end if
               k = k-1
            end if
         end do

*        istart2 now points to the top of the bottom right
*        unreduced Hessenberg block
         if (istart2+1 .ge. istop) then
            istop = istart2-1
         end if

      end subroutine

      subroutine d_normalize_eigenvalues(wantS,wantQ,wantZ,n,ilo,ihi,A,
     $   ldA,B,ldB,Q,ldQ,Z,ldZ,alphar,alphai,beta,ulp,smlnum)

*        Arguments
         logical,intent(in) :: wantS,wantQ,wantZ
         integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ
         double precision,intent(in) :: ulp,smlnum

         double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $      Z(ldZ,*),alphar(*),alphai(*),beta(*)

*        Local scalars
         double precision :: safmin,safmax,sr1,sr2,si,sb,temp,v(3),c1,
     $      s1,c2,s2,H(2,3),eshift
         integer :: istartm,istopm,istart,istop,iiter,maxit,istart2,k,
     $      ld,i
         logical :: bulge

*        External Subroutines
         external            :: dlartg

         k = ilo
         do while (k .le. ihi)
            bulge = .false.
            if (k .lt. ihi) then
               if (A(k+1,k) .ne. zero .or. B(k+1,k) .ne. zero) then
                  bulge = .true.
               end if
            end if
            if (.not. bulge) then
*              1x1 eigenvalue block
               alphar(k) = A(k,k)
               alphai(k) = zero
               beta(k) = B(k,k)
               k = k+1
            else
*              2x2 eigenvalue block
               call eig22(A(k:k+1,k:k+1),B(k:k+1,k:k+1),alphar(k),
     $            alphar(k+1),alphai(k),beta(k))
               alphai(k+1) =-alphai(k)
               beta(k+1) = beta(k)

               if (wantS) then
                  istartm = 1
                  istopm = n
               else
                  istartm = k
                  istopm = k+1
               end if

               if (alphai(k) .eq. zero) then
*                 2x2 block has real eigenvalues

                  do i = 1,3

*                    Recalculate eigenvalues
                     call eig22(A(k:k+1,k:k+1),B(k:k+1,k:k+1),alphar(k),
     $                  alphar(k+1),alphai(k),beta(k))
                     alphai(k+1) =-alphai(k)
                     beta(k+1) = beta(k)

                     H(1:2,1:2) = beta(k)*A(k:k+1,k:k+1)-alphar(k)*
     $                  B(k:k+1,k:k+1)
                     if (H(2,2) .ge. H(1,2)) then
                        call dlartg(H(2,2),H(2,1),c1,s1,temp)
                     else
                        call dlartg(H(1,2),H(1,1),c1,s1,temp)
                     end if
                     call drot(2,H(1,2),1,H(1,1),1,c1,s1)

                     call drot(k+1-istartm+1,A(istartm,k+1),1,A(istartm,
     $                  k),1,c1,s1)
                     call drot(k+1-istartm+1,B(istartm,k+1),1,B(istartm,
     $                  k),1,c1,s1)
                     if (wantZ) then
                        call drot(n,Z(1,k+1),1,Z(1,k),1,c1,s1)
                     end if

                     if (A(k+1,k) .ge. B(k+1,k)) then
                        call dlartg(A(k,k),A(k+1,k),c1,s1,temp)
                     else
                        call dlartg(B(k,k),B(k+1,k),c1,s1,temp)
                     end if

                     call drot(istopm-k+1,A(k,k),ldA,A(k+1,k),ldA,c1,s1)
                     call drot(istopm-k+1,B(k,k),ldB,B(k+1,k),ldB,c1,s1)
                     if (wantQ) then
                        call drot(n,Q(1,k),1,Q(1,k+1),1,c1,s1)
                     end if

                     if (abs(A(k+1,k)) .lt. ulp*(abs(A(k,k))+abs(A(k+1,
     $                  k+1))) .and. abs(B(k+1,k)) .lt. ulp*(abs(B(k,
     $                  k))+abs(B(k+1,k+1)))) then
                        A(k+1,k) = zero
                        B(k+1,k) = zero
                        exit
                     end if

                  end do
               end if
               if(B(k+1,k).ne.zero)then
*                 Normalize the block
                  call dlartg(B(k,k),B(k+1,k),c1,s1,temp)
                  B(k,k) = temp
                  B(k+1,k) = zero

                  call drot(istopm-k+1,A(k,k),ldA,A(k+1,k),ldA,c1,s1)
                  call drot(istopm-(k+1)+1,B(k,k+1),ldB,B(k+1,k+1),ldB,
     $               c1,s1)
                  if (wantQ) then
                     call drot(n,Q(1,k),1,Q(1,k+1),1,c1,s1)
                  end if
               end if
               k = k+2
            end if
         end do
      end subroutine

      end module