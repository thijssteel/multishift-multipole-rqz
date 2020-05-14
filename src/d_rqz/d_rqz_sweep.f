      module d_rqz_sweep
         use common
         use d_swaps
         use d_manipulate_poles
         use d_swaps2
         implicit none

      contains

      subroutine rqz_sweep_real(wantS,wantQ,wantZ,n,ilo,ihi,nshifts,
     $   nblock_desired,sr,si,ss,A,ldA,B,ldB,Q,ldQ,Z,ldZ,Qc,ldQc,Zc,
     $   ldZc,work,lwork,info)

*     Function arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,lwork,nshifts,
     $   nblock_desired,ldQc,ldZc

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),Qc(ldQc,*),Zc(ldZc,*),work(*),sr(*),si(*),ss(*)

      integer,intent(out) :: info

*     Local scalars
      integer :: i,j,ns,istartm,istopm,sheight,swidth,k,np,istartb,
     $   istopb,ishift,nblock,npos,n1,n2,swap_info
      double precision :: temp,v(3),c1,s1,c2,s2,H(2,3)
      logical :: bulge

      if (nblock_desired < nshifts+1) then
         info =-1
         return
      end if
      if (lwork .eq.-1) then
*        workspace query, quick return
         work(1) = n*nblock_desired
         info = 0
         return
      else if (lwork .lt. n*nblock_desired) then
         info =-1
         return
      end if

*     Executable statements

      if (ilo .ge. ihi) then
         return
      end if

      if (wantS) then
         istartm = 1
         istopm = n
      else
         istartm = ilo
         istopm = ihi
      end if

*     Shuffle shifts into pairs of real shifts and pairs
*     of complex conjugate shifts assuming complex
*     conjugate shifts are already adjacent to one
*     another
*     TODO

      ns = nshifts
      npos = max(nblock_desired-ns,1)
      npos = min(npos,ihi-ilo-ns)

*     The following block introduces the shifts and chases
*     them down one by one just enough to make space for
*     the other shifts. The near-the-diagonal block is
*     of size (ns+1) x ns.

      np = max(npos-1,1)
      if(A(ilo+ns+np+2-1,ilo+ns+np-1) .ne. zero) then
         np = np-1
      end if

      call dlaset('Full',ns+np+1,ns+np+1,zero,one,Qc,ldQc)
      call dlaset('Full',ns+np,ns+np,zero,one,Zc,ldZc)

      i = 1
      do while(i .le. ns)

         if( i .lt. ns .and. si(i) .ne. zero ) then
*           Introduce double shift
            n1 = 2
            bulge = .false.
            if ( ilo+3 .le. ihi ) then
               bulge = A(ilo+3,ilo+1) .ne. zero
            end if
            if (bulge) then
               call d_swap(ilo,ilo+ns+np-1,.true.,.true.,ilo,1,2,1,A,
     $            ldA,B,ldB,ns+np+1,ilo,Qc,ldQc,ns+np,ilo,Zc,ldZc,work,
     $            lwork,swap_info)
            end if
            call d_remove_double_pole_start(.true.,.true.,ilo,A,ldA,B,
     $         ldB,ns+np+1,ilo,Qc,ldQc,ns+np,ilo,Zc,ldZc,ilo,
     $         ilo+ns+np-1)

            call shiftcolumn(A(ilo,ilo),ldA,B(ilo,ilo),ldB,sr(i),
     $         sr(i+1),si(i),ss(i),v)

            temp = v(2)
            call dlartg(temp,v(3),c1,s1,v(2))
            call dlartg(v(1),v(2),c2,s2,temp)

            call drot(ns+np,A(ilo+1,ilo),ldA,A(ilo+2,ilo),ldA,c1,s1)
            call drot(ns+np,A(ilo,ilo),ldA,A(ilo+1,ilo),ldA,c2,s2)
            call drot(ns+np,B(ilo+1,ilo),ldB,B(ilo+2,ilo),ldB,c1,s1)
            call drot(ns+np,B(ilo,ilo),ldB,B(ilo+1,ilo),ldB,c2,s2)
            call drot(ns+np+1,Qc(1,2),1,Qc(1,3),1,c1,s1)
            call drot(ns+np+1,Qc(1,1),1,Qc(1,2),1,c2,s2)

         else
*           Introduce single shift
            n1 = 1
            bulge = .false.
            if ( ilo+2 .le. ihi ) then
               bulge = A(ilo+2,ilo) .ne. zero
            end if
            if (bulge) then
*              There is a 2x2 block in the way, clear it first
               call d_remove_double_pole_start(.true.,.true.,ilo,A,ldA,
     $            B,ldB,ns+np+1,ilo,Qc,ldQc,ns+np,ilo,Zc,ldZc,ilo,
     $            ilo+ns+np-1)
            end if
            call d_setpole_start(ilo,ilo+ns+np-1,.true.,.true.,ilo,
     $         sr(i),ss(i),A,ldA,B,ldB,ns+np+1,ilo,Qc,ldQc,ns+np,ilo,Zc,
     $         ldZc)

         end if

         j = 1
         do while(j .le. ns+np-i+1-n1)
            n2 = 1
            if(A(ilo+j+3-1+n1-1,ilo+j+1-1+n1-1) .ne. zero) then
               n2 = 2
            end if
            call d_swap(ilo,ilo+ns+np-1,.true.,.true.,ilo+j-1,n1,n2,1,A,
     $         ldA,B,ldB,ns+np+1,ilo,Qc,ldQc,ns+np,ilo,Zc,ldZc,work,
     $         lwork,swap_info)
            j = j+n2
         end do
         i = i+n1
         

      end do

*     Update the rest of the pencil

*     Update A(ilo:ilo+ns+np,ilo+ns+np:istopm) and B(ilo:ilo+ns+np,ilo+np:istopm)
*     from the left with Qc(1:ns+1,1:ns+1)'
      sheight = ns+np+1
      swidth = istopm-(ilo+ns+np)+1
      if (swidth > 0) then
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,A(ilo,
     $      ilo+ns+np),ldA,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,A(ilo,ilo+ns+np),
     $      ldA)
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,B(ilo,
     $      ilo+ns+np),ldB,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,B(ilo,ilo+ns+np),
     $      ldB)
      end if
      if (wantQ) then
        call dgemm('N','N',n,sheight,sheight,one,Q(1,ilo),ldQ,Qc,ldQc,
     $     zero,work,n)
         call dlacpy('ALL',n,sheight,work,n,Q(1,ilo),ldQ)
      end if

*     Update A(istartm:ilo-1,ilo:ilo+ns+np-1) and B(istartm:ilo-1,ilo:ilo+ns+np-1)
*     from the right with Zc(1:ns,1:ns)
      sheight = ilo-1-istartm+1
      swidth = ns+np
      if (sheight > 0) then
         call dgemm('N','N',sheight,swidth,swidth,one,A(istartm,ilo),
     $      ldA,Zc,ldZc,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,A(istartm,ilo),
     $      ldA)
         call dgemm('N','N',sheight,swidth,swidth,one,B(istartm,ilo),
     $      ldB,Zc,ldZc,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,B(istartm,ilo),
     $      ldB)
      end if
      if (wantZ) then
         call dgemm('N','N',n,swidth,swidth,one,Z(1,ilo),ldZ,Zc,ldZc,
     $      zero,work,n)
         call dlacpy('ALL',n,swidth,work,n,Z(1,ilo),ldZ)
      end if

*     The following block chases the shifts down to the bottom
*     right block. If possible, a shift is moved down npos
*     positions at a time

      k = ilo+np
      do while (k < ihi-ns)
         np = min(ihi-ns-k,npos)
         if( k+ns+np+1 .le. ihi ) then
            if (A(k+ns+np+1,k+ns+np-1) .ne. zero) then
               np = np-1
            end if
         end if
*        Size of the near-the-diagonal block
         nblock = ns+np
*        istartb points to the first row we will be updating
         istartb = k+1
*        istopb points to the last column we will be updating
         istopb = k+nblock-1

         call dlaset('Full',ns+np,ns+np,zero,one,Qc,ldQc)
         call dlaset('Full',ns+np,ns+np,zero,one,Zc,ldZc)

*        Near the diagonal shift chase
         i = ns
         do while (i .ge. 1)
            n1 = 1
            if (i .gt. 1) then
               if (A(k+i,k+i-2) .ne. zero) then
                  n1 = 2
               end if
            end if
            j = 0
            do while (j .lt. np)
               n2 = 1
               if (k+i+j+2 .le. ihi) then
                  if (A(k+j+i+2,k+i+j) .ne. zero) then
                     n2 = 2
                  end if
               end if
               call d_swap(istartb,istopb,.true.,.true.,k+i-n1+j,n1,n2,
     $            1,A,ldA,B,ldB,ns+np,k+1,Qc,ldQc,ns+np,k,Zc,ldZc,work,
     $            lwork,swap_info)
               j = j+n2
            end do
            i = i-n1
         end do

*        Update rest of the pencil

*        Update A(k+1:k+ns+np, k+ns+np:istopm) and
*        B(k+1:k+ns+np, k+ns+np:istopm)
*        from the left with Qc(1:ns+np,1:ns+np)'
         sheight = ns+np
         swidth = istopm-(k+ns+np)+1
         if (swidth > 0) then
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,A(k+1,
     $      k+ns+np),ldA,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,A(k+1,
     $         k+ns+np),ldA)
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,B(k+1,
     $      k+ns+np),ldB,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,B(k+1,
     $         k+ns+np),ldB)
         end if
         if (wantQ) then
        call dgemm('N','N',n,nblock,nblock,one,Q(1,k+1),ldQ,Qc,ldQc,
     $     zero,work,n)
            call dlacpy('ALL',n,nblock,work,n,Q(1,k+1),ldQ)
         end if

*        Update A(istartm:k,k:k+ns+np-1) and B(istartm:k,k:k+ns+np-1)
*        from the right with Zc(1:ns+np,1:ns+np)
         sheight = k-istartm+1
         swidth = nblock
         if (sheight > 0) then
            call dgemm('N','N',sheight,swidth,swidth,one,A(istartm,k),
     $         ldA,Zc,ldZc,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,A(istartm,k),
     $         ldA)
            call dgemm('N','N',sheight,swidth,swidth,one,B(istartm,k),
     $         ldB,Zc,ldZc,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,B(istartm,k),
     $         ldB)
         end if
         if (wantZ) then
            call dgemm('N','N',n,nblock,nblock,one,Z(1,k),ldZ,Zc,ldZc,
     $         zero,work,n)
            call dlacpy('ALL',n,nblock,work,n,Z(1,k),ldZ)
         end if

         k = k+np

      end do

      info = 0

      end subroutine

      end module
