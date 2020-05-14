      module z_rqz_sweep
         use common
         use z_swaps
         implicit none

      contains

      subroutine rqz_sweep_complex(wantS,wantQ,wantZ,n,ilo,ihi,nshifts,
     $   nblock_desired,alphas,betas,A,ldA,B,ldB,Q,ldQ,Z,ldZ,Qc,ldQc,Zc,
     $   ldZc,work,lwork,info)

*     Function arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,lwork,nshifts,
     $   nblock_desired,ldQc,ldZc

      double complex,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,
     $   *),Qc(ldQc,*),Zc(ldZc,*),work(*),alphas(*),betas(*)

      integer,intent(out) :: info

*     Local scalars
      integer :: i,j,ns,istartm,istopm,sheight,swidth,k,np,istartb,
     $   istopb,ishift,nblock,npos,lwork_req

      if (nblock_desired < nshifts+1) then
         info =-1
         return
      end if
      lwork_req = n*nblock_desired
      if (lwork .eq.-1) then
*     workspace query, quick return
         work(1) = lwork_req
         info = 0
         return
      else if (lwork .lt. lwork_req) then
         info =-1
         return
      end if

*     Executable statements

      if (nshifts .lt. 1) then
         return
      end if

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

      ns = nshifts
      npos = max(nblock_desired-ns,1)

*     The following block introduces the shifts and chases
*     them down one by one just enough to make space for
*     the other shifts. The near-the-diagonal block is
*     of size (ns+1) x ns.

      call zlaset('Full',ns+1,ns+1,czero,cone,Qc,ldQc)
      call zlaset('Full',ns,ns,czero,cone,Zc,ldZc)

*     istartb points to the first row we will be updating
      istartb = ilo
*     istopb points to the last column we will be updating
      istopb = ilo+ns-1
      do i = 1,ns
*     Introduce the shift
        call z_setpole_start(istopb,.true.,ilo,alphas(i),betas(i),A,ldA,
     $     B,ldB,ns+1,ilo,Qc,ldQc)

*     Chase the shift down
         do ishift = ilo,ilo+ns-i-1

            call z_swap11(istartb,istopb,.true.,.true.,ishift,1,A,ldA,B,
     $         ldB,ns+1,ilo,Qc,ldQc,ns,ilo,Zc,ldZc)

         end do

      end do

*     Update the rest of the pencil

*     Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
*     from the left with Qc(1:ns+1,1:ns+1)'
      sheight = ns+1
      swidth = istopm-(ilo+ns)+1
      if (swidth > 0) then
        call zgemm('C','N',sheight,swidth,sheight,cone,Qc,ldQc,A(ilo,
     $     ilo+ns),ldA,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,A(ilo,ilo+ns),
     $      ldA)
        call zgemm('C','N',sheight,swidth,sheight,cone,Qc,ldQc,B(ilo,
     $     ilo+ns),ldB,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,B(ilo,ilo+ns),
     $      ldB)
      end if
      if (wantQ) then
       call zgemm('N','N',n,sheight,sheight,cone,Q(1,ilo),ldQ,Qc,ldQc,
     $    czero,work,n)
         call zlacpy('ALL',n,sheight,work,n,Q(1,ilo),ldQ)
      end if

*     Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
*     from the right with Zc(1:ns,1:ns)
      sheight = ilo-1-istartm+1
      swidth = ns
      if (sheight > 0) then
         call zgemm('N','N',sheight,swidth,swidth,cone,A(istartm,ilo),
     $      ldA,Zc,ldZc,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,A(istartm,ilo),
     $      ldA)
         call zgemm('N','N',sheight,swidth,swidth,cone,B(istartm,ilo),
     $      ldB,Zc,ldZc,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,B(istartm,ilo),
     $      ldB)
      end if
      if (wantZ) then
         call zgemm('N','N',n,swidth,swidth,cone,Z(1,ilo),ldZ,Zc,ldZc,
     $      czero,work,n)
         call zlacpy('ALL',n,swidth,work,n,Z(1,ilo),ldZ)
      end if

*     The following block chases the shifts down to the bottom
*     right block. If possible, a shift is moved down npos
*     positions at a time

      k = ilo
      do while (k < ihi-ns)
*     do while (k < ilo+1)
         np = min(ihi-ns-k,npos)
*     Size of the near-the-diagonal block
         nblock = ns+np
*     istartb points to the first row we will be updating
         istartb = k+1
*     istopb points to the last column we will be updating
         istopb = k+nblock-1

         call zlaset('Full',ns+np,ns+np,czero,cone,Qc,ldQc)
         call zlaset('Full',ns+np,ns+np,czero,cone,Zc,ldZc)

*     Near the diagonal shift chase
         do i = ns,1,-1
            do j = 0,np-1
*     Move down the shift with index k+i+j, updating
*     the (ns+np x ns+np) block:
*     (k:k+ns+np,k:k+ns+np-1)
               ishift = k+i+j-1

             call z_swap11(istartb,istopb,.true.,.true.,ishift,1,A,ldA,
     $          B,ldB,nblock,k+1,Qc,ldQc,nblock,k,Zc,ldZc)

            end do
         end do

*     Update rest of the pencil

*     Update A(k+1:k+ns+np, k+ns+np:istopm) and
*     B(k+1:k+ns+np, k+ns+np:istopm)
*     from the left with Qc(1:ns+np,1:ns+np)'
         sheight = ns+np
         swidth = istopm-(k+ns+np)+1
         if (swidth > 0) then
        call zgemm('C','N',sheight,swidth,sheight,cone,Qc,ldQc,A(k+1,
     $     k+ns+np),ldA,czero,work,sheight)
            call zlacpy('ALL',sheight,swidth,work,sheight,A(k+1,
     $         k+ns+np),ldA)
        call zgemm('C','N',sheight,swidth,sheight,cone,Qc,ldQc,B(k+1,
     $     k+ns+np),ldB,czero,work,sheight)
            call zlacpy('ALL',sheight,swidth,work,sheight,B(k+1,
     $         k+ns+np),ldB)
         end if
         if (wantQ) then
       call zgemm('N','N',n,nblock,nblock,cone,Q(1,k+1),ldQ,Qc,ldQc,
     $    czero,work,n)
            call zlacpy('ALL',n,nblock,work,n,Q(1,k+1),ldQ)
         end if

*     Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
*     from the right with Zc(1:ns+np,1:ns+np)
         sheight = k-istartm+1
         swidth = nblock
         if (sheight > 0) then
            call zgemm('N','N',sheight,swidth,swidth,cone,A(istartm,k),
     $         ldA,Zc,ldZc,czero,work,sheight)
            call zlacpy('ALL',sheight,swidth,work,sheight,A(istartm,k),
     $         ldA)
            call zgemm('N','N',sheight,swidth,swidth,cone,B(istartm,k),
     $         ldB,Zc,ldZc,czero,work,sheight)
            call zlacpy('ALL',sheight,swidth,work,sheight,B(istartm,k),
     $         ldB)
         end if
         if (wantZ) then
           call zgemm('N','N',n,nblock,nblock,cone,Z(1,k),ldZ,Zc,ldZc,
     $        czero,work,n)
            call zlacpy('ALL',n,nblock,work,n,Z(1,k),ldZ)
         end if

         k = k+np

      end do

*     The following block removes the shifts from the bottom right corner
*     one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

      call zlaset('Full',ns,ns,czero,cone,Qc,ldQc)
      call zlaset('Full',ns+1,ns+1,czero,cone,Zc,ldZc)

*     istartb points to the first row we will be updating
      istartb = ihi-ns+1
*     istopb points to the last column we will be updating
      istopb = ihi

      do i = 1,ns,1
*     Chase the shift down to the bottom right corner

         do ishift = ihi-i,ihi-2

            call z_swap11(istartb,istopb,.true.,.true.,ishift,1,A,ldA,B,
     $         ldB,ns,ihi-ns+1,Qc,ldQc,ns+1,ihi-ns,Zc,ldZc)

         end do

*     Remove the shift
         call z_setpole_end(istartb,.true.,ihi,cone,czero,A,ldA,B,ldB,
     $      ns+1,ihi-ns,Zc,ldZc)

      end do

*     Update rest of the pencil

*     Update A(ihi-ns+1:ihi, ihi+1:istopm)
*     from the left with Qc(1:ns,1:ns)'
      sheight = ns
      swidth = istopm-(ihi+1)+1
      if (swidth > 0) then
        call zgemm('C','N',sheight,swidth,sheight,cone,Qc,ldQc,
     $     A(ihi-ns+1,ihi+1),ldA,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,A(ihi-ns+1,
     $      ihi+1),ldA)
        call zgemm('C','N',sheight,swidth,sheight,cone,Qc,ldQc,
     $     B(ihi-ns+1,ihi+1),ldB,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,B(ihi-ns+1,
     $      ihi+1),ldB)
      end if
      if (wantQ) then
        call zgemm('N','N',n,ns,ns,cone,Q(1,ihi-ns+1),ldQ,Qc,ldQc,czero,
     $     work,n)
         call zlacpy('ALL',n,ns,work,n,Q(1,ihi-ns+1),ldQ)
      end if

*     Update A(istartm:ihi-ns,ihi-ns:ihi)
*     from the right with Zc(1:ns+1,1:ns+1)
      sheight = ihi-ns-istartm+1
      swidth = ns+1
      if (sheight > 0) then
         call zgemm('N','N',sheight,swidth,swidth,cone,A(istartm,
     $      ihi-ns),ldA,Zc,ldZc,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,A(istartm,
     $      ihi-ns),ldA)
         call zgemm('N','N',sheight,swidth,swidth,cone,B(istartm,
     $      ihi-ns),ldB,Zc,ldZc,czero,work,sheight)
         call zlacpy('ALL',sheight,swidth,work,sheight,B(istartm,
     $      ihi-ns),ldB)
      end if
      if (wantZ) then
         call zgemm('N','N',n,ns+1,ns+1,cone,Z(1,ihi-ns),ldZ,Zc,ldZc,
     $      czero,work,n)
         call zlacpy('ALL',n,ns+1,work,n,Z(1,ihi-ns),ldZ)
      end if

      info = 0

      end subroutine rqz_sweep_complex

      end module z_rqz_sweep
