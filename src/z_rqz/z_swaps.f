      module z_swaps
         use common
         implicit none

      contains

*     Subroutine: z_swap11
*    
*    ------------------------------------------------------------------------------
*     DESCRIPTION:
*    >  Swaps either two eigenvalues in the Schur form, or two poles in a
*       Hessenberg-Triangular pair.
*    
*    ------------------------------------------------------------------------------
*     ARGUMENTS
*    >  TODO
*    
*    ------------------------------------------------------------------------------
      subroutine z_swap11(istartm,istopm,wantQ,wantZ,k,off,A,ldA,B,ldB,
     $   qn,qStart,Q,ldQ,zn,zStart,Z,ldZ)

*        Arguments
         integer,intent(in) :: istartm,istopm,k,ldA,ldB,ldQ,ldZ,qStart,
     $      zStart,qn,zn,off
         logical,intent(in) :: wantQ,wantZ
         double complex,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $      Z(ldZ,*)

*        Local variables
         logical :: ab
         double complex :: temp,temp2,temp3,s
         double precision :: c

*        Executable statements
         ab = abs(A(k+off,k))*abs(B(k+off+1,k+1)) .gt.abs(A(k+off+1,
     $      k+1))*abs(B(k+off,k))

         temp2 = B(k+off+1,k+1)*A(k+off,k+1)-A(k+off+1,k+1)*B(k+off,k+1)
         temp3 = B(k+off+1,k+1)*A(k+off,k)-A(k+off+1,k+1)*B(k+off,k)

         call zlartg(temp2,temp3,c,s,temp)
         call zrot(k+off+1-istartm+1,A(istartm,k+1),1,A(istartm,k),1,c,
     $      s)
         call zrot(k+off+1-istartm+1,B(istartm,k+1),1,B(istartm,k),1,c,
     $      s)
         if (wantZ) then
            call zrot(zn,Z(1,k+1-zStart+1),1,Z(1,k-zStart+1),1,c,s)
         end if

         if (ab) then
            call zlartg(B(k+off,k),B(k+off+1,k),c,s,temp)
         else
            call zlartg(A(k+off,k),A(k+off+1,k),c,s,temp)
         end if

         call zrot(istopm-k+1,A(k+off,k),ldA,A(k+off+1,k),ldA,c,s)
         call zrot(istopm-k+1,B(k+off,k),ldB,B(k+off+1,k),ldB,c,s)
         if (wantQ) then
            call zrot(qn,Q(1,k+off-qStart+1),1,Q(1,k+off+1-qStart+1),1,
     $         c,dconjg(s))
         end if

         A(k+1+off,k) = czero
         B(k+1+off,k) = czero

      end subroutine z_swap11

      subroutine z_setpole_end(istartm,wantZ,ihi,alpha,beta,A,ldA,B,ldB,
     $   nz,zStart,Z,ldZ)
*        Arguments
         integer,intent(in) :: istartm,ihi,ldA,ldB,ldZ,zStart,nz
         logical,intent(in) :: wantZ
         double complex,intent(in) :: alpha,beta
         double complex,intent(inout) :: A(ldA,*),B(ldB,*),Z(ldZ,*)

*        Local variables
         double complex :: temp,temp2,temp3,s
         double precision :: c

         temp2 = beta*A(ihi,ihi-1)-alpha*B(ihi,ihi-1)
         temp3 = beta*A(ihi,ihi)-alpha*B(ihi,ihi)

         call zlartg(temp3,temp2,c,s,temp)
         call zrot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c,
     $      s)
         call zrot(ihi-istartm+1,B(istartm,ihi),1,B(istartm,ihi-1),1,c,
     $      s)
         if (wantZ) then
            call zrot(nz,Z(1,ihi-zStart+1),1,Z(1,ihi-1-zStart+1),1,c,s)
         end if

      end subroutine z_setpole_end

      subroutine z_setpole_start(istopm,wantQ,ilo,alpha,beta,A,ldA,B,
     $   ldB,nq,qStart,Q,ldQ)
*        Arguments
         integer,intent(in) :: istopm,ilo,ldA,ldB,ldQ,qStart,nq
         logical,intent(in) :: wantQ
         double complex,intent(in) :: alpha,beta
         double complex,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*)

*        Local variables
         double complex :: temp,temp2,temp3,s
         double precision :: c

         temp2 = beta*A(ilo,ilo)-alpha*B(ilo,ilo)
         temp3 = beta*A(ilo+1,ilo)-alpha*B(ilo+1,ilo)

         call zlartg(temp2,temp3,c,s,temp)
         call zrot(istopm-ilo+1,A(ilo,ilo),ldA,A(ilo+1,ilo),ldA,c,s)
         call zrot(istopm-ilo+1,B(ilo,ilo),ldB,B(ilo+1,ilo),ldB,c,s)
         if (wantQ) then
            call zrot(nq,Q(1,ilo-qStart+1),1,Q(1,ilo+1-qStart+1),1,c,
     $         dconjg(s))
         end if

      end subroutine z_setpole_start

      end module z_swaps
