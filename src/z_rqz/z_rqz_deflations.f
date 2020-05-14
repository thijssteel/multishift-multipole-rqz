      module z_deflations
         implicit none

      contains

      subroutine check_end_deflation(wantS,wantQ,wantZ,ulp,ihi,A,ldA,B,
     $   ldB,Q,ldQ,Z,ldZ)
*     Argument definitions
      logical,intent(in):: wantS,wantQ,wantZ
      double complex,intent(inout):: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,
     $   *)
      integer,intent(in):: ldA,ldB,ldQ,ldZ
      integer,intent(inout) :: ihi
      double precision,intent(in) :: ulp

*     Local variables
      double complex :: At(2,2),s,temp
      double precision :: c

      At = A(ihi-1:ihi,ihi-1:ihi)
      call zlartg(B(ihi,ihi),B(ihi,ihi-1),c,s,temp)
      call zrot(2,At(1,2),1,At(1,1),1,c,s)

      if (abs(At(2,1)) .lt. ulp*(abs(At(1,1))+abs(At(2,2)))) then
*           The last eigenvalue has converged, apply the rotation and deflate it
      end if

      end subroutine

      end module z_deflations
