      module d_swaps
         use common
         implicit none

      contains

      subroutine d_swap_21_inf(wantQ,wantZ,k,A,ldA,B,ldB,nQ,Q,ldQ,nZ,Z,
     $   ldZ,istartm,istopm)
*     Arguments
      logical,intent(in) :: wantQ,wantZ
         integer,intent(in) :: k,ldA,ldB,ldQ,ldZ,istartm,istopm,nQ,nZ
      double precision :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,*)

*     Local variables
      double precision :: H(2,3),c1,s1,c2,s2,temp

      H = B(k+1:k+2,k:k+2)

*     Make H upper triangular
      call dlartg(H(1,1),H(2,1),c1,s1,temp)
      H(2,1) = zero
      H(1,1) = temp
      call drot(2,H(1,2),2,H(2,2),2,c1,s1)

*     Calculate Z1 and Z2
      call dlartg(H(2,3),H(2,2),c1,s1,temp)
      call drot(1,H(1,3),1,H(1,2),1,c1,s1)
      call dlartg(H(1,2),H(1,1),c2,s2,temp)

*     Apply transformations from the right
      call drot(k+3-istartm+1,A(istartm,k+2),1,A(istartm,k+1),1,c1,s1)
      call drot(k+3-istartm+1,A(istartm,k+1),1,A(istartm,k),1,c2,s2)
      call drot(k+3-istartm+1,B(istartm,k+2),1,B(istartm,k+1),1,c1,s1)
      call drot(k+3-istartm+1,B(istartm,k+1),1,B(istartm,k),1,c2,s2)
      if (wantZ) then
         call drot(nZ,Z(1,k+2),1,Z(1,k+1),1,c1,s1)
         call drot(nZ,Z(1,k+1),1,Z(1,k),1,c2,s2)
      end if
      B(k+1,k) = zero
      B(k+2,k) = zero
      B(k+3,k) = zero

*     Calculate Q1 and Q2
      call dlartg(A(k+2,k),A(k+3,k),c1,s1,temp)
      A(k+2,k) = temp
      A(k+3,k) = zero
      call dlartg(A(k+1,k),A(k+2,k),c2,s2,temp)
      A(k+1,k) = temp
      A(k+2,k) = zero

*     Apply transformations from the left
      call drot(istopm-k,A(k+2,k+1),ldA,A(k+3,k+1),ldA,c1,s1)
      call drot(istopm-k,A(k+1,k+1),ldA,A(k+2,k+1),ldA,c2,s2)

      call drot(istopm-k,B(k+2,k+1),ldB,B(k+3,k+1),ldB,c1,s1)
      call drot(istopm-k,B(k+1,k+1),ldB,B(k+2,k+1),ldB,c2,s2)
      if (wantQ) then
         call drot(nQ,Q(1,k+2),1,Q(1,k+3),1,c1,s1)
         call drot(nQ,Q(1,k+1),1,Q(1,k+2),1,c2,s2)
      end if

      end subroutine d_swap_21_inf

      subroutine d_swap_11_inf(wantQ,wantZ,n,k,A,ldA,B,ldB,Q,ldQ,Z,ldZ,
     $   istartm,istopm)
*     Arguments
      logical,intent(in) :: wantQ,wantZ
      integer,intent(in) :: k,ldA,ldB,ldQ,ldZ,istartm,istopm,n
      double precision :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,*)

*     Local variables
      double precision :: c,s,temp

*     Apply transformation from the right
      call dlartg(B(k+1,k+1),B(k+1,k),c,s,temp)
      B(k+1,k+1) = temp
      B(k+1,k) = zero
      call drot(k+2-istartm+1,A(istartm,k+1),1,A(istartm,k),1,c,s)
      call drot(k-istartm+1,B(istartm,k+1),1,B(istartm,k),1,c,s)
      if (wantZ) then
         call drot(n,Z(1,k+1),1,Z(1,k),1,c,s)
      end if

*     Apply transformation from the left
      call dlartg(A(k+1,k),A(k+2,k),c,s,temp)
      A(k+1,k) = temp
      A(k+2,k) = zero
      call drot(istopm-k,A(k+1,k+1),ldA,A(k+2,k+1),ldA,c,s)
      call drot(istopm-k,B(k+1,k+1),ldB,B(k+2,k+1),ldB,c,s)
      if (wantQ) then
         call drot(n,Q(1,k+1),1,Q(1,k+2),1,c,s)
      end if

      end subroutine d_swap_11_inf

      subroutine d_remove_double_shift(wantQ,wantZ,n,ihi,A,ldA,B,ldB,Q,
     $   ldQ,Z,ldZ,istartm,istopm)
*     Arguments
      logical,intent(in) :: wantQ,wantZ
      integer,intent(in) :: ihi,ldA,ldB,ldQ,ldZ,istartm,istopm,n
      double precision :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,*)

*     Local variables
      double precision :: H(2,3),c1,s1,c2,s2,temp

      H = B(ihi-1:ihi,ihi-2:ihi)
*     Make H upper triangular
      call dlartg(H(1,1),H(2,1),c1,s1,temp)
      H(2,1) = zero
      H(1,1) = temp
      call drot(2,H(1,2),2,H(2,2),2,c1,s1)

      call dlartg(H(2,3),H(2,2),c1,s1,temp)
      call drot(1,H(1,3),1,H(1,2),1,c1,s1)
      call dlartg(H(1,2),H(1,1),c2,s2,temp)

      call drot(ihi-istartm+1,B(istartm,ihi),1,B(istartm,ihi-1),1,c1,
     $   s1)
      call drot(ihi-istartm+1,B(istartm,ihi-1),1,B(istartm,ihi-2),1,c2,
     $   s2)
      B(ihi-1,ihi-2) = zero
      B(ihi,ihi-2) = zero
      call drot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c1,
     $   s1)
      call drot(ihi-istartm+1,A(istartm,ihi-1),1,A(istartm,ihi-2),1,c2,
     $   s2)
      if (wantZ) then
         call drot(n,Z(1,ihi),1,Z(1,ihi-1),1,c1,s1)
         call drot(n,Z(1,ihi-1),1,Z(1,ihi-2),1,c2,s2)
      end if

      call dlartg(A(ihi-1,ihi-2),A(ihi,ihi-2),c1,s1,temp)
      A(ihi-1,ihi-2) = temp
      A(ihi,ihi-2) = zero
      call drot(istopm-ihi+2,A(ihi-1,ihi-1),ldA,A(ihi,ihi-1),ldA,c1,s1)
      call drot(istopm-ihi+2,B(ihi-1,ihi-1),ldB,B(ihi,ihi-1),ldB,c1,s1)
      if (wantQ) then
         call drot(n,Q(1,ihi-1),1,Q(1,ihi),1,c1,s1)
      end if

      call dlartg(B(ihi,ihi),B(ihi,ihi-1),c1,s1,temp)
      B(ihi,ihi) = temp
      B(ihi,ihi-1) = zero
      call drot(ihi-istartm,B(istartm,ihi),1,B(istartm,ihi-1),1,c1,s1)
      call drot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c1,
     $   s1)
      if (wantZ) then
         call drot(n,Z(1,ihi),1,Z(1,ihi-1),1,c1,s1)
      end if

      end subroutine d_remove_double_shift

      subroutine d_remove_single_shift(wantZ,n,ihi,A,ldA,B,ldB,Z,ldZ,
     $   istartm)
     
*     Arguments
      logical,intent(in) :: wantZ
      integer,intent(in) :: ihi,ldA,ldB,ldZ,istartm,n
      double precision :: A(ldA,*),B(ldB,*),Z(ldZ,*)

*     Local variables
      double precision :: c,s,temp

      call dlartg(B(ihi,ihi),B(ihi,ihi-1),c,s,temp)
      B(ihi,ihi) = temp
      B(ihi,ihi-1) = zero
      call drot(ihi-istartm,B(istartm,ihi),1,B(istartm,ihi-1),1,c,s)
      call drot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c,s)
      if (wantZ) then
         call drot(n,Z(1,ihi),1,Z(1,ihi-1),1,c,s)
      end if

      end subroutine d_remove_single_shift

      end module d_swaps
