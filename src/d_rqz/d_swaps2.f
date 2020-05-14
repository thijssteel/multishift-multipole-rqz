      module d_swaps2
         use common
         implicit none

      contains
      
      subroutine d_swap(istartm,istopm,wantQ,wantZ,k,n1,n2,off,A,ldA,B,
     $   ldB,qn,qStart,Q,ldQ,zn,zStart,Z,ldZ,work,lwork,info)

*        Arguments
         integer,intent(in) :: istartm,istopm,k,ldA,ldB,ldQ,ldZ,qStart,
     $      zStart,qn,zn,off,n1,n2,lwork
         logical,intent(in) :: wantQ,wantZ
         double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $      Z(ldZ,*),work(*)
         integer,intent(out) :: info

*        Local variables
         logical :: ab
         double precision :: temp,temp2,temp3,c1,s1,c2,s2,c3,s3,c4,s4,
     $      H(2,3),H2(3,2),X(2,2),Y(2,2),scale,rdsum,rdscal,U(2,2),Vt(2,
     $      2),sval(2),Qc(4,4),Zc(4,4),Ac(4,4),Bc(4,4),ulp,M(8,8),XY(8),
     $      Qc2(4,4),Zc2(4,4)
         integer :: iwork(6),pq,sylinfo,i,ipiv(8),iiter

*        External Functions
         double precision,external :: dlamch
         double precision,external :: dlange

*        Executable statements

         info = 0

         if(n1 .eq. 1 .and. n2 .eq. 1) then

            ab = abs(A(k+off,k))*abs(B(k+off+1,k+1)) .gt.abs(A(k+off+1,
     $         k+1))*abs(B(k+off,k))

            temp2 = B(k+off+1,k+1)*A(k+off,k+1)-A(k+off+1,k+1)*B(k+off,
     $         k+1)
            temp3 = B(k+off+1,k+1)*A(k+off,k)-A(k+off+1,k+1)*B(k+off,k)

            call dlartg(temp2,temp3,c1,s1,temp)
            call drot(k+1+off-istartm+1,A(istartm,k+1),1,A(istartm,k),1,
     $         c1,s1)
            call drot(k+1+off-istartm+1,B(istartm,k+1),1,B(istartm,k),1,
     $         c1,s1)
            if (wantZ) then
               call drot(zn,Z(1,k+1-zStart+1),1,Z(1,k-zStart+1),1,c1,s1)
            end if

            if (ab) then
               call dlartg(B(k+off,k),B(k+1+off,k),c1,s1,temp)
            else
               call dlartg(A(k+off,k),A(k+off+1,k),c1,s1,temp)
            end if

            call drot(istopm-k+1,A(k+off,k),ldA,A(k+1+off,k),ldA,c1,s1)
            call drot(istopm-k+1,B(k+off,k),ldB,B(k+1+off,k),ldB,c1,s1)
            if (wantQ) then
               call drot(qn,Q(1,k+off-qStart+1),1,Q(1,k+1+off-qStart+1),
     $            1,c1,s1)
            end if

            A(k+1+off,k) = zero
            B(k+1+off,k) = zero

         else if (n1 .eq. 2 .and. n2 .eq. 1) then

*           Original criterion, replacing with eigenvalue check could
*           lead to better results
            ab = abs(A(k+2+off,k+2)) .ge. abs(B(k+2+off,k+2))

*           If required, rotate to make A upper triangular
            if(ab) then
               call dlartg(A(k+off,k),A(k+1+off,k),c1,s1,temp)
               A(k+off,k) = temp
               A(k+1+off,k) = zero

               call drot(istopm-(k+1)+1,A(k+off,k+1),ldA,A(k+1+off,k+1),
     $            ldA,c1,s1)
               call drot(istopm-k+1,B(k+off,k),ldB,B(k+1+off,k),ldB,c1,
     $            s1)
               if(wantQ) then
                  call drot(qn,Q(1,k+off-qStart+1),1,Q(1,k+1+off-qStart+
     $               1),1,c1,s1)
               end if
            end if


*           Calculate H and make it upper triangular
            H = B(k+2+off,k+2)*A(k+off:k+1+off,k:k+2)-A(k+2+off,
     $         k+2)*B(k+off:k+1+off,k:k+2)
            call dlartg(H(1,1),H(2,1),c1,s1,temp)
            H(2,1) = zero
            H(1,1) = temp
            call drot(2,H(1,2),2,H(2,2),2,c1,s1)

*           Calculate Z1 and Z2
            call dlartg(H(2,3),H(2,2),c1,s1,temp)
            call drot(1,H(1,3),1,H(1,2),1,c1,s1)
            call dlartg(H(1,2),H(1,1),c2,s2,temp)

*           Apply Z1
            call drot(k+2+off-istartm+1,A(istartm,k+2),1,A(istartm,k+1),
     $         1,c1,s1)
            call drot(k+2+off-istartm+1,B(istartm,k+2),1,B(istartm,k+1),
     $         1,c1,s1)
            if (wantZ) then
               call drot(zn,Z(1,k+2-zStart+1),1,Z(1,k+1-zStart+1),1,c1,
     $            s1)
            end if

*           Calculate and apply Q1
            if(ab) then
               call dlartg(A(k+1+off,k+1),A(k+2+off,k+1),c1,s1,temp)
            else
               call dlartg(B(k+1+off,k+1),B(k+2+off,k+1),c1,s1,temp)
            end if
            call drot(istopm-k+1,A(k+1+off,k),ldA,A(k+2+off,k),ldA,c1,
     $         s1)
            call drot(istopm-k+1,B(k+1+off,k),ldB,B(k+2+off,k),ldB,c1,
     $         s1)
            if(wantQ) then
               call drot(qn,Q(1,k+1+off-qStart+1),1,Q(1,k+2+off-qStart+
     $            1),1,c1,s1)
            end if

*           Apply Z2
            call drot(k+2+off-istartm+1,B(istartm,k+1),1,B(istartm,k),1,
     $         c2,s2)
            call drot(k+2+off-istartm+1,A(istartm,k+1),1,A(istartm,k),1,
     $         c2,s2)
            if (wantZ) then
               call drot(zn,Z(1,k+1-zStart+1),1,Z(1,k-zStart+1),1,c2,s2)
            end if

*           Calculate and apply Q2
            if(ab) then
               call dlartg(A(k+off,k),A(k+1+off,k),c1,s1,temp)
            else
               call dlartg(B(k+off,k),B(k+1+off,k),c1,s1,temp)
            end if
            call drot(istopm-k+1,A(k+off,k),ldA,A(k+1+off,k),ldA,c1,s1)
            call drot(istopm-k+1,B(k+off,k),ldB,B(k+1+off,k),ldB,c1,s1)
            if(wantQ) then
               call drot(qn,Q(1,k+off-qStart+1),1,Q(1,k+1+off-qStart+1),
     $            1,c1,s1)
            end if

*           If required, return B to standard form    
            if(ab) then
               call dlartg(B(k+1+off,k+1),B(k+2+off,k+1),c1,s1,temp)
               B(k+1+off,k+1) = temp
               B(k+2+off,k+1) = zero

               call drot(istopm-(k+1)+1,A(k+1+off,k+1),ldA,A(k+2+off,
     $            k+1),ldA,c1,s1)
               call drot(istopm-(k+2)+1,B(k+1+off,k+2),ldB,B(k+2+off,
     $            k+2),ldB,c1,s1)
               if(wantQ) then
                  call drot(qn,Q(1,k+1+off-qStart+1),1,Q(1,
     $               k+2+off-qStart+1),1,c1,s1)
               end if
            end if

            A(k+1+off,k) = zero
            A(k+2+off,k) = zero
            B(k+1+off,k) = zero
            B(k+2+off,k) = zero
            B(k+2+off,k+1) = zero

         else if (n1 .eq. 1 .and. n2 .eq. 2) then

*           Original criterion, replacing with eigenvalue check could
*           lead to better results
            ab = abs(A(k+off,k)) .ge. abs(B(k+off,k))

*           If required, make A upper triangular    
            if(ab) then
               call dlartg(A(k+1+off,k+1),A(k+2+off,k+1),c1,s1,temp)
               A(k+1+off,k+1) = temp
               A(k+2+off,k+1) = zero

               call drot(istopm-(k+2)+1,A(k+1+off,k+2),ldA,A(k+2+off,
     $            k+2),ldA,c1,s1)
               call drot(istopm-(k+1)+1,B(k+1+off,k+1),ldB,B(k+2+off,
     $            k+1),ldB,c1,s1)
               if(wantQ) then
                  call drot(qn,Q(1,k+1+off-qStart+1),1,Q(1,
     $               k+2+off-qStart+1),1,c1,s1)
               end if
            end if
            
*           Calculate H and make it upper triangular
            H2 = B(k+off,k)*A(k+off:k+2+off,k+1:k+2)-A(k+off,
     $         k)*B(k+off:k+2+off,k+1:k+2)
            call dlartg(H2(3,2),H2(3,1),c1,s1,temp)
            H2(3,1) = zero
            H2(3,2) = temp
            call drot(2,H2(1,2),1,H2(1,1),1,c1,s1)

            
*           Calculate Q1 and Q2
            call dlartg(H2(1,1),H2(2,1),c1,s1,temp)
            call drot(1,H2(1,2),3,H2(2,2),3,c1,s1)
            call dlartg(H2(2,2),H2(3,2),c2,s2,temp)

*           Apply Q1
            call drot(istopm-k+1,A(k+off,k),ldA,A(k+1+off,k),ldA,c1,s1)
            call drot(istopm-k+1,B(k+off,k),ldB,B(k+1+off,k),ldB,c1,s1)
            if (wantQ) then
               call drot(qn,Q(1,k+off-qStart+1),1,Q(1,k+1+off-qStart+1),
     $            1,c1,s1)
            end if

*           Calculate and apply Z1
            if(ab) then
               call dlartg(A(k+1+off,k+1),A(k+1+off,k),c1,s1,temp)
            else
               call dlartg(B(k+1+off,k+1),B(k+1+off,k),c1,s1,temp)
            end if
            call drot(k+2+off-istartm+1,B(istartm,k+1),1,B(istartm,k),1,
     $         c1,s1)
            call drot(k+2+off-istartm+1,A(istartm,k+1),1,A(istartm,k),1,
     $         c1,s1)
            if (wantZ) then
               call drot(zn,Z(1,k+1-zStart+1),1,Z(1,k-zStart+1),1,c1,s1)
            end if
            
*           Apply Q2
            call drot(istopm-k+1,A(k+1+off,k),ldA,A(k+2+off,k),ldA,c2,
     $         s2)
            call drot(istopm-k+1,B(k+1+off,k),ldB,B(k+2+off,k),ldB,c2,
     $         s2)
            if(wantQ) then
               call drot(qn,Q(1,k+1+off-qStart+1),1,Q(1,k+2+off-qStart+
     $            1),1,c2,s2)
            end if

*           Calculate and apply Z2
            if(ab) then
               call dlartg(A(k+2+off,k+2),A(k+2+off,k+1),c1,s1,temp)
            else
               call dlartg(B(k+2+off,k+2),B(k+2+off,k+1),c1,s1,temp)
            end if
            call drot(k+2+off-istartm+1,A(istartm,k+2),1,A(istartm,k+1),
     $         1,c1,s1)
            call drot(k+2+off-istartm+1,B(istartm,k+2),1,B(istartm,k+1),
     $         1,c1,s1)
            if (wantZ) then
               call drot(zn,Z(1,k+2-zStart+1),1,Z(1,k+1-zStart+1),1,c1,
     $            s1)
            end if

*           If required, return B to standard form  
            if(ab) then
               call dlartg(B(k+off,k),B(k+1+off,k),c1,s1,temp)
               B(k+off,k) = temp
               B(k+1+off,k) = zero

               call drot(istopm-(k+1)+1,B(k+off,k+1),ldB,B(k+1+off,k+1),
     $            ldB,c1,s1)
               call drot(istopm-k+1,A(k+off,k),ldA,A(k+1+off,k),ldA,c1,
     $            s1)
               if(wantQ) then
                  call drot(qn,Q(1,k+off-qStart+1),1,Q(1,k+1+off-qStart+
     $               1),1,c1,s1)
               end if
            end if

            A(k+off+2,k) = zero
            A(k+off+2,k+1) = zero
            B(k+off+2,k) = zero
            B(k+off+2,k+1) = zero
            B(k+off+1,k) = zero

         else if (n1 .eq. 2 .and. n2 .eq. 2) then

            ulp = dlamch('precision')

            Ac = A(k+off:k+off+3,k:k+3)
            Bc = B(k+off:k+off+3,k:k+3)

            call dlaset('All',4,4,zero,one,Qc,4)
            call dlaset('All',4,4,zero,one,Zc,4)

            M = zero
            ! I (x) A11
            M(1:2,1:2) = Ac(1:2,1:2)
            M(3:4,3:4) = Ac(1:2,1:2)
            ! I (x) B11
            M(5:6,1:2) = Bc(1:2,1:2)
            M(7:8,3:4) = Bc(1:2,1:2)
            ! A22T (x) I
            M(1,5) =-Ac(3,3)
            M(1,6) =-Ac(4,3)
            M(2,7) =-Ac(3,3)
            M(2,8) =-Ac(4,3)
            M(3,5) =-Ac(3,4)
            M(3,6) =-Ac(4,4)
            M(4,7) =-Ac(3,4)
            M(4,8) =-Ac(4,4)
            ! B22T (x) I
            M(5,5) =-Bc(3,3)
            M(5,6) =-Bc(4,3)
            M(6,7) =-Bc(3,3)
            M(6,8) =-Bc(4,3)
            M(7,5) =-Bc(3,4)
            M(7,6) =-Bc(4,4)
            M(8,7) =-Bc(3,4)
            M(8,8) =-Bc(4,4)

            ! Create RHS
            XY(1:2) =Ac(1:2,3)
            XY(3:4) =Ac(1:2,4)
            XY(5:6) =Bc(1:2,3)
            XY(7:8) =Bc(1:2,4)

            call dgesv(8,1,M,8,ipiv,XY,8,sylinfo)

            if (sylinfo .ne. 0) then
               info = 1
               ! write(*,*) "swap failed"
               return
            end if

            X(1,1) = XY(1)
            X(2,1) = XY(2)
            X(1,2) = XY(3)
            X(2,2) = XY(4)
            Y(1,1) = XY(5)
            Y(1,2) = XY(6)
            Y(2,1) = XY(7)
            Y(2,2) = XY(8)

            if (sylinfo .ne. 0) then
               info = 1
               ! write(*,*) "swap failed,solving sylvester system"
               return
            end if
            Y = transpose(Y)

*           Calculate and apply Z
            call dgesvd('A','A',2,2,X,2,sval,U,2,Vt,2,work,lwork,
     $         sylinfo)

            if(sval(1) .lt. one) then
               s1 = one/sqrt(one+sval(1)**2)
               c1 = s1*sval(1)
            else
               c1 = one/sqrt(one+one/sval(1)**2)
               s1 = c1/sval(1)
            end if
            if(sval(2) .lt. one) then
               s2 = one/sqrt(one+sval(2)**2)
               c2 = s2*sval(2)
            else
               c2 = one/sqrt(one+one/sval(2)**2)
               s2 = c2/sval(2)
            end if

            Zc(1:2,1:2) = U
            Zc(3:4,3:4) = transpose(Vt)
            call drot(4,Zc(1,1),1,Zc(1,3),1,c1,-s1)
            call drot(4,Zc(1,2),1,Zc(1,4),1,c2,-s2)

            call dgemm('N','N',4,4,4,one,Ac,4,Zc,4,zero,work,4)
            call dlacpy('ALL',4,4,work,4,Ac,4)
            call dgemm('N','N',4,4,4,one,Bc,4,Zc,4,zero,work,4)
            call dlacpy('ALL',4,4,work,4,Bc,4)

*           Calculate and apply Q
            call dgesvd('A','A',2,2,Y,2,sval,U,2,Vt,2,work,lwork,
     $         sylinfo)

            if(sval(1) .lt. one) then
               s1 = one/sqrt(one+sval(1)**2)
               c1 = s1*sval(1)
            else
               c1 = one/sqrt(one+one/sval(1)**2)
               s1 = c1/sval(1)
            end if
            if(sval(2) .lt. one) then
               s2 = one/sqrt(one+sval(2)**2)
               c2 = s2*sval(2)
            else
               c2 = one/sqrt(one+one/sval(2)**2)
               s2 = c2/sval(2)
            end if

            Qc(1:2,1:2) = transpose(Vt)
            Qc(3:4,3:4) = U
            call drot(4,Qc(1,1),1,Qc(1,3),1,c1,-s1)
            call drot(4,Qc(1,2),1,Qc(1,4),1,c2,-s2)

            call dgemm('T','N',4,4,4,one,Qc,4,Ac,4,zero,work,4)
            call dlacpy('ALL',4,4,work,4,Ac,4)
            call dgemm('T','N',4,4,4,one,Qc,4,Bc,4,zero,work,4)
            call dlacpy('ALL',4,4,work,4,Bc,4)

            iiter = 0
            do while(dlange('F',2,2,Ac(3,1),4,work) .gt. ulp*dlange('F',
     $         4,4,A(k+off,k),ldA,work) .or. dlange('F',2,2,Bc(3,1),4,
     $         work).gt.ulp*dlange('F',4,4,B(k+off,k),ldB,work))
               if(iiter .ge. 3) then
                  info = 1
   !                write(*,*) "swap failed",dlange('F',2,2,Ac(3,1),4,
   !   $               work)/dlange('F',4,4,A(k+off,k),ldA,work),
   !   $               dlange('F',2,2,Bc(3,1),4,work)/dlange('F',4,4,
   !   $               B(k+off,k),ldB,work)
                  return
               end if

               call dlaset('All',4,4,zero,one,Qc2,4)
               call dlaset('All',4,4,zero,one,Zc2,4)

               M = zero
               M(1:2,1:2) = Ac(3:4,3:4)
               M(3:4,3:4) = Ac(3:4,3:4)
               M(1,5) =-Ac(1,1)
               M(1,6) =-Ac(2,1)
               M(2,7) =-Ac(1,1)
               M(2,8) =-Ac(2,1)
               M(3,5) =-Ac(1,2)
               M(3,6) =-Ac(2,2)
               M(4,7) =-Ac(1,2)
               M(4,8) =-Ac(2,2)
               M(5:6,1:2) = Bc(3:4,3:4)
               M(7:8,3:4) = Bc(3:4,3:4)
               M(5,5) =-Bc(1,1)
               M(5,6) =-Bc(2,1)
               M(6,7) =-Bc(1,1)
               M(6,8) =-Bc(2,1)
               M(7,5) =-Bc(1,2)
               M(7,6) =-Bc(2,2)
               M(8,7) =-Bc(1,2)
               M(8,8) =-Bc(2,2)

               ! Create RHS
               XY(1:2) =Ac(3:4,1)
               XY(3:4) =Ac(3:4,2)
               XY(5:6) =Bc(3:4,1)
               XY(7:8) =Bc(3:4,2)

               call dgesv(8,1,M,8,ipiv,XY,8,sylinfo)
               if (sylinfo .ne. 0) then
                  info = 1
   !                write(*,*) "swap failed",dlange('F',2,2,Ac(3,1),4,
   !   $               work)/dlange('F',4,4,A(k+off,k),ldA,work),
   !   $               dlange('F',2,2,Bc(3,1),4,work)/dlange('F',4,4,
   !   $               B(k+off,k),ldB,work)
                  return
               end if

               X(1,1) = XY(1)
               X(2,1) = XY(2)
               X(1,2) = XY(3)
               X(2,2) = XY(4)
               Y(1,1) = XY(5)
               Y(1,2) = XY(6)
               Y(2,1) = XY(7)
               Y(2,2) = XY(8)

*              Calculate and apply Z
               call dgesvd('A','A',2,2,X,2,sval,U,2,Vt,2,work,lwork,
     $            sylinfo)

               call dlartg(one,-sval(1),c1,s1,temp)
               call dlartg(one,-sval(2),c2,s2,temp)

               Zc2(1:2,1:2) = transpose(Vt)
               Zc2(3:4,3:4) = U
               call drot(4,Zc2(1,1),1,Zc2(1,3),1,c1,s1)
               call drot(4,Zc2(1,2),1,Zc2(1,4),1,c2,s2)

               call dgemm('N','N',4,4,4,one,Ac,4,Zc2,4,zero,work,4)
               call dlacpy('ALL',4,4,work,4,Ac,4)
               call dgemm('N','N',4,4,4,one,Bc,4,Zc2,4,zero,work,4)
               call dlacpy('ALL',4,4,work,4,Bc,4)
               call dgemm('N','N',4,4,4,one,Zc,4,Zc2,4,zero,work,4)
               call dlacpy('ALL',4,4,work,4,Zc,4)

*              Calculate and apply Q
               call dgesvd('A','A',2,2,Y,2,sval,U,2,Vt,2,work,lwork,
     $            sylinfo)

               call dlartg(one,-sval(1),c1,s1,temp)
               call dlartg(one,-sval(2),c2,s2,temp)

               Qc2(1:2,1:2) = transpose(Vt)
               Qc2(3:4,3:4) = U
               call drot(4,Qc2(1,1),1,Qc2(1,3),1,c1,s1)
               call drot(4,Qc2(1,2),1,Qc2(1,4),1,c2,s2)

               call dgemm('T','N',4,4,4,one,Qc2,4,Ac,4,zero,work,4)
               call dlacpy('ALL',4,4,work,4,Ac,4)
               call dgemm('T','N',4,4,4,one,Qc2,4,Bc,4,zero,work,4)
               call dlacpy('ALL',4,4,work,4,Bc,4)
               call dgemm('N','N',4,4,4,one,Qc,4,Qc2,4,zero,work,4)
               call dlacpy('ALL',4,4,work,4,Qc,4)

               iiter = iiter+1
            end do

            call dlartg(Bc(1,1),Bc(2,1),c1,s1,temp)
            call dlartg(Bc(3,3),Bc(4,3),c2,s2,temp)
            call drot(4,Qc(1,1),1,Qc(1,2),1,c1,s1)
            call drot(4,Qc(1,3),1,Qc(1,4),1,c2,s2)

*           Apply updates to the rest of the pencil

            call dgemm('T','N',4,istopm-k+1,4,one,Qc,4,A(k+off,k),ldA,
     $         zero,work,4)
            call dlacpy('ALL',4,istopm-k+1,work,4,A(k+off,k),ldA)
            call dgemm('T','N',4,istopm-k+1,4,one,Qc,4,B(k+off,k),ldB,
     $         zero,work,4)
            call dlacpy('ALL',4,istopm-k+1,work,4,B(k+off,k),ldB)
            if (wantQ) then
               call dgemm('N','N',qn,4,4,one,Q(1,k+off-qStart+1),ldQ,Qc,
     $            4,zero,work,qn)
               call dlacpy('ALL',qn,4,work,qn,Q(1,k+off-qStart+1),ldQ)
            end if

            call dgemm('N','N',k+off+3-istartm+1,4,4,one,A(istartm,k),
     $         ldA,Zc,4,zero,work,k+off+3-istartm+1)
            call dlacpy('ALL',k+off+3-istartm+1,4,work,k+off+3-istartm+
     $         1,A(istartm,k),ldA)
            call dgemm('N','N',k+off+3-istartm+1,4,4,one,B(istartm,k),
     $         ldB,Zc,4,zero,work,k+off+3-istartm+1)
            call dlacpy('ALL',k+off+3-istartm+1,4,work,k+off+3-istartm+
     $         1,B(istartm,k),ldB)
            if (wantZ) then
               call dgemm('N','N',zn,4,4,one,Z(1,k-zStart+1),ldZ,Zc,4,
     $            zero,work,zn)
               call dlacpy('ALL',zn,4,work,zn,Z(1,k-zStart+1),ldZ)
            end if

            A(k+off+2:k+off+3,k:k+1) = zero
            B(k+off+2:k+off+3,k:k+1) = zero
            B(k+off+1,k) = zero
            B(k+off+3,k+2) = zero
            
            
         end if

      end subroutine d_swap

      end module