      module common
         implicit none

         double precision :: zero,one,half
         parameter(zero=0.0d0,one=1.0d0,half=0.5d0)

         double complex :: czero,cone
         parameter(czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0))

*        use_double_pole = 0 => complex poles at the bottom replaced with infinity, real poles at top introduced in double pole block
*        use_double_pole = 1 => use double poles for complex block
*        use_double_pole = 2 => always use double pole blocks
         integer,parameter :: use_double_pole = 1

      contains

      subroutine find_errors(A)
         double precision,intent(in) :: A(:,:)
         integer m,n,i,j

         m = size(A,1)
         n = size(A,2)

         do i = 1,m
            do j = 1,n
               if(A(i,j).gt.1.0E-15) write(*,*) "error",i,j,A(i,j)
            end do
         enddo
         write (*,*)

      end subroutine

      subroutine pretty_print(A)
         double precision,intent(in) :: A(:,:)
         integer m,n,i,j

         m = size(A,1)
         n = size(A,2)

         do i = 1,m
            write(*,"(16ES15.5)") ( A(i,j),j=1,n )
            ! write (*,"(16F6.2)") (A(i,j),j=1,n)
            ! do j=1,n
            !    write (*,"(I8,I8,F8.4)") i,j,A(i,j)
            ! end do
         enddo
         write (*,*)

      end subroutine pretty_print

      subroutine z_pretty_print(A)
         double complex,intent(in) :: A(:,:)
         integer m,n,i,j

         m = size(A,1)
         n = size(A,2)

         do i = 1,m
       write(*,"(16ES8.4)") ( A(i,j),j=1,n )
            ! write (*,"(20F6.3)") (A(i,j),j=1,n)
         enddo
         write (*,*)

      end subroutine z_pretty_print

      subroutine my_drot(N,DX,INCX,DY,INCY,C,S)

         double precision C,S
         integer INCX,INCY,N

         double precision DX(*),DY(*)

         double precision DTEMP
         integer I,IX,IY

         if (n .LE. 0) RETURN
         if (incx .EQ. 1 .AND. incy .EQ. 1) then
            do i = 1,n
               dtemp = c*dx(i)+s*dy(i)
               dy(i) = c*dy(i)-s*dx(i)
               dx(i) = dtemp
            end do
         ELSE
            ix = 1
            iy = 1
            if (incx .LT. 0) ix = (-n+1)*incx+1
            if (incy .LT. 0) iy = (-n+1)*incy+1
            do i = 1,n
               dtemp = c*dx(ix)+s*dy(iy)
               dy(iy) = c*dy(iy)-s*dx(ix)
               dx(ix) = dtemp
               ix = ix+incx
               iy = iy+incy
            end do
         end if
         RETURN

      end subroutine my_drot

      subroutine shiftcolumn(A,ldA,B,ldB,sr1,sr2,si,beta,v)

!        Arguments
         integer,intent(in) :: ldA,ldB
         double precision,intent(in) :: A(ldA,*),B(ldB,*)
         double precision,intent(out) :: sr1,sr2,si,beta,v(*)

!        Local scalars
         double precision :: w(2)

!        Calculate first shifted vector
         w(1) = beta*A(1,1)-sr1*B(1,1)
         w(2) = beta*A(2,1)-sr1*B(2,1)

!        Solve linear system
         w(2) = w(2)/B(2,2)
         w(1) = (w(1)-B(1,2)*w(2))/B(1,1)

!        Apply second shift
         v(1) = beta*(A(1,1)*w(1)+A(1,2)*w(2))-sr2*(B(1,1)*w(1)+B(1,
     $      2)*w(2))
         v(2) = beta*(A(2,1)*w(1)+A(2,2)*w(2))-sr2*(B(2,1)*w(1)+B(2,
     $      2)*w(2))
         v(3) = beta*(A(3,1)*w(1)+A(3,2)*w(2))-sr2*(B(3,1)*w(1)+B(3,
     $      2)*w(2))

!        Account for imaginary part
         v(1) = v(1)+si*si*B(1,1)

      end subroutine shiftcolumn

      subroutine polecolumn(A,ldA,B,ldB,sr1,sr2,si,beta,v)

!        Arguments
         integer,intent(in) :: ldA,ldB
         double precision,intent(in) :: A(ldA,*),B(ldB,*),sr1,sr2,si,
     $      beta
         double precision,intent(out) :: v(*)

!        Local scalars
         double precision :: w(2)

!        Calculate first shifted vector
         w(1) = beta*A(3,2)-sr1*B(3,2)
         w(2) = beta*A(3,3)-sr1*B(3,3)

!        Solve linear system
         w(1) = w(1)/B(2,2)
         w(2) = (w(2)-B(2,3)*w(1))/B(3,3)

!        Apply second shift
         v(1) = beta*A(2,1)*w(1)
         v(2) = beta*(A(2,2)*w(1)+A(3,2)*w(2))-sr2*(B(2,2)*w(1))
         v(3) = beta*(A(2,3)*w(1)+A(3,3)*w(2))-sr2*(B(2,3)*w(1)+B(3,
     $      3)*w(2))

!        Account for imaginary part
         v(3) = v(3)+si*si*B(3,3)

         if(isnan(v(1)) .or. isnan(v(2)) .or. isnan(v(3))) then
            v(1) = one
            v(2) = zero
            v(3) = zero
         end if

      end subroutine polecolumn

      subroutine eig22(A,B,alphar1,alphar2,alphai,beta)

!        Arguments
         double precision,intent(in) :: A(2,2),B(2,2)
         double precision,intent(out) :: alphar1,alphar2,alphai,beta

!        Local scalars
         double precision :: detb,ab11,ab12,ab21,ab22,t,d

         detb = B(1,1)*B(2,2)-B(1,2)*B(2,1)

         ab11 = A(1,1)*B(2,2)-A(2,1)*B(1,2)
         ab12 = A(1,2)*B(2,2)-A(2,2)*B(1,2)
         ab21 = A(2,1)*B(1,1)-A(1,1)*B(2,1)
         ab22 = A(2,2)*B(1,1)-A(1,2)*B(2,1)

         t = half*(ab11+ab22)
         d = t**2+ab12*ab21-ab11*ab22

         beta = detb
         if (d .ge. zero) then
            d = sqrt(d)
            alphar1 = t+d
            alphar2 = t-d
            alphai = zero
         else
            d = sqrt(-d)
            alphar1 = t
            alphar2 = t
            alphai = d
         end if

      end subroutine eig22

      subroutine z_eig22(A,ldA,B,ldB,alpha1,alpha2,beta)

!        Arguments
         integer,intent(in) :: ldA,ldB
         double complex,intent(in) :: A(ldA,*),B(ldB,*)
         double complex,intent(out) :: alpha1,alpha2,beta

!        Local scalars
         double complex :: detb,ab11,ab12,ab21,ab22,t,d

         detb = B(1,1)*B(2,2)-B(1,2)*B(2,1)

         ab11 = A(1,1)*B(2,2)-A(2,1)*B(1,2)
         ab12 = A(1,2)*B(2,2)-A(2,2)*B(1,2)
         ab21 = A(2,1)*B(1,1)-A(1,1)*B(2,1)
         ab22 = A(2,2)*B(1,1)-A(1,2)*B(2,1)

         t = (0.5d0,0.0d0)*(ab11+ab22)
         d = t**2+ab12*ab21-ab11*ab22

         beta = detb
         d = sqrt(d)
         alpha1 = t+d
         alpha2 = t-d

      end subroutine z_eig22

      subroutine d_apply_ct_l(v1,v2,c,s)
! DESCRIPTION
!___________________________________________________________________________
! Applies an elementary rotation from the left to two row
! vectors v1 and v2 of equal lenght :
! [ c  s ]    [ v1  ]
! [      ] *  [     ]
! [-s  c ]    [ v2  ]
!
! ARGUMENTS
!___________________________________________________________________________
! v1        double array [INOUT]
!             First row vector, on output equal to c * v1 + s * v2
! v2        double array [INOUT]
!             Second row vector, on output equal to -s * v1 + c * v2
! c         double [IN]
!             Cosine of rotation
! s         double [IN]
!             Sine of rotation
!___________________________________________________________________________
! last edit: September 5, 2018
         double precision,intent(inout)  :: v1(:),v2(:)
         double precision,intent(in)     :: c,s

         integer                       :: i,n
         double precision                 :: t

         n = size(v1)
         do i = 1,n
            t = c*v1(i)+s*v2(i)
            v2(i) =-s*v1(i)+c*v2(i)
            v1(i) = t
         end do
      end subroutine

      subroutine d_apply_ct_r(v1,v2,c,s)
! DESCRIPTION
!___________________________________________________________________________
! Applies an elementary rotation from the left to two column
! vectors v1 and v2 of equal lenght :
! [     ]   [ c  s ]
! [v1 v2] * [      ]
! [     ]   [-s  c ]
!
! ARGUMENTS
!___________________________________________________________________________
! v1        double array [INOUT]
!             First column vector, on output equal to c * v1 - s * v2
! v2        double array [INOUT]
!             Second row vector, on output equal to s * v1 + c * v2
! c         double [IN]
!             Cosine of rotation
! s         double [IN]
!             Sine of rotation
!___________________________________________________________________________
! last edit: September 5, 2018
         double precision,intent(inout)  :: v1(:),v2(:)
         double precision,intent(in)     :: c,s

         integer                       :: i,n
         double precision                 :: t

         n = size(v1)
         do i = 1,n
            t = c*v1(i)-s*v2(i)
            v2(i) = s*v1(i)+c*v2(i)
            v1(i) = t
         end do
         end subroutine

      end module common
