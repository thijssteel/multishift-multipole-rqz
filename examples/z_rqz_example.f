      program test_rqz_full
         use common
         use z_rqz_full
         implicit none

         integer,parameter :: n = 4000
         integer :: i,j,ilo,ihi,lwork,n_aed,n_sweep,n_shifts
         integer(8) :: starttime,endtime,count_rate

         double complex,allocatable :: A(:,:),B(:,:),Q(:,:),Z(:,:),
     $      A_copy(:,:),B_copy(:,:),C(:,:),alpha(:),beta(:),work(:)
         double precision,allocatable :: rwork(:)

         double complex :: dummywork(1)
         double precision :: r1,r2,normA,normB,errorA,errorB

         double precision,external :: zlange

         write(*,*) "Example program for complex RQZ"
         write(*,*)

         allocate (A(n,n))
         allocate (B(n,n))
         allocate (Q(n,n))
         allocate (Z(n,n))
         allocate (A_copy(n,n))
         allocate (B_copy(n,n))
         allocate (C(n,n))
         allocate (alpha(n))
         allocate (beta(n))

         call zlaset('Full',n,n,czero,czero,A,n)
         call zlaset('Full',n,n,czero,czero,B,n)

         write(*,*) "Generating example pencil"

         do j = 1,n
            do i = 1,min(j+1,n)
               call random_number(r1)
               A(i,j) = r1
            end do
         end do
         do j = 1,n
            do i = 1,j
               call random_number(r1)
               B(i,j) = r1
            end do
         end do

         ilo = 1
         ihi = n

         A_copy = A
         B_copy = B

         allocate (rwork(n))
         normA = zlange('F',n,n,A_copy,n,rwork)
         normB = zlange('F',n,n,B_copy,n,rwork)
         deallocate (rwork)

         write(*,*) "Running z_rqz_f"

         call zlaset('Full',n,n,czero,cone,Q,n)
         call zlaset('Full',n,n,czero,cone,Z,n)

         call system_clock(starttime)
*        Workspace query
         call z_rqz_f(.true.,.true.,.true.,n,ilo,ihi,A,n,B,n,Q,n,Z,n,
     $      alpha,beta,dummywork,-1,n_aed,n_sweep,n_shifts)
         lwork = int(dummywork(1))
         allocate (work(lwork))
*        Run the algorithm
         call z_rqz_f(.true.,.true.,.true.,n,ilo,ihi,A,n,B,n,Q,n,Z,n,
     $      alpha,beta,work,lwork,n_aed,n_sweep,n_shifts)
         deallocate (work)
         call system_clock(endtime,count_rate)

         write (*,*) "z_rqz_f finished,elapsed time (s):  ",
     $      (endtime-starttime)/real(count_rate,8)

         write (*,*) "calculating backward error"

*        Calculate backward error
         call zgemm('C','N',n,n,n,cone,Q,n,A_copy,n,czero,C,n)
         A_copy = C
         call zgemm('N','N',n,n,n,cone,A_copy,n,Z,n,czero,C,n)
         A_copy = C
         A_copy = A_copy-A

         call zgemm('C','N',n,n,n,cone,Q,n,B_copy,n,czero,C,n)
         B_copy = C
         call zgemm('N','N',n,n,n,cone,B_copy,n,Z,n,czero,C,n)
         B_copy = C
         B_copy = B_copy-B

         allocate (rwork(n))
         errorA = zlange('F',n,n,A_copy,n,rwork)
         errorB = zlange('F',n,n,B_copy,n,rwork)
         deallocate (rwork)

         write (*,*) "||Q*AZ-H||/||A||: ",errorA/normA
         write (*,*) "||Q*BZ-T||/||B||: ",errorB/normB

         deallocate (A)
         deallocate (B)
         deallocate (Q)
         deallocate (Z)
         deallocate (A_copy)
         deallocate (B_copy)
         deallocate (C)
         deallocate (alpha)
         deallocate (beta)

      end program test_rqz_full
