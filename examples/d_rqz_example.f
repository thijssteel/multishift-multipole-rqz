      program test_rqz_full
         use common
         use d_rqz_full
         implicit none

         integer,parameter :: n = 4000
         integer :: i,j,ilo,ihi,lwork,n_aed,n_sweep,n_shifts
         integer(8) :: starttime,endtime,count_rate

         double precision,allocatable :: A(:,:),B(:,:),Q(:,:),Z(:,:),
     $      A_copy(:,:),B_copy(:,:),C(:,:),alphar(:),alphai(:),beta(:),
     $      work(:)
         double precision :: dummywork(1),normA,normB,errorA,errorB

         integer,allocatable :: seed(:)
         integer :: seed_size

         double precision,external :: dlange

         write(*,*) "Example program for complex RQZ"
         write(*,*)

         allocate (A(n,n))
         allocate (B(n,n))
         allocate (Q(n,n))
         allocate (Z(n,n))
         allocate (A_copy(n,n))
         allocate (B_copy(n,n))
         allocate (C(n,n))
         allocate (alphar(n))
         allocate (alphai(n))
         allocate (beta(n))

         call random_seed(size=seed_size)
         allocate (seed(seed_size))

         do i = 1,seed_size
            seed(i) = 0
         end do
         seed(1) = 1
         call random_seed(put=seed)

         write(*,*) "Generating example pencil"

         do j = 1,n
            do i = 1,min(j+1,n)
               call random_number(A(i,j))
            end do
            do i = j+2,n
               A(i,j) = zero
            end do
         end do
         do j = 1,n
            do i = 1,j
               call random_number(B(i,j))
            end do
            do i = j+1,n
               B(i,j) = zero
            end do
         end do

         ilo = 1
         ihi = n

         A_copy = A
         B_copy = B

         allocate (work(n))
         normA = dlange('F',n,n,A_copy,n,work)
         normB = dlange('F',n,n,B_copy,n,work)
         deallocate (work)

         call dlaset('Full',n,n,zero,one,Q,n)
         call dlaset('Full',n,n,zero,one,Z,n)

         write(*,*) "Running d_rqz_f"

         call system_clock(starttime)
         call d_rqz_f(.true.,.true.,.true.,n,ilo,ihi,A,n,B,n,Q,n,Z,n,
     $      alphar,alphai,beta,dummywork,-1,n_aed,n_sweep,n_shifts)
         lwork = int(dummywork(1))
         allocate (work(lwork))
         call d_rqz_f(.true.,.true.,.true.,n,ilo,ihi,A,n,B,n,Q,n,Z,n,
     $      alphar,alphai,beta,work,lwork,n_aed,n_sweep,n_shifts)
         deallocate (work)
         call system_clock(endtime,count_rate)

         write (*,*) "d_rqz_f finished,elapsed time (s):  ",
     $      (endtime-starttime)/real(count_rate,8)

         write (*,*) "calculating backward error"

         call dgemm('C','N',n,n,n,one,Q,n,A_copy,n,zero,C,n)
         A_copy = C
         call dgemm('N','N',n,n,n,one,A_copy,n,Z,n,zero,C,n)
         A_copy = C
         A_copy = A_copy-A

         call dgemm('C','N',n,n,n,one,Q,n,B_copy,n,zero,C,n)
         B_copy = C
         call dgemm('N','N',n,n,n,one,B_copy,n,Z,n,zero,C,n)
         B_copy = C
         B_copy = B_copy-B

         allocate (work(n))
         errorA = dlange('F',n,n,A_copy,n,work)
         errorB = dlange('F',n,n,B_copy,n,work)
         deallocate (work)

         write (*,*) "||Q*AZ-H||/||A||: ",errorA/normA
         write (*,*) "||Q*BZ-T||/||B||: ",errorB/normB

         deallocate (A)
         deallocate (B)
         deallocate (Q)
         deallocate (Z)
         deallocate (A_copy)
         deallocate (B_copy)
         deallocate (C)
         deallocate (alphar)
         deallocate (alphai)
         deallocate (beta)

      end program test_rqz_full
