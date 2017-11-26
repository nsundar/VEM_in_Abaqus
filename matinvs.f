	 contains 
! -------------------------------------------------
!                inv
! -------------------------------------------------
! Matrix Inversion
! Ax=b
! PA = LU or A=P'LU
! P'LUx=b
! LUx=Pb
! Solve Ld=Pb using Forward  sub where d=Ux
! Solve Ux=d using Backward sub interface inv
          function inv(A) result(Ainv)

  	   real*8, intent(in), dimension(:,:) :: A
  	   real*8, dimension(size(A,1),size(A,1)) :: Ainv

   	  integer :: n
  	  integer :: i,j,k,bb ! Running variables

! 	 Variables for calculating Permutation matrix
  	  real*8 :: max_elem
 	   real*8, dimension(size(A,1),size(A,1)) :: A_dummy
  	  real*8, dimension(size(A,1)) :: P_swap, A_swap

!	 Variables for LU Decomposition
  	  real*8, dimension(size(A,1),size(A,1)) :: L,U,P
  	  real*8, dimension(size(A,1)) :: Pb,d,x,bvec
   	  real*8 :: sumu, suml, det,eps

  	  n=size(A,1)
  	  A_dummy=A

!       Find Permutation Matrix for partial pivoting
!       Creating P as Identity Matrix
  	   P=0.d0
  	   do i=1,n
 	     P(i,i)=1.
  	   enddo

 	   do j=1,n
   	   max_elem=maxval(A_dummy(j:n,j))
   	   do i=j,n
    	    if (A(i,j) .eq. max_elem) then
           P_swap=P(i,:)
           P(i,:)=P(j,:)
           P(j,:)=P_swap

           A_swap=A_dummy(i,:)
           A_dummy(i,:)=A_dummy(j,:)
           A_dummy(j,:)=A_swap
           exit
   	     endif
   	   enddo
   	  enddo

!	 LU decomposition using Doolittle algorithm on PA
!         A_dummy is now P*A
 	   U(1,:) = A_dummy(1,:)
  	  L(:,1) = A_dummy(:,1)/A_dummy(1,1)
	
 	   sumu=0.d0
  	  suml=0.d0
	
  	  do i=2,n
  	    do j=2,n
         sumu=0.d0
         suml=0.d0
         do k=1,i-1
          sumu=sumu+L(i,k)*U(k,j)
          suml=suml+L(j,k)*U(k,i)
         enddo
         U(i,j)= A_dummy(i,j)-sumu
!       if (U(i,i) .eq. 0.d0) stop 'Error: Division by zero!! In inversion subroutine '
        ! Possible bug
          L(j,i)=(A_dummy(j,i)-suml)/U(i,i)
   	   enddo
  	  enddo

!	  Assigning all zero elements in triangular matrices
   	  det=1.d0
   	  do i=1,n
    	  det=det*U(i,i)
    	  do j=1,n
     	   if (i>j) then
          U(i,j)=0.d0
          elseif(j>i) then
          L(i,j)=0.d0
          endif
    	  enddo
  	  enddo

    ! Checking Determinant for singularity
	  eps = 10e-8
   	 if (abs(det)<eps) then
!     	 print*,'det',det
!     	 print*,'ERROR: Matrix is Singular or Ill-conditioned!!'
!    	  call print_mat(A)
!    	  print*,'Determinant was found to be:'
!    	  print*,det
!     	 call print_mat(U)
!    	  stop 404
   	 endif

!	 Changing RHS loop
  	  do bb=1,n
     	 bvec=0.d0
    	  bvec(bb)=1.d0
	
    	  Pb = matmul(P,bvec)
    	  d=0.d0
    	  x=0.d0

!	 Forward Substitution
    	  d(1) = Pb(1)
     	 do i=2,n
        suml=0.d0
     	   do k=1,i-1
          suml=suml+L(i,k)*d(k)
    	    enddo
    	    d(i)=Pb(i)-suml
     	 enddo

!	 Backward Substitution
   	   x(n)=d(n)/U(n,n)
   	   do i=n-1,1,-1
         sumu=0.d0
         do k=i+1,n
          sumu=sumu+U(i,k)*x(k)
     	   enddo
     	   x(i)=(d(i)-sumu)/U(i,i)
     	 enddo

    	  Ainv(:,bb) = x
        	enddo

 	 end function inv
