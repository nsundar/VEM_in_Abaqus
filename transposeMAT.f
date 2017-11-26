	! contains 
	
	 function trans(A) result(Atrns)
  	  real*8, intent(in), dimension(:,:) :: A
  	  real*8, dimension(size(A,2),size(A,1)) :: Atrns
	  integer::i,j
          do i =1,size(A,1)  !rows of A
	  do j =1,size(A,2)  ! colms of A
          Atrns(j,i) = A(i,j)
 	  end do!i
	  end do!j

 	 end function trans    


   
