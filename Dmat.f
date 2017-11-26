!=======GENERATING DMAT MATRIX
	do i = 1,NNODE
	xx = (COORDS(1,i)-xc)/dia
	yy = (COORDS(2,i)-yc)/dia	
	
	Dmat((i-1)*2+1,1)= 1.d0
	Dmat((i-1)*2+1,2)= 0.d0
	Dmat((i-1)*2+1,3)= -yy
	Dmat((i-1)*2+1,4)= xx
	Dmat((i-1)*2+1,5)= yy
	Dmat((i-1)*2+1,6)= xx
	Dmat((i-1)*2+1,7)= (xx*yy)
	Dmat((i-1)*2+1,8)= 0.d0
	Dmat((i-1)*2+1,9)= (xx*xx)
	Dmat((i-1)*2+1,10)= 0.d0
	Dmat((i-1)*2+1,11)= yy*yy
	Dmat((i-1)*2+1,12)= 0.d0
!==================================
	Dmat((i-1)*2+2,1)= 0.d0
	Dmat((i-1)*2+2,2)= 1.d0
	Dmat((i-1)*2+2,3)= xx
	Dmat((i-1)*2+2,4)= -yy
	Dmat((i-1)*2+2,5)= xx
	Dmat((i-1)*2+2,6)= yy
	Dmat((i-1)*2+2,7)= 0.d0
	Dmat((i-1)*2+2,8)= xx*yy
	Dmat((i-1)*2+2,9)= 0.d0
	Dmat((i-1)*2+2,10)= xx*xx
	Dmat((i-1)*2+2,11)= 0.d0
	Dmat((i-1)*2+2,12)= yy*yy
	end do !2*NNODE
!======= INTEGRATION FOR THE LAST TWO ROWS OF Dmat,:'(
!=======the do loop for the inegration is:-
	s12 =0.d0
	s13 =0.d0
	s15 =0.d0
	s16 =0.d0
	s17 =0.d0
	
	do i = 1,verts
	if(i == verts) then
	ip1 = 1
	else
	ip1 = i + 1
	end if 
	mn = i + verts  !
	DELY = (COORDS(2,ip1) -COORDS(2,i))/2.d0
	do j = 1,3 ! i need to declare j and  coeff 
	if (j==1) then
	w = i
	coeff = (1.d0/3.d0)
	elseif(j ==2) then
	w = mn
	coeff = (4.d0/3.d0)
	else
	w = ip1
	coeff = (1.d0/3.d0)
	end if
	xx = (COORDS(1,w)-xc)/dia
	yy = (COORDS(2,w)-yc)/dia


	s12 = S12 + (dia/2.d0)*(xx**2)*coeff*DELY 
	s13 = S13 + dia*xx*yy*coeff*DELY
	s15 = S15 + (dia/2.d0)*(xx**2)*yy*coeff*DELY
	s16 = S16 + (dia/3.d0)*(xx**3)*coeff*DELY
	s17 = S17 +(dia*xx)*yy*yy*coeff*DELY
	
	 end do !j 1,3
	end do ! for intgration of d 
	
	Dmat(NNODE*2+1,1) = 1.d0
	Dmat(NNODE*2+1,3) = -s13/area
	Dmat(NNODE*2+1,4) = s12/area
	Dmat(NNODE*2+1,5) = s13/area
	Dmat(NNODE*2+1,6) = s12/area
	Dmat(NNODE*2+1,7) = s15/area
	Dmat(NNODE*2+1,9) = s16/area
	Dmat(NNODE*2+1,11) = s17/area
	!======================================
	Dmat(NNODE*2+2,2) = 1.d0
	Dmat(NNODE*2+2,3) = s12/area
	Dmat(NNODE*2+2,4) = -s13/area
	Dmat(NNODE*2+2,5) = s12/area
	Dmat(NNODE*2+2,6) = s13/area
	Dmat(NNODE*2+2,8) = s15/area
	Dmat(NNODE*2+2,10) = s16/area
	Dmat(NNODE*2+2,12) = s17/area
!====================================================================	
