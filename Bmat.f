!==============FOR b1 MATRIX======================================

	  do i = 1,NNODE
          if(i.le.verts) then
        	if(i.gt.1)then
	        im1=i-1
		else
		im1=verts
		endif
       		if(i.lt.verts)then
		ip1=i+1
		else
		ip1=1
		endif

        	xx = (COORDS(1,i)-xc)/dia
        	yy = (COORDS(2,i)-yc)/dia
       	
	 call normal(COORDS(1,i),COORDS(1,im1),nx,ny,length)
	 
	iv1(1)=0.d0
	iv1(2)=0.d0
	iv1(3)=0.d0
	iv1(4) =(1.d0/dia)*nx*(ca1-cb1)
	iv1(5) =(1.d0/dia)*ny*2*cc1
	iv1(6) =(1.d0/dia)*nx*(ca1+cb1)
	iv1(7) =(1.d0/dia)*(nx*yy*ca1+ny*xx*cc1)
	iv1(8) =(1.d0/dia)*(cb1*nx*xx+cc1*ny*yy)
	iv1(9) =(1.d0/dia)*2.d0*nx*xx*ca1
	iv1(10) =(1.d0/dia)*2.d0*ny*xx*cc1
	iv1(11) =(1.d0/dia)*2.d0*ny*yy*cc1      
        iv1(12) =(1.d0/dia)*2.d0*nx*yy*cb1

	
	
	iv7 = (length/6.d0)*iv1  ! 

	call normal(COORDS(1,ip1),COORDS(1,i),nx,ny,length)

	
	iv1(1)=0.d0
        iv1(2)=0.d0
        iv1(3)=0.d0
        iv1(4) =(1.d0/dia)*nx*(ca1-cb1)
        iv1(5) =(1.d0/dia)*ny*2*cc1
        iv1(6) =(1.d0/dia)*nx*(ca1+cb1)
        iv1(7) =(1.d0/dia)*(nx*yy*ca1+ny*xx*cc1)
        iv1(8) =(1.d0/dia)*(cb1*nx*xx+cc1*ny*yy)
        iv1(9) =(1.d0/dia)*2.d0*nx*xx*ca1
        iv1(10) =(1.d0/dia)*2.d0*ny*xx*cc1
        iv1(11) =(1.d0/dia)*2.d0*ny*yy*cc1
        iv1(12) =(1.d0/dia)*2.d0*nx*yy*cb1 

        iv2 = (length/6.d0)*iv1  !	
	
	iv = izero
	iv = iv7+iv2	
	b1(:,i) = iv
!         print*,'nx,ny',nx,ny
        else
       
        im1 = i -verts
       
        if(i == NNODE)then 
        ip1 = 1
        else
        ip1 = im1+1
        end if 
       
	xx = (COORDS(1,i)-xc)/dia
	yy = (COORDS(2,i)-yc)/dia
	
	call normal(COORDS(1,ip1),COORDS(1,im1),nx,ny,length)
!	 print*,'IMd',coords(1,i),coords(2,i),nx,ny,ip1,im1
        iv1(1)=0.d0
        iv1(2)=0.d0
        iv1(3)=0.d0
        iv1(4) =(1.d0/dia)*nx*(ca1-cb1)
        iv1(5) =(1.d0/dia)*ny*2*cc1
        iv1(6) =(1.d0/dia)*nx*(ca1+cb1)
        iv1(7) =(1.d0/dia)*(nx*yy*ca1+ny*xx*cc1)
        iv1(8) =(1.d0/dia)*(cb1*nx*xx+cc1*ny*yy)
        iv1(9) =(1.d0/dia)*2.d0*nx*xx*ca1
        iv1(10) =(1.d0/dia)*2.d0*ny*xx*cc1
        iv1(11) =(1.d0/dia)*2.d0*ny*yy*cc1
        iv1(12) =(1.d0/dia)*2.d0*nx*yy*cb1

	iv3 = (4.d0*length/6.d0)*iv1   
!            print*,'nx,ny',nx,ny
	 b1(:,i) = iv3
         endif
         end do!NNODE
!==================================================================

! =================FOR b2 MATRIX===================================
	  
	do i = 1,NNODE
        if(i.le.verts) then
	if(i.gt.1)then
	im1=i-1
	else
	im1=verts
	endif
       	if(i.lt.verts)then
	ip1=i+1
	else
	ip1=1
	endif

	xx = (COORDS(1,i)-xc)/dia
	yy = (COORDS(2,i)-yc)/dia

	call normal(COORDS(1,i),COORDS(1,im1),nx,ny,length)
	
	iv1(1)=0.d0
        iv1(2)=0.d0
        iv1(3)=0.d0
	iv1(4) =(1.d0/dia)*ny*(cb1-ca1)
	iv1(5) =(1.d0/dia)*2.d0*nx*cc1
	iv1(6) =(1.d0/dia)*ny*(cb1+ca1)
	iv1(7) =(1.d0/dia)*(cc1*xx*nx+cb1*ny*yy)
	iv1(8) =(1.d0/dia)*(nx*yy*cc1+ny*xx*ca1)
        iv1(9) =(1.d0/dia)*2.d0*ny*xx*cb1
	iv1(10) =(1.d0/dia)*2.d0*nx*xx*cc1
	iv1(11) =(1.d0/dia)*2.d0*nx*yy*cc1
	iv1(12) =(1.d0/dia)*2.d0*ny*yy*ca1
	
	iv7 = (length/6.d0)*iv1  ! 
	
	
	call normal(COORDS(1,ip1),COORDS(1,i),nx,ny,length)
	
	
        iv1(1)=0.d0
        iv1(2)=0.d0
        iv1(3)=0.d0
        iv1(4) =(1.d0/dia)*ny*(cb1-ca1)
        iv1(5) =(1.d0/dia)*2.d0*nx*cc1
        iv1(6) =(1.d0/dia)*ny*(cb1+ca1)
        iv1(7) =(1.d0/dia)*(cc1*xx*nx+cb1*ny*yy)
        iv1(8) =(1.d0/dia)*(nx*yy*cc1+ny*xx*ca1)
        iv1(9) =(1.d0/dia)*2.d0*ny*xx*cb1
        iv1(10) =(1.d0/dia)*2.d0*nx*xx*cc1
        iv1(11) =(1.d0/dia)*2.d0*nx*yy*cc1
        iv1(12) =(1.d0/dia)*2.d0*ny*yy*ca1

	iv2 = (length/6.d0)*iv1  ! 
	iv =izero
	iv = iv7+iv2
	
	b2(:,i) = iv
	
        else
       
        im1 = i-verts
       
        if(i == NNODE)then 
        ip1 = 1
        else
        ip1 = im1+1
        end if 
	xx = (COORDS(1,i)-xc)/dia
	yy = (COORDS(2,i)-yc)/dia
	
	call normal(COORDS(1,ip1),COORDS(1,im1),nx,ny,length)
	
        iv1(1)=0.d0
        iv1(2)=0.d0
        iv1(3)=0.d0
        iv1(4) =(1.d0/dia)*ny*(cb1-ca1)
        iv1(5) =(1.d0/dia)*2.d0*nx*cc1
        iv1(6) =(1.d0/dia)*ny*(cb1+ca1)
        iv1(7) =(1.d0/dia)*(cc1*xx*nx+cb1*ny*yy)
        iv1(8) =(1.d0/dia)*(nx*yy*cc1+ny*xx*ca1)
        iv1(9) =(1.d0/dia)*2.d0*ny*xx*cb1
        iv1(10) =(1.d0/dia)*2.d0*nx*xx*cc1
        iv1(11) =(1.d0/dia)*2.d0*nx*yy*cc1
        iv1(12) =(1.d0/dia)*2.d0*ny*yy*ca1
	
	iv3 = (4.d0*length/6.d0)*iv1  !
	
	b2(:,i) = iv3
       endif
       
      end do!NNODE
!====================================================================
! for adding the two internal dof 
	iv4(1) = 0.d0
	iv4(2) = 0.d0
	iv4(3) = 0.d0
	iv4(4) = 0.d0
	iv4(5) = 0.d0
	iv4(6) = 0.d0
	iv4(7) = 0.d0
	iv4(8) = (1.d0/dia**2.d0)*(cc1+cb1)
	iv4(9) =(1.d0/dia**2.d0)*2.d0*ca1
	iv4(10)=0.d0
	iv4(11) =(1.d0/dia**2.d0)*2.d0*cc1
	iv4(12)=0.d0


	b1(:,NNODE+1) = -area*iv4

	iv5(1) = 0.d0
	iv5(2) = 0.d0
	iv5(3) = 0.d0
	iv5(4) = 0.d0
	iv5(5) = 0.d0
	iv5(6) = 0.d0
	iv5(7) = (1.d0/dia**2.d0)*(cc1+cb1)
	iv5(8) = 0.d0
	iv5(9) =0.d0
	iv5(10) =(1.d0/dia**2.d0)*2.d0*cc1
	iv5(11) = 0.d0
	iv5(12) =(1.d0/dia**2.d0)*2.d0*ca1
	
	b2(:,NNODE+1) = -area*iv5
!====================================================================
!=======GENERATING B MATRIX
	do i = 1,NNODE+1
	indx = i*2-1
	B(:,indx+1) = b2(:,i)
	B(:,indx) = b1(:,i)
	end do !i NNODE+1

!======= FOR FIRST THREE ROWS		
	B(1,2*NNODE+1) = 1.d0
	B(2,2*NNODE+2) = 1.d0

	do i =1,NNODE
	indx = i*2-1
	B(3,indx) = -(1.d0/NNODE)*((COORDS(2,i)-yc)/dia)
	B(3,indx+1) = (1.d0/NNODE)*((COORDS(1,i)-xc)/dia)
	
	end do
!====================================================================
