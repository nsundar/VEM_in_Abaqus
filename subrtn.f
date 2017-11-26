
!=======SUBROUTINS
             subroutine area_polygon(COORDS,area,xc,yc,NNODE)	  
	     INCLUDE 'ABA_PARAM.INC'  
              real*8 COORDS(2,*)
	      real*8 x1,y1,s,area,xc,yc
	      integer i,j,NNODE
	      s=0.d0
	      xc=0.d0
	      yc=0.d0
	      area=0.d0  
	   
	      do j=1,int(NNODE/2.d0)

	      if(j.eq.int(NNODE/2.d0))then   
	       x1= COORDS(1,1)
	      y1= COORDS(2,1)
	      else      
	      x1 = COORDS(1,j+1)
	      y1 = COORDS(2,j+1)
	      endif

	      s = COORDS(1,j)*y1 -COORDS(2,j)*x1 
!	      print*,'in area polygon S',J,s
	      xc=xc+s*(COORDS(1,j)+x1)
	      yc=yc+s*(COORDS(2,j)+y1)
	      area=area+s
	      end do
	      area = 0.5*abs(area)
!	      print*,'in area polygon S',area
! 		centroid of the polygon xc and yc 
	      xc = xc/(6*area)
	      yc = yc/(6*area)
		  RETURN
                END
	!    end subroutine area_polygon  
!==================================================================
             subroutine diameterpoly(COORDS,NNODE,hD)    
	      INCLUDE 'ABA_PARAM.INC'  
	      real*8 COORDS(2,*)
              real*8 temp,hD
	      integer j
!	      print*,'in dia',NNODE
              hD=0.d0
              temp=0.d0
              do i=1,NNODE
              do j=1,NNODE-i
              temp=dsqrt((COORDS(1,i)-COORDS(1,j))**2+
     &          (COORDS(2,i)-COORDS(2,j))**2)
	      if(temp.gt.hD) then
               hD=temp
               endif 
	       enddo
	       enddo 
	
	       RETURN
               END
C===================================================================
	  subroutine normal(x1,x2,nx,ny,length)
      	     INCLUDE 'ABA_PARAM.INC' 
  
          real*8 x1(2),x2(2),dx,dy,length
          real*8 nx,ny
        
          dx = x1(1)-x2(1)
          dy = x1(2)-x2(2)
          
          length = sqrt(dx**2+dy**2)
          nx = (1.d0/length)*(dy)
          ny = (1.d0/length)*(-dx)

   !      end subroutine normal 
	       RETURN
               END
C===================================================================
