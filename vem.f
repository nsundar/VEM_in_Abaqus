      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     
c      implicit none
      INCLUDE 'ABA_PARAM.INC'

       parameter(zero=0.d0, half=0.5, one=1.d0, two=2.d0, three=3.d0, 
     1 four=4.d0, six=6.d0, eight=8.d0, twelve=12.d0,NGPT=4)
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
! NNODE= NUMBER OF NODES PER ELEMENT(8-FOR Q8)
! NDOF = NUM OF DEOF PER ELEMENT (16-FOR Q8) 

C=================================================================
C======= variables for VEM
C=================================================================
        real*8 dia,hD,coeff,E1,NU,ktrace
        real*8 xc,yc ,area,nx,ny,xx,yy,detg    
        real*8 s11,s12,s13,s15,s16,s17,length,DELY,mu,lamda
        real*8 iv1(12),iv2(12),iv3(12),iv4(12),iv5(12),iv(12),iv7(12)
	real*8 izero(12),DivH(1:12,1:12),KIINV(1:2,1:2),CMAT(1:3,1:3)
        real*8 GMAT(1:12,1:12),Ginv(1:12,1:12),Gtilda(1:12,1:12)

	real*8,allocatable::b1(:,:),b2(:,:),B(:,:)
	real*8,allocatable::Dmat(:,:),shp(:,:),strans(:,:)
        real*8,allocatable::kmu(:,:),klamda(:,:),stiff(:,:)
        real*8,allocatable::kmu_temp(:,:),klamda_temp(:,:),Kstiff(:,:) 
!=================================
       ! declare for the stability part 
        real*8,allocatable::mat1(:,:),mat2(:,:),mat2t(:,:)
        real*8,allocatable::kstability(:,:),ID(:,:)
!=================================
	real*8::KII(1:2,1:2)
        real*8,allocatable::KBB(:,:),KBI(:,:)
        real*8,allocatable::KIB(:,:),MBI(:,:),MBB(:,:)
        real*8,allocatable::Rtemp(:)
!=================================
!=======INTEGER DECLARATION 
        integer NNODE,w,mn,i1,j1,gs
        integer k,l,im,pt,e,nseg,i,k1,k4
        integer ip1,im1,status,j,verts,indx  
!=================================================================
        allocate(b1(1:12,1:NNODE+1),stat=status)         
	allocate(b2(1:12,1:NNODE+1),stat=status)
	allocate(B(1:12,1:(2*NNODE+2)),stat=status)
	allocate(Dmat(1:(2*NNODE+2),1:12),stat=status)
	allocate(shp(1:12,1:(2*NNODE+2)),stat=status)
	allocate(strans(1:(2*NNODE+2),1:12),stat=status)   
        allocate(kmu(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status)  
        allocate(kmu_temp(1:(2*NNODE+2),1:12),stat=status) 
        allocate(klamda_temp(1:(2*NNODE+2),1:12),stat=status)  
        allocate(klamda(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status) 
        allocate(stiff(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status)
        allocate(Kstiff(1:2*NNODE,1:2*NNODE),stat=status)
C=================================================================
        allocate(mat1(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status) 
        allocate(mat2(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status)
        allocate(mat2t(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status)
        allocate(kstability(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status)
        allocate(ID(1:(2*NNODE+2),1:(2*NNODE+2)),stat=status) 
C=================================================================
        allocate(KBB(1:2*NNODE,1:2*NNODE),stat=status) 
        allocate(KBI(1:2*NNODE,1:2),stat=status)  
        allocate(KIB(1:2,1:2*NNODE),stat=status) 
        allocate(MBI(1:2*NNODE,1:2),stat=status)  
        allocate(MBB(1:2*NNODE,1:2*NNODE),stat=status)
        allocate(Rtemp(1:NDOFEL),stat=status)
C=================================================================
C======= VEM MAIN PGM STARTS HERE,:)
!======= DUMMY ARRY FOR CREATING ZERO VECTOR,:'(
	do i =1,12
	izero(i) = 0.d0
	end do  
!======= Dmat 
	do i =1,2*NNODE+2
	do j = 1,12
	Dmat(i,j) = 0.d0
	end do 
	end do
	do i =1,12
	do j = 1,12
	GMAT(i,j) = 0.d0
!        Ginv(i,j) = 0.d0
	DivH(i,j) = 0.d0
	end do 
	end do
	
	do i=1,12
	do j=1,2*NNODE+2
	B(i,j)=0.d0
	Dmat(j,i)=0.d0
	end do
	end do


        do i=1,12
        do j=1,NNODE+1
        b1(i,j)=0.d0
        b2(i,j)=0.d0
        end do
        end do
!=================================================================
!======= for geting the diameter of the polygon 
        call diameterpoly(COORDS,NNODE,hD)
        dia=hD
!=================================================================
        E1 = PROPS(2)
        NU = PROPS(3)
        CMAT(1,1:3) = (/1.d0,NU,0.d0 /)
        CMAT(2,1:3) = (/NU,1.d0,0.d0 /)
        CMAT(3,1:3) = (/0.d0,0.d0,((1.d0-NU)/2.d0)/)
            CMAT=(E1/(1.d0-(NU**2.d0)))*CMAT
        ca1 =CMAT(1,1)
        Cb1 =CMAT(1,2)
        cc1 =CMAT(3,3)

!=================================================================
!=================================================================
c for geting the polygon area and centroid
       call area_polygon(COORDS,area,xc,yc,NNODE)
       verts=int(NNODE*0.5d0)    ! FOR THE CASE OF Q8 VERTS -4
C=================================================================
	include 'Bmat.f'
	include 'Dmat.f'
!	include 'DivH.f'	
!=======GMATRIX MATRIX MULT,

        GMAT = matmul(B, DMAT)
	Ginv = inv(GMAT)
 	Shp = matmul(Ginv,B)
        strans = transpose(Shp) 
        ! need to take the transpose of S matrix, Strans 

!====================================================================
         Gtilda = GMAT
         do i=1,3
         do j=1,12
         Gtilda(i,j) = 0.d0
         enddo
         enddo 
!====================================================================
         do i=1,2*NNODE+2
         do j=1,2*NNODE+2
          if (i.eq.j) then
          ID(i,j) =1.d0
          else
          ID(i,j) = 0.d0
          end if 
	 end do
	 end do

        mat1 = matmul(Dmat,shp)
        mat2 = ID -mat1         
        mat2t = transpose(mat2)
        kstability  =matmul(mat2t,mat2)

         Kmu_temp = matmul(Strans,Gtilda)
         kmu = matmul(kmu_temp,Shp)
!==================================================================== 
         ktrace = 0.D0
         do i=1,2*NNODE+2

          ktrace = ktrace + Kmu(i,i)

          end do

          stiff = Kmu + ktrace*kstability
!================================================================

	 KBB =stiff(1:2*NNODE,1:2*NNODE)
	 KBI =stiff(1:2*NNODE,2*NNODE+1:2*NNODE+2)
	 KIB =stiff(2*NNODE+1:2*NNODE+2,1:2*NNODE)
	 KII =stiff(2*NNODE+1:2*NNODE+2,2*NNODE+1:2*NNODE+2)
	 KIINV = inv(KII)
	 MBI = matmul(KBI,KIINV)
	 MBB = matmul(MBI,KIB)
	 Kstiff = KBB-MBB
	 
!================================================================
	amatrx=kstiff
!        print*,'sz',size(u,1)!,size(rhs,2)
              Rtemp = -matmul(amatrx,u)
         do i=1,NDOFEL
             RHS(i,1)=RHS(i,1)+Rtemp(i)
         enddo
!	  RETURN
!=====================================================================
!======= DEALLOCATING SECTION
         deallocate(b1,stat=status)
	 deallocate(b2,stat=status)
 	 deallocate(B,stat=status)
	 deallocate(Dmat,stat=status)
	 deallocate(shp,stat=status)
	 deallocate(strans,stat=status)
!================================
         deallocate(kmu,stat=status)
         deallocate(kmu_temp,stat=status)
         deallocate(klamda,stat=status)
         deallocate(klamda_temp,stat=status) 
         deallocate(stiff,stat=status)
         deallocate(Kstiff,stat=status)
!================================
         deallocate(mat1,stat=status)
         deallocate(mat2,stat=status)
         deallocate(mat2t,stat=status)
         deallocate(kstability,stat=status)
         deallocate(ID,stat=status)
!===============================
	 deallocate(KBB,stat=status)
         deallocate(KBI,stat=status)
         deallocate(KIB,stat=status)

         deallocate(MBI,stat=status)
	 deallocate(MBB,stat=status)
         deallocate(Rtemp,stat=status)
!=====================================================================
!            CONTAINS
            include 'matinvs.f'

          END

C=======END OF THE MAIN PROGRAM 
C=====================================================================
            include 'subrtn.f'
C=====================================================================
