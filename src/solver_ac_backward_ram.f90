!! iprop=1 Forward propagation
!! iprop=2 Adjoint propagation

subroutine solver_ac_backward_ram(rank,NumSou,NumRec,nxSou,nySou,nxRec,nyRec,&
Sadj,dmodel,dt,ny,nx,ctmp,&
store_snap,nt_snap,nt,Fields,kernel)

USE mod_parfile, ONLY: dPML, fs

implicit none

	!!!! Simulation Geometry
	integer :: op,op2
	integer :: unit_fields,ierr,shotid
	integer :: rank,iprop,i,nt,ny,nx,j,k,NumSou
	integer :: snap,store_snap,nt_snap,icont,nn,dn
	integer :: NumRec,l,m,jj,acc,free_surface
	integer :: dny,dnx,dny2,nxt,nyt

	integer :: nxSou(NumSou),nySou(NumSou)
	integer :: nxRec(NumRec),nyRec(NumRec)

	!!!! Simulation parameters
	real :: dt,dmodel
	real :: Sadj(nt,NumRec)
	real :: ctmp(ny,nx),kernel(ny,nx)

	integer, allocatable :: nxs(:),nys(:)
        integer, allocatable :: nxreco(:),nyreco(:)

	real, allocatable :: coef1(:),coef2(:)
        real :: Fields(ny,nx,nt_snap)
	
	!!!! Variables ecuacion acustica
	real, allocatable :: source(:,:)
	real, allocatable :: p(:,:),p_1(:,:),p_3(:,:)
	real, allocatable :: c(:,:),c2(:,:)
	real, allocatable :: deriv2_p_x(:,:),deriv2_p_y(:,:)
	real, allocatable :: deriv_p_x(:,:),deriv_p_y(:,:)

	!!!! Variables PML
	integer, allocatable :: vyt(:),vxt(:),my(:),mx(:)
	integer, allocatable :: yu(:),yd(:),xl(:),xr(:) 

	real, allocatable :: factorA(:),factorB(:)
	real, allocatable :: ld_factorA(:,:),ld_factorB(:,:)
	real, allocatable :: rd_factorA(:,:),rd_factorB(:,:)
	real, allocatable :: lu_factorA(:,:),lu_factorB(:,:)
	real, allocatable :: ru_factorA(:,:),ru_factorB(:,:)

	real, allocatable :: x_prop_pml(:,:),y_prop_pml(:,:),corner_prop_pml(:,:)

	!!!! Variables PML y
	real, allocatable :: yu_pml_p_1(:,:),yu_pml_p(:,:),yu_pml_p_3(:,:)
	real, allocatable :: yu_phi_p_y_1(:,:),yu_phi_p_y(:,:),yu_chi_p_y_1(:,:),yu_chi_p_y(:,:)

	real, allocatable :: yd_pml_p_1(:,:),yd_pml_p(:,:),yd_pml_p_3(:,:)
	real, allocatable :: yd_phi_p_y_1(:,:),yd_phi_p_y(:,:),yd_chi_p_y_1(:,:),yd_chi_p_y(:,:)

	!!!! Variables PML x direction
	real, allocatable :: xl_pml_p_1(:,:),xl_pml_p(:,:),xl_pml_p_3(:,:)
	real, allocatable :: xl_phi_p_x_1(:,:),xl_phi_p_x(:,:),xl_chi_p_x_1(:,:),xl_chi_p_x(:,:)

	real, allocatable :: xr_pml_p_1(:,:),xr_pml_p(:,:),xr_pml_p_3(:,:)
	real, allocatable :: xr_phi_p_x_1(:,:),xr_phi_p_x(:,:),xr_chi_p_x_1(:,:),xr_chi_p_x(:,:)
	
	!!!! Variables PML corner up direction
	real, allocatable :: lu_pml_p_1(:,:),lu_pml_p(:,:),lu_pml_p_3(:,:)
	real, allocatable :: lu_phi_p_x_1(:,:),lu_phi_p_x(:,:),lu_chi_p_x_1(:,:),lu_chi_p_x(:,:)
	real, allocatable :: lu_phi_p_y_1(:,:),lu_phi_p_y(:,:),lu_chi_p_y_1(:,:),lu_chi_p_y(:,:)

	real, allocatable :: ru_pml_p_1(:,:),ru_pml_p(:,:),ru_pml_p_3(:,:)
	real, allocatable :: ru_phi_p_x_1(:,:),ru_phi_p_x(:,:),ru_chi_p_x_1(:,:),ru_chi_p_x(:,:)
	real, allocatable :: ru_phi_p_y_1(:,:),ru_phi_p_y(:,:),ru_chi_p_y_1(:,:),ru_chi_p_y(:,:)
	
	!!!! Variables PML corner down direction
	real, allocatable :: ld_pml_p_1(:,:),ld_pml_p(:,:),ld_pml_p_3(:,:)
	real, allocatable :: ld_phi_p_x_1(:,:),ld_phi_p_x(:,:),ld_chi_p_x_1(:,:),ld_chi_p_x(:,:)
	real, allocatable :: ld_phi_p_y_1(:,:),ld_phi_p_y(:,:),ld_chi_p_y_1(:,:),ld_chi_p_y(:,:)

	real, allocatable :: rd_pml_p_1(:,:),rd_pml_p(:,:),rd_pml_p_3(:,:)
	real, allocatable :: rd_phi_p_x_1(:,:),rd_phi_p_x(:,:),rd_chi_p_x_1(:,:),rd_chi_p_x(:,:)
	real, allocatable :: rd_phi_p_y_1(:,:),rd_phi_p_y(:,:),rd_chi_p_y_1(:,:),rd_chi_p_y(:,:)

	character(len=500) :: file_name,file_name2,file_name3,file_fields_snap,folder_fields
	character(len=500) :: Str_snap,Str_shot,access,form,command

	logical :: input_exists

        write(Str_shot,*) shotid

	iprop=2	!! backward propagation
	!dPML=20
	acc=3
	free_surface=fs	!! si =0, PML absorbente

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!! PML layers and points discretization
	
	dn=dPML+2*acc
	dnx=dn
	dny=dn
	dny2=dn

!	if(free_surface.ne.0)dny2=0

	allocate(nxreco(NumRec),nyreco(NumRec))

	if(iprop.eq.1)  then
        	allocate(nxs(NumSou),nys(NumSou))
		nxs=dn+nxSou
		nys=dny2+nySou
		nxreco=dn+nxRec
		nyreco=dny2+nyRec	
	endif
	
	if(iprop.eq.2)	then	
		allocate(nxs(NumRec),nys(NumRec))
		nxs=dn+nxRec
		nys=dny2+nyRec	
	endif

	nxt=nx+2*dn
	nyt=dny2+ny+dn

	allocate(vxt(nxt),vyt(nyt))
	allocate(mx(nx),my(ny))
	allocate(yu(dny2))
	allocate(xl(dnx))
	allocate(xr(dnx))
	allocate(yd(dny))

	vxt=0;vyt=0;mx=0;my=0;yu=0;yd=0;xl=0;xr=0;

	allocate(coef1(acc+1),coef2(acc+1))
	coef1=0.;coef2=0.;

	allocate(source(nyt,nxt))

	allocate(p_1(nyt,nxt),p(nyt,nxt),p_3(nyt,nxt))
	source=0.;p_1=0.;p=0.;p_3=0.;

	allocate(c(nyt,nxt),c2(nyt,nxt))
	c=0.;c2=0.; !!6

	allocate(deriv2_p_y(nyt,nxt),deriv2_p_x(nyt,nxt))
	allocate(deriv_p_y(nyt,nxt),deriv_p_x(nyt,nxt))		
	deriv2_p_y=0.;deriv2_p_x=0.;
	deriv_p_y=0.;deriv_p_x=0.;	!!16	
		
	!!! Variables PML
	allocate(x_prop_pml(ny,dn),y_prop_pml(dn,nx),corner_prop_pml(dn,dn))
	x_prop_pml=0.;y_prop_pml=0.;corner_prop_pml=0.;
	
	allocate(factorA(dn),factorB(dn))
	allocate(lu_factorA(dn,dn),lu_factorB(dn,dn),ru_factorA(dn,dn),ru_factorB(dn,dn))
	allocate(ld_factorA(dn,dn),ld_factorB(dn,dn),rd_factorA(dn,dn),rd_factorB(dn,dn))
	factorA=0.;factorB=0.;
	lu_factorA=0;lu_factorB=0;ru_factorA=0;ru_factorB=0;!!8
	ld_factorA=0;ld_factorB=0;rd_factorA=0;rd_factorB=0;!!8
	
	!!! Variables PML y	
	allocate(yu_pml_p_1(dn,nx),yu_pml_p(dn,nx),yu_pml_p_3(dn,nx))
	allocate(yu_phi_p_y_1(dn,nx),yu_phi_p_y(dn,nx),yu_chi_p_y_1(dn,nx),yu_chi_p_y(dn,nx))
	yu_pml_p_1=0.;yu_pml_p=0.;yu_pml_p_3=0.;
	yu_phi_p_y_1=0.;yu_phi_p_y=0.;yu_chi_p_y_1=0.;yu_chi_p_y=0.;
	
	!!! Variables PML corner up
	allocate(lu_pml_p_1(dn,dn),lu_pml_p(dn,dn),lu_pml_p_3(dn,dn))
	allocate(lu_phi_p_x_1(dn,dn),lu_phi_p_x(dn,dn),lu_chi_p_x_1(dn,dn),lu_chi_p_x(dn,dn))
	allocate(lu_phi_p_y_1(dn,dn),lu_phi_p_y(dn,dn),lu_chi_p_y_1(dn,dn),lu_chi_p_y(dn,dn))
	lu_pml_p_1=0.;lu_pml_p=0.;lu_pml_p_3=0.;
	lu_phi_p_x_1=0.;lu_phi_p_x=0.;lu_chi_p_x_1=0.;lu_chi_p_x=0.;
	lu_phi_p_y_1=0.;lu_phi_p_y=0.;lu_chi_p_y_1=0.;lu_chi_p_y=0.;
	
	allocate(ru_pml_p_1(dn,dn),ru_pml_p(dn,dn),ru_pml_p_3(dn,dn))
	allocate(ru_phi_p_x_1(dn,dn),ru_phi_p_x(dn,dn),ru_chi_p_x_1(dn,dn),ru_chi_p_x(dn,dn))
	allocate(ru_phi_p_y_1(dn,dn),ru_phi_p_y(dn,dn),ru_chi_p_y_1(dn,dn),ru_chi_p_y(dn,dn))
	ru_pml_p_1=0.;ru_pml_p=0.;ru_pml_p_3=0.;
	ru_phi_p_x_1=0.;ru_phi_p_x=0.;ru_chi_p_x_1=0.;ru_chi_p_x=0.;
	ru_phi_p_y_1=0.;ru_phi_p_y=0.;ru_chi_p_y_1=0.;ru_chi_p_y=0.;

	allocate(yd_pml_p_1(dn,nx),yd_pml_p(dn,nx),yd_pml_p_3(dn,nx))
	allocate(yd_phi_p_y_1(dn,nx),yd_phi_p_y(dn,nx),yd_chi_p_y_1(dn,nx),yd_chi_p_y(dn,nx))
	yd_pml_p_1=0.;yd_pml_p=0.;yd_pml_p_3=0.;
	yd_phi_p_y_1=0.;yd_phi_p_y=0.;yd_chi_p_y_1=0.;yd_chi_p_y=0.;

	!!! Variables PML x 
	allocate(xl_pml_p_1(ny,dn),xl_pml_p(ny,dn),xl_pml_p_3(ny,dn))
	allocate(xl_phi_p_x_1(ny,dn),xl_phi_p_x(ny,dn),xl_chi_p_x_1(ny,dn),xl_chi_p_x(ny,dn))
	xl_pml_p_1=0.;xl_pml_p=0.;xl_pml_p_3=0.;
	xl_phi_p_x_1=0.;xl_phi_p_x=0.;xl_chi_p_x_1=0.;xl_chi_p_x=0.;
	
	allocate(xr_pml_p_1(ny,dn),xr_pml_p(ny,dn),xr_pml_p_3(ny,dn))
	allocate(xr_phi_p_x_1(ny,dn),xr_phi_p_x(ny,dn),xr_chi_p_x_1(ny,dn),xr_chi_p_x(ny,dn))
	xr_pml_p_1=0.;xr_pml_p=0.;xr_pml_p_3=0.;
	xr_phi_p_x_1=0.;xr_phi_p_x=0.;xr_chi_p_x_1=0.;xr_chi_p_x=0.;
	
	!!! Variables PML corner down
	allocate(ld_pml_p_1(dn,dn),ld_pml_p(dn,dn),ld_pml_p_3(dn,dn))
	allocate(ld_phi_p_x_1(dn,dn),ld_phi_p_x(dn,dn),ld_chi_p_x_1(dn,dn),ld_chi_p_x(dn,dn))
	allocate(ld_phi_p_y_1(dn,dn),ld_phi_p_y(dn,dn),ld_chi_p_y_1(dn,dn),ld_chi_p_y(dn,dn))
	ld_pml_p_1=0.;ld_pml_p=0.;ld_pml_p_3=0.;
	ld_phi_p_x_1=0.;ld_phi_p_x=0.;ld_chi_p_x_1=0.;ld_chi_p_x=0.;
	ld_phi_p_y_1=0.;ld_phi_p_y=0.;ld_chi_p_y_1=0.;ld_chi_p_y=0.;
	
	allocate(rd_pml_p_1(dn,dn),rd_pml_p(dn,dn),rd_pml_p_3(dn,dn))
	allocate(rd_phi_p_x_1(dn,dn),rd_phi_p_x(dn,dn),rd_chi_p_x_1(dn,dn),rd_chi_p_x(dn,dn))
	allocate(rd_phi_p_y_1(dn,dn),rd_phi_p_y(dn,dn),rd_chi_p_y_1(dn,dn),rd_chi_p_y(dn,dn))
	rd_pml_p_1=0.;rd_pml_p=0.;rd_pml_p_3=0.;
	rd_phi_p_x_1=0.;rd_phi_p_x=0.;rd_chi_p_x_1=0.;rd_chi_p_x=0.;
	rd_phi_p_y_1=0.;rd_phi_p_y=0.;rd_chi_p_y_1=0.;rd_chi_p_y=0.;


	!!!!!!
	!!!!!! Define matrix pointers
	!!!!!!

	m=0
	do i=1,nxt
		m=m+1
		vxt(m)=i
	enddo
	m=0
	do i=1,nyt
		m=m+1
		vyt(m)=i
	enddo
	
	!!!! Pointers bottom PML 

	m=0
	do i=dnx+1,dnx+nx
		m=m+1
		mx(m)=i
	enddo

	m=0
	do i=dny2+1,dny2+ny
		m=m+1
		my(m)=i
	enddo
	
	!!!! Pointers PML right-y corner 
	m=0
	do i=1,dny2
		m=m+1
		yu(m)=i
	enddo
	
	m=0
	do i=dny2+ny+1,dny2+ny+dny
		m=m+1
		yd(m)=i
	enddo

	!!!! Pointers PML x-y corner

	m=0
	do i=1,dnx
		m=m+1
		xl(m)=i
	enddo

	m=0
	do i=dnx+nx+1,dnx+nx+dnx
		m=m+1
		xr(m)=i
	enddo
	
	!!!!!	Define Acquisition Geometry
	!!!!!

        do j=1,nx
        do i=1,ny 
                c(dny2+i,dnx+j)=ctmp(i,j);
        enddo
        enddo

        do j=1,dnx
        do i=1,dny2
                c(i,j)=ctmp(1,1)                !! corner left up
                c(i,dnx+j)=ctmp(1,j)            !! y up
                c(i,dnx+nx+j)=ctmp(1,nx)        !! corner right up
        enddo
        enddo

        do j=1,dnx
        do i=1,dny
                c(i+dny2+ny,dnx+j)=ctmp(ny,j);          !! 4
                c(i+dny2+ny,j)=ctmp(ny,1);              !! 5
                c(i+dny2+ny,dnx+nx+j)=ctmp(ny,nx);      !! 6
        enddo
        enddo
        
        do j=1,dnx
        do i=1,ny
                c(dny2+i,j)=ctmp(i,1);          !!capa xl
                c(dny2+i,dnx+nx+j)=ctmp(i,nx);  !!capa xr
        enddo
        enddo

	c2=c**2*dt**2/dmodel**2;	!!ojo, que pasa si dx es distinto a dz?

	deallocate(c)
	
!!!!!!!!!!!!!!   PML'S ABSORTION COEFFICIENTS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call PML_coeffs(acc,dPML,dt,dmodel,dn,factorA,factorB,&
	ld_factorA,rd_factorA,lu_factorA,ru_factorA,ld_factorB,rd_factorB,lu_factorB,ru_factorB)
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!! Start Computation !!!!!

	l=0
	if(iprop.eq.1)snap=0
	if(iprop.eq.2)snap=nt_snap

	call finite_difference_coef(1,acc,coef1)			
	call finite_difference_coef(2,acc,coef2)			

	nn=nt/4
	icont=0

	op=0
	op2=0
	do k=2,nt

		l=l+1

		if(iprop.eq.2)	then
			do j=1,NumRec
				source(nys(j),nxs(j))=Sadj(k-1,j)
			enddo
		endif


	        call space_deriv2x(acc,coef2,nyt,nxt,p,deriv2_p_x)
	        call space_deriv2y(acc,coef2,nyt,nxt,p,deriv2_p_y)  


!!!		START PML'S CALCULATION
	        
		call space_derivx(acc,coef1,nyt,nxt,p,deriv_p_x)
		call space_derivy(acc,coef1,nyt,nxt,p,deriv_p_y)
		
!!!		START PML'S CALCULATION

!!!		PML x- CALCULATION

		call pml_matrices_x(1,acc,coef1,ny,dn,nyt,nxt,my,xr,factorA,factorB,&
		deriv_p_x,deriv2_p_x,xr_phi_p_x,xr_phi_p_x_1,xr_chi_p_x,xr_chi_p_x_1,x_prop_pml)		
		xr_pml_p_3=-xr_pml_p_1+2*xr_pml_p+c2(my,xr)*&
		(deriv2_p_x(my,xr)+deriv2_p_y(my,xr)+x_prop_pml)
		
		call pml_matrices_x(2,acc,coef1,ny,dn,nyt,nxt,my,xl,factorA,factorB,&
		deriv_p_x,deriv2_p_x,xl_phi_p_x,xl_phi_p_x_1,xl_chi_p_x,xl_chi_p_x_1,x_prop_pml)		
		xl_pml_p_3=-xl_pml_p_1+2*xl_pml_p+c2(my,xl)*&
		(deriv2_p_x(my,xl)+deriv2_p_y(my,xl)+x_prop_pml)		

!!!		PML y- CALCULATION
	
		call pml_matrices_y(1,acc,coef1,dn,nx,nyt,nxt,yd,mx,factorA,factorB,&
		deriv_p_y,deriv2_p_y,yd_phi_p_y,yd_phi_p_y_1,yd_chi_p_y,yd_chi_p_y_1,y_prop_pml)		
		yd_pml_p_3=-yd_pml_p_1+2*yd_pml_p+c2(yd,mx)*&
		(deriv2_p_x(yd,mx)+deriv2_p_y(yd,mx)+y_prop_pml)
		
		if(free_surface.eq.0)	then

		call pml_matrices_y(2,acc,coef1,dn,nx,nyt,nxt,yu,mx,factorA,factorB,&
		deriv_p_y,deriv2_p_y,yu_phi_p_y,yu_phi_p_y_1,yu_chi_p_y,yu_chi_p_y_1,y_prop_pml)
		yu_pml_p_3=-yu_pml_p_1+2*yu_pml_p+c2(yu,mx)*&
		(deriv2_p_x(yu,mx)+deriv2_p_y(yu,mx)+y_prop_pml)	

!!!		PML corner up CALCULATION

		call pml_matrices_xycorner(1,acc,coef1,dn,dn,nyt,nxt,yu,xl,lu_factorA,lu_factorB,&
		deriv_p_x,deriv2_p_x,lu_phi_p_x,lu_phi_p_x_1,lu_chi_p_x,lu_chi_p_x_1,corner_prop_pml)
		lu_pml_p_3=-lu_pml_p_1+2*lu_pml_p+c2(yu,xl)*&
		(deriv2_p_x(yu,xl)+deriv2_p_y(yu,xl)+corner_prop_pml)
		
		call pml_matrices_xycorner(2,acc,coef1,dn,dn,nyt,nxt,yu,xl,lu_factorA,lu_factorB,&
		deriv_p_y,deriv2_p_y,lu_phi_p_y,lu_phi_p_y_1,lu_chi_p_y,lu_chi_p_y_1,corner_prop_pml)
		lu_pml_p_3=lu_pml_p_3+c2(yu,xl)*corner_prop_pml

		call pml_matrices_xycorner(1,acc,coef1,dn,dn,nyt,nxt,yu,xr,ru_factorA,ru_factorB,&
		deriv_p_x,deriv2_p_x,ru_phi_p_x,ru_phi_p_x_1,ru_chi_p_x,ru_chi_p_x_1,corner_prop_pml)
		ru_pml_p_3=-ru_pml_p_1+2*ru_pml_p+c2(yu,xr)*&
		(deriv2_p_x(yu,xr)+deriv2_p_y(yu,xr)+corner_prop_pml)
				
		call pml_matrices_xycorner(2,acc,coef1,dn,dn,nyt,nxt,yu,xr,ru_factorA,ru_factorB,&
		deriv_p_y,deriv2_p_y,ru_phi_p_y,ru_phi_p_y_1,ru_chi_p_y,ru_chi_p_y_1,corner_prop_pml)
		ru_pml_p_3=ru_pml_p_3+c2(yu,xr)*corner_prop_pml
				
		endif
		
!!!		PML corner down CALCULATION
			
		call pml_matrices_xycorner(1,acc,coef1,dn,dn,nyt,nxt,yd,xl,ld_factorA,ld_factorB,&
		deriv_p_x,deriv2_p_x,ld_phi_p_x,ld_phi_p_x_1,ld_chi_p_x,ld_chi_p_x_1,corner_prop_pml)
		ld_pml_p_3=-ld_pml_p_1+2*ld_pml_p+c2(yd,xl)*&
		(deriv2_p_x(yd,xl)+deriv2_p_y(yd,xl)+corner_prop_pml)
		
		call pml_matrices_xycorner(2,acc,coef1,dn,dn,nyt,nxt,yd,xl,ld_factorA,ld_factorB,&
		deriv_p_y,deriv2_p_y,ld_phi_p_y,ld_phi_p_y_1,ld_chi_p_y,ld_chi_p_y_1,corner_prop_pml)
		ld_pml_p_3=ld_pml_p_3+c2(yd,xl)*corner_prop_pml				

		call pml_matrices_xycorner(1,acc,coef1,dn,dn,nyt,nxt,yd,xr,rd_factorA,rd_factorB,&
		deriv_p_x,deriv2_p_x,rd_phi_p_x,rd_phi_p_x_1,rd_chi_p_x,rd_chi_p_x_1,corner_prop_pml)
		rd_pml_p_3=-rd_pml_p_1+2*rd_pml_p+c2(yd,xr)*&
		(deriv2_p_x(yd,xr)+deriv2_p_y(yd,xr)+corner_prop_pml)
		
		call pml_matrices_xycorner(2,acc,coef1,dn,dn,nyt,nxt,yd,xr,rd_factorA,rd_factorB,&
		deriv_p_y,deriv2_p_y,rd_phi_p_y,rd_phi_p_y_1,rd_chi_p_y,rd_chi_p_y_1,corner_prop_pml)
		rd_pml_p_3=rd_pml_p_3+c2(yd,xr)*corner_prop_pml


!!!		 Propagation through the medium

	        p_3(my,mx)=-p_1(my,mx)+2.*p(my,mx)+(deriv2_p_x(my,mx)+&
		deriv2_p_y(my,mx)+source(my,mx))*c2(my,mx)	!!check if this have to be before or after free surface

		!!!!! Add PML layers in propagator model

!		if(free_surface.eq.0)	then
		        p_3(yu,mx)=yu_pml_p_3                   !!capa "y" arriba
		        p_3(yu,xl)=lu_pml_p_3			!!esquina "x-y" izda arriba
		        p_3(yu,xr)=ru_pml_p_3			!!esquina "x-y" derecha arriba
!		endif

                p_3(yd,mx)=yd_pml_p_3                   !!capa "y" abajo
                p_3(yd,xl)=ld_pml_p_3            	!!esquina "x-y" izda abajo
                p_3(yd,xr)=rd_pml_p_3            	!!esquina "x-y" derecha abajo

                p_3(my,xl)=xl_pml_p_3                 	!!capa "x" izda
                p_3(my,xr)=xr_pml_p_3                 	!!capa "x" derecha

!!!             FREE SURFACE CALCULATION

                if(free_surface.ne.0)   then

	                p_3(1+dny2,:)=0.

                        !do j=1,dny2
                        !        yu_pml_p_3(j,:)=p_3(dny2+1+(1-j),mx)
                        !        lu_pml_p_3(j,:)=p_3(dny2+1+(1-j),xl)
                        !        ru_pml_p_3(j,:)=p_3(dny2+1+(1-j),xr)
                        !enddo

                endif


		!! Imaging Condition

		if(iprop.eq.2)	then

			if(l.eq.store_snap) then
			snap=snap-1
			if(snap.gt.nt_snap)write(*,*)'snap overflow 2',snap

			do i=1,nx
			do j=1,ny

				kernel(j,i)=kernel(j,i)+0.5*(p_3(dny2+j,dnx+i)-p_1(dny2+j,dnx+i))*Fields(j,i,snap)*store_snap
			enddo
			enddo

			l=0

                                !if(MOD(snap,50).eq.0) then
                                !op2=op2+1
                                !do j=1,ny
                                !        write(200+op2,'(20000(e12.5,2x))') (kernel(j,i),i=1,nx)
                                !enddo
                                !endif

			endif

		endif !iprop


		!!!!! Time actualization
	
		p_1=p
		p=p_3

		yu_pml_p_1=yu_pml_p
		yu_pml_p=yu_pml_p_3
		
		yd_pml_p_1=yd_pml_p
		yd_pml_p=yd_pml_p_3

		xl_pml_p_1=xl_pml_p
		xl_pml_p=xl_pml_p_3
		
		xr_pml_p_1=xr_pml_p
		xr_pml_p=xr_pml_p_3

		lu_pml_p_1=lu_pml_p
		lu_pml_p=lu_pml_p_3
	
		ru_pml_p_1=ru_pml_p
		ru_pml_p=ru_pml_p_3

		ld_pml_p_1=ld_pml_p
		ld_pml_p=ld_pml_p_3
						
		rd_pml_p_1=rd_pml_p
		rd_pml_p=rd_pml_p_3
			
		yu_phi_p_y_1=yu_phi_p_y
		yu_chi_p_y_1=yu_chi_p_y
		
		yd_phi_p_y_1=yd_phi_p_y
		yd_chi_p_y_1=yd_chi_p_y

		xl_phi_p_x_1=xl_phi_p_x		
		xl_chi_p_x_1=xl_chi_p_x
		xr_phi_p_x_1=xr_phi_p_x		
		xr_chi_p_x_1=xr_chi_p_x

		lu_phi_p_x_1=lu_phi_p_x		
		lu_phi_p_y_1=lu_phi_p_y
		lu_chi_p_x_1=lu_chi_p_x
		lu_chi_p_y_1=lu_chi_p_y
		
		ru_phi_p_x_1=ru_phi_p_x		
		ru_phi_p_y_1=ru_phi_p_y
		ru_chi_p_x_1=ru_chi_p_x
		ru_chi_p_y_1=ru_chi_p_y
						
		ld_phi_p_x_1=ld_phi_p_x		
		ld_phi_p_y_1=ld_phi_p_y
		ld_chi_p_x_1=ld_chi_p_x
		ld_chi_p_y_1=ld_chi_p_y
		
		rd_phi_p_x_1=rd_phi_p_x		
		rd_phi_p_y_1=rd_phi_p_y
		rd_chi_p_x_1=rd_chi_p_x
		rd_chi_p_y_1=rd_chi_p_y

	
	enddo
	

	!!!! Variables propagador

	deallocate(nxs,nys,nxreco,nyreco)
	deallocate(coef1,coef2)
	deallocate(vyt,vxt,my,mx,yd,xl,xr)
	deallocate(c2)
	deallocate(deriv2_p_x,deriv2_p_y)
	deallocate(deriv_p_x,deriv_p_y)
  	deallocate(source)
	deallocate(p,p_1,p_3)

	!!!! Variables PML 
	deallocate(factorA,factorB)
	deallocate(lu_factorA,lu_factorB,ru_factorA,ru_factorB)
	deallocate(ld_factorA,ld_factorB,rd_factorA,rd_factorB)
	deallocate(x_prop_pml,y_prop_pml,corner_prop_pml)
	
	!!!! Variables PML y	
	deallocate(yu_pml_p_1,yu_pml_p,yu_pml_p_3)
	deallocate(yu_phi_p_y_1,yu_phi_p_y,yu_chi_p_y_1,yu_chi_p_y)
	
	deallocate(yd_pml_p_1,yd_pml_p,yd_pml_p_3)
	deallocate(yd_phi_p_y_1,yd_phi_p_y,yd_chi_p_y_1,yd_chi_p_y)

	!!!! Variables PML x direction
	deallocate(xl_pml_p_1,xl_pml_p,xl_pml_p_3)
	deallocate(xl_phi_p_x_1,xl_phi_p_x,xl_chi_p_x_1,xl_chi_p_x)

	deallocate(xr_pml_p_1,xr_pml_p,xr_pml_p_3)
	deallocate(xr_phi_p_x_1,xr_phi_p_x,xr_chi_p_x_1,xr_chi_p_x)
	
	!!!! Variables PML corner up direction
	deallocate(lu_pml_p_1,lu_pml_p,lu_pml_p_3)
	deallocate(lu_phi_p_x_1,lu_phi_p_x,lu_chi_p_x_1,lu_chi_p_x)
	deallocate(lu_phi_p_y_1,lu_phi_p_y,lu_chi_p_y_1,lu_chi_p_y)

	deallocate(ru_pml_p_1,ru_pml_p,ru_pml_p_3)
	deallocate(ru_phi_p_x_1,ru_phi_p_x,ru_chi_p_x_1,ru_chi_p_x)
	deallocate(ru_phi_p_y_1,ru_phi_p_y,ru_chi_p_y_1,ru_chi_p_y)
	
	!!!! Variables PML corner down direction
	deallocate(ld_pml_p_1,ld_pml_p,ld_pml_p_3)
	deallocate(ld_phi_p_x_1,ld_phi_p_x,ld_chi_p_x_1,ld_chi_p_x)
	deallocate(ld_phi_p_y_1,ld_phi_p_y,ld_chi_p_y_1,ld_chi_p_y)

	deallocate(rd_pml_p_1,rd_pml_p,rd_pml_p_3)
	deallocate(rd_phi_p_x_1,rd_phi_p_x,rd_chi_p_x_1,rd_chi_p_x)
	deallocate(rd_phi_p_y_1,rd_phi_p_y,rd_chi_p_y_1,rd_chi_p_y)
	
end subroutine

