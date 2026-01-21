	subroutine pml_matrices_x(idir,acc,coef1,nvy,nvx,nyt,nxt,vy,vx,factorA,factorB,&
	deriv,deriv2,phi,phi_1,chi,chi_1,prop_pml)

	implicit none

	integer :: idir,acc,coef1,nvy,nvx,nxt,nyt,i,j,vx(nvx),vy(nvy)
	real :: deriv(nyt,nxt),deriv2(nyt,nxt),factorA(nvx),factorB(nvx)
	real :: prop_pml(nvy,nvx),phi_1(nvy,nvx),phi(nvy,nvx),chi_1(nvy,nvx),chi(nvy,nvx)
	real, allocatable :: deriv_phi(:,:)
	allocate(deriv_phi(nvy,nvx))	

	prop_pml=0.;deriv_phi=0.;

	call space_derivx(acc,coef1,nvy,nvx,phi,deriv_phi)
           
        if(idir.eq.1)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(i)*phi_1(j,i)+factorB(i)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(i)*chi_1(j,i)+factorB(i)*(deriv2(vy(j),vx(i))+deriv_phi(j,i))	
	enddo
	enddo
	endif
	
        if(idir.eq.2)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(nvx-i+1)*phi_1(j,i)+factorB(nvx-i+1)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(nvx-i+1)*chi_1(j,i)+factorB(nvx-i+1)*(deriv2(vy(j),vx(i))+deriv_phi(j,i))	
	enddo
	enddo
	endif

	prop_pml=deriv_phi+chi
		
	deallocate(deriv_phi)	
				
	return
	end
	
	subroutine pml_matrices_y(idir,acc,coef1,nvy,nvx,nyt,nxt,vy,vx,factorA,factorB,&
	deriv,deriv2,phi,phi_1,chi,chi_1,prop_pml)

	implicit none

	integer :: idir,acc,coef1,nvy,nvx,nxt,nyt,i,j,vx(nvx),vy(nvy)
	real :: deriv(nyt,nxt),deriv2(nyt,nxt),factorA(nvy),factorB(nvy)
	real :: prop_pml(nvy,nvx),phi_1(nvy,nvx),phi(nvy,nvx),chi_1(nvy,nvx),chi(nvy,nvx)
	
	real, allocatable :: deriv_phi(:,:)
	allocate(deriv_phi(nvy,nvx))	

	prop_pml=0.;deriv_phi=0.;

	call space_derivy(acc,coef1,nvy,nvx,phi,deriv_phi)
                       
        if(idir.eq.1)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(j)*phi_1(j,i)+factorB(j)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(j)*chi_1(j,i)+factorB(j)*(deriv2(vy(j),vx(i))+deriv_phi(j,i))	
	enddo
	enddo
	endif
	
        if(idir.eq.2)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(nvy-j+1)*phi_1(j,i)+factorB(nvy-j+1)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(nvy-j+1)*chi_1(j,i)+factorB(nvy-j+1)*(deriv2(vy(j),vx(i))+deriv_phi(j,i))	
	enddo
	enddo
	endif
	
	prop_pml=deriv_phi+chi
		
	deallocate(deriv_phi)	
				
	return
	end
	
	subroutine pml_matrices_xycorner(id,acc,coef1,nvy,nvx,nyt,nxt,vy,vx,factorA,factorB,&
	deriv,deriv2,phi,phi_1,chi,chi_1,prop_pml)

	implicit none

	integer :: id,acc,coef1,nvy,nvx,nxt,nyt,i,j,vx(nvx),vy(nvy)
	real :: deriv(nyt,nxt),deriv2(nyt,nxt),factorA(nvy,nvx),factorB(nvy,nvx)
	real :: prop_pml(nvy,nvx),phi_1(nvy,nvx),phi(nvy,nvx),chi_1(nvy,nvx),chi(nvy,nvx)
	
	real, allocatable :: deriv_phi(:,:)
	allocate(deriv_phi(nvy,nvx))	

	prop_pml=0.;deriv_phi=0.;
	
	if(id.eq.1)call space_derivx(acc,coef1,nvy,nvx,phi,deriv_phi)
	if(id.eq.2)call space_derivy(acc,coef1,nvy,nvx,phi,deriv_phi)
	               
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(j,i)*phi_1(j,i)+factorB(j,i)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(j,i)*chi_1(j,i)+factorB(j,i)*(deriv2(vy(j),vx(i))+deriv_phi(j,i))	
	enddo
	enddo
	
	prop_pml=deriv_phi+chi
		
	deallocate(deriv_phi)	
				
	return
	end
	
