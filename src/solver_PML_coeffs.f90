	subroutine PML_coeffs(acc,dPML,dt,d,dn,factorA,factorB,&
	ld_factorA,rd_factorA,lu_factorA,ru_factorA,ld_factorB,rd_factorB,lu_factorB,ru_factorB)


	implicit none

	integer dn,dPML,acc,i,j,ix
	real pi,dt,d,Refl,wPML
	real factorA(dn),factorB(dn)
	real lu_factorA(dn,dn),lu_factorB(dn,dn),ru_factorA(dn,dn),ru_factorB(dn,dn)
	real ld_factorA(dn,dn),ld_factorB(dn,dn),rd_factorA(dn,dn),rd_factorB(dn,dn)
	
	real, allocatable :: yPML(:),alpha(:),sigma(:)
	allocate(yPML(dn),alpha(dn),sigma(dn))

	pi=3.14159265			
	Refl=0.00001
	wPML=(dPML-1)*d

	do i=1+acc,dPML+acc	!!from interior to exterior
		yPML(i)=(i-(acc+1))*d
		alpha(i)=30.*pi*(1.-(yPML(i)/wPML)**1)
		sigma(i)=(-(2+1)*2000*log(Refl)/(2.*wPML))*(yPML(i)/wPML)**2
	enddo	

	do i=1+acc,dPML+acc	!!from interior to exterior
		factorA(i)=exp(-(sigma(i)+alpha(i))*dt)
		factorB(i)=(factorA(i)-1)*sigma(i)/(sigma(i)+alpha(i))
	enddo

	do i=dPML+acc+1,dn	!!from interior to exterior
		factorA(i)=factorA(dPML+acc)
		factorB(i)=factorB(dPML+acc)	
	enddo
		
	do i=1,dn
		do j=1,i
			rd_factorA(j,i)=factorA(i)
			rd_factorB(j,i)=factorB(i)
		enddo
	enddo
	do i=1,dn
		if(i+1.le.dn)then
			do j=i+1,dn
				rd_factorA(j,i)=rd_factorA(i,j)	
				rd_factorB(j,i)=rd_factorB(i,j)	
			enddo
		endif
	enddo
		
	do i=1,dn
		do j=1,dn-i+1
			ld_factorA(j,i)=factorA(dn-i+1)
			ld_factorB(j,i)=factorB(dn-i+1)
		enddo
	enddo
	do ix=1,dn-1
		j=2+(ix-1)
		do i=dn-(ix-1),dn
			ld_factorA(j,i)=factorA(j)
			ld_factorB(j,i)=factorB(j)
		enddo
	enddo

	do j=1,dn
		lu_factorA(j,:)=ld_factorA(dn-j+1,:)		
		lu_factorB(j,:)=ld_factorB(dn-j+1,:)	
		ru_factorA(j,:)=rd_factorA(dn-j+1,:)		
		ru_factorB(j,:)=rd_factorB(dn-j+1,:)	
	enddo
	
	deallocate(yPML,alpha,sigma)
	
	return
	end
