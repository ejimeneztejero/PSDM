subroutine space_deriv2x(acc,coef,nyt,nxt,p,deriv2_px)

	integer nthreads
	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2,k
	real p(nyt,nxt),deriv2_px(nyt,nxt)
	real coef(acc+1)

	i1=1+acc
	i2=nxt-acc

        deriv2_px=0.

        do k=1,acc+1

	do i=i1,i2
	do j=1,nyt
	        deriv2_px(j,i)=deriv2_px(j,i)+(p(j,i+(k-1))+p(j,i-(k-1)))*coef(k)
        enddo
        enddo

	enddo

        do i=1,i1-1
	do j=1,nyt
                deriv2_px(j,i)=deriv2_px(j,acc+1)
        enddo
        enddo

	do i=i2+1,nxt
	do j=1,nyt
                deriv2_px(j,i)=deriv2_px(j,nxt-acc)
        enddo
        enddo

return
end
	
subroutine space_deriv2y(acc,coef,nyt,nxt,p,deriv2_py)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2,k
	real p(nyt,nxt),deriv2_py(nyt,nxt)
	real coef(acc+1)

	j1=1+acc
	j2=nyt-acc

	deriv2_py=0.

	do i=1,nxt

		do j=j1,j2
		do k=1,acc+1
			deriv2_py(j,i)=deriv2_py(j,i)+(p(j+(k-1),i)+p(j-(k-1),i))*coef(k)
		enddo
		enddo

		do j=1,j1-1
			deriv2_py(j,i)=deriv2_py(acc+1,i)
		enddo

		do j=1+j2,nyt
			deriv2_py(j,i)=deriv2_py(nyt-acc,i)
		enddo
	enddo

return
end

subroutine space_derivx(acc,coef,nyt,nxt,p,deriv_px)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2
	real p(nyt,nxt),deriv_px(nyt,nxt)
	real coef(acc+1)

	i1=1+acc
	i2=nxt-acc

	deriv_px=0.

	do k=1,acc+1

	do i=i1,i2
	do j=1,nyt
		deriv_px(j,i)=deriv_px(j,i)+(p(j,i+(k-1))-p(j,i-(k-1)))*coef(k)
	enddo
	enddo

	enddo				    

	do i=1,i1-1
	do j=1,nyt
		deriv_px(j,i)=deriv_px(j,acc+1)
	enddo
	enddo				    

	do i=1+i2,nxt
	do j=1,nyt
		deriv_px(j,i)=deriv_px(j,nxt-acc)
	enddo				    
	enddo				    


return
end
		
subroutine space_derivy(acc,coef,nyt,nxt,p,deriv_py)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2
	real p(nyt,nxt),deriv_py(nyt,nxt)
	real coef(acc+1)
	
	j1=1+acc
	j2=nyt-acc

	deriv_py=0.

	do i=1,nxt

		do j=j1,j2
		do k=1,acc+1
			deriv_py(j,i)=deriv_py(j,i)+(p(j+(k-1),i)-p(j-(k-1),i))*coef(k)
		enddo
		enddo

		do j=1,j1-1
			deriv_py(j,i)=deriv_py(acc+1,i)
		enddo

		do j=1+j2,nyt
			deriv_py(j,i)=deriv_py(nyt-acc,i)
		enddo

	enddo
	
return
end
