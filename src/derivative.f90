	subroutine derivative(n,x,y,dy)
  
	integer :: n
	real :: x(n),y(n),dy(n)

! compute derivation
	do i = 1, n-1
		dy(i) = (y(i+1) - y(i)) / (x(i+1) - x(i))
	enddo

! computing last derivation using linear extrapolation
	dy(n) = dy(n-1) + (dy(n-1)-dy(n-2))/(x(n-1)-x(n-2))*(x(n)-x(n-1))

	return
	end
