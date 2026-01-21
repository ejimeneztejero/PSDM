subroutine finite_difference_coef(num,accuracy,coef)
implicit none

integer	num,accuracy
real coef(accuracy+1)

	if(num.eq.1)	then

        if(accuracy.eq.1)  then
		coef(1)=0.
		coef(2)=1/2
	endif

        if(accuracy.eq.3)  then
                coef(1)=0.
                coef(2)=3./4.
                coef(3)=-3./20.
                coef(4)=1./60.  
        endif

	if(accuracy.eq.4)  then
                coef(1)=0.
                coef(2)=4./5.
                coef(3)=-1./5.
                coef(4)=4./105.
                coef(5)=-1./280.
        endif

	if(accuracy.eq.5)  then
                coef(1)=0.
                coef(2)=0.
                coef(3)=0.
                coef(4)=0.
                coef(5)=0.
                coef(6)=0.
        endif

	endif


	if(num.eq.2)	then

        if(accuracy.eq.3)  then
!                coef(1)=-3.+3./10.-1./40. !!modificacion Jan
                coef(1)=-49./18.
                coef(2)=3./2.
                coef(3)=-3./20.
!                coef(4)=1./80.		!!modificacion Jan
                coef(4)=1./90.                  
        endif

	if(accuracy.eq.4)  then
                coef(1)=-205./72.     
                coef(2)=8./5.
                coef(3)=-1./5.
                coef(4)=8./315.     
                coef(5)=-1./560.
        endif

	if(accuracy.eq.5)  then
                coef(1)=-5269./1800.    
                coef(2)=5./3.
                coef(3)=-5./21.
                coef(4)=5./126.     
                coef(5)=-5./1008.
                coef(6)=1./3150.
        endif
	coef(1)=coef(1)/2.
	endif


return
end
