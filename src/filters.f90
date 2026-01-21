subroutine taper_pad_data(field,NumRec,nt_prop,N_taper,N_zeros,nt_new)
  implicit none

  integer :: i,j,k
  integer	:: NumRec
  integer	:: nt_prop
  integer	:: nt_new
  integer	:: N_taper
  integer	:: N_zeros
  real		:: field(nt_new,NumRec)

  do i=1,NumRec

	call taper_pad_trace(field(:,i),nt_prop,N_taper,N_zeros,nt_new)

  enddo

end subroutine taper_pad_data

subroutine taper_pad_trace(trace, nt_prop,  N_taper, N_zeros, nt_new)
  implicit none

  integer :: i0, n
  real :: pi
  real :: w
  integer	:: nt_prop
  integer	:: nt_new
  integer	:: N_taper
  integer	:: N_zeros
  real		:: trace(nt_new)

  pi = 4.0 * atan(1.0)

  ! apply half-cosine taper over last N_taper points of trace_in
  i0 = max(1, nt_prop - N_taper + 1)
  do n = 0, N_taper-1
     if (i0 + n <= nt_prop) then
        w = 0.5 * (1.0 + cos(pi * real(n) / real(N_taper)))
        trace(i0 + n) = trace(i0 + n) * w
     end if
  end do

  ! pad zeros after nt_prop
  if (N_zeros > 0) then
     trace(nt_prop+1:nt_new) = 0.0
  end if

end subroutine taper_pad_trace

function filters_function1D(S,ts,fe,Num,typef,f1,fc)
implicit none

	real :: f0, Q, fc, fs, f1, f2, G0, fc1, fc2, fs1, fs2, fe
	real :: f_1, f_2, f_1prima, f_2prima, f, ripple
	integer :: i, j, ts, n, typef, lensav, IER, lenwrk, filter, parity,ind,Num

	real, parameter :: pi = 3.14159265
	complex, parameter :: imag=(0.,1.)
	complex, allocatable :: wsave(:), work(:), H(:), Y(:) 
	real :: S(ts,Num),yf(ts),Stmp(ts)
	complex, allocatable :: Sc(:)
	CHARACTER(LEN=30) :: Format1
	logical :: filters_function1D
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!
	!!!!!!! Parameters for all filters (except do Ormsby, look the note bedoe the end of this section)
	!!
	!filter=3! Filter selected, 1 for ideal filter, 2 for ormsby, 3 for Butterworth, 4 for Chebyshev and 5 for Bessel
	n=6! Filter order
	!!!typef entra desde Inversion.f90 ! typef of filter, typef=1 for low-pass, typef=2 for high-pass, typef=3 for pass-band
	parity=0 ! parity=0 if n is even and 1 if n is odd, just useful in Butterworth filter

	!fc : Cut-off frecuency (in Hz) do low-pass
	fs=40!40! Start frecuency for high-pass (in Hz)
	!f1 Start frecuency for pass-band filter (in Hz)
	f2=fc  ! Cut-off frecuency for pass-band filter (in Hz)
	G0=1 ! Maximum gain for all typef of filter

	!!!!!------- Other parameters, calculus
	
	f0=sqrt(f1*f2)
	Q=f0/(f2-f1)

	!!!!!------- FFT caculation do data input (data with noise, S)
	
	!! Calling the filter selected
        call butterworth(n,parity,fc,fs,f1,f2,f,fe,f0,Q,ts,typef,H)	!!!solo se usa este
	
	lensav=2*ts+ceiling(log(ts/1.))+4
	lenwrk=2*ts
	allocate(wsave(lensav),work(lenwrk),Sc(ts))
	allocate(Y(1:ts))
	do ind=1,Num
		
		do i=1,ts
			Stmp(i)=S(i,ind)
			Sc(i)=Stmp(i)
		end do
	
		!!transformada Fourier
		call cfft1i(ts,wsave,lensav,IER)
		call cfft1f(ts,1,Sc,ts,wsave,lensav,work,lenwrk,IER)
		
		do i=1,ts
			Y(i)=Sc(i)*H(i)
		end do

		call cfft1b(ts,1,Y,ts,wsave,lensav,work,lenwrk,IER)	!!inversa

		do i=1,ts
			S(i,ind)=(Y(i))
		end do
	
	enddo
	
	!! Saving the results

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	filters_function1D = .true. !logical output to make the compiler stop

	CONTAINS	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine butterworth(n,parity,fc,fs,f1,f2,f,fe,f0,Q,ts,typef,H)
	real :: f, fc, fs, f1, f2, fe, f0, Q 
	complex, allocatable :: H(:), Bn(:)
	integer :: typef, j, ts, n, parity ,k
	real, parameter :: pi = 3.14159265
	complex, parameter :: imag=(0.,1.)
	allocate(H(1:ts), Bn(1:ts))
	f=0.
	do i=1,ts
	Bn(i)=1.
	end do
	if(typef.eq.1) then
		if(parity.eq.0) then
		    do j=1,ts/2+1
		        do k=1,n/2
		            Bn(j)=Bn(j)*((imag*f/fc)**2-2*(imag*f/fc)*cos((2*k+n-1)*pi/(2*n))+1)
		        end do
		    f=f+fe/ts
		    end do
		else
		    if(n.le.1) then
		        do j=1,ts/2+1
		            Bn(j)=(imag*f/fc)+1
		            f=f+fe/ts
		        end do
		    else 
		        do j=1,ts/2+1
		            do k=1,(n-1)/2
		                Bn(j)=Bn(j)*((imag*f/fc)+1)*((imag*f/fc)**2-2*(imag*f/fc)*cos((2*k+n-1)*pi/(2*n))+1)
		            end do
		        f=f+fe/ts
		        end do
		    end if
		end if
	else                
		! We change in low-pass filter i*w/wc do wc'/(i*w), where wc'=ws,
		! in order to get a high-pass filter
		if(typef.eq.2) then ! High-pass filter
		    if (parity.eq.0) then
		        do j=1,ts/2+1
		            do k=1,n/2
		                Bn(j)=Bn(j)*((fs/(imag*f))**2-2*(fs/(imag*f))*cos((2*k+n-1)*pi/(2*n))+1)
		            end do
		        f=f+fe/ts
		        end do
		    else
		        if(n.eq.1) then
		            do j=1,ts/2+1
		                Bn(j)=(fs/(imag*f))+1
		                f=f+fe/ts
		            end do
		        else
		            do j=1,ts/2+1
		                do k=1,(n-1)/2
		                    Bn(j)=Bn(j)*((fs/(imag*f))+1)*((fs/(imag*f))**2-2*(fs/(imag*f))*cos((2*k+n-1)*pi/(2*n))+1)
		                end do
		            f=f+fe/ts
		            end do
		        end if
		    end if      
		else ! Band-pass filter
		! We change in low-pass filter i*w/wc do Q(iw/w0+w0/(iw)), where 
		!w0=square(w3w4), incw=w4-w3 and Q=w0/incw, in order to get a band-pass filter
		    if (parity.eq.0) then
		        do j=1,ts/2+1
		            do k=1,n/2
		                Bn(j)=Bn(j)*((Q*((imag*f)/f0+f0/(imag*f)))**2-2*(Q*((imag*f)/f0+f0/(imag*f)))&
				*cos((2*k+n-1)*pi/(2*n))+1)
		            end do
		        f=f+fe/ts
		        end do
		    else
		        if(n.eq.1) then
		            do j=1,ts/2+1
		                Bn(j)=(Q*((imag*f)/f0+f0/(imag*f)))+1
		                f=f+fe/ts                  
		            end do
		        else
		            do j=1,ts/2+1
		                do k=1,(n-1)/2
		                    Bn(j)=Bn(j)*((Q*((imag*f)/f0+f0/(imag*f)))+1)*((Q*((imag*f)/f0+f0/(imag*f)))**2&
					-2*(Q*((imag*f)/f0+f0/(imag*f)))*cos((2*k+n-1)*pi/(2*n))+1)
		                end do
		            f=f+fe/ts
		            end do
		        end if
		    end if
		end if
	end if
	do j=1,ts/2+1
	  H(j)=(1/Bn(j))
	end do

	if (typef .eq. 1) then
		H(1)=1.
	else
		H(1)=0.
	endif
	
	do j=2,ts/2
	  H(ts/2+j)=(real(H(ts/2+1-(j-1)))-imag*aimag(H(ts/2+1-(j-1))))
	end do

end subroutine butterworth

end function filters_function1D
