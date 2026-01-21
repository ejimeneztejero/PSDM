subroutine first_arrivals(Data_trace,nt,ntt)
implicit none		

	integer :: nt,FAj
	integer :: j,k,ntt,l
	real :: Data_trace(nt)
        real, allocatable :: T(:) 

	allocate(T(nt))

!!!	PRIMERAS LLEGADAS

	T=Data_trace/maxval(abs(Data_trace(:)))

	call F_A(T,nt,ntt)!FAT, el tiempo donde la señal comienza a no ser cero

	deallocate(T)


return
end

subroutine F_A(S,nt,FA)
implicit none

        integer :: nt,k,FA
        real :: umbral,S(nt)            !!! tiene que estar normalizada.

        umbral=0.001
        do k=3,nt-2

                if(abs(S(k)).ge.umbral) then
                if(abs(S(k+1)).ge.umbral.and.abs(S(k-1)).ge.umbral)     then
                if(abs(S(k+2)).ge.umbral.and.abs(S(k-2)).ge.umbral)     then
                        FA=k
                        goto 22
                endif
                endif
                endif

        enddo
        22 continue


return
end

