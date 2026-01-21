function interpolation_akima(Data_tmp,dt_Data,nt_Data,Data_trace,NumRec,dt,nt)
implicit none

	logical interpolation_akima
	integer nt_Data,nt,NumRec,i,j,k,ERR
	real :: dt_Data,dt,Time_Data(nt_Data),time(nt)
	real :: Data_trace(nt,NumRec),Data_tmp(nt_Data,NumRec)
	real :: Data_1(nt_Data),Data_2(nt)

	ERR=0

	do i=1,nt_Data
		Time_Data(i)=(i-1)*dt_Data
	end do

	do i=1,nt
		time(i)=(i-1)*dt
	end do

	do k=1,NumRec

		Data_1(:)=Data_tmp(:,k)
		call  INTRPL(nt_Data,Time_Data,Data_1,nt,time,Data_2,ERR)
		Data_trace(:,k)=Data_2(:)

	enddo

	interpolation_akima=1
	
end function interpolation_akima

