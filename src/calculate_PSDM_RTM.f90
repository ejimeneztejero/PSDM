subroutine calculate_PSDM_RTM(Vpmodel,Batmodel,new_comm,PSDM_depth)

use mod_parfile
use mod_data_arrays
use mpi

implicit none

integer :: numtasks,rank, ierr, TAG
integer :: new_comm
integer :: itimes,icount,ntimes
integer :: nsamples,last,trid
integer :: i,j,k,kk,kkk,l,it
integer :: NumRec
integer :: nt_data
integer :: nt_taper,N_zeros,N_taper,nt_prop2

integer :: SourceNum
integer :: minS,maxS,minR,maxR
integer :: x_i,x_f,y_i,y_f,nx_w,nx_w_max

integer :: mid_x
real    :: offset,D_src,D_rec,leg
real	:: twt_sf,twt_multiple

real :: tt_tot,twt_bat
real :: ddx,f1,f2
real :: tfin,dt_data,max_data

real :: Vpmodel(ny_prop,nx_prop)
real :: PSDM_depth(ny_prop,nx_prop)
real :: Batmodel(nx_prop)

integer :: nxSou(1),nySou(1)

integer, allocatable :: nxRec(:),nyRec(:)
integer, allocatable :: nt_cut(:)
integer, allocatable :: nt_sf(:)

real, allocatable :: Cov(:,:),Cov_sum(:,:)
real, allocatable :: Cov_w(:,:),Cov_thread(:,:),Cov_itimes(:,:)

real, allocatable :: PSDM_sum(:,:)
real, allocatable :: PSDM_w(:,:),PSDM_thread(:,:),PSDM_itimes(:,:)

real, allocatable :: Data_source(:)
real, allocatable :: Data_synth(:,:)
real, allocatable :: Fields(:,:,:)
real, allocatable :: Data_target(:,:)
real, allocatable :: Field_Data(:,:),Raw_Data(:,:),Delay_Data(:,:)
real, allocatable :: BackSource(:,:)
real, allocatable :: Vp_w(:,:)

character(len=1000) :: file_name,str_snap,file_V_twt,file_P_twt

call MPI_COMM_SIZE(new_comm,numtasks,ierr)
call MPI_COMM_RANK(new_comm,rank,ierr)
TAG=0

allocate(PSDM_sum(ny_prop,nx_prop)) 
allocate(Cov_sum(ny_prop,nx_prop)) 
Cov_sum=0;PSDM_sum=0;

nt_prop2 = nt_prop
nt_taper = nt_prop2

if(rank.eq.0)   then
        write(*,*)
        write(*,*)"oooooooooooooooooooooooooooooooooooooooooooooooooooooo"
        write(*,*)'CALCULATE PSDM'
        write(*,*)"oooooooooooooooooooooooooooooooooooooooooooooooooooooo"
        write(*,*)
endif

if(rtm_mcs.ne.0)	then

	if(shot_i.eq.0)then
		shot_i=shotID_nav_mcs(1)
	endif

	if(shot_f.eq.0)then
		shot_f=shotID_nav_mcs(NumShots_mcs)
	endif

	nt_data=nt_mcs
	dt_data=dt_mcs
	tfin=tfin_mcs
	NumRec=NumRec_MCS

endif

if(rtm_was.ne.0)	then

	if(obs_i.eq.0)then
		obs_i=obsID_nav(1)
	endif
	if(obs_f.eq.0)then
		obs_f=obsID_nav(NumOBS)
	endif

	nt_data=nt_was
	dt_data=dt_was
	tfin=tfin_was
	NumRec=NumRec_WAS

endif

nsamples=line_f-line_i+1

if(rank.eq.0)write(*,*) "Allocating RTM dimensions"

allocate(nt_cut(NumRec))
allocate(nt_sf(NumRec))
allocate(nxRec(NumRec),nyRec(NumRec))

allocate(Data_source(nt_prop)) 
allocate(Data_synth(nt_prop,NumRec))

allocate(Field_Data(nt_data,NumRec))
allocate(Raw_Data(nt_data,NumRec))
allocate(Delay_Data(nt_data,NumRec))

if(cut_multiple.eq.0)nt_cut=nt_prop;
if(cut_sf.eq.0)nt_sf=1;

nxRec=0;nyRec=0;
nxSou=-10;nySou=-10;
nx_w=1;x_i=1;x_f=1;

Data_Source=0
Data_synth=0
Field_Data=0
Raw_Data=0
Delay_Data=0

!!!!!-----------------------------------
!!!!!-----------------------------------
!!!!!-----------------------------------

!! Fuente sintetica a la resolucion nt_prop, dt_prop
if(rtm_was.ne.0)call Ricker_was(rank,Data_Source)
if(rtm_mcs.ne.0)call Ricker_mcs(rank,Data_Source)

if(f_init.gt.0.and.f_final.gt.0) then

	if(rank.eq.0)write(*,*) "Filtering Source"
	call freq_filter(Data_Source,nt_prop,dt_prop,1,typef,f_init,f_final)

endif

ntimes=nsamples
if(numtasks.gt.1)    then
        if(nsamples.gt.numtasks)        then
                ntimes=ceiling(1.*nsamples/numtasks)
        endif
        if(nsamples.le.numtasks)        then
               ntimes=1
        endif
endif

if(rank.eq.0)write(*,*)"Init and final lines ... : ",line_i," to ",line_f
if(rank.eq.0)write(*,*)"Init and final shots ... : ",shotID_nav_mcs(line_i)," to ",shotID_nav_mcs(line_f)

l=0

do itimes=1,ntimes

icount=(itimes-1)*numtasks+rank+line_i !! paralelizacion

!!write(*,*)"kk ",rank,icount,line_i,line_f,shotID_nav_mcs(icount),shot_i,shot_f,NumShots_mcs

if(rank.eq.0.and.ntimes.gt.1)   then
        write(*,*)
        write(*,*)"oooooooooooooooooooooooooooooooooooooooooooooooooooooo"
        write(*,*)'ROUND ',itimes, 'OUT OF', ntimes
        write(*,*)"oooooooooooooooooooooooooooooooooooooooooooooooooooooo"
        write(*,*)
endif

if(icount.ge.line_i.and.icount.le.line_f)  then

	if(itimes.lt.ntimes) last=icount+numtasks-1
	if(itimes.eq.ntimes) last=line_f

!!	write(*,*)"OO ",rank,icount,shotID_nav_mcs(icount),shotID_su_mcs(icount),shotID_nav_mcs(last),shot_f

	if(rank.eq.0) then
		write(*,*)"Shots from: ",shotID_nav_mcs(icount)," to ",shotID_nav_mcs(last)
		write(*,*)"Reading lines from: ",icount," to ",last
	endif

	if(rtm_WAS.ne.0)    then ! WAS data

	    nxSou=nxSou_WAS(icount)
	    nySou=nySou_WAS(icount)
	    nxRec=nxRec_WAS(:,icount)
	    nyRec=nyRec_WAS(:,icount)

	    call open_obs(rank,icount)		    !! opens unit file for WAS data
	    call get_obs_data(icount,Field_data)    !! raw_data SU-WAS files
	    call close_was_files(icount)
		
	endif

	if(rtm_MCS.ne.0)    then ! MCS data

	    nxSou=nxSou_MCS(icount)
	    nySou=nySou_MCS(icount)
	    nxRec=nxRec_MCS(:,icount)
	    nyRec=nyRec_MCS(:,icount)

	    call open_shot_mcs(rank,icount)       	    !! opens unit file for MCS data
	    call get_shot_data(icount,Raw_Data,Delay_Data)  !! raw_data SU-MCS files		
	    call close_mcs_files(icount)

	    if(delay_time_data.ne.0) then
		Field_Data=Delay_Data
	    else	
		Field_Data=Raw_Data
	    endif
	
!!           write(23,*)icount,nxSou
!!           write(23,*)icount,nxRec(1),nxRec(NumRec)

	endif

!!!!!-----------------------------------
!!!!!-----------------------------------
!!!!!-----------------------------------

if(cut_multiple.gt.0)	then

!	!! calculo nt maximo para cada traza, cortar antes del multiple usando la info de Batmodel

	do j=1,NumRec

                ! Offset absoluto (m)
                offset = abs( nxSou(1)*dmodel_prop - nxRec(j)*dmodel_prop )

                ! Midpoint horizontal (puntos)
                mid_x = ceiling(0.5*(nxSou(1) + nxRec(j))) !!adimensional

                ! Distancias verticales efectivas al fondo (metros)
                D_src =  Batmodel(mid_x)  - nySou(1)*dmodel_prop
                D_rec =  Batmodel(mid_x)  - nyRec(j)*dmodel_prop

                ! Aproximacion simetrica del trayecto del multiple
                leg = sqrt( (offset*0.5)**2 + (0.5*(D_src + D_rec))**2 )

                twt_multiple = 4.0 * leg / 1500.
                nt_cut(j) = ceiling( 0.9*twt_multiple/dt_prop)	!!lo corto un poco antes por si acaso

	enddo

endif

if(cut_sf.gt.0)	then

!	!! calculo nt maximo para cada traza, cortar antes del multiple usando la info de Batmodel

	do j=1,NumRec

                ! Offset absoluto (m)
                offset = abs( nxSou(1)*dmodel_prop - nxRec(j)*dmodel_prop )

                ! Midpoint horizontal (puntos)
                mid_x = ceiling(0.5*(nxSou(1) + nxRec(j))) !!adimensional

                ! Distancias verticales efectivas al fondo (metros)
                D_src =  Batmodel(mid_x)  - nySou(1)*dmodel_prop
                D_rec =  Batmodel(mid_x)  - nyRec(j)*dmodel_prop

                ! Aproximacion simetrica del trayecto del multiple
                leg = sqrt( (offset*0.5)**2 + (0.5*(D_src + D_rec))**2 )

                twt_sf = 2.0 * leg / 1500.
                nt_sf(j) = ceiling(1.1* twt_sf/dt_prop)       !!lo corto un poco antes por si acaso

	enddo

endif

!! suavizar el corte brusco a cero para cada receptor (i)
!! taper and padding, time limits per trace

if(taper_option.ne.0.or.cut_multiple.ne.0)	then
	N_taper = max(1, int(T_taper / dt_prop))
	N_zeros = max(0, int(T_zeros / dt_prop))
	nt_prop2= min(nt_prop, maxval(nt_cut))
	nt_taper = nt_prop2 + N_zeros
	!write(77,*)i,nt_prop,nt_prop2,n_zeros,nt_taper
endif

allocate(Data_target(nt_taper,NumRec))
allocate(BackSource(nt_taper,NumRec))
Data_target=0
BackSource=0

!!!!!-----------------------------------
!!!!!------ Frequency filter of source and data
!!!!!-----------------------------------

	if(f_init.gt.0.and.f_final.gt.0) then

	    if(rank.eq.0)       then
	        write(*,*) "Filtering Field Data"
	        if(typef.eq.1)write(*,*) 'Lowpass at FREQ = ',f_init
	        if(typef.eq.3)write(*,*) 'Bandpass between FREQ = ',f_init,' and ',f_final
	    endif

	    call freq_filter(Field_Data,nt_data,dt_data,NumRec,typef,f_init,f_final)

	endif

!	write(*,*) "ojo 1"
!!!!!-----------------------------------
!!!!!-----------------------------------
!!!!!-----------------------------------

	if(dt_data.ne.dt_prop)	then
	        call interpolation_akima(Field_Data,dt_data,nt_data,Data_target(1:nt_prop2,:),NumRec,dt_prop,nt_prop2)
	else
		Data_target(1:nt_prop2,:)=Field_Data(1:nt_prop2,:)
	endif

!	write(*,*) "ojo 2"

!!!!!-----------------------------------
!!!!!-----------------------------------
!!!!!-----------------------------------

        open(unit=12,file="Target.dat",status='unknown')
        do j=1,nt_prop2
	        write(12,'(20000(e12.5,2x))') (Data_target(j,i),i=1,NumRec)
        enddo

	if(cut_sf.gt.0) then
		do i=1,NumRec
                        Data_target(1:nt_sf(i),i)=0.    !! para sacar el ruido del agua
		enddo
	endif

	if(cut_multiple.eq.0) then

		if(taper_option.ne.0)call taper_pad_data(Data_target,NumRec,nt_prop2,N_taper,N_zeros,nt_taper)

	else if(cut_multiple.gt.0) then

		!! cortar datos justo antes del multiple
		do i=1,NumRec

			Data_target(nt_cut(i):nt_taper,i)=0.			

			call taper_pad_trace(Data_target(:,i),nt_cut(i),N_taper,N_zeros,nt_taper)

		enddo

	endif

        open(unit=12,file="Target_window.dat",status='unknown')
        do j=1,nt_taper
	        write(12,'(20000(e12.5,2x))') (Data_target(j,i),i=1,NumRec)
        enddo


!write(*,*)"PARO YA"
!stop

!	write(*,*) "ojo 3"

!!!!!-----------------------------------
!!!!!-----------------------------------
!!!!!-----------------------------------

endif	!!icount

if(nxSou(1).ge.0) then
	call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
	dmodel_prop,nx_prop,added_space,nx_w,x_i,x_f)
endif

allocate(PSDM_w(ny_prop,nx_w)) 
allocate(Cov_w(ny_prop,nx_w)) 
PSDM_w=0;Cov_w=0;

if(icount.ge.line_i.and.icount.le.line_f)  then

	if(rank.eq.0)	then
		write(*,*) "Dynamic window in X (ptos) (nx_w,x_i,x_f): ",nx_w,x_i,x_f
	endif

	if(rank.eq.0)write(*,*) "Allocating Field dimensions: ",ny_prop,nx_w,nt_taper," : "
	if(rank.eq.0)write(*,*) "YOU NEED PER CPU: ",ceiling(4.0d0*ny_prop*nx_w*nt_taper/1e9)," GB of RAM"

	if(strategy_num.eq.1)	then
		allocate(Fields(ny_prop,nx_w,nt_taper))
		Fields=0
	endif

        allocate(Vp_w(ny_prop,nx_w))
        Vp_w=0.;
        Vp_w(:,1:nx_w)=Vpmodel(:,x_i:x_f)

        !!!!!---------------------          
        !!!!!------ Forward solver
        !!!!!---------------------

        if(rank.eq.0)write(*,*)
        if(rank.eq.0)write(*,*)"Forward solver ..."

	if(strategy_num.eq.1)	then

		call dynamic_window_solver_forward_ram(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
		Data_source,dmodel_prop,dt_prop,ny_prop,nx_prop,nx_w,Vp_w,&
		added_space,store_snap,nt_snap,nt_prop2,Fields(:,:,1:nt_prop2),Data_synth)
		
!		do k=1,nt_taper
!			write(33,*) (k-1)*dt_prop,Fields(ny_prop/2,nx_w/2,k)
!		enddo

		do i=1,nx_w
		do j=1,ny_prop

		if(taper_option.ne.0) call taper_pad_trace(Fields(j,i,:),nt_prop2, N_taper, N_zeros,nt_taper)

		enddo
		enddo

!		do k=1,nt_taper
!			write(34,*) (k-1)*dt_prop,Fields(ny_prop/2,nx_w/2,k)
!		enddo

	endif

	if(strategy_num.eq.2)	then
		call dynamic_window_solver_forward_disk(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
		Data_source,dmodel_prop,dt_prop,ny_prop,nx_prop,nx_w,Vp_w,&
		icount,unit_fields,folder_fields,&
		added_space,store_snap,nt_snap,nt_prop2,Data_synth)
	endif

!!!!!----------------------
!!!!!------ Backward propagation
!!!!!----------------------

!!	this needs to me called with nt_new

        do j=1,NumRec
		BackSource(:,j) = Data_target(nt_taper:1:-1,j)
        enddo

       	if(rank.eq.0) write(*,*)
       	if(rank.eq.0) write(*,*)"Propagation of MCS data back in time ..."

	if(strategy_num.eq.1)	then
		call dynamic_window_solver_backward_ram(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
		BackSource,dmodel_prop,dt_prop,ny_prop,nx_prop,nx_w,Vp_w,&
		added_space,store_snap,nt_taper,nt_taper,Fields,PSDM_w)
	endif

	if(strategy_num.eq.2)	then
		call dynamic_window_solver_backward_disk(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
		BackSource,dmodel_prop,dt_prop,ny_prop,nx_prop,nx_w,Vp_w,&
		icount,unit_fields,folder_fields,&
		added_space,store_snap,nt_taper,nt_taper,PSDM_w)
	endif

	if(strategy_num.eq.1)deallocate(Fields)
	deallocate(Vp_w)

	do i = 1, nx_w
	   do j = 1, ny_prop
	      Cov_w(j,i) = 1.
	   enddo
	enddo


deallocate(BackSource,Data_target)

endif !! icount

	allocate(Cov_thread(ny_prop,nx_prop)) 
	allocate(PSDM_thread(ny_prop,nx_prop)) 
	Cov_thread=0;PSDM_thread=0;

	PSDM_thread(:,x_i:x_f)=PSDM_w(:,:)
	Cov_thread(:,x_i:x_f)=Cov_w(:,:)

	allocate(Cov_itimes(ny_prop,nx_prop)) 
	allocate(PSDM_itimes(ny_prop,nx_prop)) 
	Cov_itimes=0.;PSDM_itimes=0.

	call MPI_ALLREDUCE(Cov_thread,Cov_itimes,nx_prop*ny_prop,MPI_REAL,MPI_SUM,new_comm,ierr)
	call MPI_ALLREDUCE(PSDM_thread,PSDM_itimes,nx_prop*ny_prop,MPI_REAL,MPI_SUM,new_comm,ierr)

	PSDM_sum = PSDM_sum + PSDM_itimes
	Cov_sum  = Cov_sum  + Cov_itimes

	deallocate(Cov_w,PSDM_w)
	deallocate(Cov_thread,PSDM_thread)
	deallocate(Cov_itimes,PSDM_itimes)

	if(rank.eq.0) then

		file_name=trim(folder_output) // trim("PSDM_RTM")
                file_V_twt=trim(folder_output) // trim("Velocity_twt")
                file_P_twt=trim(folder_output) // trim("PSDM_RTM_twt")

		if(write_interval.gt.0.and.mod(itimes,write_interval) == 0) then !! escribir cada N loops

                        call normalize_psdm(nx_prop,ny_prop,dmodel_prop,dmodel_prop,PSDM_sum,Cov_sum,PSDM_depth)        
			call write_psdm(nx_prop,ny_prop,dmodel_prop,dmodel_prop,file_name,PSDM_depth,200,0)

		        write(*,*)
			write(*,*)"------------------------------------------------------------------------"
			write(*,*)"     PARTIAL PSDM UNTIL ROUND ",itimes," WRITTEN IN FOLDER ", trim(folder_output)
			write(*,*)"------------------------------------------------------------------------"
			write(*,*)

		endif !write interval

		if(itimes.eq.ntimes) then !! escritura final

		        write(*,*)
			write(*,*)"------------------------------------------------------------------------"
			write(*,*)"     FINAL PSDM WRITTEN IN FOLDER ", trim(folder_output)
			write(*,*)"------------------------------------------------------------------------"
			write(*,*)

		        call normalize_psdm(nx_prop,ny_prop,dmodel_prop,dmodel_prop,PSDM_sum,Cov_sum,PSDM_depth)
			call write_psdm(nx_prop,ny_prop,dmodel_prop,dmodel_prop,file_name,PSDM_depth,200,1)

			!! convierte PSDM a TWT
                        if(run_psdm_twt.ne.0)call psdm2twt(nx_prop,ny_prop,dmodel_prop,dt_data,&
			file_V_twt,file_P_twt,Vpmodel,PSDM_depth)

		endif ! ntimes

	endif ! rank

enddo   !itimes 

!!write(*,*) "LLEGA AL FINAL, rank: ",rank

deallocate(nt_cut)
deallocate(PSDM_sum,Cov_sum)
deallocate(nxRec,nyRec)
deallocate(Data_synth,Data_source)
deallocate(Field_Data)
deallocate(Raw_Data,Delay_Data)

end subroutine calculate_PSDM_RTM
