subroutine calculate_PSDM_kirchoff(Vpmodel,Batmodel,new_comm,PSDM_depth)

use mod_parfile
use mod_data_arrays
use mpi

implicit none

real, parameter :: tmin        = 0.02      ! s  (evita singularidades)
real, parameter :: eps         = 1.0e-6
real, parameter :: pi          = 3.14159265

integer :: numtasks,rank, ierr, TAG
integer :: new_comm
integer :: itimes,icount,ntimes
integer :: nsamples,last,trid
integer :: i,j,k,kk,kkk,l
integer :: NumRec
integer :: nt_data
integer :: nt_taper,N_zeros,N_taper,nt_prop2

integer :: SourceNum
integer :: minS,maxS,minR,maxR
integer :: x_i,x_f,y_i,y_f,nx_w,nx_w_max
integer :: it

integer :: mid_x
real    :: offset,D_src,D_rec,leg
real    :: twt_sf,twt_multiple
real :: weight,rr,rs,amp

real :: tt_tot,twt_bat,Veff
real :: ddx,f1,f2
real :: tfin,dt_data,max_data
real :: Misfit,Misfit_,Misfit_thread

real :: x_rec,x_img,z_img,z_water
real :: theta,w_theta,cos_theta,pxs,pzs,pxr,pzr

real :: PSDM_depth(ny_prop,nx_prop)
real :: Vpmodel(ny_prop,nx_prop)
real :: Batmodel(nx_prop)

integer :: nxSou(1),nySou(1)

integer, allocatable :: nxRec(:),nyRec(:)
integer, allocatable :: nt_cut(:)
integer, allocatable :: nt_sf(:)

real, allocatable :: Cov(:,:),Cov_sum(:,:)
real, allocatable :: Cov_w(:,:),Cov_thread(:,:),Cov_itimes(:,:)

real, allocatable :: PSDM_sum(:,:)
real, allocatable :: PSDM_w(:,:),PSDM_thread(:,:),PSDM_itimes(:,:)

real, allocatable :: Data_target(:,:)
real, allocatable :: Field_Data(:,:)
real, allocatable :: Raw_Data(:,:)
real, allocatable :: Delay_Data(:,:)
real, allocatable :: Vp_w(:,:)
real, allocatable :: TT_s(:,:),TT_r(:,:,:)
real, allocatable :: Vprom(:)

character(len=1000) :: file_name,str_snap,file_V_twt,file_P_twt

call MPI_COMM_SIZE(new_comm,numtasks,ierr)
call MPI_COMM_RANK(new_comm,rank,ierr)
TAG=0


allocate(Vprom(ny_prop))


nt_prop2 = nt_prop
nt_taper = nt_prop2

!if(taper_option.ne.0)	then
!	N_taper = max(1, int(T_taper / dt_prop))
!	N_zeros = max(0, int(T_zeros / dt_prop))
!	nt_taper = nt_prop + N_zeros
!endif

!write(*,*) "dt_data,dt_prop: ",dt_mcs,dt_prop
!write(*,*) "nt_data,nt_prop,nt_taper: ",nt_prop,nt_taper

allocate(PSDM_sum(ny_prop,nx_prop)) 
allocate(Cov_sum(ny_prop,nx_prop)) 

Cov_sum=0;PSDM_sum=0;

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


if(rank.eq.0)write(*,*) "Allocating kirchoff dimensions"

allocate(nt_cut(NumRec))
allocate(nt_sf(NumRec))
allocate(nxRec(NumRec),nyRec(NumRec))
allocate(Field_Data(nt_data,NumRec))
allocate(Raw_Data(nt_data,NumRec))
allocate(Delay_Data(nt_data,NumRec))
allocate(Data_target(nt_taper,NumRec))

nxRec=0;nyRec=0;
nxSou=-10;nySou=-10;
nx_w=1;x_i=1;x_f=1;

Field_Data=0
Delay_Data=0
Raw_Data=0
Data_target=0

if(cut_multiple.eq.0)nt_cut=nt_prop;
if(cut_sf.eq.0)nt_sf=1;

ntimes=nsamples
if(numtasks.gt.1)    then
        if(nsamples.gt.numtasks)        then
                ntimes=ceiling(1.*nsamples/numtasks)
        endif
        if(nsamples.le.numtasks)        then
               ntimes=1
        endif
endif

!if(rank.eq.0)write(*,*)"Init and final lines ... : ",line_i," to ",line_f

l=0

do itimes=1,ntimes

icount=(itimes-1)*numtasks+rank+line_i        !!aquí sucede la paralelización

!write(*,*)"icount: ",icount

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

	    call open_shot_mcs(rank,icount)       !!opens unit file for MCS data
	    call get_shot_data(icount,Raw_Data,Delay_Data)  !! raw_data SU-MCS files		
	    call close_mcs_files(icount)

            if(delay_time_data.ne.0) then
                Field_Data=Delay_Data
            else
                Field_Data=Raw_Data
            endif

	endif

if(cut_multiple.gt.0)   then

!       !! calculo nt maximo para cada traza, cortar antes del multiple usando la info de Batmodel

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
                nt_cut(j) = ceiling( 0.9*twt_multiple/dt_prop)      !!lo corto un poco antes por si acaso

        enddo

endif


if(cut_sf.gt.0) then

!       !! calculo nt maximo para cada traza, cortar antes del multiple usando la info de Batmodel

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
                nt_sf(j) = ceiling( 1.1*twt_sf/dt_prop)       !!lo corto un poco antes por si acaso

        enddo

endif

!! suavizar el corte brusco a cero para cada receptor (i)
!! taper and padding, time limits per trace

if(taper_option.ne.0.or.cut_multiple.ne.0)      then
        N_taper = max(1, int(T_taper / dt_prop))
        N_zeros = max(0, int(T_zeros / dt_prop))
        nt_prop2= min(nt_prop, maxval(nt_cut))
        nt_taper = nt_prop2 + N_zeros
        !write(77,*)i,nt_prop,nt_prop2,n_zeros,nt_taper
endif

!!!!!-----------------------------------
!!!!!------ Frequency filter of data
!!!!!-----------------------------------

if(f_init.gt.0.and.f_final.gt.0) then

    if(rank.eq.0)       then
        write(*,*) "Filtering Field Data"
        if(typef.eq.1)write(*,*) 'Lowpass at FREQ = ',f_init
        if(typef.eq.3)write(*,*) 'Bandpass between FREQ = ',f_init,' and ',f_final
    endif

    call freq_filter(Field_Data,nt_data,dt_data,NumRec,typef,f_init,f_final)

endif

	if(dt_data.ne.dt_prop)	then
	        call interpolation_akima(Field_Data,dt_data,nt_data,Data_target(1:nt_prop2,:),NumRec,dt_prop,nt_prop2)
	else
		Data_target(1:nt_prop2,:)=Field_Data(1:nt_prop2,:)
	endif

!        open(unit=12,file="Target.dat",status='unknown')
!        do j=1,nt_prop2
!	        write(12,'(20000(e12.5,2x))') (Data_target(j,i),i=1,NumRec)
!        enddo

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

!	if(rank.eq.0)	then
!        open(unit=12,file="Target_window.dat",status='unknown')
!        do j=1,nt_taper
!               write(12,'(20000(e12.5,2x))') (Data_target(j,i),i=1,NumRec)
!        enddo
!	endif

endif	!!icount

if(nxSou(1).ge.0) then
	call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
	dmodel_prop,nx_prop,added_space,nx_w,x_i,x_f)
endif

allocate(PSDM_w(ny_prop,nx_w)) 
allocate(Cov_w(ny_prop,nx_w)) 
PSDM_w=0;Cov_w=0;

allocate(TT_s(ny_prop,nx_w))
allocate(TT_r(ny_prop,nx_w,NumRec))
allocate(Vp_w(ny_prop,nx_w))

if(icount.ge.line_i.and.icount.le.line_f)  then

	if(rank.eq.0)write(*,*) "Dynamic window in X (ptos) (nx_w,x_i,x_f): ",nx_w,x_i,x_f

        Vp_w(:,1:nx_w)=Vpmodel(:,x_i:x_f)

        !!!!!---------------------          
        !!!!!------ Eikonal equation
        !!!!!---------------------

        if(rank.eq.0)write(*,*)
        if(rank.eq.0)write(*,*)"Solve equation ..."
        !if(rank.eq.0)write(*,*)"values: ",NumRec,nxSou(1),nySou(1),nxRec(1),nyRec(1),&
        !dmodel_prop,ny_prop,nx_prop,nx_w,Vp_w(1,1),added_space

	call dynamic_window_tt_forward(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        dmodel_prop,ny_prop,nx_prop,nx_w,Vp_w,added_space,TT_s)

	call dynamic_window_tt_backward(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        dmodel_prop,ny_prop,nx_prop,nx_w,Vp_w,added_space,TT_r)

!        open(unit=12,file="TT_s.dat",status='unknown')
!        do j=1,ny_prop
!	        write(12,'(20000(e12.5,2x))') (TT_s(j,i),i=1,nx_w)
!        enddo

!!	Imaging condition:
!	do j = 1, ny_prop
!	   Vprom(j) = 0.0
!	   do i = 1, nx_w
!	      Vprom(j) = Vprom(j) + Vpmodel(j,i)
!	   end do
!	   Vprom(j) = Vprom(j) / nx_w
!	   Vprom(j) = min(max(Vprom(j), 1500.0), 5000.0)
!	end do

	PSDM_w=0.	!!ojo
	do k=1,NumRec
           x_rec = nxRec(k) * dmodel_prop

	   do i=1,nx_w
	     z_water = Batmodel(i)
	     x_img = (i-1) * dmodel_prop

	      do j=1,ny_prop
	        z_img = (j-1) * dmodel_prop

        	! Tiempo total
        	tt_tot = TT_s(j,i) + TT_r(j,i,k)
        	it = 1 + floor(tt_tot / dt_prop)
        	if (it <= 1 .or. it >= nt_prop) cycle

!!		call v_eff(z_img,z_water,Veff)
		Veff=vp0

!!		if(z_img.le.z_water)Veff=1500.
!!		if(z_img.gt.z_water)Veff=Vprom(j)

		rs = Veff * max(TT_s(j,i), tmin)
		rr = Veff * max(TT_r(j,i,k), tmin)

		!!if (z_img.lt.z_shallow) then
		!!	theta_max = 80.0; theta_taper = 15.0; wmax = 1.0
		!!else if (z_img.ge.z_shallow.and.z_img.lt.z_deep) then
		!!	theta_max = 60.0; theta_taper = 10.0; wmax = 0.7		
		!!else
		!!	theta_max = 45.0; theta_taper = 8.0;  wmax = 0.4
		!!endif		

		! Rays and angles
		cos_theta = z_img / (rs + rr + 1.e-6)
		cos_theta = max(min(cos_theta,1.0),0.0)	!!new
		theta = acos(cos_theta) * 180.0 / pi

		if (theta <= theta_max) then
			w_theta = 1.0
		else
			w_theta = exp(-((theta-theta_max)/theta_taper)**2)
		endif

		! Kirchhoff weight
		weight = w_theta * min( cos_theta / sqrt(rs*rr + eps), wmax ) !new
		!!weight = w_theta *  cos_theta / sqrt(rs*rr + eps)

		 !! Imaging condition
		 amp = (Data_target(it+1,k) - Data_target(it-1,k)) / (2.0*dt_prop)
	
	         PSDM_w(j,i) = PSDM_w(j,i) + weight * amp

	      end do
	   end do
	end do


	do i = 1, nx_w
	   do j = 1, ny_prop
	      Cov_w(j,i) = 1.
	   enddo
	enddo


endif !!line_i,line_f

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

	deallocate(TT_s,TT_r)
	deallocate(Vp_w)

        if(rank.eq.0) then

                file_name=trim(folder_output) // trim("PSDM_kirchoff")

                file_V_twt=trim(folder_output) // trim("Velocity_twt")
                file_P_twt=trim(folder_output) // trim("PSDM_kirchoff_twt")

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
deallocate(nt_sf)
deallocate(PSDM_sum,Cov_sum)
deallocate(nxRec,nyRec)
deallocate(Data_target)
deallocate(Field_Data)
deallocate(Raw_Data)
deallocate(Delay_Data)

end subroutine calculate_PSDM_kirchoff
