program main_RTM

use mod_parfile
use mod_data_arrays

implicit none
include 'mpif.h'

integer :: i,j,k,kkk
integer :: numtasks,rank,ierr
integer :: new_comm,new_rank,color

real, allocatable :: PSDM_depth(:,:)
real, allocatable :: Vpmodel(:,:)
real, allocatable :: Batmodel(:)
real, allocatable :: Batmodel_shots(:)

character(len=500) :: command
character(len=500) :: file_name
character(len=500) :: file_vel
character(len=500) :: file_psdm

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!-------------------
!!!!! Getting input data parameters
!!!!!-------------------

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*******************************'
        write(*,*)'READING PARAMETERS FROM PARFILE'
        write(*,*)'*******************************'
endif

call read_parfile(rank,numtasks)

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*** PARAMETER CHECK OK ***'
endif

call allocate_data_arrays()

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*******************************'
        write(*,*)'READING AND CHECKING INPUT FILES '
        write(*,*)'*******************************'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call allocate_model_arrays(rank)    !! From was and/or mcs

!! FOR SHOTS-MCS NAV

if(rtm_WAS.ne.0) then
	call allocate_model_arrays_was(rank)
	call check_nav_shot_was(rank)
	call read_nav_shot_was(rank)
	call check_nav_obs(rank)
	call read_nav_obs(rank)
	if(nav_was_cols.eq.1111) then
		call get_bathymetry_model_was(rank)
		if(rank.eq.0.and.SUM(bat_model_was).eq.0)write(*,*)" *** WARNING: Bathymetry not given in WAS"
	endif
endif

if(rtm_MCS.ne.0) then
	call allocate_model_arrays_mcs(rank)
	call check_nav_shot_mcs(rank)
	call read_nav_shot_mcs(rank)
	if(nav_mcs_cols.eq.1111) then
		call get_bathymetry_model_mcs(rank)
		if(rank.eq.0.and.SUM(bat_model_mcs).eq.0)write(*,*)" *** WARNING: Bathymetry not given in MCS"
	endif
endif

if(rank.eq.0)   then
        write(*,*)
        write(*,*)'*********************************'
        write(*,*)'GET INITIAL MODEL '
        write(*,*)'*********************************'
endif

    call get_vel(rank)
    call optimization(rank)

    if(rank.eq.0)write(*,*)"params out: nx_prop: ",nx_prop
    if(rank.eq.0)write(*,*)"params out: ny_prop: ",ny_prop
    if(rank.eq.0)write(*,*)"params out: dmodel_prop: ",dmodel_prop
    if(rank.eq.0)write(*,*)"params out: dt_mcs: ",dt_mcs
    if(rank.eq.0)write(*,*)"params out: nt_data: ",nt_mcs
    if(rank.eq.0)write(*,*)"params out: dt_prop: ",dt_prop
    if(rank.eq.0)write(*,*)"params out: nt_prop: ",nt_prop

!    if(rank.eq.0)write(*,*)"params out: n_cores_real: ",n_cores_real
!    if(rank.eq.0)write(*,*)"params out: strategy_num: ",strategy_num
!    if(rank.eq.0)write(*,*)"params out: strategy: ",trim(adjustl(strategy))
!    if(rank.eq.0)write(*,*)

	allocate(PSDM_depth(ny_prop,nx_prop))
        allocate(Vpmodel(ny_prop,nx_prop))
	allocate(Batmodel(nx_prop))

	PSDM_depth=0.
	Vpmodel=0.
	Batmodel=0.

	if(strategy_num.eq.2.and.rank.eq.0)	then
		command="mkdir " // trim(adjustl(folder_fields))
		call system(command)
		write(*,*) trim(adjustl(command))
	endif
	
	if (rank < n_cores_real) then
		color = 1
	else
		color = MPI_UNDEFINED
	endif
	
	call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, rank, new_comm, ierr)
	
	if (color /= MPI_UNDEFINED) then
	 
		!!! Interpolation velocity model
	
		if(dmodel.eq.dmodel_prop)	then
			Vpmodel=vel	
		else
			call linear_2D(nxmodel,nymodel,dmodel,vel,nx_prop,ny_prop,dmodel_prop,Vpmodel)
		endif
	
		if(rank.eq.0)   then  !ojo
			write(*,*)"Writting Vel model interpolated to a grid of ",dmodel_prop," meters"
			write(*,*)"with nx, ny: ",nx_prop,ny_prop	!!ojo
			file_name=trim(folder_output) // trim("Velocity")
			call write_vel(nx_prop,ny_prop,dmodel_prop,dmodel_prop,200,file_name,Vpmodel)        
		endif
	
		!!! Interpolation of bathymetry
	
		do i=1,nx_prop
		kkk=0
		do j=1,ny_prop

			if(Vpmodel(j,i).gt.water_velocity.and.kkk.eq.0) then
	
				Batmodel(i)=(j-1-1)*dmodel_prop		
				kkk=1
			endif

		enddo
		enddo		

		if(nav_mcs_cols.eq.1111)	then

			allocate(Batmodel_shots(nx_prop))
			Batmodel_shots=0.

			if(dmodel.eq.dmodel_prop)	then
				if(RTM_was.ne.0)Batmodel_shots=bat_model_was
				if(RTM_mcs.ne.0)Batmodel_shots=bat_model_mcs
				!write(*,*)"ENTRA",bat_model_mcs(NumShots_mcs)
			else
				if(RTM_was.ne.0)call linear_1D(nxmodel,dmodel,bat_model_was,nx_prop,dmodel_prop,Batmodel_shots)
				if(RTM_mcs.ne.0)call linear_1D(nxmodel,dmodel,bat_model_mcs,nx_prop,dmodel_prop,Batmodel_shots)
			endif

			!do i=1,nx_prop
			!	write(1111,*)i,Batmodel(i),Batmodel_shots(i)
			!enddo

			deallocate(Batmodel_shots)

		endif
 

		if(rank.eq.0)   then !ojo
			write(*,*)"Writting Bat model interpolated to a grid of ",dmodel_prop ," meters"
			file_name=trim(folder_output) // trim("Bat_model_meters.txt")
		        open(12,file=file_name,status='unknown')
		        do i=1,nx_prop
			        write(12,*)(i-1)*dmodel_prop,Batmodel(i)
			enddo
			close(12)

		endif

		if(RTM_was.ne.0)then

			do i=1,NumOBS

				call open_obs(rank,i)     !!opens unit file for WAS data
				call get_geom_was_data(rank,i)
				call close_was_files(i)
			enddo

		endif

!		if(RTM_mcs.eq.11110)then	!ojo
		if(RTM_mcs.ne.0)then	
	
			file_name=trim(folder_output) // trim("Shot_position_model.txt")
			open(unit=12,file=file_name,status='unknown')
	
			do i=1,NumShots_mcs
		

				if(rank.eq.0.and.i.eq.NumShots_mcs/4)write(*,*)"25% geom"
				if(rank.eq.0.and.i.eq.NumShots_mcs/2)write(*,*)"50% geom"
				if(rank.eq.0.and.i.eq.3*NumShots_mcs/4)write(*,*)"75% geom"

				call open_shot_mcs(rank,i)       !!opens unit file for MCS data
				call get_geom_mcs_data(rank,i)
				call close_mcs_files(i)
	
				write(12,*)i,shotID_nav_mcs(i),nxSou_mcs(i),pos_shot_mcs(i)

				if(rank.eq.0) then !ojo

				if ( min( nxSou_MCS(i), minval(nxRec_MCS(:,i)) ) < 1 ) then
					if(rank.eq.0)write(*,*)"Warning geometry is out of lower boundary from Model at shot: ",shotID_nav_MCS(i)
					if(rank.eq.0)write(*,*)"nxSou min,nxRec min ",nxSou_MCS(i),minval(nxRec_MCS(:,i))
					if(rank.eq.0)write(*,*)
				endif
	
				if ( max( nxSou_MCS(i), maxval(nxRec_MCS(:,i)) ) > nx_prop ) then
					if(rank.eq.0)write(*,*)"Warning geometry is out of upper boundary from Model at shot: ",shotID_nav_MCS(i)
					if(rank.eq.0)write(*,*)"nxSou max,nxRec max",nxSou_MCS(i),maxval(nxRec_MCS(:,i))
					if(rank.eq.0)write(*,*)
				endif

				endif
			enddo

			close(12)

		endif	!mcs

		if(run_psdm.ne.0)	then

			!!if(psdm_ttwave.ne.0) call calculate_PSDM_hibrid(Vpmodel,batmodel,new_comm,PSDM_depth)

			if(psdm_rtm.ne.0) call calculate_PSDM_RTM(Vpmodel,Batmodel,new_comm,PSDM_depth)
			if(psdm_eikonal.ne.0) call calculate_PSDM_eikonal(Vpmodel,Batmodel,new_comm,PSDM_depth)
			if(psdm_kirchoff.ne.0) call calculate_PSDM_kirchoff(Vpmodel,Batmodel,new_comm,PSDM_depth)

		endif

		if(rank.eq.0)	then

		if(run_psdm.eq.0.and.run_psdm_twt.ne.0)	then

			!! PSDM and Vel to TWT

			call get_psdm(nx_prop,ny_prop,PSDM_depth)

			file_vel=trim(folder_output) // trim("Velocity_twt")
			file_psdm=trim(folder_output) // trim("PSDM_twt")

			if(RTM_mcs.gt.0) call psdm2twt(nx_prop,ny_prop,dmodel_prop,dt_mcs,file_vel,file_psdm,Vpmodel,PSDM_depth)
			if(RTM_was.gt.0) call psdm2twt(nx_prop,ny_prop,dmodel_prop,dt_was,file_vel,file_psdm,Vpmodel,PSDM_depth)

		endif

		endif

	endif !if color

	deallocate(PSDM_depth)
	deallocate(Vpmodel)
	deallocate(Batmodel)

	call deallocate_arrays()

!    if(strategy_num.eq.2.and.rank.eq.0)	then
!	    command="rm -r " // trim(adjustl(folder_fields))
!	    call system(command)
!	    write(*,*) trim(adjustl(command))
!    endif

	write(*,*)'end rank',rank
	call MPI_FINALIZE(ierr)

end program main_RTM
