!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains the modules:
!!	(mod_parfile, mod_data_arrays)
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_parfile

!! INFORMATION ABOUT: WAS PARAMETERS, MCS PARAMETERS, MODEL PARAMETERS, FWI PARAMETERS

implicit none

  integer, parameter :: unit_was=300,unit_mcs=400, unit_fields=500
!  integer, parameter :: unit_FWI=50
 
  INTEGER, parameter :: size_su_header = 60
  INTEGER(4), parameter :: byte_fldr = 9
  INTEGER(4), parameter :: byte_tracl = 1
  INTEGER(4), parameter :: byte_tracr = 5
  INTEGER(4), parameter :: byte_offset = 37
  INTEGER(4), parameter :: byte_scalco = 71
  INTEGER(4), parameter :: byte_sx = 73
  INTEGER(4), parameter :: byte_sy = 77
  
  INTEGER(4) :: byte_shotnumber_was
  INTEGER(4) :: byte_shotnumber_mcs

  integer :: run_psdm,run_psdm_twt
  integer :: psdm_rtm,psdm_eikonal,psdm_kirchoff,psdm_ttwave
  integer :: automatic_optimization,write_interval,taper_option
  integer :: mute_forward
  integer :: post_filt
  integer :: fs,dPML
  integer :: bat_clean
  integer :: dir_boat
  integer :: shot_i,shot_f
  integer :: obs_i,obs_f
  integer :: rtm_was,rtm_mcs
  integer :: NumShots_was,NumShots_mcs
  integer :: NumOBS,NumChannels	!!fuentes reales y OBS reales
  integer :: NumSou_MCS,NumSou_WAS,NumRec_MCS,NumRec_WAS
  integer :: seed_option
  integer :: niter
  integer :: nPML
  integer :: cut_multiple
  integer :: cut_sf
  integer :: line_i,line_f
  !! inputs given by user:
  real :: delay_ricker,delay_time_data
  real :: T_zeros, T_taper
  real :: RAM_total_GB, Disk_total_GB	!!input computer
  integer :: n_cores_total,n_cores_real
  integer :: nx_prop,ny_prop	!!input model
  integer :: nt_prop 	!! input data	
  integer :: nt_snap,store_snap,strategy_num	!!outputs for modelling
  real    :: v_min, v_max, x_max, y_max 	!!input from model
  real    :: t_max			!!input from data
  real    :: dt_prop  !!outputs for good modelling
  integer :: nxmodel,nymodel,nxmodel0,nymodel0
  integer :: nxmodel_was,nxmodel_mcs

  integer :: nt_was,nt_was_max,nt_mcs,nt_mcs_max
  integer :: read_vel    !!ojo, un modelo de velocidad para was, y otro para mcs? o uno comun?

  integer :: endianness_data,endianness_machine
  integer :: save_txt,save_xzv,save_v,save_su
  integer :: vel_v,vel_xzv,vel_su
  integer :: psdm_v,psdm_xzv,psdm_su
  integer :: offset_header
  integer :: maxbytes

  integer :: BP_type
  integer :: method,typef
  integer :: nav_was_cols,nav_mcs_cols

  real :: vp0
  real :: xmodel_min,xmodel_max
  real :: ymodel_max
  real :: theta_max,theta_taper,wmax,z_shallow,z_deep
  real :: dvel
  real :: offset_unit
  real :: drec
  real :: dshots_was,dshots_mcs
  real :: near_offset,dmodel,dmodel_prop
  real :: dt_was,dt_mcs
  real :: shot_depth_was,shot_depth_mcs
  real :: streamer_depth
  real :: streamer_length
  real :: added_space,water_velocity
  real :: d_shot_model
  real :: tfin_was,tfin_mcs
  real :: water_depth
  real :: freq_ricker
  real :: f_init,f_final
  real :: vel_i,vel_f
  real :: length_model  !!ojo
  real :: length_model_was  !!ojo
  real :: length_model_mcs  !!ojo

  character(len=500) :: folder_input_was,folder_input_mcs,folder_input_model,psdm_file
  character(len=500) :: folder_output,folder_fields
  character(len=500) :: folder_VEL,folder_GRAD,folder_DATA,folder_SOU
  character(len=500) :: su_file_FWI
  character(len=500) :: was_data,mcs_data
  character(len=500) :: nav_obs !! navegacion (posiciones OBS-WAS)
  character(len=500) :: nav_shot_was,nav_shot_mcs
  character(len=500) :: vel_file
  character(len=500) :: par_file

  character(len=50)  :: strategy

  character(len=100), allocatable :: su_file_was(:)
  character(len=100), allocatable :: su_file_mcs(:)

  contains

  subroutine read_parfile(rank,numtasks)

  implicit none
  include 'mpif.h'

  integer :: numtasks,rank,ierr,errcode,status(MPI_STATUS_SIZE)

  ! Input related variables
  character(len=200) :: buffer,label
  integer :: ps,icount,i,ifile,interval,nlines
  integer, parameter :: fh = 10
  integer :: ios = 0
  integer :: line = 0

  character(len=500) :: command0,command
  character(len=50) :: Str,access,form,num_split
  character(len=500) :: file_name,file_name2

  logical :: su_exist,nav_exist,vel_exist,data_exist
  logical :: input_exist, output_exist

  call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  ierr=0;

  access = 'STREAM'
  form = 'UNFORMATTED'

  icount = iargc()  ! The number does not includes the executable name, so if user passed one argument, 1 is returned.
  if ( icount.eq.1 ) then
	call getarg(1, par_file)	! The file name of the executable.
	if(rank.eq.0)write(*,*)'name par_file: ',trim(adjustl(par_file))
	file_name = trim(adjustl(par_file))
	INQUIRE(file=file_name,EXIST=su_exist)
	if(.NOT. su_exist)     then
  		if(rank.eq.0)write(*,*)'ERROR: Par_file named: ', trim(adjustl(par_file)),' does not exist'
        	if(rank.eq.0)call ascii_art(2)
		call MPI_barrier(MPI_COMM_WORLD,ierr)
		call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
		stop
  	endif
  endif

  run_psdm=1
  run_psdm_twt=1
  n_cores_total=numtasks
  strategy_num=1
  su_file_FWI = 'null'

  nav_shot_was = 'null'
  nav_shot_mcs = 'null'
  nav_obs = 'null'
  vel_file = 'null'

  nav_mcs_cols=11
  nav_was_cols=11
  psdm_ttwave=0
  psdm_rtm=0
  psdm_eikonal=0
  psdm_kirchoff=0

  psdm_su=0;psdm_v=0;psdm_xzv=0;

  byte_shotnumber_was= byte_tracr
  byte_shotnumber_mcs= byte_fldr
  endianness_data=1;endianness_machine=0;
  offset_header=0;offset_unit=1;

  cut_multiple=0;
  cut_sf=0;
  taper_option=0;
  T_zeros=0.5		! second units

  d_shot_model=0
  dPML=20

  vp0=1500.

  !! ojo, automatizar para los tres marcos de profundidad
  theta_max=60
  theta_taper=10
  wmax=1.
  z_shallow=3000.;
  z_deep=8000.;

  dvel=0.;
  delay_ricker=0.;
  delay_time_data=0.;
  fs=1;
  bat_clean=0;
  post_filt=0;
  write_interval=0;
  shot_i=0;shot_f=0;
  obs_i=0;obs_f=0;
  automatic_optimization=0
  store_snap=1
  rtm_WAS=0;rtm_MCS=0;
  seed_option=0;
  added_space=500!!en metros

  water_velocity=1500.
  dmodel=0;
  dmodel_prop=0;
  shot_depth_was=0;
  shot_depth_mcs=0;
  streamer_depth=0;
  streamer_length=0;
  f_init=-10;f_final=-10;
  save_txt=0;save_v=0;save_xzv=0;save_su=1;
  vel_v=0;vel_xzv=0;vel_su=0;

  method=1
  typef=3  
  BP_type=1

  tfin_was=-100.;
  tfin_mcs=-100.;
  nt_was=-100;dt_was=-100;
  nt_mcs=-100;dt_mcs=-100;
  nt_mcs_max=-100;  
  nt_was_max=-100;  
  NumOBS=0;NumChannels=0;NumShots_was=0;
  NumShots_mcs=0;
  NumSou_WAS=0;NumSou_MCS=0;
  NumRec_WAS=0;NumRec_MCS=0;
  read_vel=0
  vel_i=1500.;vel_f=4000;
  nxmodel=1
  nymodel=1
  nxmodel0=1
  nymodel0=1
  near_offset=0  
  dir_boat=0
  open(fh, file=par_file)

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected. It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)

     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

	! Verificar si la línea es un comentario
        if (buffer(1:1) .ne. '#') then

        ! Find the first instance of whitespace.  Split label and data.
        ps = scan(buffer,' ')
        label = buffer(1:ps)
        buffer = buffer(ps+1:)

        select case (label)

        case ('run_psdm:')
           read(buffer, *, iostat=ios) run_psdm
        case ('run_psdm_twt:')
           read(buffer, *, iostat=ios) run_psdm_twt
        case ('psdm_rtm:')
           read(buffer, *, iostat=ios) psdm_rtm
        case ('psdm_eikonal:')
           read(buffer, *, iostat=ios) psdm_eikonal
        case ('psdm_kirchoff:')
           read(buffer, *, iostat=ios) psdm_kirchoff
        case ('psdm_ttwave:')
           read(buffer, *, iostat=ios) psdm_ttwave
        case ('post_filt:')
           read(buffer, *, iostat=ios) post_filt
        case ('mute_forward:')
           read(buffer, *, iostat=ios) mute_forward
        case ('vp0:')
           read(buffer, *, iostat=ios) vp0
        case ('theta_max:')
           read(buffer, *, iostat=ios) theta_max
        case ('theta_taper:')
           read(buffer, *, iostat=ios) theta_taper
        case ('z_shallow:')
           read(buffer, *, iostat=ios) z_shallow
        case ('z_deep:')
           read(buffer, *, iostat=ios) z_deep
        case ('dvel:')
           read(buffer, *, iostat=ios) dvel
        case ('delay_time_data:')
           read(buffer, *, iostat=ios) delay_time_data
        case ('delay_ricker:')
           read(buffer, *, iostat=ios) delay_ricker
        case ('dPML:')
           read(buffer, *, iostat=ios) dPML
        case ('free_surface:')
           read(buffer, *, iostat=ios) fs
        case ('T_taper:')
           read(buffer, *, iostat=ios) T_taper
        case ('T_zeros:')
           read(buffer, *, iostat=ios) T_zeros
        case ('cut_multiple:')
           read(buffer, *, iostat=ios) cut_multiple
        case ('cut_sf:')
           read(buffer, *, iostat=ios) cut_sf
        case ('bat_clean:')
           read(buffer, *, iostat=ios) bat_clean
        case ('fi:')
           read(buffer, *, iostat=ios) f_init
        case ('ff:')
           read(buffer, *, iostat=ios) f_final
        case ('freq_ricker:')
           read(buffer, *, iostat=ios) freq_ricker
        case ('endianness_data:')
           read(buffer, *, iostat=ios) endianness_data
        case ('endianness_machine:')
           read(buffer, *, iostat=ios) endianness_machine
        case ('psdm_xzv:')
           read(buffer, *, iostat=ios) psdm_xzv
        case ('psdm_v:')
           read(buffer, *, iostat=ios) psdm_v
        case ('psdm_su:')
           read(buffer, *, iostat=ios) psdm_su
        case ('save_su:')
           read(buffer, *, iostat=ios) save_su
        case ('save_v:')
           read(buffer, *, iostat=ios) save_v
        case ('save_xzv:')
           read(buffer, *, iostat=ios) save_xzv
        case ('method:')
           read(buffer, *, iostat=ios) method
        case ('typef:')
           read(buffer, *, iostat=ios) typef
        case ('BP_type:')
           read(buffer, *, iostat=ios) BP_type
        case ('vel_file:')
           read(buffer, *, iostat=ios) vel_file
        case ('vel_su:')
           read(buffer, *, iostat=ios) vel_su
        case ('vel_v:')
           read(buffer, *, iostat=ios) vel_v
        case ('vel_xzv:')
           read(buffer, *, iostat=ios) vel_xzv
        case ('vel_i:')
           read(buffer, *, iostat=ios) vel_i
        case ('vel_f:')
           read(buffer, *, iostat=ios) vel_f
        case ('psdm_file:')
           read(buffer, *, iostat=ios) psdm_file
        case ('output_folder:')
           read(buffer, *, iostat=ios) folder_output
        case ('nxmodel:') 
	   read(buffer, *, iostat=ios) nxmodel
        case ('nymodel_max:')
           read(buffer, *, iostat=ios) nymodel
        case ('nymodel:')
           read(buffer, *, iostat=ios) nymodel0
        case ('xmodel_min:')
           read(buffer, *, iostat=ios) xmodel_min
        case ('xmodel_max:')
           read(buffer, *, iostat=ios) xmodel_max
        case ('ymodel_max:')
           read(buffer, *, iostat=ios) ymodel_max
        case ('dmodel:')
           read(buffer, *, iostat=ios) dmodel
        case ('water_velocity:')
           read(buffer, *, iostat=ios) water_velocity
        case ('seed_option:')
           read(buffer, *, iostat=ios) seed_option 
        case ('input_folder_model:')
	    read(buffer, *, iostat=ios) folder_input_model
        case ('rtm_WAS:')
           read(buffer, *, iostat=ios) rtm_WAS
        case ('rtm_MCS:')
           read(buffer, *, iostat=ios) rtm_MCS

        case ('RAM_total_GB:')
           read(buffer, *, iostat=ios) RAM_total_GB
        case ('Disk_total_GB:')
           read(buffer, *, iostat=ios) Disk_total_GB

        case ('automatic_optimization:')
           read(buffer, *, iostat=ios) automatic_optimization
        case ('dmodel_prop:')
           read(buffer, *, iostat=ios) dmodel_prop
        case ('dt_prop:')
           read(buffer, *, iostat=ios) dt_prop
        case ('store_snap:')
           read(buffer, *, iostat=ios) store_snap
        case ('added_space:')
           read(buffer, *, iostat=ios) added_space
        case ('d_shot_model:')
           read(buffer, *, iostat=ios) d_shot_model
        case ('strategy_num:')
           read(buffer, *, iostat=ios) strategy_num
        case ('write_interval:')
           read(buffer, *, iostat=ios) write_interval
        case ('dir_boat:')
           read(buffer, *, iostat=ios) dir_boat

!        case default
!           if(rank.eq.0)print *, 'Skipping invalid label at line', line

        end select

        end if  ! Fin del if para verificar comentarios
     end if
  end do

close(fh)


if(rtm_WAS.ne.0)	then

ios=0
open(fh, file=par_file)

do while (ios == 0)
read(fh, '(A)', iostat=ios) buffer
if (ios == 0) then
        line = line + 1
	! Verificar si la línea es un comentario
        if (buffer(1:1) .ne. '#') then
        ! Find the first instance of whitespace.  Split label and data.
        ps = scan(buffer,' ')
        label = buffer(1:ps)
        buffer = buffer(ps+1:)

        select case (label)

	        case ('byte_shotnumber_was:')
	           read(buffer, *, iostat=ios) byte_shotnumber_was
	        case ('was_file_list:')
	           read(buffer, *, iostat=ios) was_data
	        case ('nav_was_file:')
	           read(buffer, *, iostat=ios) nav_shot_was
	        case ('nav_obs:')
	           read(buffer, *, iostat=ios) nav_obs
	        case ('input_folder_was:')
	           read(buffer, *, iostat=ios) folder_input_was
	        case ('dt_was:')
		read(buffer, *, iostat=ios) dt_was
		case ('nt_was:')
	           read(buffer, *, iostat=ios) nt_was
		case ('nt_was_max:')
	           read(buffer, *, iostat=ios) nt_was_max
	        case ('NumShots_was:')
	           read(buffer, *, iostat=ios) NumShots_was
	        case ('dshots_was:')
	           read(buffer, *, iostat=ios) dshots_was
	        case ('NumOBS:')
	           read(buffer, *, iostat=ios) NumOBS
	        case ('shot_depth_was:')
	           read(buffer, *, iostat=ios) shot_depth_was
        	case ('obs_i:')
        	   read(buffer, *, iostat=ios) obs_i
        	case ('obs_f:')
        	   read(buffer, *, iostat=ios) obs_f

end select
end if  ! Fin del if para verificar comentarios
end if
end do
close(fh)

endif

if(rtm_MCS.ne.0)	then

ios=0
open(fh, file=par_file)

do while (ios == 0)

read(fh, '(A)', iostat=ios) buffer
if (ios == 0) then
        line = line + 1

	! Verificar si la línea es un comentario
        if (buffer(1:1) .ne. '#') then

        ! Find the first instance of whitespace.  Split label and data.
        ps = scan(buffer,' ')
        label = buffer(1:ps)
        buffer = buffer(ps+1:)

        select case (label)

	        case ('byte_shotnumber_mcs:')
	            read(buffer, *, iostat=ios) byte_shotnumber_mcs
	        case ('mcs_file_list:')
	           read(buffer, *, iostat=ios) mcs_data
	        case ('nav_mcs_file:')
	           read(buffer, *, iostat=ios) nav_shot_mcs
	        case ('nav_mcs_cols:')
	           read(buffer, *, iostat=ios) nav_mcs_cols
	        case ('offset_header_mcs:')
	            read(buffer, *, iostat=ios) offset_header
	        case ('offset_unit_mcs:')
		    read(buffer, *, iostat=ios) offset_unit
        	case ('input_folder_mcs:')
        	    read(buffer, *, iostat=ios) folder_input_mcs
        	case ('dt_mcs:')
        	   read(buffer, *, iostat=ios) dt_mcs
        	case ('nt_mcs:')
        	   read(buffer, *, iostat=ios) nt_mcs_max
        	case ('nt_mcs_cut:')
        	   read(buffer, *, iostat=ios) nt_mcs
        	case ('NumShots_mcs:')
        	   read(buffer, *, iostat=ios) NumShots_mcs
	        case ('dshots_mcs:')
	           read(buffer, *, iostat=ios) dshots_mcs
        	case ('shot_i:')
        	   read(buffer, *, iostat=ios) shot_i
        	case ('shot_f:')
        	   read(buffer, *, iostat=ios) shot_f
        	case ('near_offset:')
        	   read(buffer, *, iostat=ios) near_offset
        	case ('NumRec:')
        	   read(buffer, *, iostat=ios) NumChannels
        	case ('drec:')
        	   read(buffer, *, iostat=ios) drec
        	case ('shot_depth_mcs:')
        	    read(buffer, *, iostat=ios) shot_depth_mcs
        	case ('streamer_depth:')
        	   read(buffer, *, iostat=ios) streamer_depth

	case default

!if(rank.eq.0)print *, 'WARNING in file ',trim(adjustl(par_file)),': skipping invalid label at line', line

end select

end if  ! Fin del if para verificar comentarios
end if
end do
close(fh)

end if


if(T_taper.eq.0)T_taper=T_zeros/4.

if(rtm_WAS.ne.0)folder_input_was = trim(adjustl(folder_input_was)) // '/'
if(rtm_MCS.ne.0)folder_input_mcs = trim(adjustl(folder_input_mcs)) // '/'
folder_input_model = trim(adjustl(folder_input_model)) // '/'

folder_output = trim(adjustl(folder_output)) // '/'
folder_fields = trim(adjustl(folder_output)) // "FIELDS" // '/'

call MPI_barrier(MPI_COMM_WORLD,ierr)
if(rank.eq.0)	then
	command="mkdir " // trim(adjustl(folder_output))
	call system(command)
	command="cp " // trim(adjustl(par_file)) // " " // trim(adjustl(folder_output))
	call system(command)
endif
call MPI_barrier(MPI_COMM_WORLD,ierr)

if(rank.eq.0)call warnings_errors(numtasks)

RAM_total_GB=RAM_total_GB*0.9

if(nt_mcs.le.0)nt_mcs=nt_mcs_max
if(nt_was.le.0)nt_was=nt_was_max
if(dmodel.eq.0)dmodel=drec
if(dmodel_prop.eq.0)dmodel_prop=dmodel
length_model=(nxmodel-1)*dmodel

if(nymodel.eq.1)nymodel=nymodel0

NumRec_WAS=NumShots_was	!! IN WAS DATA: RECIPROCITY. SHOTS ACT AS RECEIVERS:
NumSou_WAS=NumOBS	!! sources in WAS are OBSs

NumSou_MCS=NumShots_mcs	!! source in MCS are shotgathers
NumRec_MCS=NumChannels	!! Receivers are channels

tfin_was=(nt_was-1)*dt_was
tfin_mcs=(nt_mcs-1)*dt_mcs

streamer_length=drec*(NumRec_MCS-1)

if(rtm_WAS.ne.0) then
    allocate(su_file_was(NumOBS))
    file_name= trim(adjustl(folder_input_was)) // trim(adjustl(was_data))
    open(unit=10,file=file_name,status='old')
    do i=1,NumOBS
        read(10,*)file_name2
        su_file_was(i) = trim(adjustl(file_name2))
    enddo
    close(10)
endif

if(rtm_MCS.ne.0) then
    allocate(su_file_mcs(NumShots_mcs))
    file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(mcs_data))
    open(unit=10,file=file_name,status='old')
    do i=1,NumShots_mcs
        read(10,*)file_name2
        su_file_mcs(i) = trim(adjustl(file_name2))
	!write(10+i,*) trim(adjustl(su_file_mcs(i)))
    enddo
    close(10)
endif

maxbytes=1900000000
	
if(save_v.ne.0.or.save_xzv.ne.0)save_txt=1

file_name = trim(adjustl(folder_input_model)) // trim(adjustl(vel_file))
INQUIRE(FILE=file_name, EXIST=vel_exist)
if(vel_exist)read_vel=1

if(rank.eq.0)	then

	write(*,*)
	write(*,*)'*******************************************'
	write(*,*)'RUNNING PROGRAM FOR NEXT SET OF PARAMETERS:'
	write(*,*)'*******************************************'
        write(*,*)'Output folder: ',adjustl(trim(folder_output))
	write(*,*)'dmodel (m): ',dmodel
	write(*,*)'nxmodel: ',nxmodel
	write(*,*)'nymodel: ',nymodel
	write(*,*)'method: ',method
	write(*,*)'typef: ',typef
	write(*,*)'BP_type: ',BP_type
	if(read_vel.eq.1)    then
        	write(*,*)'Input folder model: ',adjustl(trim(folder_input_model))
        	write(*,*)'Vp file: ',adjustl(trim(vel_file))
	endif
        if(rtm_WAS.ne.0)    then
        	write(*,*)'Inversion parameters for WAS: '
        	write(*,*)'Input WAS folder: ',adjustl(trim(folder_input_was))
        	write(*,*)'NumOBS: ',NumOBS
        	write(*,*)'NumShots WAS: ',NumShots_was
        	write(*,*)'name of navigation obs file: ',adjustl(trim(nav_obs))
        	write(*,*)'name of navigation shots WAS file: ',adjustl(trim(nav_shot_was))
        	write(*,*)'name of file with the name of the WAS-SU-data files: ',adjustl(trim(was_data))
	        write(*,*)'dt_was: ',dt_was
		write(*,*)'nt_was_cut: ',nt_was
		write(*,*)'nt_was: ',nt_was_max
        endif
        if(rtm_MCS.ne.0) then
        	write(*,*)'Inversion parameters for MCS: '
        	write(*,*)'Input MCS folder: ',adjustl(trim(folder_input_mcs))
        	write(*,*)'NumChannels: ',NumChannels
        	write(*,*)'NumShots MCS: ',NumShots_mcs
        	write(*,*)'**name of navigation shots MCS file: ',adjustl(trim(nav_shot_mcs))
        	write(*,*)'**name of file with the name of the MCS (shotgathers) SU-data files: ',adjustl(trim(mcs_data))
		write(*,*)'Offset header for mcs: ',offset_header
	        write(*,*)'dt_mcs: ',dt_mcs
		write(*,*)'nt_mcs_cut: ',nt_mcs
		write(*,*)'nt_mcs: ',nt_mcs_max
        endif

endif

end subroutine read_parfile

subroutine warnings_errors(numtasks)

implicit none
integer :: rank,numtasks

  ! Input related variables
  character(len=200) :: buffer,label
  integer :: ps,icount,i,ifile,interval,nlines
  integer, parameter :: fh = 10
  integer :: ios = 0
  integer :: line = 0
  integer :: error = 0
  character(len=500) :: command0,command
  character(len=50) :: Str,access,form,num_split
  character(len=500) :: file_name,file_name2

  logical :: su_exist,nav_exist,vel_exist,data_exist
  logical :: input_exist, output_exist

file_name = trim(adjustl(folder_input_model)) // trim(adjustl(vel_file))
INQUIRE(FILE=file_name, EXIST=vel_exist)
if(.NOT. vel_exist)	then
	write(*,*)"WARNING: vel_file not given or not found: using only water velocity"
	!error=1
endif

if(endianness_data.ne.0.and.endianness_data.ne.1)	then
	write(*,*)'ERROR: endianness_data should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
	error=1
endif

if(endianness_machine.ne.0.and.endianness_machine.ne.1)	then
	write(*,*)'ERROR: endianness_machine should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
	error=1
endif

if(rtm_was.ne.0)    then

    file_name = trim(adjustl(folder_input_was)) // trim(adjustl(was_data))
    INQUIRE(FILE=file_name, EXIST=data_exist)
    if(.NOT. data_exist)     then
        write(*,*)'ERROR: "file_data_list: " file not found'
        error=1
    else
        !!	OBS DATA FILE
        file_name= trim(adjustl(folder_input_was)) // trim(adjustl(was_data))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=10)
            nlines = nlines + 1
        enddo
        10 close (10)

        if(nlines.lt.NumOBS)    then
            write(*,*)'ERROR: ', abs(nlines-NumOBS),' OBS data files are missing in was_file_list'
            error=1
        endif
    endif

    file_name = trim(adjustl(folder_input_was)) // trim(adjustl(nav_shot_was))
    INQUIRE(FILE=file_name, EXIST=nav_exist)
    if(.NOT. nav_exist)	then
        write(*,*)'ERROR: nav_shot_was not found'
        error=1
    endif
    file_name = trim(adjustl(folder_input_was)) // trim(adjustl(nav_obs))
    INQUIRE(FILE=file_name, EXIST=nav_exist)
    if(.NOT. nav_exist)	then
        write(*,*)'ERROR: nav_obs not found'
        error=1
    endif

if(dt_was.eq.0) then
	write(*,*) 'ERROR: Please give a value to dt (seconds) in ',trim(adjustl(par_file))
	error=1
endif

if(nt_was.eq.0) then
	write(*,*)'ERROR: Please give a value to nt in ',trim(adjustl(par_file))
	error=1
endif

if(dir_boat.eq.0) then
	write(*,*)'ERROR: Please give a value to dir_boat in ',trim(adjustl(par_file))
	error=1
endif


endif


if(rtm_MCS.ne.0)    then

    file_name = trim(adjustl(folder_input_mcs)) // trim(adjustl(mcs_data))
    INQUIRE(FILE=file_name, EXIST=data_exist)
    if(.NOT. data_exist)     then
        write(*,*)'ERROR: "file_data_list: " file not found'
        error=1
    else
        file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(mcs_data))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=11)
            nlines = nlines + 1
        enddo
        11 close (10)

        if(nlines.ne.NumShots_mcs)    then
            write(*,*)'ERROR: ', abs(nlines-NumShots_mcs),' MCS data files differ in navigation and shot data list'
            error=1
        endif
    endif

    file_name = trim(adjustl(folder_input_mcs)) // trim(adjustl(nav_shot_mcs))
    INQUIRE(FILE=file_name, EXIST=nav_exist)
    if(.NOT. nav_exist)    then
        write(*,*)'ERROR: nav_shot_mcs not found: ',trim(adjustl(file_name))
        write(*,*)'ERROR: nav_shot_mcs not found: ',trim(adjustl(folder_input_mcs))
        write(*,*)'ERROR: nav_shot_mcs not found: ',trim(adjustl(nav_shot_mcs))
        error=1
    endif

    if(dt_mcs.eq.0) then
	write(*,*) 'ERROR: Please give a value to dt (seconds) in ',trim(adjustl(par_file))
	error=1
    endif

    if(nt_mcs.eq.0) then
	write(*,*)'ERROR: Please give a value to nt in ',trim(adjustl(par_file))
	error=1
    endif

endif

if(rtm_WAS.ne.0.and.shot_depth_was.eq.0) then
	write(*,*)'ERROR: Please give a value to shot_depth_was (meters) in ',trim(adjustl(par_file))
	error=1
endif

if(rtm_MCS.ne.0.and.shot_depth_mcs.eq.0) then
    write(*,*)'ERROR: Please give a value to shot_depth_mcs (meters) in ',trim(adjustl(par_file))
    error=1
endif

if(error.eq.1) then
        call ascii_art(2)
        stop
endif

end subroutine warnings_errors

end module mod_parfile

module mod_data_arrays

implicit none

!!      WAS DATA
        integer :: shotID_was_1, shotID_was_n
        integer, allocatable :: shotID_nav_was(:),shotID_su_was(:,:),shotID_was_(:)
        real, allocatable :: pos_bat_was(:) !!position shot
        real, allocatable :: pos_shot_was(:) !!position shot

        integer, allocatable :: obsID_nav(:)
        integer, allocatable :: pos_xobs_grid(:) !!position obs x
        integer, allocatable :: pos_zobs_grid(:) !!position obs z
        real, allocatable :: pos_xobs(:) !!position obs x
        real, allocatable :: pos_zobs(:) !!position obs z
        real, allocatable    :: bat_model_was(:)
        integer, allocatable :: nxSou_WAS(:),nySou_WAS(:)
        integer, allocatable :: nxRec_WAS(:,:),nyRec_WAS(:,:)
        integer, allocatable :: sizeof_was(:)

!!      MCS DATA
        integer :: shotID_mcs_1, shotID_mcs_n
        integer, allocatable :: shotID_nav_mcs(:),shotID_su_mcs(:),shotID_mcs_(:)
        real, allocatable :: pos_bat_mcs(:) !!position shot
        real, allocatable :: pos_shot_mcs(:) !!position shot
        real, allocatable :: pos_rec(:,:) !!position mcs x
        real, allocatable :: MCS_raw_Source(:,:),MCS_raw_Data(:,:)
        real, allocatable    :: bat_model_mcs(:)
        real, allocatable    :: offset_su(:,:)
        integer, allocatable :: nxSou_MCS(:),nySou_MCS(:)
        integer, allocatable :: nxRec_MCS(:,:),nyRec_MCS(:,:)

!!      GENERAL
        real, allocatable :: vel(:,:)

contains

subroutine allocate_data_arrays()
use mod_parfile

implicit none

    if(rtm_WAS.ne.0)    then

	allocate(sizeof_was(NumShots_was))
        allocate(shotID_nav_was(NumShots_was))
        allocate(shotID_su_was(NumShots_was,NumOBS),shotID_was_(NumShots_was))
        allocate(pos_shot_was(NumShots_was))
        allocate(pos_bat_was(NumShots_was))

        shotID_was_1=0;shotID_was_n=0;
        shotID_nav_was=0;
        shotID_su_was=0;shotID_was_=0;
        pos_shot_was=0;
        pos_bat_was=0;

	allocate(obsID_nav(NumOBS))
        allocate(pos_xobs(NumOBS),pos_zobs(NumOBS))
        allocate(pos_xobs_grid(NumOBS),pos_zobs_grid(NumOBS))
	obsID_nav=0;
        pos_xobs=0;pos_xobs_grid=0;
        pos_zobs=0;pos_zobs_grid=0;

!!      Geometry en WAS
        allocate(nxSou_WAS(NumSou_WAS),nySou_WAS(NumSou_WAS))
        allocate(nxRec_WAS(NumRec_WAS,NumSou_WAS),nyRec_WAS(NumRec_WAS,NumSou_WAS))
        nxSou_WAS=0;nySou_WAS=0;nxRec_WAS=0;nyRec_WAS=0;

    endif

        if(rtm_MCS.ne.0)    then

            allocate(shotID_nav_mcs(NumShots_mcs))
            allocate(shotID_su_mcs(NumShots_mcs),shotID_mcs_(NumShots_mcs))
            allocate(pos_shot_mcs(NumShots_mcs))
            allocate(pos_bat_mcs(NumShots_mcs))

            shotID_mcs_1=0;shotID_mcs_n=0;
            shotID_nav_mcs=0;
            shotID_su_mcs=0;shotID_mcs_=0;
            pos_shot_mcs=0;
            pos_bat_mcs=0;

            allocate(offset_su(NumChannels,NumShots_mcs))
            allocate(pos_rec(NumChannels,NumShots_mcs))
            offset_su=0;pos_rec=0;

!!		Geometry en MCS
            allocate(nxSou_MCS(NumSou_MCS),nySou_MCS(NumSou_MCS))
            allocate(nxRec_MCS(NumRec_MCS,NumSou_MCS),nyRec_MCS(NumRec_MCS,NumSou_MCS))
            allocate(MCS_raw_Source(nt_mcs,NumSou_MCS),MCS_raw_Data(nt_mcs,NumRec_MCS))
            nxSou_MCS=0;nySou_MCS=0;nxRec_MCS=0;nyRec_MCS=0;
            MCS_raw_Source=0.;MCS_raw_Data=0.;

        endif


end subroutine allocate_data_arrays

subroutine allocate_model_arrays_was(rank)
use mod_parfile

    implicit none
	integer	:: rank

    allocate(bat_model_was(nxmodel))
    bat_model_was=0.;

end subroutine allocate_model_arrays_was

subroutine allocate_model_arrays_mcs(rank)
use mod_parfile

    implicit none
    integer :: rank
    allocate(bat_model_mcs(nxmodel))
    bat_model_mcs=0.;

end subroutine allocate_model_arrays_mcs

subroutine allocate_model_arrays(rank)
use mod_parfile

    implicit none
    integer :: rank

	allocate(vel(nymodel,nxmodel))
        vel=0.;

end subroutine allocate_model_arrays

subroutine deallocate_arrays()
use mod_parfile


        if(rtm_WAS.ne.0)    then
	    deallocate(obsID_nav,pos_xobs_grid,pos_zobs_grid,pos_xobs,pos_zobs)
            deallocate(shotID_nav_was)
            deallocate(shotID_su_was,shotID_was_)
            deallocate(pos_shot_was)
            deallocate(pos_bat_was)
            deallocate(nxSou_WAS,nySou_WAS,nxRec_WAS,nyRec_WAS)
            deallocate(bat_model_was)
	    deallocate(sizeof_was)
        endif

        if(rtm_MCS.ne.0)    then
            deallocate(shotID_nav_mcs)
            deallocate(shotID_su_mcs,shotID_mcs_)
            deallocate(pos_shot_mcs)
            deallocate(offset_su)
            deallocate(pos_rec)
            deallocate(pos_bat_mcs)
            deallocate(nxSou_MCS,nySou_MCS,nxRec_MCS,nyRec_MCS)
            deallocate(MCS_raw_Source,MCS_raw_Data)
            deallocate(bat_model_mcs)
        endif

        deallocate(vel)

end subroutine deallocate_arrays


end module mod_data_arrays
