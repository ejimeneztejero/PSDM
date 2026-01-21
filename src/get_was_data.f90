!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines for extracting input data features
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_geom_was_data(rank,iOBS)

USE mod_parfile
USE mod_data_arrays

implicit none

integer :: i,j,k,dif
integer :: rank,iOBS
real	:: dd
character(len=500) :: file_name
logical :: file_exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call check_obs(rank,iOBS)    !!check shotnumber content in SU file
call get_geometry_was(rank,iOBS,dmodel_prop)
call relocation_zOBS(rank)	!!ojo, arreglar para iOBS

END subroutine get_geom_was_data

subroutine relocation_zOBS(rank)

USE mod_parfile
USE mod_data_arrays

implicit none
integer	:: i,rank
real	:: dif

	dif=0

	if(rank.eq.0)	then
		do i=1,NumOBS
			dif=pos_zobs(i)-bat_model_was(pos_xobs_grid(i))
			write(*,*)i,". Depth difference (m) between OBS ",i," and the corresponding bathymetry: ", dif
			write(*,*)i,". OBS, Z Position (m): ",pos_zobs(i),", bathymetry: ",bat_model_was(pos_xobs_grid(i))
			!!write(*,*)i,"pos_xobs_grid(i)= ", pos_xobs_grid(i)
		enddo
	endif

	do i=1,NumOBS
		pos_zobs(i)=bat_model_was(pos_xobs_grid(i))
		pos_zobs_grid(i)=1+ceiling(pos_zobs(i)/dmodel)
	enddo

	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)"All the OBSs are now relocated at the exact model bathymetry projection"

	dif=0
	if(rank.eq.0)	then
		do i=1,NumOBS
			dif=pos_zobs(i)-bat_model_was(pos_xobs_grid(i))
			write(*,*)i,". OBS, Z Position (m): ",pos_zobs(i)
			write(*,*)i,". OBS, Z Position (model dim): ",pos_zobs_grid(i)
!			write(*,*)"Depth difference (m) between OBS ",i," and the corresponding bathymetry: ", dif
		enddo
	endif

end subroutine relocation_zOBS

subroutine read_nav_shot_was(rank)

USE mod_parfile
USE mod_data_arrays

implicit none

integer :: i,j,k,rank
integer :: iline,itr,ERR,nlines
integer :: shotID
real :: xs,ys,bat,x1,xn,y1,yn,add1
character (len=500) :: file_name

add1=d_shot_model

!!	LECTURA DATOS SHOTGATHERS
file_name= trim(adjustl(folder_input_was)) // trim(adjustl(nav_shot_was))
open(unit=10,file=file_name,status='old')

do i=1,NumShots_was !!number of lines nav_shot

    shotID=0;xs=0;ys=0;bat=0.;!ojo

    if(nav_was_cols.eq.11)read(10,*)shotID,xs
    if(nav_was_cols.eq.111)read(10,*)shotID,xs,ys
    if(nav_was_cols.eq.1111)read(10,*)shotID,xs,ys,bat

    if(i.eq.1)     then
        x1=xs;y1=ys;
    endif
    if(i.eq.NumShots_was)     then
        xn=xs;yn=ys;
    endif

    shotID_nav_was(i)=shotID
    pos_shot_was(i) = add1+sqrt((xs-x1)**2.+(ys-y1)**2.)
    pos_bat_was(i) = abs(bat)

    if(shot_i.eq.shotID)line_i=i
    if(shot_f.eq.shotID)line_f=i

enddo

   dshots_was=abs(xn-x1)/(NumShots_was-1)
   length_model_was=xs+2*added_space
   nxmodel_was=1+ceiling(length_model_was/dmodel)

   if(rank.eq.0)	then
       write(*,*)"dshots_was: ",dshots_was
	   if(length_model_was.gt.length_model)write(*,*)"WARNING: length_model_was>length_model"
	   if(nxmodel_was.gt.nxmodel)write(*,*)"WARNING: nxmodel_was>nxmodel"
   endif

CLOSE(10)

end subroutine  read_nav_shot_was

subroutine  read_nav_obs(rank)

USE mod_parfile
USE mod_data_arrays

implicit none

integer :: i,j,k,io,rank
real :: xs,xz
character (len=500) :: file_name

!!	LECTURA DATOS obs

file_name= trim(adjustl(folder_input_was)) // trim(adjustl(nav_obs))
open(unit=10,file=file_name,status='old')

do i=1,NumOBS !!number of lines nav_obs

    io=0;xs=0;xz=0
    read(10,*)io,xs,xz	!num obs, posX, posZ

    obsid_nav(i)=io

    pos_xobs(i) = abs(xs)
    pos_xobs_grid(i)=1+ceiling(pos_xobs(i)/dmodel)

    pos_zobs(i) = abs(xz)
    pos_zobs_grid(i)=1+ceiling(pos_zobs(i)/dmodel)

enddo

CLOSE(10)

end subroutine  read_nav_OBS

subroutine get_geometry_WAS(rank,iOBS,dd)

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: i,j,k,iOBS,irec,rank
real :: dd

!! SOURCES IN THE INVERSION ARE THE OBSs

!do iOBS=1,NumOBS

nxSou_WAS(iOBS)=pos_xobs_grid(iOBS)  !!cambiarlo por "pos_shot_grid()"
nySou_WAS(iOBS)=pos_zobs_grid(iOBS)  !!cambiarlo por "pos_shot_grid()"

do irec=1,NumRec_WAS	!!receivers are the shots
	nxRec_WAS(irec,iOBS)=1+ceiling(pos_shot_was(irec)/dd)
	nyRec_WAS(irec,iOBS)=1+ceiling(shot_depth_was/dd)
enddo

!do irec=1,NumRec_WAS
!	write(22,*)irec,nxSou_WAS(iOBS),nySou_WAS(iOBS),nxRec_WAS(irec,iOBS),nyRec_WAS(irec,iOBS)
!enddo

!enddo
   
end subroutine get_geometry_WAS

subroutine get_obs_data(iOBS,raw_data)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE

    INTEGER :: nh,rank
    INTEGER :: i,j,k,iOBS,icount
    INTEGER :: iSS,NumSou_SS,NumRec_SS
    integer :: SourceNum

    integer(4) :: pos_r
    real    :: raw_data(nt_was,NumRec_WAS)
    real(4), allocatable :: sudata(:,:)

    nh = size_su_header

    allocate(sudata(nt_was_max+nh,NumRec_SS))
    sudata=0.;
    raw_Data=0.

    do j=1,NumRec_SS
	pos_r = 1 + (j-1)*(nh+nt_was_max)*4
	READ(unit_was+iOBS, pos=pos_r)  sudata(1:nt_was_max+nh,j)
    enddo

    do j = 1,NumRec_SS
        do i = 1, nt_was
                raw_data(i,j) = sudata(i+nh,j)
        enddo
    enddo


    deallocate(sudata)

END SUBROUTINE get_obs_data

subroutine check_nav_obs(rank)

use mod_parfile
USE mod_data_arrays

implicit none

integer :: i,rank
integer :: nlines,nobs
integer :: obsi,obsf

CHARACTER(len=50) :: access,form,Str
CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     NAVIGATION FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

file_name= trim(adjustl(folder_input_was)) // trim(adjustl(nav_obs))
open(unit=10,file=file_name,status='old')
nlines=0
do
    read(10,*, END=10)
    nlines = nlines + 1
enddo
10 close (10)

if(rank.eq.0.and.nlines.lt.NumOBS)      then
    write(*,*)'ERROR: navigation, some OBS are are missing in nav_obs_file'
    call ascii_art(2)
    stop
endif

open(unit=10,file=file_name,status='old')
do i=1,NumOBS
    read(10,*)nobs
    if(i.eq.1)obsi=nobs
    if(i.eq.NumOBS)obsf=nobs
enddo
close(10)

if(rank.eq.0)   then
    write(*,*)'Navigation file contains OBS from ',obsi,' to ',obsf
    write(*,*)'Num of OBS:', obsf-obsi+1,', Num Lines nav_obs:',nlines!
endif

end subroutine check_nav_obs

subroutine check_nav_shot_was(rank)

use mod_parfile
USE mod_data_arrays

implicit none

        integer :: i,nlines,nums,rank
        integer :: shoti,shotf,nx_user
        real    :: xi,xf,xs,bat
        CHARACTER(len=50) :: access,form,Str
        CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     NAVIGATION FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        file_name= trim(adjustl(folder_input_was)) // trim(adjustl(nav_shot_was))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=10)
            nlines = nlines + 1
        enddo
        10 close (10)

        open(unit=10,file=file_name,status='old')
        do i=1,nlines

                read(10,*)nums,xs,bat

                if(i.eq.1)      then
                        shoti=nums
                        xi=xs
                endif

                if(i.eq.nlines) then
                        shotf=nums
                        xf=xs
                endif

        enddo
        close(10)

        if(rank.eq.0)   then
                write(*,*)'Navigation file contains shots from ',shoti,' to ',shotf
                write(*,*)'NumShots total:', shotf-shoti+1,', Num Lines Nav_file:',nlines
        endif


        if(nlines.ne.NumShots_was)  then

                if(rank.eq.0)then
                        write(*,*)'ERROR: navigation, some shots are missing in nav_shot_was'
                        call ascii_art(2)
                        stop
                endif
        endif

        nx_user=1+ceiling(xf/dmodel)

        if(rank.eq.0)   then
                write(*,*)'FOR WAS: nxmodel should be >', nx_user
        endif

        if(nx_user.gt.nxmodel)    then

            if(rank.eq.0)then
                write(*,*)
                write(*,*)"ERROR: fix you nxmodel. It should be bigger or equal to:"
                write(*,*)"1+distance(shotN-shot1)/dmodel = ", nx_user
                write(*,*)", where: distance(shotN-shot1) = ",xf-xi," and dmodel= ",dmodel
                call ascii_art(2)
                stop
            endif

    endif

end subroutine check_nav_shot_was

SUBROUTINE open_obs(rank,iOBS)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE
integer     :: i,j,k,iOBS,rank
CHARACTER(len=50) :: access,form,OBS_num
CHARACTER(len=500) :: file_name

access = 'stream'
form = 'unformatted'

!do iOBS=1,NumOBS

file_name = trim(adjustl(folder_input_was)) // su_file_was(iOBS)

if(endianness_machine.eq.0)     then
        if(endianness_data.eq.0)open(unit_was+iOBS,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')
        if(endianness_data.eq.1)open(unit_was+iOBS,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')
endif

if(endianness_machine.eq.1)     then
        if(endianness_data.eq.0)open(unit_was+iOBS,FILE=file_name,ACCESS=access,FORM=form,CONVERT='LITTLE_ENDIAN',STATUS='old')
        if(endianness_data.eq.1)open(unit_was+iOBS,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')
endif

INQUIRE(FILE=file_name, SIZE=sizeof_was(iOBS))

!enddo

END SUBROUTINE open_obs

subroutine check_obs(rank,iOBS)

USE mod_parfile
USE mod_data_arrays

implicit none

INTEGER :: i,j,k,icount,iOBS,rank
INTEGER :: nh

INTEGER*4 :: shotID,pos_read
INTEGER*4 :: pos_byte,pos_byte2

CHARACTER(len=50) :: Str
CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    .SU FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! TO READ ShotID:

pos_read=byte_shotnumber_was
nh = size_su_header    !Size su header = 60*4 bytes

!do iOBS=1,NumOBS

! READ SU FILES
file_name = trim(adjustl(folder_input_was)) // su_file_was(iOBS)

READ(unit_was+iOBS,pos=pos_read) shotID_was_1
READ(unit_was+iOBS,pos=sizeof_was(iOBS)-(nh+nt_was_max)*4+pos_read) shotID_was_n

if(rank.eq.0)	then
	write(*,*) "OBS file: ",trim(adjustl(file_name))," contains shots from ",shotID_was_1,' to ',shotID_was_n
	write(*,*)'A total number of ',1+shotID_was_n-shotID_was_1

	if(1+shotID_was_n-shotID_was_1.ne.NumShots_was)    then
	    write(*,*)"Warning, the number of shotgathers at OBS num: ",iOBS," not coincident with variable NumShots"
	endif
endif

!------------------------------------------------------------------------
! Loop shot-by-shot

pos_byte=1
icount=1

READ(unit_was+iOBS,pos=pos_read) shotID

do while(icount.le.NumShots_was)

    shotID_su_was(icount,iOBS)=shotID
    pos_byte=pos_byte+(nh+nt_was_max)*4

        if(pos_byte.ge.sizeof_was(iOBS)) then
                exit
        endif
        if(pos_byte.gt.0.and.pos_byte.lt.sizeof_was(iOBS))     then
                READ(unit_was+iOBS,pos=pos_byte+pos_read-1) shotID
                icount=icount+1
        endif

enddo

do j=1,NumShots_was
    if(shotID_su_was(j,iOBS).eq.0)    then
        write(*,*)'WARNING, not recorded in .su file, shot number: ',j,' with shotID: ',shotID_su_was(icount,iOBS)
    endif
enddo

!enddo

end subroutine check_obs

subroutine Ricker_was(rank,S0)

use mod_parfile
use mod_data_arrays

implicit none

    integer :: i,j,k,rank,iOBS
!    integer :: SourceNum
    real :: pi,f0,wc,t0,tmp,S0(nt_prop),t(nt_prop)
    real :: n

        !!!!!------ Using a Ricker wavelet as impulse

        pi=3.14159265
        f0=freq_ricker
        wc=f0*2*pi
        t0=1.5/f0
        tmp=10./f0

        do i=1,nt_prop
                t(i)=(i-1)*dt_prop
        enddo

	S0=(1.-2.*(pi*f0*(t-t0))**2)*exp(-1.*(pi*f0*(t-t0))**2)

        do i=1,nt_prop
                if(t(i).ge.tmp)then
                        S0(i)=0.
                endif
        enddo


end subroutine Ricker_was

subroutine get_bathymetry_model_was(rank)

USE mod_parfile
USE mod_data_arrays

implicit none

        integer :: i,j,ERR,rank

        real, allocatable :: xmodel(:),xmodel_(:),xbat_(:)
        character (len=500) :: file_name

        allocate(xmodel(nxmodel))
        allocate(xmodel_(NumShots_was+2),xbat_(NumShots_was+2))

        do i=2,NumShots_was+1
                xbat_(i)=pos_bat_was(i-1)
                xmodel_(i)=pos_shot_was(i-1)
        enddo

        xmodel_(1)=0;
        xmodel_(NumShots_was+2)=length_model;!!x

        xbat_(1)=xbat_(2);
        xbat_(NumShots_was+2)=xbat_(NumShots_was+1);!!y

!       do i=1,NumShots_was+2
!		    write(17,*)i,xmodel_(i),xbat_(i)
!	    enddo

        do i=1,nxmodel
            bat_model_was(i)=0
            xmodel(i)=(i-1)*dmodel
!           write(18,*)i,xmodel(i),bat_model_was(i)
        enddo

        ERR=0
        call INTRPL(NumShots_was+2,xmodel_,xbat_,nxmodel,xmodel,bat_model_was,ERR)

        file_name=trim(adjustl(folder_output))// "bathymetry_model_was_meters.txt"
        open(12,file=file_name,status='unknown')
        do i=1,nxmodel
                write(12,*)xmodel(i),-bat_model_was(i)
        enddo
        close(12)

        deallocate(xmodel)
        deallocate(xmodel_,xbat_)

end subroutine get_bathymetry_model_was

subroutine close_was_files(iOBS)

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: k,iOBS

!do iOBS=1,NumSou_WAS
!       write(*,*)"close obs, unit: ",unit_was+iOBS
        close(unit_was+iOBS)
!enddo

end subroutine close_was_files
