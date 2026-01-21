!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines for extracting input data features
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_geom_mcs_data(rank,icount)

USE mod_parfile
USE mod_data_arrays

implicit none

integer :: i,j,k,dif
integer :: rank,icount
character(len=500) :: file_name
logical :: file_exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! FOR SHOTS-MCS

call get_pos_channels(rank,icount)    !!ojo, repasar, añadir lectura en cabecera: sx
call get_geometry_mcs(rank,icount,dmodel_prop)

END subroutine get_geom_mcs_data

subroutine read_nav_shot_mcs(rank)	

USE mod_parfile
USE mod_data_arrays

implicit none

integer :: i,j,k,icount,rank
integer :: iline,itr,ERR,nlines
integer :: shotID
integer :: ios

real :: xs,ys,bat,x1,xn,y1,yn,add1

character (len=500) :: file_name,line

add1=d_shot_model

!!	LECTURA DATOS SHOTGATHERS

file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(nav_shot_mcs))
open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=10)
            nlines = nlines + 1
        enddo
        10 close (10)

open(unit=10,file=file_name,status='old')
do icount=1,nlines

    shotID=0;xs=0;ys=0.;bat=0.;

    if(nav_mcs_cols.eq.11)read(10,*)shotID,xs
    if(nav_mcs_cols.eq.111)read(10,*)shotID,xs,ys
    if(nav_mcs_cols.eq.1111)read(10,*)shotID,xs,ys,bat

    if(dir_boat.eq.1)	then
    if(icount.eq.1)	then
	x1=xs;y1=ys;
    endif
    endif

    if(dir_boat.eq.-1)	then
    if(icount.eq.nlines)	then
	x1=xs;y1=ys;
    endif
    endif

    if(icount.eq.1.and.shot_i.eq.0)shot_i=shotID
    if(icount.eq.nlines.and.shot_f.eq.0)shot_f=shotID

    if(shot_i.gt.0.and.shot_i.eq.shotID)line_i=icount
    if(shot_f.gt.0.and.shot_f.eq.shotID)line_f=icount

enddo

CLOSE(10)
open(unit=10,file=file_name,status='old')

do icount=1,nlines

    shotID=0;xs=0;ys=0.;bat=0.;

    if(nav_mcs_cols.eq.11)read(10,*)shotID,xs
    if(nav_mcs_cols.eq.111)read(10,*)shotID,xs,ys
    if(nav_mcs_cols.eq.1111)read(10,*)shotID,xs,ys,bat

    shotID_nav_mcs(icount)=shotID
    pos_shot_mcs(icount) = add1+sqrt((xs-x1)**2.+(ys-y1)**2.) 
    pos_bat_mcs(icount) = abs(bat)

!    write(13,*)shotID,xs,bat,pos_bat_mcs(icount)

enddo
CLOSE(10)

if(rank.eq.0)	then

write(*,*)"shot init, shot final, nlines",shot_i,shot_f,nlines,line_i,line_f
write(*,*)"pos_shot1,pos_shotn: ",pos_shot_mcs(1),pos_shot_mcs(nlines)

if(line_i.eq.0)	then
	write(*,*)"ERROR shot_i not found ", shot_i,line_i
	stop
endif

if(line_f.eq.0)	then
	write(*,*)"ERROR shot_f not found ",shot_f,line_f
	stop
endif


endif

end subroutine  read_nav_shot_mcs

subroutine get_geometry_mcs(rank,icount,dd)

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: i,j,k,icount,irec,rank
real	:: dd

nxSou_MCS(icount)=1+floor(pos_shot_mcs(icount)/dd)
nySou_MCS(icount)=1+floor(shot_depth_mcs/dd)

!write(22,*)icount,nxSou_MCS(icount),pos_shot_mcs(icount)

do irec=1,NumRec_MCS
	nxRec_MCS(irec,icount)=1+floor(pos_rec(irec,icount)/dd)
	nyRec_MCS(irec,icount)=1+floor(streamer_depth/dd)
!	write(22,*)icount,irec,nxRec_MCS(irec,icount),pos_rec(irec,icount)
enddo

end subroutine get_geometry_mcs

SUBROUTINE get_shot_data(icount,raw_data,delay_data)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE

INTEGER :: nh
INTEGER :: i,j,k,icount,id,idp
INTEGER(4) :: pos_r
real    :: raw_data(nt_mcs,NumRec_MCS)
real    :: delay_data(nt_mcs,NumRec_MCS)

REAL(4), ALLOCATABLE :: sudata(:,:)	

nh = size_su_header

allocate(sudata(nt_mcs_max+nh,NumRec_MCS))
sudata=0.;

!!! Loop trace by trace
do j=1,NumRec_MCS
        pos_r = 1 + (j-1)*(nh+nt_mcs_max)*4
        READ(unit_mcs+icount, pos=pos_r)  sudata(1:nt_mcs_max+nh,j)
enddo

!write(112,*)id,idp,delay_data

do j = 1,NumRec_MCS
do i = 1,nt_mcs
	raw_Data(i,j) = sudata(i+nh,j)
enddo
enddo

if(id.ne.0)	then

	id=delay_time_data/dt_mcs
	idp=abs(id)

	if(id.ge.0)	then
		do j = 1,NumRec_MCS
	        do i = 1,nt_mcs-id
	                delay_data(i+id,j) = raw_data(i,j)
	        enddo
		enddo
	endif

	if(id.lt.0)	then
		do j = 1,NumRec_MCS
	        do i = idp+1,nt_mcs
	                delay_data(i-idp,j) = raw_data(i,j)
	        enddo
		enddo
	endif

endif

!do j=1,nt_mcs
!	write(1001,'(20000(e12.5,2x))') (raw_Data(j,i),i=1,NumRec_MCS)
!	write(3159,'(20000(e12.5,2x))') (raw_Data(j,i),i=1,NumRec_MCS)
!enddo

deallocate(sudata)

END SUBROUTINE get_shot_data

subroutine check_nav_shot_mcs(rank)

use mod_parfile
USE mod_data_arrays

implicit none

        integer :: i,nlines,rank
        integer :: shoti,shotf,shotID
        real    :: xi,xf,xs,ys,bat
        CHARACTER(len=50) :: access,form,Str
        CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     NAVIGATION FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(nav_shot_mcs))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=10)
            nlines = nlines + 1
        enddo
        10 close (10)

        open(unit=10,file=file_name,status='old')
        do i=1,nlines

    		shotID=0;xs=0;ys=0.;bat=0.;

    		if(nav_mcs_cols.eq.11)read(10,*)shotID,xs
    		if(nav_mcs_cols.eq.111)read(10,*)shotID,xs,ys
    		if(nav_mcs_cols.eq.1111)read(10,*)shotID,xs,ys,bat

                if(i.eq.1)      then
                        shoti=shotID
                        xi=xs
                endif

                if(i.eq.nlines) then
                        shotf=shotID
                        xf=xs
                endif

        enddo
        close(10)

        if(rank.eq.0)   then
                write(*,*)'Navigation file in mcs contains shots from ',shoti,' to ',shotf
                write(*,*)'Num Lines Nav_file:',nlines
        endif

        if(nlines.lt.NumShots_mcs)  then

                if(rank.eq.0)then
                        write(*,*)'ERROR: navigation, some shots are missing in nav_shot_mcs: ',nlines,NumShots_mcs
                        call ascii_art(2)
                        stop
                endif
        endif
        if(nlines.gt.NumShots_mcs)  then

                if(rank.eq.0)then
                        write(*,*)'WARNING: navigation contains more shots than wanted to be read: ',nlines,NumShots_mcs
                endif
        endif


end subroutine check_nav_shot_mcs

SUBROUTINE open_shot_mcs(rank,icount)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE
integer     :: i,j,k,icount,rank
CHARACTER(len=50) :: OBS_num
CHARACTER(len=500) :: file_name

CHARACTER(len=50) :: access,form

access = 'stream'
form = 'unformatted'

file_name = trim(adjustl(folder_input_mcs)) // trim(adjustl(su_file_mcs(icount)))

if(endianness_machine.eq.0)     then
        if(endianness_data.eq.0)open(unit_mcs+icount,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')
        if(endianness_data.eq.1)open(unit_mcs+icount,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')
endif

if(endianness_machine.eq.1)     then
        if(endianness_data.eq.0)open(unit_mcs+icount,FILE=file_name,ACCESS=access,FORM=form,CONVERT='LITTLE_ENDIAN',STATUS='old')
        if(endianness_data.eq.1)open(unit_mcs+icount,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')
endif


END SUBROUTINE open_shot_mcs


subroutine Ricker_mcs(rank,S0)

use mod_parfile
use mod_data_arrays

implicit none

    integer :: i,j,k,rank,NumSou_SS
    integer :: SourceNum
    real :: pi,f0,wc,t0,tmp,S0(nt_prop),t(nt_prop)
    real :: n,dt

	!!!!!------ Using a Ricker wavelet as impulse

        pi=3.14159265
        f0=freq_ricker
        wc=f0*2*pi
        t0=delay_ricker	
        tmp=10./f0
	!tmp=0.

        do i=1,nt_prop
                t(i)=(i-1)*dt_prop
        enddo

	S0=(1.-2.*(pi*f0*(t-t0))**2)*exp(-1.*(pi*f0*(t-t0))**2)

        do i=1,nt_prop
                if(t(i).ge.tmp)then
                        S0(i)=0.
                endif
        enddo

end subroutine Ricker_mcs

subroutine get_bathymetry_model_mcs(rank)

USE mod_parfile
USE mod_data_arrays

implicit none

        integer :: i,j,ERR,rank

        real, allocatable :: xmodel(:),xmodel_(:),xbat_(:)
        character (len=500) :: file_name

        allocate(xmodel(nxmodel))
        allocate(xmodel_(NumShots_mcs+2),xbat_(NumShots_mcs+2))

        do i=2,NumShots_mcs+1
                xbat_(i)=pos_bat_mcs(i-1)
                xmodel_(i)=pos_shot_mcs(i-1)
        enddo

        xmodel_(1)=0;
        xmodel_(NumShots_mcs+2)=length_model;!!x

        xbat_(1)=xbat_(2);
        xbat_(NumShots_mcs+2)=xbat_(NumShots_mcs+1);!!y

        do i=1,nxmodel
            bat_model_mcs(i)=0
            xmodel(i)=(i-1)*dmodel
        enddo

        do i=1,NumShots_mcs+2
		write(14,*) i,xmodel_(i),xbat_(i)
	enddo

        do i=1,nxmodel
		write(15,*) i-1,xmodel(i)
	enddo


        ERR=0
        call INTRPL(NumShots_mcs+2,xmodel_,xbat_,nxmodel,xmodel,bat_model_mcs,ERR)

!        do i=1,NumShots_mcs+2
!		write(16,*) i,xmodel_(i),xbat_(i)
!	enddo

!        do i=1,nxmodel
!		write(17,*) i-1,xmodel(i),bat_model_mcs(i)
!	enddo


        file_name=trim(adjustl(folder_output))// "bat_mcs.txt"
        open(12,file=file_name,status='unknown')
        do i=1,nxmodel
                write(12,*)xmodel(i),-bat_model_mcs(i)
        enddo
        close(12)

        deallocate(xmodel)
        deallocate(xmodel_,xbat_)

end subroutine get_bathymetry_model_mcs


subroutine get_pos_channels(rank,icount)

USE mod_parfile
USE mod_data_arrays

implicit none
include 'mpif.h'

    integer :: numtasks,rank,ierr,errcode,status(MPI_STATUS_SIZE)
    integer :: i,j,k,icount,irec
    INTEGER :: nh

    INTEGER*2 :: pos_scalco,scalco

    INTEGER*4 :: shotID,pos_read
    INTEGER*4 :: pos_byte_trace

    INTEGER*4 :: pos_offset,pos_sx,pos_sy
    INTEGER*4 :: sx,sy,offset

    call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    ierr=0;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    .SU FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! TO READ ShotID:

nh=size_su_header ! Size su header = 60*4 bytes
pos_read=byte_shotnumber_mcs
pos_offset = byte_offset

!------------------------------------------------------------------------

!do icount=1,NumShots_MCS

	READ(unit_mcs+icount,pos=pos_read) shotID
	
	shotID_su_mcs(icount)=shotID

        if(offset_header.eq.1.and.offset_unit.gt.0)    then

	!    Loop trace by trace
	    do irec=1,NumRec_MCS

            pos_byte_trace=(irec-1)*(nh+nt_mcs_max)*4
            READ(unit_mcs+icount,pos=pos_byte_trace+pos_offset) offset

            offset_su(irec,icount)=abs(offset)*offset_unit

	    pos_rec(irec,icount)=pos_shot_mcs(icount)&
	    -dir_boat*offset_su(irec,icount)	! (+1 hacia la derecha, -1 hacia la izda)

	    enddo   !irec

        endif !offset_header

        if(offset_header.eq.1.and.offset_unit.lt.0)    then

	!    Loop trace by trace
	    do irec=1,NumRec_MCS

            pos_byte_trace=(irec-1)*(nh+nt_mcs_max)*4
            READ(unit_mcs+icount,pos=pos_byte_trace+pos_offset) offset

            offset_su(irec,icount)=abs(offset)/abs(offset_unit)

	    pos_rec(irec,icount)=pos_shot_mcs(icount)&
	    -dir_boat*offset_su(irec,icount)	! (+1 hacia la derecha, -1 hacia la izda)

	    enddo   !irec

        endif !offset_header

        if(offset_header.eq.0)    then

	!    Loop trace by trace
	    do irec=1,NumRec_MCS

		pos_rec(irec,icount)=pos_shot_mcs(icount)&
		-dir_boat*(near_offset+(irec-1)*drec)

	    enddo   !irec

        endif !offset_header


!enddo	 !icount

end subroutine get_pos_channels

subroutine close_mcs_files(ishot)

USE mod_parfile
USE mod_data_arrays

implicit none
integer    :: k,ishot

!do ishot=1,NumSou_MCS
!    write(*,*)"close obs, unit: ",unit_was+iOBS
    close(unit_mcs+ishot)
!enddo

end subroutine close_mcs_files
