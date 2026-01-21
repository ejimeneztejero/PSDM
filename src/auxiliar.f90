!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains axiliary subroutines, mostly for reading/writting data
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_su_data(ni,ns,ns_max,unit_su,vp_data)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE

INTEGER :: nh
INTEGER :: i,j,k,ni,ns,ns_max,unit_su
INTEGER(4) :: pos_r
REAL    :: vp_data(ns_max,ni)
REAL(4), ALLOCATABLE :: Data_tmp(:,:),sudata(:,:)

nh = size_su_header

allocate(Data_tmp(ns_max,ni))
allocate(sudata(ns+nh,ni))
Data_tmp=0.;sudata=0.;

do i=1,ni
        pos_r = 1 + (i-1)*(nh+ns)*4
        READ(unit_su, pos=pos_r)  sudata(1:ns+nh,i)
enddo

do i=1,ni
        do j=1,ns_max
                Data_tmp(j,i) = sudata(j+nh,i)
        enddo
enddo

vp_data=Data_tmp

deallocate(Data_tmp)
deallocate(sudata)

END SUBROUTINE get_su_data

subroutine get_vel(rank)

  USE mod_parfile
  USE mod_data_arrays

  implicit none

  integer :: i,j,jwater,rank
  real	:: xi,xj
  real, allocatable :: vel0(:,:)
  character(len=500) :: file_name
  character(len=1000) :: file_init
  CHARACTER(len=50) :: access,form

  access = 'stream'
  form = 'unformatted'

  allocate(vel0(nymodel0,nxmodel))

!!! Vp

    vel=water_velocity

    if(read_vel.eq.1)    then


!!	ARREGLAR

        file_name=trim(adjustl(folder_input_model)) // trim(adjustl(vel_file))

	if(vel_v.eq.1)	then
	        open(unit=10,file=file_name,status='old')
	        do j=1,nymodel
			read(10,*) (vel(j,i),i=1,nxmodel)
		enddo
	endif
	if(vel_xzv.eq.1)	then
	        open(unit=10,file=file_name,status='old')
                do i=1,nxmodel
                        do j=1,nymodel0
                                read(10,*) xi,xj,vel0(j,i)
                        enddo
			vel(1:nymodel,i)=vel0(1:nymodel,i)
                enddo
	endif
	if(vel_su.eq.1)	then
	        open(10,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')
		call get_su_data(nxmodel,nymodel0,nymodel,10,vel)
	endif

	close(10)

	if(dvel.gt.0)	then
        do i=1,nxmodel
        do j=1,nymodel
		if(vel(j,i).gt.1500)vel(j,i)=vel(j,i)+dvel
        enddo
        enddo
	endif

    endif

  deallocate(vel0)

end subroutine get_vel


subroutine is_not_cero(nt,data1d,isnotCero)
implicit none

    integer :: i,nt
    real :: data1d(nt)
    logical :: isnotCero

    isnotCero=.FALSE.

    ! Verificar si todas las componentes son cero
    do i = 1, nt
        if (data1d(i) /= 0.0) then
            isnotCero = .TRUE.
            exit ! Puedes salir del bucle en cuanto encuentres un valor no cero
        end if
    end do

end subroutine is_not_cero

subroutine write_vel(ni,nj,di,dj,trid,file_name,model)	!!paralelizar

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: i,j,k,ni,nj,trid
real	:: xi,xj
real	:: di,dj
real	:: model(nj,ni)

character(len=500) :: Str,Str_nj,Str_di,Str_dj,Str_trid,command
character(len=500) :: file_name,file_name2
character(len=500) :: file_temp,file_temp_filt,file_temp_filt2
character(len=500) :: file_bin,file_su,file_su_filt

	if(save_v.ne.0)	then

		file_name2= trim(adjustl(file_name)) // ".v"

		open(unit=12,file=file_name2,status='unknown')

		do j=1,nj
			write(12,'(*(f0.0,2x))') (model(j,i),i=1,ni)
		enddo
		close(12)

	endif


	if(save_xzv.ne.0)       then

		file_name2= trim(adjustl(file_name)) // ".xzv"
		open(unit=12,file=file_name2,status='unknown')

	        do i=1,ni
			xi=(i-1)*dmodel_prop
        	        do j=1,nj
				xj=(j-1)*dmodel_prop
        	                write(12,*) xi,xj,model(j,i)
        	        enddo
        	enddo

		close(12)

	endif   


	if(save_su.ne.0)       then

                file_temp=trim(adjustl(file_name)) // ".temp"
                open(unit=12,file=file_temp,status='unknown')
                do i=1,ni
                do j=1,nj
                        write(12,*) model(j,i)
                enddo
                enddo
                close(12)

                file_bin= trim(adjustl(file_name)) // ".bin"
                open(unit=14,file=file_bin,status='unknown')

                file_su= trim(adjustl(file_name)) // ".su"
                open(unit=16,file=file_su,status='unknown')

                write(Str_nj,*) nj
                write(Str_di,*) di
                write(Str_dj,*) dj
                write(Str_trid,*) trid

                command= "a2b < " // trim(adjustl(file_temp)) // " n1=1 > " // trim(adjustl(file_bin))
                call system(command)

		command = "suaddhead < " // trim(adjustl(file_bin)) // " ns=" // trim(adjustl(Str_nj)) // " | " // &
        	   "sushw key=d1 a=" // trim(adjustl(Str_dj)) // " | " // &
        	   "sushw key=f1 a=0 | " // &
        	   "sushw key=d2 a=" // trim(adjustl(Str_di)) // " | " // &
        	   "sushw key=f2 a=0 j=" // trim(adjustl(Str_di)) // " | " // &
        	   "suchw key1=cdp key2=tracl a=-" // trim(adjustl(Str_di)) // " b=" // trim(adjustl(Str_di)) // " | " // &
        	   "sushw key=trid a=" // trim(adjustl(Str_trid)) // " > " // trim(adjustl(file_su))

                call system(command)

                command="rm " // trim(adjustl(file_bin)) // " " // trim(adjustl(file_temp))
                call system(command)
		
                close(14)
                close(16)

	endif

!	write(*,*)"FINISHED WRITE"

end subroutine write_vel

subroutine write_psdm(ni,nj,di,dj,file_name,model,trid,last)	!!paralelizar

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: i,j,k,ni,nj,trid,last
real	:: xi,xj
real	:: di,dj
real	:: model(nj,ni)

character(len=500) :: Str,Str_nj,Str_di,Str_dj,Str_trid,command
character(len=500) :: file_name,file_name2
character(len=500) :: file_temp,file_temp_filt,file_temp_filt2
character(len=500) :: file_bin,file_su,file_su_filt

	if(save_v.ne.0)	then

		file_name2= trim(adjustl(file_name)) // ".v"

		open(unit=12,file=file_name2,status='unknown')

		do j=1,nj
			write(12,'(20000(e12.5,2x))') (model(j,i),i=1,ni)
		enddo
		close(12)

	endif


	if(save_xzv.ne.0)       then

		file_name2= trim(adjustl(file_name)) // ".xzv"
		open(unit=12,file=file_name2,status='unknown')

	        do i=1,ni
			xi=(i-1)*dmodel_prop
        	        do j=1,nj
				xj=(j-1)*dmodel_prop
        	                write(12,*) xi,xj,model(j,i)
        	        enddo
        	enddo

		close(12)

	endif   


	if(save_su.ne.0)       then

                file_temp=trim(adjustl(file_name)) // ".temp"
                open(unit=12,file=file_temp,status='unknown')
                do i=1,ni
                do j=1,nj
                        write(12,*) model(j,i)
                enddo
                enddo
                close(12)

                file_bin= trim(adjustl(file_name)) // ".bin"
                open(unit=14,file=file_bin,status='unknown')

                file_su= trim(adjustl(file_name)) // ".su"
                open(unit=16,file=file_su,status='unknown')

                write(Str_nj,*) nj
                write(Str_di,*) di
                write(Str_dj,*) dj
                write(Str_trid,*) trid

                command= "a2b < " // trim(adjustl(file_temp)) // " n1=1 > " // trim(adjustl(file_bin))
                call system(command)

		command = "suaddhead < " // trim(adjustl(file_bin)) // " ns=" // trim(adjustl(Str_nj)) // " | " // &
        	   "sushw key=d1 a=" // trim(adjustl(Str_dj)) // " | " // &
        	   "sushw key=f1 a=0 | " // &
        	   "sushw key=d2 a=" // trim(adjustl(Str_di)) // " | " // &
        	   "sushw key=f2 a=0 j=" // trim(adjustl(Str_di)) // " | " // &
        	   "suchw key1=cdp key2=tracl a=-" // trim(adjustl(Str_di)) // " b=" // trim(adjustl(Str_di)) // " | " // &
        	   "sushw key=trid a=" // trim(adjustl(Str_trid)) // " > " // trim(adjustl(file_su))

                call system(command)

                command="rm " // trim(adjustl(file_bin)) // " " // trim(adjustl(file_temp))
                call system(command)
		
		if(post_filt.ne.0)	then

		if(last.ne.0)	then

	                file_temp_filt=trim(adjustl(file_name)) // "filt.temp"
	                open(unit=18,file=file_temp_filt,status='unknown')

	                file_temp_filt2=trim(adjustl(file_name)) // "filt2.temp"
	                open(unit=19,file=file_temp_filt2,status='unknown')
			
	                file_su_filt= trim(adjustl(file_name)) // "_filt.su"
	                open(unit=20,file=file_su_filt,status='unknown')
		
			command = "sushw < " // trim(adjustl(file_su)) // &
			" key=trid a=1 > " // trim(adjustl(file_temp_filt))
		        call system(command)

			command = "sufilter < " // trim(adjustl(file_temp_filt)) // &
			" f=3,6,60,80 > " // trim(adjustl(file_temp_filt2))
	                call system(command)

			command = "sushw < " // trim(adjustl(file_temp_filt2)) // &
			" key=trid a=" // trim(adjustl(Str_trid)) // " > " // trim(adjustl(file_su_filt))
	                call system(command)
		
	                command="rm " // trim(adjustl(file_temp_filt)) 
	                call system(command)
	                command="rm " // trim(adjustl(file_temp_filt2)) 
	                call system(command)

			close(18)
			close(19)
			close(20)

		endif
		endif

                close(14)
                close(16)

	endif

!	write(*,*)"FINISHED WRITE"

end subroutine write_psdm

SUBROUTINE ascii_art(i)

USE mod_parfile

implicit none


integer :: i

if(i.eq.1)	then

write(*,*)
write(*,*)
write(*,*)                                                                      
write(*,*)"		Release (2025)			"
write(*,*)"		Author: Clara Estela Jiménez Tejero	"
write(*,*)"		email: ejimenez@icm.csic.es 		"
write(*,*)"		Barcelona Center for Subsurface Imaging "
write(*,*)"		Instituto de Ciencias Marinas (ICM-CSIC)"
write(*,*)
write(*,*)
write(*,*)"	                     |				"
write(*,*)"	                     |				"
write(*,*)"	            |        |				"
write(*,*)"	          |-|-|      |				"
write(*,*)"	            |        |				"
write(*,*)"	            | {O}    |				"
write(*,*)"	            '--|     |				"
write(*,*)"	              .|]_   |				"
write(*,*)"	        _.-=.' |     |				"	
write(*,*)"	       |    |  |]_   |				"
write(*,*)"	       |_.-='  |   __|__				"
write(*,*)"	        _.-='  |\   /|\				"
write(*,*)"	       |    |  |-'-'-'-'-.				"
write(*,*)"	       |_.-='  '========='				"
write(*,*)"	            `   |     |				"
write(*,*)"	             `. |    / \				"
write(*,*)"	               ||   /   \____.--=''''==--.._		"
write(*,*)"	               ||_.'--=='    |__  __  __  _.'	"
write(*,*)"	               ||  |    |    |\ ||  ||  || |                        ___	"	
write(*,*)"	  ____         ||__|____|____| \||__||__||_/    __________________/|   |	"
write(*,*)"	 |    |______  |===.---. .---.========''''=-._ |     |     |     / |   |	"
write(*,*)"	 |    ||     |\| |||   | |   |      '===' ||  \|_____|_____|____/__|___|	"
write(*,*)"	 |-.._||_____|_\___'---' '---'______....---===''======//=//////========|	"
write(*,*)"	 |--------------\------------------/-----------------//-//////---------/	"
write(*,*)"	 |               \                /                 // //////         /	"
write(*,*)"	 |                \______________/                 // //////         /	"
write(*,*)"	 |                                        _____===//=//////=========/	"
write(*,*)"	 |=================================================================/		"
write(*,*)"	  -----------------------------------------------------------------		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)
write(*,*)	
write(*,*)"				         ______		"
write(*,*)"				        /      \ 	"	
write(*,*)"				       /        \ 	"
write(*,*)"				       |        |	"
write(*,*)"				    )  o        o   (	"
write(*,*)"				   (    \      /    )	"
write(*,*)"				  _ \___/||||||\___/ _	"
write(*,*)"				   \____/ |||| \____/ `	"
write(*,*)"				   ,-.___/ || \__,-._	"
write(*,*)"				  /    ___/  \__	"
write(*,*)"				     _/         `--	"
write(*,*)
write(*,*)
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)
write(*,*)


endif

if(i.eq.2)	then

write(*,*)
write(*,*)
write(*,*)'Fix your inputs please'
write(*,*)'                                |     |'
write(*,*)'                                \\_V_//'
write(*,*)'                                \/=|=\/'
write(*,*)'                                 [=v=]'
write(*,*)'                               __\___/_____'
write(*,*)'                              /..[  _____  ]'
write(*,*)'                             /_  [ [  M /] ]'
write(*,*)'                            /../.[ [ M /@] ]'
write(*,*)'                           <-->[_[ [M /@/] ]'
write(*,*)'                          /../ [.[ [ /@/ ] ]'
write(*,*)'     _________________]\ /__/  [_[ [/@/ C] ]'
write(*,*)'    <_________________>>0---]  [=\ \@/ C / /'
write(*,*)'       ___      ___   ]/000o   /__\ \ C / /'
write(*,*)'          \    /              /....\ \_/ /'
write(*,*)'       ....\||/....           [___/=\___/'
write(*,*)'      .    .  .    .          [...] [...]'
write(*,*)'     .      ..      .         [___/ \___]'
write(*,*)'     .    0 .. 0    .         <---> <--->'
write(*,*)'  /\/\.    .  .    ./\/\      [..]   [..]'
write(*,*)' / / / .../|  |\... \ \ \    _[__]   [__]_'
write(*,*)'/ / /       \/       \ \ \  [____>   <____]'
write(*,*)
write(*,*)

endif

if(i.eq.3)	then

	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)"	 _______  ___   __    _  ___   _______  __   __  _______  ______  	"
	write(*,*)"	|       ||   | |  |  | ||   | |       ||  | |  ||       ||      | 	"
	write(*,*)"	|    ___||   | |   |_| ||   | |  _____||  |_|  ||    ___||  _    |	"
	write(*,*)"	|   |___ |   | |       ||   | | |_____ |       ||   |___ | | |   |	"
	write(*,*)"	|    ___||   | |  _    ||   | |_____  ||       ||    ___|| |_|   |	"
	write(*,*)"	|   |    |   | | | |   ||   |  _____| ||   _   ||   |___ |       |	"
	write(*,*)"	|___|    |___| |_|  |__||___| |_______||__| |__||_______||______|	"
	write(*,*)
	write(*,*)

endif

if(i.eq.4)	then

write(*,*)""
write(*,*)""
write(*,*)"*** WARNING: numtasks do not need to be greater than variable NumOBS"
write(*,*)"***          You are wasting energy. Please, next time you run a job, set in you MPI execution line, numtasks=",NumOBS
write(*,*)""

endif

END SUBROUTINE ascii_art


subroutine freq_filter(Data_trace,nt,dt,Num,typef,f1,f2)
implicit none
        
        integer :: nt,Num,nt2,zero_pad,ind,typef
        real :: dt,fe,var,f1,f2
        real :: Data_trace(nt,Num)
        real,allocatable :: tmp(:)
        integer :: j,i,k 
        logical :: isnotCero            ! Variable para almacenar si el vector es cero
        
        zero_pad=24
        nt2=zero_pad*nt
        allocate(tmp(nt2))
        fe=1/dt
        
        do i=1,Num
        
                tmp=0.
                tmp((zero_pad-1)*nt+1:nt2)=Data_trace(:,i)
        
                call is_not_cero(nt,Data_trace(:,i),isnotCero)
        
                if(isnotCero)   then
                        call filters_function1D(tmp,nt2,fe,1,typef,f1,f2)
                        Data_trace(:,i)=tmp(nt2-nt+1:nt2)
                endif
        
        enddo
        
        
        deallocate(tmp)  
        
return          
end
