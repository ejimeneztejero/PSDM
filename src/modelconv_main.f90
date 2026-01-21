program modelconv_main

implicit none

integer		:: i,j,k,icount
integer		:: nx,nz,nt,n2	! model size
integer		:: format_in, format_out
integer		:: in_su,in_v,in_xzv
integer		:: out_su,out_v,out_xzv
integer		:: type_model
real		:: d2
real		:: dx,dz		! depth sampling (m)
real		:: dt		! TWT sampling (s)
real, allocatable	:: model(:,:)

character(len=1000)     :: arg,flag,str
character(len=1000)	:: file_model
character(len=1000)	:: out_model

logical                 :: model_exist

nx=0
nz=0
nt=0

dx=0
dz=0
dt=0

format_in=0
format_out=0
in_su=0
in_v=0
in_xzv=0
out_su=0
out_v=0
out_xzv=0

file_model='null'
out_model='null'

  icount=command_argument_count()

  do i = 1, icount      !command_argument_count()

    call get_command_argument(i, value=arg)
    call get_command_argument(i + 1, value=flag)

      if (arg == '-nx') then
         read(flag, *) nx
    endif
      if (arg == '-nz') then
        read(flag, *) nz
    endif
      if (arg == '-dx') then
        read(flag, *) dx
    endif
      if (arg == '-dz') then
        read(flag, *) dz
    endif
      if (arg == '-nt') then
        read(flag, *) nt
    endif
      if (arg == '-dt') then
        read(flag, *) dt
    endif
      if (arg == '-model.in') then
        if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
            flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
        endif
        file_model = ADJUSTL(flag)
    endif
      if (arg == '-model.out') then
        if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
            flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
        endif
        out_model = ADJUSTL(flag)
    endif
      if (arg == '-format.in') then
         read(flag, *) format_in
    endif
      if (arg == '-format.out') then
         read(flag, *) format_out
    endif

enddo

  INQUIRE(file=trim(adjustl(file_model)),EXIST=model_exist)
  if(.NOT.model_exist)     then
        write(*,*)'ERROR: -model.in flag is empty or model file is not found'
        stop
  endif

  if(out_model.eq.'null')     then
        write(*,*)'WARNING: -model.out flag is empty. Model_out will be use as name'
	out_model="Model_out"
  endif

  if(nx.eq.0)        then
        write(*,*)'ERROR: -nx flag must not be empty'
        stop
  endif
  if(dx.eq.0)        then
        write(*,*)'ERROR: -dx flag must not be empty'
        stop
  endif
  if(nz.eq.0)        then
	type_model=2
  endif

  if(nz.ne.0.and.dz.eq.0)        then
        write(*,*)'ERROR: -dz flag is empty'
	stop
  endif

  if(nt.eq.0)        then
	type_model=1
  endif

  if(nt.ne.0.and.dt.eq.0)        then
        write(*,*)'ERROR: -dt is empty'
	stop
  endif

  if(format_in.eq.0)        then
        write(*,*)'ERROR: -format.in flag must not be empty'
        stop
  endif
  if(format_out.eq.0)        then
        write(*,*)'ERROR: -format.out flag must not be empty'
        stop
  endif

if(format_in.eq.1) in_su=1
if(format_in.eq.2) in_xzv=1
if(format_in.eq.3) in_v=1

if(format_out.eq.1)out_su=1
if(format_out.eq.2)out_xzv=1
if(format_out.eq.3)out_v=1
if(format_out.eq.123)	then
	out_su=1
	out_xzv=1
	out_v=1
endif

if(type_model.eq.1)	then
	n2=nz
	d2=dz
endif
if(type_model.eq.2)	then
	n2=nt
	d2=dt
endif

if(model_exist)	then

	allocate(model(n2,nx))

	call get_model_main(nx,n2,file_model,in_v,in_xzv,in_su,model)
	call write_model_main(nx,n2,dx,d2,out_model,model,out_su,out_v,out_xzv,type_model)

	deallocate(model)

endif


end program modelconv_main

subroutine write_model_main(ni,nj,di,dj,file_name,model,save_su,save_v,save_xzv,type)	!!paralelizar

implicit none
integer :: i,j,k,ni,nj,type
integer	:: save_su,save_v,save_xzv
real    :: xi,xj
real    :: di,dj
real    :: model(nj,ni)

character(len=1000) :: Str,Str_nj,Str_di,Str_dj,command
character(len=1000) :: file_name,file_name2
character(len=1000) :: file_temp,file_temp_filt,file_temp_filt2
character(len=1000) :: file_bin,file_su,file_su_filt

        if(save_v.ne.0) then

                file_name2= trim(adjustl(file_name)) // ".f"

                open(unit=12,file=file_name2,status='unknown')

                do j=1,nj
                        write(12,'(*(f0.0,2x))') (model(j,i),i=1,ni)
                enddo
                close(12)

        endif


        if(save_xzv.ne.0)       then

                file_name2= trim(adjustl(file_name)) // ".xzf"
                open(unit=12,file=file_name2,status='unknown')

                do i=1,ni
                        xi=(i-1)*di
                        do j=1,nj
                                xj=(j-1)*dj
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

                command= "a2b < " // trim(adjustl(file_temp)) // " n1=1 > " // trim(adjustl(file_bin))
                call system(command)

		if(type.eq.1)	then
                command = "suaddhead < " // trim(adjustl(file_bin)) // " ns=" // trim(adjustl(Str_nj)) // " | " // &
                   "sushw key=d1 a=" // trim(adjustl(Str_dj)) // " | " // &
                   "sushw key=f1 a=0 | " // &
                   "sushw key=d2 a=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=f2 a=0 j=" // trim(adjustl(Str_di)) // " | " // &
                   "suchw key1=cdp key2=tracl a=-" // trim(adjustl(Str_di)) // " b=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=trid a=200 > " // trim(adjustl(file_su))
		endif

		if(type.eq.2)	then
                command = "suaddhead < " // trim(adjustl(file_bin)) // " ns=" // trim(adjustl(Str_nj)) // " | " // &
                   "sushw key=d1 a=" // trim(adjustl(Str_dj)) // " | " // &
                   "sushw key=f1 a=0 | " // &
                   "sushw key=d2 a=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=f2 a=0 j=" // trim(adjustl(Str_di)) // " | " // &
                   "suchw key1=cdp key2=tracl a=-" // trim(adjustl(Str_di)) // " b=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=trid a=1 > " // trim(adjustl(file_su))
		endif

                call system(command)

                command="rm " // trim(adjustl(file_bin)) // " " // trim(adjustl(file_temp))
                call system(command)

                close(14)
                close(16)

        endif

!       write(*,*)"FINISHED WRITE"

end subroutine write_model_main

subroutine get_model_main(nx,nz,file_name,model_v,model_xzv,model_su,model)

  implicit none

  integer :: i,j,nx,nz
  integer :: model_v,model_xzv,model_su
  real  :: xi,xj
  real  :: model(nz,nx)
  character(len=500) :: file_name
  CHARACTER(len=50) :: access,form

  access = 'stream'
  form = 'unformatted'

!!      ARREGLAR

        if(model_v.eq.1)  then
                open(unit=10,file=file_name,status='old')
                do j=1,nz
                        read(10,*) (model(j,i),i=1,nx)
                enddo
        endif
        if(model_xzv.eq.1)        then
                open(unit=10,file=file_name,status='old')
                do i=1,nx
                        do j=1,nz
                                read(10,*) xi,xj,model(j,i)
                        enddo
                enddo
        endif
        if(model_su.eq.1) then
                open(10,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')
                call get_su_data_main(nx,nz,10,model)
        endif

        close(10)


end subroutine get_model_main

SUBROUTINE get_su_data_main(ni,ns,unit_su,mdata)

IMPLICIT NONE

INTEGER :: nh
INTEGER :: i,j,k,ni,ns,unit_su
INTEGER(4) :: pos_r
REAL    :: mdata(ns,ni)
REAL(4), ALLOCATABLE :: Data_tmp(:,:),sudata(:,:)

nh = 60

allocate(Data_tmp(ns,ni))
allocate(sudata(ns+nh,ni))
Data_tmp=0.;sudata=0.;

do i=1,ni
	pos_r = 1 + (i-1)*(nh+ns)*4
        READ(unit_su, pos=pos_r)  sudata(1:ns+nh,i)
enddo

do i=1,ni
	do j=1,ns
                Data_tmp(j,i) = sudata(j+nh,i)
        enddo
enddo

mdata=Data_tmp

deallocate(Data_tmp)
deallocate(sudata)

END SUBROUTINE get_su_data_main
