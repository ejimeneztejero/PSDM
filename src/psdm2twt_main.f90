!###########################################################
!###########################################################
!###########################################################
!##  psdm2twt -v.in -v.out -psdm.in -psdm.out -nx 3000 -nz 1000 -dz 12.5 -dt 0.004 -format.in 1 -format.out 1 
!###########################################################
!###########################################################
!###########################################################

program psdm2twt_main

implicit none

integer		:: i,j,k,icount
integer		:: nx, nz	! model size
integer		:: format_in, format_out
integer		:: in_su,in_v,in_xzv
integer		:: out_su,out_v,out_xzv
real		:: dz		! depth sampling (m)
real		:: dt		! TWT sampling (s)
real, allocatable	:: vel_model(:,:)
real, allocatable	:: psdm_model(:,:)

character(len=1000)     :: arg,flag,str
character(len=1000)	:: file_vel,file_psdm
character(len=1000)	:: out_vel,out_psdm

logical                 :: vel_exist,psdm_exist

format_in=0
format_out=0
in_su=0
in_v=0
in_xzv=0
out_su=0
out_v=0
out_xzv=0

file_vel='null'
file_psdm='null'
out_vel='null'
out_psdm='null'

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
      if (arg == '-dz') then
        read(flag, *) dz
    endif
      if (arg == '-dt') then
        read(flag, *) dt
    endif

      if (arg == '-v.in') then
        if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
            flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
        endif
        file_vel = ADJUSTL(flag)
    endif

      if (arg == '-psdm.in') then
        if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
            flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
        endif
        file_psdm = ADJUSTL(flag)
    endif

      if (arg == '-v.out') then
        if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
            flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
        endif
        out_vel = ADJUSTL(flag)
    endif

      if (arg == '-psdm.out') then
        if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
            flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
        endif
        out_psdm = ADJUSTL(flag)
    endif

      if (arg == '-format.in') then
         read(flag, *) format_in
    endif

      if (arg == '-format.out') then
         read(flag, *) format_out
    endif

enddo

  INQUIRE(file=trim(adjustl(file_vel)),EXIST=vel_exist)
  if(.NOT.vel_exist)     then
        write(*,*)'ERROR: -v.in flag is empty or velocity file is not found'
        stop
  endif

  if(out_vel.eq.'null')     then
        write(*,*)'ERROR: -v.out flag is empty'
        stop
  endif

  INQUIRE(file=trim(adjustl(file_psdm)),EXIST=psdm_exist)
  if(.NOT. psdm_exist)     then
        write(*,*)'WARNING: -psdm.in flag is empty or velocity file is not found'
  endif

  if(nx.eq.0)        then
        write(*,*)'ERROR: -nx flag must not be empty'
        stop
  endif
  if(nz.eq.0)        then
        write(*,*)'ERROR: -nz flag must not be empty'
        stop
  endif
  if(dz.eq.0)        then
        write(*,*)'ERROR: -dz flag must not be empty'
        stop
  endif
  if(dt.eq.0)        then
        write(*,*)'ERROR: -dt flag must not be empty'
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
if(format_in.eq.123)	then
	in_su=1
	in_xzv=1
	in_v=1
endif

if(format_out.eq.1)out_su=1
if(format_out.eq.2)out_xzv=1
if(format_out.eq.3)out_v=1
if(format_out.eq.123)	then
	out_su=1
	out_xzv=1
	out_v=1
endif

allocate(vel_model(nz,nx))
allocate(psdm_model(nz,nx))

if(vel_exist)	then
	call get_model_main(nx,nz,file_vel,in_v,in_xzv,in_su,vel_model)
endif

if(psdm_exist)	then
	call get_model_main(nx,nz,file_psdm,in_v,in_xzv,in_su,psdm_model)
endif

call psdm2twt_sub(nx,nz,dz,dt,out_vel,out_psdm,vel_model,psdm_exist,psdm_model,out_su,out_v,out_xzv)

deallocate(vel_model)
deallocate(psdm_model)

end program psdm2twt_main

subroutine psdm2twt_sub(nx,nz,dd,dt,file_vel,file_psdm,vel_total,psdm_exist,psdm_total,save_su,save_v,save_xzv)
    implicit none

    !-------------------------
    ! INPUTS
    !-------------------------
    integer, intent(in) :: nx, nz     ! model size
    integer, intent(in)	:: save_su,save_v,save_xzv
    real, intent(in) :: dd	! depth sampling (m)
    real, intent(in) :: dt	! desired TWT sampling (s)
    real, intent(in) :: vel_total(nz,nx)
    real, intent(in) :: psdm_total(nz,nx)
    character(len=1000), intent(in) :: file_vel,file_psdm
    logical, intent(in)	:: psdm_exist

    !-------------------------
    ! LOCAL ARRAYS
    !-------------------------
    real, allocatable :: time(:)
    real, allocatable :: time_total(:,:)
    real, allocatable :: vel_twt(:,:)
    real, allocatable :: psdm_twt(:,:)
    real :: time_max
    integer :: i, k, nt

    !**********************************************************************
    ! 1. Compute cumulative TWT at each depth: time_total(k,i)
    !**********************************************************************
    allocate(time_total(nz,nx))

    do i = 1, nx
        time_total(1,i) = 0.0
        do k = 2, nz
            if (abs(vel_total(k,i)) > 0.0) then
                time_total(k,i) = time_total(k-1,i) + 2.0*dd/vel_total(k,i)
            else
                time_total(k,i) = time_total(k-1,i)
            end if
        end do
    end do


    !**********************************************************************
    ! 2. Build uniform time axis from 0 to time_max
    !**********************************************************************
    time_max = maxval(abs(time_total))
    nt = 1 + ceiling(time_max/dt)

    allocate(time(nt))
    allocate(vel_twt(nt,nx))
    allocate(psdm_twt(nt,nx))
    do i = 1, nt
        time(i) = (i-1)*dt
    end do


    !**********************************************************************
    ! 3. Interpolate velocity from depth sampling to TWT sampling
    !**********************************************************************

    do i = 1, nx
        call linear_interpolation_1D(time_total(:,i),vel_total(:,i),nz,time,vel_twt(:,i),nt)
        call linear_interpolation_1D(time_total(:,i),psdm_total(:,i),nz,time,psdm_twt(:,i),nt)
    end do

    call write_vel_main(nx,nt,dd,dt,file_vel,vel_twt,save_su,save_v,save_xzv)
    if(psdm_exist)call write_psdm_main(nx,nt,dd,dt,file_psdm,psdm_twt,save_su,save_v,save_xzv)

    deallocate(time,time_total,vel_twt,psdm_twt)

end subroutine psdm2twt_sub

subroutine write_vel_main(ni,nj,di,dj,file_name,model,save_su,save_v,save_xzv)	!!paralelizar

implicit none
integer :: i,j,k,ni,nj
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

                command = "suaddhead < " // trim(adjustl(file_bin)) // " ns=" // trim(adjustl(Str_nj)) // " | " // &
                   "sushw key=d1 a=" // trim(adjustl(Str_dj)) // " | " // &
                   "sushw key=f1 a=0 | " // &
                   "sushw key=d2 a=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=f2 a=0 j=" // trim(adjustl(Str_di)) // " | " // &
                   "suchw key1=cdp key2=tracl a=-" // trim(adjustl(Str_di)) // " b=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=trid a=1 > " // trim(adjustl(file_su))

                call system(command)

                command="rm " // trim(adjustl(file_bin)) // " " // trim(adjustl(file_temp))
                call system(command)

                close(14)
                close(16)

        endif

!       write(*,*)"FINISHED WRITE"

end subroutine write_vel_main

subroutine write_psdm_main(ni,nj,di,dj,file_name,model,save_su,save_v,save_xzv)    !!paralelizar

implicit none
integer :: i,j,k,ni,nj,last
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
                        write(12,'(20000(e12.5,2x))') (model(j,i),i=1,ni)
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

                command = "suaddhead < " // trim(adjustl(file_bin)) // " ns=" // trim(adjustl(Str_nj)) // " | " // &
                   "sushw key=d1 a=" // trim(adjustl(Str_dj)) // " | " // &
                   "sushw key=f1 a=0 | " // &
                   "sushw key=d2 a=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=f2 a=0 j=" // trim(adjustl(Str_di)) // " | " // &
                   "suchw key1=cdp key2=tracl a=-" // trim(adjustl(Str_di)) // " b=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=trid a=1 > " // trim(adjustl(file_su))

                call system(command)

                command="rm " // trim(adjustl(file_bin)) // " " // trim(adjustl(file_temp))
                call system(command)


                close(14)
                close(16)

        endif

!       write(*,*)"FINISHED WRITE"

end subroutine write_psdm_main

subroutine linear_interpolation_1D(x1, y1, n1, x2, y2, n2)

    implicit none

    integer :: i, j, n1, n2
    real :: m
    real :: a1, a2, b1, b2, a3
    real :: x1(n1), y1(n1)
    real :: x2(n2), y2(n2)

    y2 = 0.0

    do j = 1, n2
        a3 = x2(j)

        ! Extrapolaci      n a la izquierda
        if (a3 <= x1(1)) then
            y2(j) = y1(1)

        ! Extrapolaci      n a la derecha
        elseif (a3 >= x1(n1)) then
            y2(j) = y1(n1)

        ! Interpolaci      n dentro del rango
        else
            do i = 1, n1 - 1
                a1 = x1(i+1)
                a2 = x1(i)
                b1 = y1(i+1)
                b2 = y1(i)

                if (a3 >= a2 .and. a3 <= a1 .and. abs(b1) > 0 .and. abs(b2) > 0) then
                    m = (b1 - b2) / (a1 - a2)
                    y2(j) = b2 + (a3 - a2) * m
                    exit
                endif
            enddo
        endif
    enddo
end subroutine


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
