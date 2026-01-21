!###########################################################
!###########################################################
!###########################################################
!##  tomo_grids -v Vel_AT01.tomo -in.km -out.m -dx 12.5 -dz 12.5 -out.su
!###########################################################
!###########################################################
!###########################################################

program TOMO_GRIDS_MAIN

implicit none

integer			:: i,j,k,icount,nr,art,err,k1,k2,arg_pos
integer			:: nx,nz,nx2,nz2,nt,nbat,nzmodel,nxtotal,nztotal
integer			:: numDec_x,numDec_t,numDec_z,numDec_v,numDec_d
integer			:: irefl,numrefl
integer			:: save_v,save_xzv,save_su
integer			:: units_in,units_out

real			:: xx,zz,vwater,vair
real			:: resox,resoz,resot
real			:: time_max,timer
real			:: minb,minr,res_test,misfit,av
real			:: sizefile

real, allocatable	:: x(:),z(:),xr(:),xdws(:),zdws(:,:)
real, allocatable	:: x2(:),z2(:),time(:),lim1(:),lim2(:)
real, allocatable	:: time_total(:,:),vel_total(:,:),dws_total(:,:),x_total(:),z_total(:)

real, allocatable	:: vel(:,:),velz(:,:),velxz(:,:),vel_time(:,:)
real, allocatable	:: refl(:),reflxz(:),refl_time(:)
real, allocatable	:: bat(:),batxz(:),bat_time(:)
real, allocatable	:: dws(:,:),dwsz(:,:),dwsxz(:,:),dws_time(:,:)

character(len=1000) 	:: arg,flag,str
character(len=1000)	:: vel_file,dws_file
character(len=1000)	:: dx_str,dz_str,dt_str
character(len=1000)	:: file_out
character(len=1000)	:: format_xtv,format_xtd,format_xt,format_xz,format_xzv,format_xzd,format_xr,format_xb
character(len=50)	:: format_x,format_z,format_t,format_v,format_d
character(len=1000), allocatable	:: refl_file(:)

logical 		:: vel_exist,refl_exist,dws_exist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ COMMAND LINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  err=0
  vwater=1.5
  numDec_x=5;numDec_z=3;numDec_t=3;numDec_v=3;numDec_d=2;
  resox=0.;resoz=0;resot=0;numrefl=0;irefl=0;
  save_v=0;save_xzv=0;save_su=0;art=1;units_in=1;units_out=1;

  icount=command_argument_count()

  do i = 1, icount	!command_argument_count()

    call get_command_argument(i, value=arg)
    call get_command_argument(i + 1, value=flag)

      if (arg == '-v') then
	if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
	    flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
	endif
	vel_file = ADJUSTL(flag)
    endif
      if (arg == '-d') then
	if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
	    flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
	endif
	dws_file = ADJUSTL(flag)
    endif
      if (arg == '-dx') then
     	 read(flag, *) resox
     	 dx_str=flag
    endif
      if (arg == '-dz') then
      	read(flag, *) resoz
      	dz_str=flag
    endif
      if (arg == '-dt') then
      	read(flag, *) resot
      	dt_str=flag
    endif
      if (arg == '-rx') then
      	read(flag, *) numDec_x
    endif
      if (arg == '-rz') then
      	read(flag, *) numDec_z
    endif
      if (arg == '-rt') then
      	read(flag, *) numDec_t
    endif
      if (arg == '-rv') then
      	read(flag, *) numDec_v
    endif
      if (arg == '-rd') then
      	read(flag, *) numDec_d
    endif
     if (arg == '-r') then
      	irefl=irefl+1
   endif
     if (arg == '-out.v') then
      	save_v=1
   endif
     if (arg == '-out.xzv') then
      	save_xzv=1
   endif
     if (arg == '-out.su') then
      	save_su=1
   endif
     if (arg == '-in.km') then
      	units_in=1
   endif
     if (arg == '-in.m') then
      	units_in=1000
   endif
     if (arg == '-out.km') then
      	units_out=1
   endif
     if (arg == '-out.m') then
      	units_out=1000
   endif

  end do

  if(irefl.gt.0) then
	numrefl=irefl
	allocate(refl_file(numrefl))
  endif

  irefl=0
  do i = 1, icount	!command_argument_count()
    call get_command_argument(i, value=arg)
    call get_command_argument(i + 1, value=flag)
      if (arg == '-r') then
        irefl=irefl+1
	if (INDEX(flag, '"') == 1 .and. INDEX(flag, '"', BACK=.TRUE.) == LEN_TRIM(flag)) then
	    flag = ADJUSTL(flag(2:LEN_TRIM(flag)-1))
	endif
	refl_file(irefl) = ADJUSTL(flag)
    endif
  end do

  INQUIRE(file=trim(adjustl(vel_file)),EXIST=vel_exist)
  if(.NOT. vel_exist)     then
        write(*,*)'ERROR: -v flag is empty or velocity file is not found'
        stop
  endif
  if(resox.eq.0)	then
        write(*,*)'ERROR: -dx flag must not be empty'
        stop
  endif
  if(resoz.eq.0)	then
        write(*,*)'ERROR: -dz flag must not be empty'
        stop
  endif
  if(resot.eq.0)	then
        write(*,*)'WARNING: -dt flag is empty'
!        stop
  endif

  write(*,*)"number of reflectors",numrefl
  INQUIRE(file=trim(adjustl(dws_file)),EXIST=dws_exist)

if (numrefl .gt. 0) then
	do i=1,numrefl
		INQUIRE(file=trim(adjustl(refl_file(i))),EXIST=refl_exist)
		if(refl_exist)write(*,*)"Reflector number: ",i," exist: ", trim(adjustl(refl_file(i)))
	enddo
	refl_exist = .true.
else
	refl_exist = .false.
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END READ COMMAND LINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ INPUT VALUES: X, bathymetry, Z, Velocity model, and reflector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! READ Vp FILE

  open(10, file=trim(adjustl(vel_file)), status='old')
  read(10,*) nx,nz,vwater,vair
  
  allocate(x(nx),bat(nx),z(nz),vel(nz,nx))
  read(10,*) (x(i), i=1,nx)
  read(10,*) (bat(i), i=1,nx)
  read(10,*) (z(k), k=1,nz)
  do i=1,nx
	read(10,*) (vel(k,i), k=1,nz)
  enddo
  close(10)

  if(units_in.eq.1.and.units_out.eq.1000)	then
	  bat=bat*1000
	  vel=vel*1000
	  x=x*1000
	  z=z*1000
	  vwater=vwater*1000
  endif
  if(units_in.eq.1000.and.units_out.eq.1)	then
	  bat=bat/1000
	  vel=vel/1000
	  x=x/1000
	  z=z/1000
	  vwater=vwater/1000
  endif

  nx2=1+(x(nx)-x(1))/resox
  allocate(x2(nx2))
  do i=1,nx2
	  x2(i)=(i-1)*resox
  enddo

  nz2=1+(z(nz)-z(1))/resoz   
  allocate(z2(nz2))
  do i=1,nz2
          z2(i)=(i-1)*resoz
  enddo

!! Checkeo resolucion batimetria y consistencia con dz
  minb=10000
  do i=2,nx
	res_test=abs(bat(i)-bat(i-1))
	if(res_test.gt.0.and.res_test.lt.minb)minb=res_test;
  enddo
  
!! READ DWS FILE
if(dws_exist)	then

	  allocate(dws(nz,nx),xdws(nx),zdws(nz,nx))
	  open(10, file=trim(adjustl(dws_file)), status='old')
	  do i=1,nx
	  do k=1,nz
	  	read(10,*) xdws(i),zdws(k,i),dws(k,i)
	  enddo
	  enddo
	  close(10)

endif

  write(*,*)
  write(*,*)
  write(*,*)' **********   YOUR INPUT ARGUMENTS *****************'
  write(*,*)' Name of velocity model: ',trim(adjustl(vel_file))
  if(refl_exist)write(*,*)' Number of Reflectors: ',numrefl
  if(dws_exist) write(*,*)' Name of dws model: ',trim(adjustl(dws_file))
  write(*,*)' dx value is:', real(resox)
  write(*,*)' dz value is:', real(resoz)
  write(*,*)' dt value is:', real(resot)
  if(nx2.eq.nx)write(*,*)'No interpolation in X, nx2=nx'
  write(*,*)' Bathymetry minimum resolution in Z is: ',minb!,' km'
  write(*,*)' **********************************************'
  write(*,*)
  print *, achar(27)//'[34m dz must not differ too much from the bathymetry resolution'//achar(27)//'[0m.'
  print *, achar(27)//'[34m dz is', abs(resoz/minb), 'times bigger than bathymetry resolution'//achar(27)//'[0m.'
  if(abs(resoz/minb).gt.50)	then
	print *, achar(27)//'[31m **** WARNING **** '//achar(27)//'[0m.'
	print *, achar(27)//'[31m You should consider using a smaller value for dz:'//achar(27)//'[0m.'
	print *, achar(27)//'[31m', ' For example, use dz value, from: ', minb, ' to ', 20*minb, achar(27)//'[0m.'
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   FORMATS and DECIMALS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(format_x,*) numDec_x
write(format_z,*) numDec_z
write(format_t,*) numDec_t
write(format_v,*) numDec_v
write(format_d,*) numDec_d

format_x="F10." // trim(adjustl(format_x));
format_z="F10." // trim(adjustl(format_z));
format_t="F10." // trim(adjustl(format_t)); 
format_v="F10." // trim(adjustl(format_v));
format_d="F10." // trim(adjustl(format_d));

format_xz= "(" //trim(adjustl(format_x))//","//trim(adjustl(format_z)) // ")" !!x,t
format_xt= "(" //trim(adjustl(format_x))//","//trim(adjustl(format_t)) // ")" !!x,t

format_xzv= "(" //trim(adjustl(format_x))//","//trim(adjustl(format_z)) // "," // trim(adjustl(format_v)) // ")" !!x,t,v
format_xzd= "(" //trim(adjustl(format_x))//","//trim(adjustl(format_z)) // "," // trim(adjustl(format_d)) // ")" !!x,t,dws

format_xtv= "(" //trim(adjustl(format_x))//","//trim(adjustl(format_t)) // "," // trim(adjustl(format_v)) // ")" !!x,t,v
format_xtd= "(" //trim(adjustl(format_x))//","//trim(adjustl(format_t)) // "," // trim(adjustl(format_d)) // ")" !!x,t,dws

!write(*,*)format_bat,format_refl
!write(*,*)"xb,xr: ",format_xb,format_xr
!write(*,*)"xz,xt: ",format_xz,format_xt
!write(*,*)"xzv,xtv: ",format_xzv,format_xtv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   INPUT FILES		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,*)"WRITE velocity and Bathymetry in separate files"
  write(*,*)

  file_out= "velDepth.inputfile"
  open(12, file=file_out, status='unknown')
  do i=1,nx
  do k=1,nz
        write(12, format_xzv ) x(i),bat(i)+z(k),vel(k,i)
  enddo
  enddo
  close(12)

  file_out= "batDepth.inputfile"
  open(12, file=file_out, status='unknown')
  do i=1,nx
	write(12,format_xz) x(i),bat(i)
  enddo  
  close(12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   INTRPOLATION IN Z AND X  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,*)' ****************************************************************'
  write(*,*)' 1. Interpolation with new dx and dz: (nx,dx) and (nz,dz)= ', nx2,resox, ' and ' ,nz2,resoz
  write(*,*)' ****************************************************************'  

  allocate(velz(nz2,nx))
  allocate(velxz(nz2,nx2))

  do i=1,nx
       call interpolation_lineal1D(z,vel(:,i),nz,z2,velz(:,i),nz2) 
  enddo

  if(nx2.ne.nx)	then
	  do i=1,nz2
	       call interpolation_lineal1D(x,velz(i,:),nx,x2,velxz(i,:),nx2)   
	  enddo

  else
	velxz=velz
  endif

  allocate(batxz(nx2))
  if(nx2.ne.nx)	then
	call interpolation_lineal1D(x,bat,nx,x2,batxz,nx2)
  else
	batxz=bat
  endif

  if(dws_exist)	then

	  allocate(dwsz(nz2,nx))
	  allocate(dwsxz(nz2,nx2))
	  do i=1,nx
       		call interpolation_lineal1D(z,dws(:,i),nz,z2,dwsz(:,i),nz2) 
	  enddo

	  if(nx2.ne.nx) then
		  do i=1,nz2
	       		call interpolation_lineal1D(x,dwsz(i,:),nx,x2,dwsxz(i,:),nx2)   
		  enddo

	  else
		dwsxz=dwsz
	  endif

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!    GRIDS     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Genero el modelo de velocidad con agua incluido
  nbat=1+ceiling(maxval(batxz)/resoz)
  nztotal=nz2+nbat-1
  nxtotal=nx2

  allocate(vel_total(nztotal,nxtotal))
  vel_total=0.
  allocate(x_total(nxtotal))
  x_total=x2

  allocate(z_total(nztotal))
  do k=1,nztotal
          z_total(k)=(k-1)*resoz
  enddo

  do i=1,nxtotal

	  k1=1+ceiling(batxz(i)/resoz)
	  k2=k1+nz2

	  do k=1,k1-1
		vel_total(k,i)=vwater
	  enddo

	  vel_total(k1,i)=velxz(1,i)	!!en el seafloor
	
	  do k=1+k1,k2-1	
		vel_total(k,i)=velxz(k-k1+1,i)
	  enddo

	  do k=k2,nztotal	
		vel_total(k,i)=vel_total(k2-1,i)
	  enddo

  enddo

  if(dws_exist)	then

	  allocate(dws_total(nztotal,nxtotal))
	  dws_total=0

	  do i=1,nxtotal

		  k1=1+ceiling(batxz(i)/resoz)
		  k2=k1+nz2

		  dws_total(k1,i)=dwsxz(1,i)	!!en el seafloor
	
		  do k=1+k1,k2-1			
			dws_total(k,i)=dwsxz(k-k1+1,i)
		  enddo

		  do k=k2,nztotal	
			dws_total(k,i)=dws_total(k2-1,i)
		  enddo
	
	  enddo
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 		writting outputs with depth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,*)"WRITE OUTPUTS WITH DEPTH: ",nxtotal,nztotal

  sizefile=nxtotal*nztotal*9*3/1024./1024./1024. 
  write(*,*)"Your velocity file in depth will occupy about ",sizefile," GB"
  if(dws_exist)write(*,*)"Your dws file in depth will also occupy ",sizefile," GB"
  sizefile=nxtotal*9*2/1024./1024./1024.

  file_out="Bat_dx" //trim(adjustl(dx_str))// "_dz" //trim(adjustl(dz_str))//".out"
  open(12, file=file_out, status='unknown')
  do i=1,nxtotal
	write(12, format_xz ) x_total(i),batxz(i)
  enddo  
  close(12)

  if(save_v.eq.1)	then

	  file_out= "Vel_dx" //trim(adjustl(dx_str))// "_dz" //trim(adjustl(dz_str))//".v"
	  open(12, file=file_out, status='unknown')
	  do k=1,nztotal
	        write(12,'(*(f0.0,2x))') (vel_total(k,i),i=1,nxtotal)
	  enddo
	  close(12)

	  if(dws_exist)	then
	
		  file_out= "DWS_dx" //trim(adjustl(dx_str))// "_dz" //trim(adjustl(dz_str))//".dws"
		  open(12, file=file_out, status='unknown')
		  do k=1,nztotal
		        write(12,'(*(f0.0,2x))') (dws_total(k,i),i=1,nxtotal)
		  enddo
		  close(12)
	
	  endif

  endif

  if(save_xzv.eq.1)	then

	  file_out= "Vel_dx" //trim(adjustl(dx_str))// "_dz" //trim(adjustl(dz_str))//".xzv"
	  open(12, file=file_out, status='unknown')
	  do i=1,nxtotal
	  do k=1,nztotal
	        write(12, format_xzv ) x_total(i),z_total(k),vel_total(k,i)
	  enddo
	  enddo

	  close(12)

	  if(dws_exist)	then
	
		  file_out= "DWS_dx" //trim(adjustl(dx_str))// "_dz" //trim(adjustl(dz_str))//".xzdws"
		  open(12, file=file_out, status='unknown')
		  do i=1,nxtotal
		  do k=1,nztotal
		        write(12, format_xzd ) x_total(i),z_total(k),dws_total(k,i)
		  enddo
		  enddo
		  close(12)
	
	  endif

  endif

  if(save_su.eq.1)	then

	  file_out= "Vel_dx" //trim(adjustl(dx_str))// "_dz" //trim(adjustl(dz_str))//".su"
	  call write_su(nxtotal,nztotal,resox,resoz,file_out,vel_total,1)

	  if(dws_exist)	then
		  file_out= "DWS_dx" //trim(adjustl(dx_str))// "_dz" //trim(adjustl(dz_str))//".su"
		  call write_su(nxtotal,nztotal,resox,resoz,file_out,dws_total,1)
	  endif

  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   TIME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(time_total(nztotal,nxtotal))
  do i=1,nxtotal
  	  time_total(1,i)=0.
	  do k=2,nztotal
	  	if(abs(vel_total(k,i)).gt.0) then
			time_total(k,i)=time_total(k-1,i)+abs(resoz*2./vel_total(k,i))		!!seconds
			if(i.eq.1)write(111,*)k,resoz,vel_total(k,i),time_total(k,i)
		endif

	  enddo
  enddo

if(resot.ne.0)	then

!!!!	INTERPOLACION TIME

  time_max=maxval(abs(time_total))
  nt=1+ceiling(time_max/resot)
  allocate(time(nt))
  do i=1,nt
	  time(i)=(i-1)*resot
  enddo

!!!!	INTERPOLACION VEL

  write(*,*)
  write(*,*)' ****************************************************************'
  write(*,*)' 2. Interpolation with time: (time_max,nt,dt)= ', time_max,nt,resot
  write(*,*)' ****************************************************************'

  allocate(vel_time(nt,nxtotal))
  allocate(bat_time(nxtotal))

  do i=1,nxtotal
	call interpolation_lineal1D(time_total(:,i),vel_total(:,i),nztotal,time,vel_time(:,i),nt)
  enddo

  if(dws_exist)	then

!! INTERPOLATION DWS
	  allocate(dws_time(nt,nxtotal))
	  do i=1,nx2
		call interpolation_lineal1D(time_total(:,i),dws_total(:,i),nztotal,time,dws_time(:,i),nt)
	  enddo
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   OUTPUT TWT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,*)"WRITE OUTPUTS WITH TWT"

  sizefile=nxtotal*nt*(8+1)*3/1024./1024./1024. 
  write(*,*)"Your velocity file in TWT will occupy about ",sizefile," GB"
  if(dws_exist)write(*,*)"Your dws file in TWT will also occupy ",sizefile," GB"

  file_out="Bat_twt_exact_dx" //trim(adjustl(dx_str))// ".out"
  open(12, file=file_out, status='unknown')
  do i=1,nxtotal
	write(12,format_xt) x_total(i),abs(batxz(i)*2/vwater)
  enddo  
  close(12)

  file_out="Bat_twt_dx" //trim(adjustl(dx_str))// "_dz"  //trim(adjustl(dz_str))
  file_out= trim(adjustl(file_out))//"_dt" //trim(adjustl(dt_str))//".out"
  open(12, file=file_out, status='unknown')
  do i=1,nxtotal
	 k=1+ceiling(batxz(i)/resoz)
	write(12,format_xt) x_total(i),time_total(k,i)
  enddo  
  close(12)

  if(save_v.eq.1)	then

	  file_out="Vel_dx" //trim(adjustl(dx_str))// "_dz"  //trim(adjustl(dz_str))
	  file_out= trim(adjustl(file_out))// "_dt" //trim(adjustl(dt_str))//".v"

	  open(12, file=file_out, status='unknown')
	  do k=1,nt
	        write(12,'(*(f0.0,2x))') (vel_time(k,i), i=1,nxtotal)
	  enddo

	  close(12)

	  if(dws_exist)	then
	
		  file_out="DWS_dx" //trim(adjustl(dx_str))// "_dz"  //trim(adjustl(dz_str))
		  file_out= trim(adjustl(file_out))// "_dt" //trim(adjustl(dt_str))//".dws"
		  open(12, file=file_out, status='unknown')
		  do k=1,nt
		        write(12,'(*(f0.0,2x))') (dws_time(k,i), i=1,nxtotal)
		  enddo
		  close(12)
	
	  endif

  endif

  if(save_xzv.eq.1)	then

	  file_out="Vel_dx" //trim(adjustl(dx_str))// "_dz"  //trim(adjustl(dz_str))
	  file_out= trim(adjustl(file_out))// "_dt" //trim(adjustl(dt_str))//".xtv"
	  open(12, file=file_out, status='unknown')
	  do i=1,nxtotal
	  do k=1,nt
	       write(12, format_xtv ) x_total(i),time(k),vel_time(k,i)
	  enddo
	  enddo
	  close(12)

	  if(dws_exist)	then

		  file_out="DWS_dx" //trim(adjustl(dx_str))// "_dz"  //trim(adjustl(dz_str))
		  file_out= trim(adjustl(file_out))// "_dt" //trim(adjustl(dt_str))//".xtdws"
		  open(12, file=file_out, status='unknown')
	  	  do i=1,nxtotal
		  do k=1,nt
		        write(12, format_xtd) x_total(i),time(k),dws_time(k,i)
		  enddo
		  enddo
		  close(12)
	  endif

  endif


  if(save_su.eq.1)	then

	  write(*,*)"Writting vel twt ..."
	  file_out="Vel_dx" //trim(adjustl(dx_str))// "_dz"  //trim(adjustl(dz_str))
	  file_out= trim(adjustl(file_out))// "_dt" //trim(adjustl(dt_str))//".su"
	  call write_su(nxtotal,nt,resox,resot,file_out,vel_time,2)


	  write(*,*)"Writting dws twt ..."
	  file_out="DWS_dx" //trim(adjustl(dx_str))// "_dz"  //trim(adjustl(dz_str))
	  file_out= trim(adjustl(file_out))// "_dt" //trim(adjustl(dt_str))//".su"
	  call write_su(nxtotal,nt,resox,resot,file_out,dws_time,2)


  endif


endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   REFLECTORS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(refl_exist)	then

  write(*,*)
  IF(NUMREFL.EQ.1)write(*,*)"INTERPOLATION OF REFLECTOR "
  IF(NUMREFL.GT.1)write(*,*)"INTERPOLATION OF REFLECTORS"

do irefl=1,numrefl

	  !! READ EACH REFLECTOR FILE
	  write(str,*)irefl

	  open(10, file=trim(adjustl(refl_file(irefl))), status='old')
	  nr=0
	  do
	  	read(10,*, END=20)
		nr = nr+1
	  enddo
	20 close(10)

	  allocate(xr(nr),refl(nr))
	  open(10, file=trim(adjustl(refl_file(irefl))), status='old')
	  do i=1,nr
	  	read(10,*) xr(i),refl(i)
	  enddo
	  close(10)

	  minr=10000
	  do i=2,nr
		res_test=abs(refl(i)-refl(i-1))
		if(res_test.gt.0.and.res_test.lt.minr)minr=res_test;
	  enddo

	  allocate(reflxz(nx2))
          if(nx2.ne.nr) then
		 call interpolation_lineal1D(xr,refl,nr,x2,reflxz,nx2)
	  else
		reflxz=refl
	  endif

	  file_out="reflDepth_dx" //trim(adjustl(dx_str)) // "_dz" // trim(adjustl(dz_str))
	  if(numrefl.eq.1)	then
		  file_out= trim(adjustl(file_out))// "_dt"//trim(adjustl(dt_str))//".out"
	  endif
	  if(numrefl.gt.1)	then
		  file_out= trim(adjustl(file_out))// "_dt"//trim(adjustl(dt_str))
		  file_out= trim(adjustl(file_out))// "_refl"//trim(adjustl(str))//".out"
	  endif

	  open(12, file=file_out, status='unknown')
	  do i=1,nxtotal
		write(12,format_xz) x_total(i),reflxz(i)
	  enddo  
	  close(12)

	if(resot.ne.0)	then 

!!!!	INTERPOLACION REFLECTOR
	  allocate(refl_time(nxtotal))
	  do i=1,nxtotal
	        k=1+ceiling(reflxz(i)/resoz)
	        refl_time(i)=time_total(k,i)
	  enddo

	  file_out="reflTWT_dx" //trim(adjustl(dx_str)) // "_dz" // trim(adjustl(dz_str))
	  if(numrefl.eq.1)	then
		  file_out= trim(adjustl(file_out))// "_dt"//trim(adjustl(dt_str))//".out"
	  endif
	  if(numrefl.gt.1)	then
		  file_out= trim(adjustl(file_out))// "_dt"//trim(adjustl(dt_str))
		  file_out= trim(adjustl(file_out))// "_refl"//trim(adjustl(str))//".out"
	  endif

	  open(12, file=file_out, status='unknown')
	  do i=1,nxtotal
		write(12,format_xt) x_total(i),refl_time(i)
	  enddo  
	  close(12)

	  deallocate(refl_time)

	endif

	  deallocate(xr,refl,reflxz)

enddo

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  write(*,*)
  write(*,*)" *********************************"
  write(*,*)" FINISHED INTERPOLATION FOR GRIDS:"
  write(*,*)" dx,nx: ",resox,nxtotal
  write(*,*)" dz,nz: ",resoz,nztotal
  if(resot.ne.0)write(*,*)" dt,nt: ",resot,nt
  write(*,*)" *********************************"

  if(resot.ne.0)	then

  av=0
  do i=1,nxtotal
	 k=1+ceiling(batxz(i)/resoz)
	 misfit=abs(abs(batxz(i)*2/vwater)-time_total(k,i))
	 av=av+misfit	 
  enddo
  av=av/nxtotal

  write(*,*)
  write(*,*)" The average misfit between exact bathymetry and its value in your grid is: ",av*1000," meters"
  write(*,*)
  if(av*1000.le.15)print *, achar(27)//'[32m The value dz is appropriate'//achar(27)//'[0m.'
  if(av*1000.le.15)print *, achar(27)//'[32m Your results are robust with the chosen depth and time grids'//achar(27)//'[0m.'
  art=1
  if(av*1000.gt.15)	then
	art=0
	print *, achar(27)//'[31m **** WARNING **** DO NOT TRUST YOUT GRID'//achar(27)//'[0m.'
	print *, achar(27)//'[31m The value dz is too big, use a smaller value'//achar(27)//'[0m.'
	print *, achar(27)//'[31m',' For example, use dz value, from: ', minb, ' to ', 20*minb, achar(27)//'[0m.'

  endif
  write(*,*)


  deallocate(x,z)
  deallocate(x2,z2)
  deallocate(vel_total,time_total,x_total,z_total)
  deallocate(vel,velz,velxz)
  deallocate(bat,batxz)
  if(dws_exist)  deallocate(xdws,zdws,dws,dwsz,dwsxz,dws_total)
  if(resot.ne.0) deallocate(time,vel_time,bat_time)
  if(dws_exist)  then
	if(resot.ne.0) deallocate(dws_time)
  endif
  endif

if(art.eq.0)	then

write(*,*)
write(*,*)"                     ^`.                     o"
write(*,*)"     ^_              \  \                  o  o"
write(*,*)"     \ \             {   \                 o"
write(*,*)"     {  \           /     `~~~--__"
write(*,*)"     {   \___----~~'              `~~-_     ______          _____"
write(*,*)"      \                         /// a  `~._(_||___)________/___"
write(*,*)"      / /~~~~-, ,__.    ,      ///  __,,,,)      o  ______/    \"
write(*,*)"      \/      \/    `~~~;   ,---~~-_`~= \ \------o-'            \"
write(*,*)"                       /   /            / /"
write(*,*)"                      '._.'           _/_/"
write(*,*)"                                      ';|\"
write(*,*)		

endif

if(art.eq.1)	then

write(*,*)
write(*,*)"                                   __"
write(*,*)"                               _.-~  )"
write(*,*)"                    _..--~~~~,'   ,-/     _"
write(*,*)"                 .-'. . . .'   ,-','    ,' )"
write(*,*)"               ,'. . . _   ,--~,-'__..-'  ,'"
write(*,*)"             ,'. . .  (@)' ---~~~~      ,'"
write(*,*)"            /. . . . '~~             ,-'"
write(*,*)"           /. . . . .             ,-'"
write(*,*)"          ; . . . .  - .        ,'"
write(*,*)"         : . . . .       _     /"
write(*,*)"        . . . . .          `-.:"
write(*,*)"       . . . ./  - .          )"
write(*,*)"      .  . . |  _____..---.._/ _____"
write(*,*)"~---~~~~----~~~~             ~~"
write(*,*)

endif


CONTAINS


end program TOMO_GRIDS_MAIN

subroutine interpolation_lineal1D(x1,y1,n1,x2,y2,n2)

implicit none

        integer i,j,n1,n2
	real :: m
	real :: a1,a2,b1,b2,a3
        real :: x1(n1),y1(n1)
        real :: x2(n2),y2(n2)

	y2=0;

        do i=1,n1-1

                a1=x1(i+1)
		a2=x1(i)
                b1=y1(i+1)
		b2=y1(i)

                m=(b1-b2)/(a1-a2)
		
		if(abs(b1).gt.0.and.abs(b2).gt.0)	then

                do j=1,n2
			a3=x2(j)
	                if (a3.ge.a2.and.a3.le.a1) then
                               y2(j)=b2+(a3-a2)*m
                        end if
                enddo

		endif

        enddo


end subroutine

subroutine write_su(ni,nj,di,dj,file_name,model,type)   

implicit none
integer :: i,j,k,ni,nj,type
integer :: save_su,save_v,save_xzv
real    :: xi,xj
real    :: di,dj
real    :: model(nj,ni)

character(len=1000) :: Str,Str_nj,Str_di,Str_dj,command
character(len=1000) :: file_name,file_name2
character(len=1000) :: file_temp,file_temp_filt,file_temp_filt2
character(len=1000) :: file_bin,file_su,file_su_filt


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

                file_su= trim(adjustl(file_name))
                open(unit=16,file=file_su,status='unknown')

                write(Str_nj,*) nj
                write(Str_di,*) di
                write(Str_dj,*) dj

                command= "a2b < " // trim(adjustl(file_temp)) // " n1=1 > " // trim(adjustl(file_bin))
                call system(command)

                if(type.eq.1)   then	!! profundidad
                command = "suaddhead < " // trim(adjustl(file_bin)) // " ns=" // trim(adjustl(Str_nj)) // " | " // &
                   "sushw key=d1 a=" // trim(adjustl(Str_dj)) // " | " // &
                   "sushw key=f1 a=0 | " // &
                   "sushw key=d2 a=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=f2 a=0 j=" // trim(adjustl(Str_di)) // " | " // &
                   "suchw key1=cdp key2=tracl a=-" // trim(adjustl(Str_di)) // " b=" // trim(adjustl(Str_di)) // " | " // &
                   "sushw key=trid a=200 > " // trim(adjustl(file_su))
                endif

                if(type.eq.2)   then	!! TWT
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

end subroutine write_su
