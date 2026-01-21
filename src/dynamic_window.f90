subroutine dynamic_window(rank,NumSou,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)

implicit none
integer :: i,j,isou,rank
integer :: nx,dir_boat
integer :: max_n,dmin,dmax
integer :: NumSou,NumRec
integer :: min_d,max_d,max_nx
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real    :: dd,dd_space
integer :: nxSou(NumSou),nxRec(NumRec,NumSou)

        ddx=1+floor(dd_space/dd) !!added model space in x-axis
        max_nx=-1000

        do isou=1,NumSou

                dmin=min(nxSou(isou),minval(nxRec(:,isou)))
                dmax=max(nxSou(isou),maxval(nxRec(:,isou)))

                x_i=dmin-ddx+1
                x_f=dmax+ddx

                if(x_i.le.0)x_i=1
                if(x_f.gt.nx)x_f=nx;

                nx_w=x_f-x_i+1

                if(nx_w.gt.max_nx)  then
                        max_nx=nx_w
                        min_d=dmin
                        max_d=dmax
                endif

        enddo

        x_i=min_d-ddx+1
        if(x_i.le.0)    then
                x_i=1
        endif

        x_f=max_d+ddx

        if(x_f.gt.nx)   then
                x_f=nx;
        endif

        nx_w=x_f-x_i+1

end subroutine

subroutine dynamic_window_solver_forward_ram(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        Data_source,dd,dt,ny,nx,nx_w,Vp_w,&
        dd_space,store_snap,nt_snap,nt,Fields,Data_synth)

implicit none

integer :: i,j,rank,iprop
integer :: nx,ny,nt,nt_snap,store_snap
integer :: NumRec
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real	:: dd_space
real	:: dd,dt
integer :: nxSou(1),nySou(1)
integer :: nxs(1),nys(1)
integer :: nxRec(NumRec),nyRec(NumRec)

real :: Data_source(nt)
real :: Data_synth(nt,NumRec) 
real :: Fields(ny,nx_w,nt_snap)
real :: Vp_w(ny,nx_w)

integer, allocatable :: nxr(:),nyr(:)

allocate(nxr(NumRec),nyr(NumRec))

        call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)

!        write(*,*)"dentro dynamic solver: ",nx_w,x_i,x_f

        nxs=nxSou-x_i+1; nys=nySou;
        nxr=nxRec-x_i+1; nyr=nyRec;

!!---------- PROPAGATOR

        call solver_ac_forward_ram(rank,1,NumRec,nxs,nys,nxr,nyr,&
        Data_source,dd,dt,ny,nx_w,Vp_w,&
        store_snap,nt_snap,nt,Fields,Data_synth)

    deallocate(nxr,nyr)

end subroutine

subroutine dynamic_window_solver_forward_disk(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        Data_source,dd,dt,ny,nx,nx_w,Vp_w,&
        shotid,unit_fields,folder_fields,&
        dd_space,store_snap,nt_snap,nt,Data_synth)

implicit none

integer :: unit_fields,shotid
integer :: i,j,rank,iprop
integer :: nx,ny,nt,nt_snap,store_snap
integer :: NumRec
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real	:: dd_space
real	:: dd,dt
integer :: nxSou(1),nySou(1)
integer :: nxs(1),nys(1)
integer :: nxRec(NumRec),nyRec(NumRec)

real :: Data_source(nt)
real :: Data_synth(nt,NumRec) 
real :: Vp_w(ny,nx_w)

integer, allocatable :: nxr(:),nyr(:)

character(len=500) :: folder_fields

allocate(nxr(NumRec),nyr(NumRec))

        call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)

        nxs=nxSou-(x_i-1);nys=nySou;
        nxr=nxRec-(x_i-1);nyr=nyRec;

!!---------- PROPAGATOR
        call solver_ac_forward_disk(rank,1,NumRec,nxs,nys,nxr,nyr,&
        Data_source,dd,dt,ny,nx_w,Vp_w,&
        shotid,unit_fields,folder_fields,&
        store_snap,nt_snap,nt,Data_synth)

    deallocate(nxr,nyr)

end subroutine

subroutine dynamic_window_solver_backward_ram(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        AdjSource,dd,dt,ny,nx,nx_w,Vp_w,&
        dd_space,store_snap,nt_snap,nt,Fields,PSDM_thread)

implicit none

integer :: i,j,rank,iprop
integer :: nx,ny,nt,nt_snap,store_snap
integer :: NumRec
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real	:: dd_space
real	:: dd,dt
integer :: nxSou(1),nySou(1)
integer :: nxs(1),nys(1)
integer :: nxRec(NumRec),nyRec(NumRec)

real :: AdjSource(nt,NumRec) 
real :: Fields(ny,nx_w,nt_snap)
real :: Vp_w(ny,nx_w)
real :: PSDM_thread(ny,nx_w)

integer, allocatable :: nxr(:),nyr(:)

allocate(nxr(NumRec),nyr(NumRec))

        call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)

        nxs=nxSou-x_i+1; nys=nySou;
        nxr=nxRec-x_i+1; nyr=nyRec;

!!---------- PROPAGATOR

        call solver_ac_backward_ram(rank,1,NumRec,nxs,nys,nxr,nyr,&
        AdjSource,dd,dt,ny,nx_w,Vp_w,&
        store_snap,nt_snap,nt,Fields,PSDM_thread)

    deallocate(nxr,nyr)

end subroutine

subroutine dynamic_window_solver_backward_disk(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        AdjSource,dd,dt,ny,nx,nx_w,Vp_w,&
        shotid,unit_fields,folder_fields,&
        dd_space,store_snap,nt_snap,nt,PSDM_thread)

implicit none

integer :: unit_fields,shotid
integer :: i,j,rank,iprop
integer :: nx,ny,nt,nt_snap,store_snap
integer :: NumRec
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real	:: dd_space
real	:: dd,dt
integer :: nxSou(1),nySou(1)
integer :: nxs(1),nys(1)
integer :: nxRec(NumRec),nyRec(NumRec)

real :: AdjSource(nt,NumRec) 
real :: Vp_w(ny,nx_w)
real :: PSDM_thread(ny,nx_w)

integer, allocatable :: nxr(:),nyr(:)

character(len=500) :: folder_fields

allocate(nxr(NumRec),nyr(NumRec))

        call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)

        nxs=nxSou-(x_i-1);nys=nySou;
        nxr=nxRec-(x_i-1);nyr=nyRec;

!!---------- PROPAGATOR
        call solver_ac_backward_disk(rank,1,NumRec,nxs,nys,nxr,nyr,&
        AdjSource,dd,dt,ny,nx_w,Vp_w,&
        shotid,unit_fields,folder_fields,&
        store_snap,nt_snap,nt,PSDM_thread)

    deallocate(nxr,nyr)

end subroutine

subroutine dynamic_window_tt_forward(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        dd,ny,nx,nx_w,Vp_w,dd_space,TT_s)

implicit none

integer :: i,j,rank,iprop
integer :: nx,ny,nt,nt_snap,store_snap
integer :: NumRec
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real	:: dd_space
real	:: dd,dt
integer :: nxSou(1),nySou(1)
integer :: nxs(1),nys(1)
integer :: nxRec(NumRec),nyRec(NumRec)

real :: Vp_w(ny,nx_w)
real :: TT_s(ny,nx_w)

integer, allocatable :: nxr(:),nyr(:)

 allocate(nxr(NumRec),nyr(NumRec))

        call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)
        nxs=nxSou-x_i+1;nys=nySou;
        nxr=nxRec-x_i+1;nyr=nyRec;

!	do j=1,NumRec
!		 write(222,*)j,nx_w,nxs(1),nxr(j),nxRec(j),x_i
!		 write(333,*)j,ny,nyr(j),nyRec(j)
!	enddo

!	write(*,*) 'nx,ny,dx,dy,NumRec:', nx_w, ny, NumRec
!	write(*,*) 'source idx (nxs,nys):', nxs, nys
!	write(*,*) 'min/max v:', minval(Vp_w), maxval(Vp_w)

	call solve_eikonal_FSM(nx_w,ny,dd,dd,Vp_w,nxs(1),nys(1),TT_s)

    deallocate(nxr,nyr)

end subroutine

subroutine dynamic_window_tt_backward(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        dd,ny,nx,nx_w,Vp_w,dd_space,TT_r)

implicit none

integer :: i,j,rank,iprop
integer :: nx,ny,nt,nt_snap,store_snap
integer :: NumRec
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real	:: dd_space
real	:: dd,dt
integer :: nxSou(1),nySou(1)
integer :: nxs(1),nys(1)
integer :: nxRec(NumRec),nyRec(NumRec)

real :: Vp_w(ny,nx_w)
real :: TT_r(ny,nx_w,NumRec)

integer, allocatable :: nxr(:),nyr(:)
real, allocatable :: TTr(:,:)

 allocate(nxr(NumRec),nyr(NumRec))
 allocate(TTr(ny,nx_w))

        call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)
        nxs=nxSou-x_i+1;nys=nySou;
        nxr=nxRec-x_i+1;nyr=nyRec;

	do j=1,NumRec

	   TTr=0.;

	   if (nxr(j)<1 .or. nxr(j)>nx_w) then
		write(*,*)"nxr out of bounds"
		write(*,*)j,nxr(j),nx_w
		stop 


  	   endif

	   if (nyr(j)<1 .or. nyr(j)>ny) stop "nyr out of bounds"

	   call solve_eikonal_FSM(nx_w,ny,dd,dd,Vp_w,nxr(j),nyr(j),TTr)
	   TT_r(:,:,j)=TTr(:,:)

	enddo

    deallocate(nxr,nyr)
    deallocate(TTr)

end subroutine

subroutine dynamic_window_solver_TTs(rank,NumRec,nxSou,nySou,nxRec,nyRec,&
        Data_source,dd,dt,ny,nx,nx_w,Vp_w,&
        dd_space,store_snap,nt_snap,nt,Data_synth,TT_s)

implicit none

integer :: i,j,rank,iprop
integer :: nx,ny,nt,nt_snap,store_snap
integer :: NumRec
integer :: x_i,x_f,y_i,y_f,nx_w
integer :: ddx
real    :: dd_space
real    :: dd,dt
integer :: nxSou(1),nySou(1)
integer :: nxs(1),nys(1)
integer :: nxRec(NumRec),nyRec(NumRec)

real :: Data_source(nt)
real :: Data_synth(nt,NumRec)
real :: Vp_w(ny,nx_w)
real :: TT_s(ny,nx_w)

integer, allocatable :: nxr(:),nyr(:)

allocate(nxr(NumRec),nyr(NumRec))

        call dynamic_window(rank,1,NumRec,nxSou,nxRec,&
        dd,nx,dd_space,nx_w,x_i,x_f)

!        write(*,*)"dentro dynamic solver: ",nx_w,x_i,x_f

        nxs=nxSou-x_i+1; nys=nySou;
        nxr=nxRec-x_i+1; nyr=nyRec;

!!---------- PROPAGATOR

        call solver_ac_forward_ram_tt(rank,1,NumRec,nxs,nys,nxr,nyr,&
        Data_source,dd,dt,ny,nx_w,Vp_w,&
        store_snap,nt_snap,nt,Data_synth,TT_s)

    deallocate(nxr,nyr)

end subroutine
