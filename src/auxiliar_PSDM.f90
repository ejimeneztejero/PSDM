subroutine normalize_psdm(nx,ny,dx,dy,PSDM_sum,Cov_sum,PSDM)
   integer, intent(in) :: nx,ny
   real, intent(in) :: dy,dx
   real, intent(in) :: PSDM_sum(ny,nx)
   real, intent(in) :: Cov_sum(ny,nx)
   real, intent(out)	:: PSDM(ny,nx)
   integer :: i,j,k

   PSDM = PSDM_sum
   do i = 1, nx
      do j = 1, ny
         if (Cov_sum(j,i) > 0) then
            PSDM(j,i) = PSDM(j,i) / Cov_sum(j,i)
         else
            PSDM(j,i) = 0.0
         endif
      enddo
   enddo

end subroutine normalize_psdm

subroutine psdm2twt(nx,nz,dd,dt,file_vel,file_psdm,vel_total,psdm_total)
    implicit none

    !-------------------------
    ! INPUTS
    !-------------------------
    integer, intent(in) :: nx, nz     ! model size
    real, intent(in) :: dd                   ! depth sampling (m)
    real, intent(in) :: dt                   ! desired TWT sampling (s)
    real, intent(in) :: vel_total(nz,nx)
    real, intent(in) :: psdm_total(nz,nx)
    character(len=1000), intent(in) :: file_vel,file_psdm


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

    call write_vel(nx,nt,dd,dt,1,file_vel,vel_twt)

    call write_psdm(nx,nt,dd,dt,file_psdm,psdm_twt,1,1)!!trid,last)

    deallocate(time,time_total,vel_twt,psdm_twt)

end subroutine psdm2twt


subroutine get_psdm(nx,nz,psdm)

  USE mod_parfile
  USE mod_data_arrays

  implicit none

  integer :: i,j,nx,nz
  real  :: xi,xj
  real  :: psdm(nz,nx)
  real, allocatable :: mi(:,:)
  character(len=500) :: file_name
  CHARACTER(len=50) :: access,form

  access = 'stream'
  form = 'unformatted'


!!! PSDM

!!      ARREGLAR

        file_name=trim(adjustl(psdm_file))

        if(psdm_v.eq.1)  then
                open(unit=10,file=file_name,status='old')
                do j=1,nz
                        read(10,*) (psdm(j,i),i=1,nx)
                enddo
        endif
        if(psdm_xzv.eq.1)        then
                open(unit=10,file=file_name,status='old')
                do i=1,nx
                        do j=1,nz
                                read(10,*) xi,xj,psdm(j,i)
                        enddo
                enddo
        endif
        if(psdm_su.eq.1) then

                open(10,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')
                call get_su_data(nx,nz,nz,10,psdm)

        endif

        close(10)


end subroutine get_psdm



SUBROUTINE TT_source_onfly(p_3, p_1, TT_source, amp_prev, snap, dt, dnx, dny, nx, ny)

  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: nx, ny, snap, dnx, dny
  REAL,    INTENT(IN)    :: dt
  REAL,    INTENT(IN)    :: p_3(dny+ny, dnx+nx), p_1(dny+ny, dnx+nx)
  REAL,    INTENT(INOUT) :: TT_source(ny,nx), amp_prev(ny,nx)

  INTEGER :: i, j
  REAL :: amp_now, amp_prev_local, amp_max_global, frac
  REAL, SAVE :: max_initialized = 0.0
  REAL, SAVE :: amp_max_global_saved = 0.0
  REAL, PARAMETER :: threshold_min = 0.02  ! ~2% del máximo global

!----------------------------------------------------------------------
! 1. Initialize on first step
!----------------------------------------------------------------------
  IF (snap == 1) THEN
     TT_source(:,:) = 0.0
     amp_prev(:,:)  = 0.0
     amp_max_global_saved = 0.0
     max_initialized = 0.0
  ENDIF


!----------------------------------------------------------------------
! 2. Compute max amplitude (only once or every few steps if needed)
!----------------------------------------------------------------------
  IF (max_initialized == 0.0) THEN
     amp_max_global = 0.0
     DO i = 1, nx
        DO j = 1, ny
           amp_max_global = MAX(amp_max_global, ABS(p_3(dny+j, dnx+i)))
        ENDDO
     ENDDO
     amp_max_global_saved = amp_max_global + 1.e-20
     max_initialized = 1.0
  ENDIF

!----------------------------------------------------------------------
! 3. Detect first arrival using derivative amplitude normalized
!----------------------------------------------------------------------
  DO i = 1, nx
     DO j = 1, ny
        IF (TT_source(j,i) == 0.0) THEN
           amp_now = ABS(p_3(dny+j, dnx+i) - p_1(dny+j, dnx+i)) / amp_max_global_saved
           amp_prev_local = amp_prev(j,i)

           !--- Detect threshold crossing with sub-sample interpolation
           IF (amp_prev_local < threshold_min .AND. amp_now >= threshold_min) THEN
              frac = (threshold_min - amp_prev_local) / (amp_now - amp_prev_local + 1.e-10)
              TT_source(j,i) = (snap - 1 + frac) * dt
           ENDIF

           amp_prev(j,i) = amp_now
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE TT_source_onfly
!======================================================================


SUBROUTINE TT_from_Fields(nx,ny,nt,dt,Fields,TT)

IMPLICIT NONE

integer :: i,j,k
integer :: nx,ny,nt,kk
real    :: dt
real    :: threshold_min,amp_max,norm_val
real    :: Fields(ny,nx,nt),TT(ny,nx)

        TT(:,:) = 0.0 ! j,i
        threshold_min = 0.02

        do i = 1, nx
            do j = 1, ny

                kk=0
                amp_max = maxval(abs(Fields(j,i,:)))
                if (amp_max <= 0.0) cycle

                do k = 1, nt
                    norm_val = abs(Fields(j,i,k)) / amp_max

                    if (kk.eq.0.and.norm_val.gt.threshold_min) then
                        TT(j,i) = (k-1) * dt
                        kk=1
                        !write(*,*)"sale"
                        !exit
                    endif
                enddo

            enddo
        enddo

END SUBROUTINE TT_from_Fields

subroutine v_eff(z,z_sf,Veff)

  implicit none
  real, intent(in) :: z, z_sf
  real, intent(out) :: Veff

  real, parameter :: v_water = 1500.0
  real, parameter :: v_max   = 5000.0
  real, parameter :: z_ramp  = 10000.0

  if (z.le.z_sf) then
     Veff = v_water

  else if (z.gt.z_sf.and.z.le.(z_sf+z_ramp)) then
     Veff = v_water + (v_max - v_water) * (z - z_sf) / z_ramp

  else
     Veff = v_max
  endif

end subroutine v_eff
