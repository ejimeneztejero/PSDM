subroutine optimization(rank)

use mod_parfile
use mod_data_arrays

implicit none
integer :: rank
real :: dt_a,dt_b

!!!!!-----------------------------------
!!!!!------ Optimization
!!!!!-----------------------------------

!if(rank.eq.0)write(*,*)"INPUT compute_grid_params: ", freq_ricker, nxmodel,nymodel,dmodel

if(RTM_MCS.eq.1)	then

    if(automatic_optimization.eq.0)     then

        n_cores_real=n_cores_total
        nx_prop=1+(nxmodel-1)*(dmodel/dmodel_prop)
        ny_prop=1+(nymodel-1)*(dmodel/dmodel_prop)

        dt_a=dmodel_prop/maxval(vel)/2.
	dt_b=0.5d0*dmodel_prop/maxval(vel)/sqrt(2.0d0)
	dt_prop=dt_a	!this is the original

	if(rank.eq.0)write(*,*)"cositas dt_a,dt_prop,dmodel_prop,maxval(vel): ",dt_a,dt_prop,dmodel_prop,maxval(vel)

        nt_prop=1+ceiling(tfin_mcs/dt_prop)
	nt_snap=ceiling(1.*nt_prop/store_snap)
        nt_prop=nt_snap*store_snap
        
    endif

    if(automatic_optimization.eq.1)     then
    
            v_min=abs(minval(vel))
            v_max=abs(maxval(vel))
            y_max=(nymodel-1)*dmodel    ! en metros  
            x_max=(nxmodel-1)*dmodel    ! en metros  
            !x_max=(NumRec_mcs-1)*drec+1000     ! ojo, tb para was
            t_max=tfin_mcs
                        
            call compute_grid_params(freq_ricker, v_min, v_max, drec, x_max, y_max, t_max, &
                nx_prop, ny_prop, nt_prop, dmodel_prop, dt_prop)
                
            call rtm_strategy(nx_prop, ny_prop, nt_prop, RAM_total_GB, Disk_total_GB, &
                n_cores_total, n_cores_real, strategy, strategy_num, store_snap, nt_snap)
                
!            if(rank.eq.0)write(*,*)"After compute_grid_params: ", nx_prop, ny_prop, nt_prop, dmodel_prop, dt_prop
!            if(rank.eq.0)write(*,*)"After rtm_strategy: ", n_cores_real, strategy, store_snap,nt_snap
            
    endif 

endif

!if(rank.eq.0)write(*,*)"After compute_grid_params: ", nx_prop, ny_prop, dmodel_prop, nt_prop, dt_prop,store_snap,nt_snap
!if(rank.eq.0)write(*,*)"After rtm_strategy: ", n_cores_real, strategy

end subroutine optimization


subroutine compute_grid_params(freq_ricker, v_min, v_max, drec, x_max, y_max, t_max, &
                               nx_prop, ny_prop, nt_prop, dmodel_prop, dt_prop)
  implicit none
  ! INPUT
  real :: freq_ricker   ! Central frequency of Ricker wavelet (Hz)
  real :: v_min, v_max    ! Min and max velocities (m/s)
  real :: drec          ! Receiver spacing (m)
  real :: x_max, y_max  ! Model dimensions (m)
  real :: dt_a,dt_b,t_max         ! Maximum propagation time (s)
  ! OUTPUT
  integer, intent(out) :: nx_prop, ny_prop, nt_prop
  real, intent(out) :: dmodel_prop,dt_prop

  ! Internal variables
  real :: dx_prop,lambda_min, points_per_wavelength
  real :: cfl

  ! --- Step 1: dispersion criterion ---
  points_per_wavelength = 10.0d0   ! typical: 8–12 points per λ
  lambda_min = v_min / freq_ricker
  dx_prop = lambda_min / points_per_wavelength
  dmodel_prop=min(dx_prop,drec)

!  write(*,*)"lambda_min: ",lambda_min,v_min,freq_ricker
!  write(*,*)"dmodel_prop: ",dmodel_prop

  ! Ensure we do not oversample more than necessary in horizontal
  if (dmodel_prop < drec / 2.0d0) dmodel_prop = drec / 2.0d0

  ! --- Step 2: CFL stability condition (2D acoustic) ---
  dt_a=dmodel_prop/v_max/2.
  dt_b=0.5d0*dmodel_prop/v_max/sqrt(2.0d0)
  dt_prop=dt_a	!this is the original

!  write(*,*)"dt_prop: ",dt_prop

  ! --- Step 3: grid points ---
  nx_prop = ceiling(x_max / dmodel_prop)
  ny_prop = ceiling(y_max / dmodel_prop)

!  write(*,*)"nx_prop: ",nx_prop
!  write(*,*)"ny_prop: ",ny_prop

  ! --- Step 4: time steps ---
  nt_prop = ceiling(t_max / dt_prop)

!  write(*,*)"nt_prop: ",nt_prop

end subroutine compute_grid_params

subroutine rtm_strategy(nx_prop, ny_prop, nt_prop, RAM_total_GB, Disk_total_GB, &
                        n_cores_total, n_cores_real, strategy, strategy_num,store_snap, nt_snap)

  implicit none
  ! Inputs
  integer, intent(in) :: nx_prop, ny_prop, nt_prop
  real, intent(in) :: RAM_total_GB, Disk_total_GB
  integer, intent(in) ::  n_cores_total
  integer, intent(out) :: n_cores_real
  integer :: strategy_num
  ! Outputs
  character(len=20), intent(out) :: strategy
  integer, intent(out) :: store_snap, nt_snap

  ! Locals
  real :: tmp
  real :: bytes_per_number,bytes_per_field, RAM_needed_full, RAM_needed_snap
  real :: Disk_needed_full, Disk_needed_snap
  real, parameter :: GB = 1024.0d0**3
  integer :: max_shots_in_parallel

  ! Memoria para un campo completo en RAM (double precision = 8 bytes)
    bytes_per_number = storage_size(tmp) / 8    ! en bytes
    bytes_per_field  = dble(nx_prop) * dble(ny_prop) * dble(bytes_per_number)

!  bytes_per_field = dble(nx_prop) * dble(ny_prop) * 8.0d0

  ! Memoria total para forward completo (todos los tiempos) por un shot
  RAM_needed_full = bytes_per_field * dble(nt_prop) / GB
  Disk_needed_full = RAM_needed_full  ! igual porque binario ocupa lo mismo que en RAM

  write(*,*)"RAM_needed_full: ",RAM_needed_full,bytes_per_field

  ! ====== 1) Calcular máximo número de shots en paralelo ======
  max_shots_in_parallel = floor(RAM_total_GB / RAM_needed_full)
  if (max_shots_in_parallel < 1) max_shots_in_parallel = 1
  n_cores_real = min(n_cores_total, max_shots_in_parallel)

  ! ====== 2) Decidir estrategia ======
  if (RAM_needed_full <= RAM_total_GB / n_cores_real) then
     ! Todo cabe en RAM
     strategy = 'RAM_FULL'
     strategy_num = 1
     store_snap = 1
     nt_snap = nt_prop
  else if (Disk_needed_full <= Disk_total_GB) then
     ! Guardar forward completo en disco
     strategy = 'DISK_FULL'
     strategy_num = 2
     store_snap = 1
     nt_snap = nt_prop
  else
     ! Necesario usar snapshots
     ! Calcular store_snap para que RAM/disk no se exceda
     store_snap = ceiling(RAM_needed_full / (RAM_total_GB / n_cores_real))
     if (store_snap < 1) store_snap = 1

     nt_snap = floor(real(nt_prop) / real(store_snap)) + 1

     ! Comprobar si es RAM o DISK checkpoints
     RAM_needed_snap = bytes_per_field * nt_snap / GB
     Disk_needed_snap = RAM_needed_snap

     if (RAM_needed_snap <= RAM_total_GB / n_cores_real) then
        strategy = 'RAM_CHECKPOINTS'
     strategy_num = 3
     else if (Disk_needed_snap <= Disk_total_GB) then
        strategy = 'DISK_CHECKPOINTS'
     strategy_num = 4
     else
        strategy = 'ERROR_NO_SPACE'
     end if
  end if

end subroutine rtm_strategy

