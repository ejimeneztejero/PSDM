subroutine solve_eikonal_FSM(nx_prop, ny_prop, dx_prop, dy_prop, vmodel_prop, sx, sz, Traveltime)
  implicit none
  integer, intent(in) :: nx_prop, ny_prop
  real, intent(in) :: dx_prop, dy_prop
  real, intent(in) :: vmodel_prop(ny_prop, nx_prop)
  integer, intent(in) :: sx(1), sz(1)
  real, intent(out) :: Traveltime(ny_prop, nx_prop)

  integer :: i, j, sweep
  real :: BIG
  parameter (BIG = 1.0e9)

  Traveltime = BIG
  Traveltime(sz, sx) = 0.0

  ! barridos para convergencia
  do sweep = 1,4
     call sweep_pass(nx_prop, ny_prop, dx_prop, dy_prop, vmodel_prop, Traveltime, sweep)
  enddo

contains

  subroutine sweep_pass(nx, ny, dx, dy, v, T, direction)
    integer, intent(in) :: nx, ny, direction
    real, intent(in) :: dx, dy
    real, intent(in) :: v(ny,nx)
    real, intent(inout) :: T(ny,nx)
    integer :: i, j, istart, iend, istep, jstart, jend, jstep
    real :: a, b, t_new, v_inv

    ! Selecciona la dirección del barrido
    select case(direction)
    case(1)
       istart=1; iend=nx; istep=1
       jstart=1; jend=ny; jstep=1
    case(2)
       istart=nx; iend=1; istep=-1
       jstart=1; jend=ny; jstep=1
    case(3)
       istart=nx; iend=1; istep=-1
       jstart=ny; jend=1; jstep=-1
    case(4)
       istart=1; iend=nx; istep=1
       jstart=ny; jend=1; jstep=-1
    end select

    do j=jstart,jend,jstep
       do i=istart,iend,istep
          if (T(j,i)==0.0) cycle  ! punto fuente

          ! vecinos más cercanos
          a = min(T(j, max(i-1,1)), T(j, min(i+1,nx)))
          b = min(T(max(j-1,1), i), T(min(j+1,ny), i))
          v_inv = 1.0 / v(j,i)

          ! cálculo local
          t_new = local_solver(a,b,dx,dy,v_inv)
          if (t_new < T(j,i)) T(j,i) = t_new
       enddo
    enddo
  end subroutine sweep_pass


  function local_solver(a,b,dx,dy,v_inv) result(t)
    real, intent(in) :: a,b,dx,dy,v_inv
    real :: t, aa, bb, cc, disc

    ! ecuación cuadrática local
    aa = 1.0/(dx*dx) + 1.0/(dy*dy)
    bb = -2.0*(a/(dx*dx) + b/(dy*dy))
    cc = (a*a)/(dx*dx) + (b*b)/(dy*dy) - v_inv*v_inv

    disc = bb*bb - 4.0*aa*cc
    if (disc >= 0.0) then
       t = (-bb + sqrt(disc)) / (2.0*aa)
       if (t < max(a,b)) t = max(a,b) + v_inv*min(dx,dy)
    else
       t = min(a,b) + v_inv*min(dx,dy)
    endif

  end function local_solver

end subroutine solve_eikonal_FSM
