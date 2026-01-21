subroutine linear_2D(nxmodel,nymodel,dmodel,vel,nx_prop,ny_prop,dmodel_prop,Vpmodel)

implicit none

integer :: i,j,k
integer :: nxmodel,nymodel,nx_prop,ny_prop
real	:: dmodel,dmodel_prop
real	:: xmodel(nxmodel)
real	:: ymodel(nymodel)
real	:: xmodel_prop(nx_prop)
real	:: ymodel_prop(ny_prop)
real	:: vel(nymodel,nxmodel)
real	:: vel_temp(ny_prop,nxmodel)
real	:: Vpmodel(ny_prop,nx_prop)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   INTRPOLATION FROM nxmodel to nx_prop and from nymodel to ny_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do i=1,nxmodel
        xmodel(i)=(i-1)*dmodel
  enddo
  do i=1,nymodel
        ymodel(i)=(i-1)*dmodel
  enddo
  do i=1,nx_prop
        xmodel_prop(i)=(i-1)*dmodel_prop
  enddo
  do i=1,ny_prop
        ymodel_prop(i)=(i-1)*dmodel_prop
  enddo

  do i=1,nxmodel
       call  linear_interpolation_1D(ymodel,vel(:,i),nymodel,ymodel_prop,vel_temp(:,i),ny_prop)
  enddo
  
  do i=1,ny_prop
       call  linear_interpolation_1D(xmodel,vel_temp(i,:),nxmodel,xmodel_prop,Vpmodel(i,:),nx_prop)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine 

subroutine linear_1D(nx1,dmodel1,model1,nx2,dmodel2,model2)

implicit none

integer :: i,j,k
integer :: nx1,nx2
real	:: dmodel1,dmodel2
real	:: xmodel1(nx1)
real	:: xmodel2(nx2)
real	:: model1(nx1)
real	:: model2(nx2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   INTRPOLATION FROM nxmodel to nx_prop and from nymodel to ny_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do i=1,nx1
        xmodel1(i)=(i-1)*dmodel1
  enddo
  do i=1,nx2
        xmodel2(i)=(i-1)*dmodel2
  enddo

       call  linear_interpolation_1D(xmodel1,model1,nx1,xmodel2,model2,nx2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine 


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

        ! Extrapolación a la izquierda
        if (a3 <= x1(1)) then
            y2(j) = y1(1)

        ! Extrapolación a la derecha
        elseif (a3 >= x1(n1)) then
            y2(j) = y1(n1)

        ! Interpolación dentro del rango
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

subroutine interpola_lineal1D(l,xc,yc,nl,tml,S_e2)

implicit none

        integer i,j,k,l,nl
        real :: xc(l),yc(l),tml(nl),S_e2(nl),m

        do i=1,l
        do j=1,nl

        m=(yc(i+1)-yc(i))/(xc(i+1)-xc(i))

        if (tml(j).ge.xc(i).and.tml(j).lt.xc(i+1)) then
                S_e2(j)=yc(i)+(tml(j)-xc(i))*m
        end if

        enddo
        enddo

        do j=1,nl
                S_e2(j)=S_e2(j)+0.001
        enddo

return
end
