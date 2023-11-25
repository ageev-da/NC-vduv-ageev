!************************************************************************************************      
SUBROUTINE NavierStokes(NI,NJ,dx,dy,U_c,V_c,P_c,CFL,Uref,NITER,EPS,U0,NU,k,X_Cell)
IMPLICIT NONE
integer :: NI,NJ,NITER
double precision :: dx,dy,CFL,Uref,EPS,U0,NU
double precision :: k
double precision, dimension(0:NI,0:NJ) :: U_c,V_c,P_c
double precision,dimension(0:NI,0:NJ) :: X_Cell
       
integer :: i,j,iter_num
double precision :: A,dt, u12, v12, p12, uc12, vc12 
double precision :: Mass,Imp_x,Imp_y
double precision, dimension(0:NI,0:NJ) :: Res_u,Res_v,Res_p
       
A = 1.0/(Uref*Uref)
       
iter_num = 0
       
if (dx > dy) then
  dt = CFL*dy/U0
else 
  dt = CFL*dx/U0
end if
       
open(1,file='Residuals_vduv.plt')
write(1,*) 'VARIABLES = "Iter", "P_res", "U_res", "V_res"' 
write(1,*) 'Zone i=',NITER
       
do while ((iter_num == 0) .or. (((maxval(abs(Res_u)) > EPS) .or. &
(maxval(abs(Res_v)) > EPS) .or. (maxval(abs(Res_p)) > EPS)) .and. &
(iter_num < NITER)))
         
Res_u(:,:) = 0.0
Res_v(:,:) = 0.0
Res_p(:,:) = 0.0
         
!входные граничные условия
do j = 1,NJ-1
  U_c(0,j) = U0 !2*U0-U_c(1,j)
  V_c(0,j) = 0.0 !-V_c(1,j)
  P_c(0,j) = P_c(1,j)
end do
         
!выходные граничные условия
do j = 1,NJ-1
U_c(NI,j) = U_c(NI-1,j)
V_c(NI,j) = V_c(NI-1,j)
P_c(NI,j) = 0.0
end do
         
!условия на стенке
do i = 1,NI-1
  U_c(i,0) = -U_c(i,1)
  V_c(i,0) = -V_c(i,1) + 2.0 * k * X_Cell(i,0) * U0
  P_c(i,0) = P_c(i,1)
end do
         
!условия на верхней границе
do i = 1,NI-1
V_c(i,NJ) = V_c(i,NJ-1)
if (V_c(i,NJ-1) > 0.0) then
  U_c(i,NJ) = U_c(i,NJ-1)
  P_c(i,NJ) = 0.0
else
  U_c(i,NJ) = U0
  P_c(i,NJ) = P_c(i,NJ-1)
end if
end do
         
!перерасчет переменных (параллельно х)
do j = 1,NJ-1
do i = 0,NI-1
           
  uc12 = (U_c(i+1,j) + U_c(i,j))/2.0
             
  if (uc12 >= 0.0) then
    u12 = U_c(i,j)
    v12 = V_c(i,j)
    p12 = P_c(i+1,j)
  else
    u12 = U_c(i+1,j)
    v12 = V_c(i+1,j)
    p12 = P_c(i,j)
  end if
             
  Mass = u12
  Imp_x = uc12*u12 + p12 - NU*(U_c(i+1,j)-U_c(i,j))/dx
  Imp_y = uc12*v12 - NU*(V_c(i+1,j)-V_c(i,j))/dx
             
  !i,j
  Res_p(i,j) = Res_p(i,j) + Mass/dx
  Res_u(i,j) = Res_u(i,j) + Imp_x/dx
  Res_v(i,j) = Res_v(i,j) + Imp_y/dx
             
  !i+1,j
  Res_p(i+1,j) = Res_p(i+1,j) - Mass/dx
  Res_u(i+1,j) = Res_u(i+1,j) - Imp_x/dx
  Res_v(i+1,j) = Res_v(i+1,j) - Imp_y/dx
             
end do

Res_p(0,j) = 0.0
Res_u(0,j) = 0.0
Res_v(0,j) = 0.0
Res_p(NI,j) = 0.0
Res_u(NI,j) = 0.0
Res_v(NI,j) = 0.0

end do
         
!перерасчет переменных (паралелльно y)
do i = 1,NI-1
do j = 0,NJ-1
           
  vc12 = (V_c(i,j+1) + V_c(i,j))/2.0
             
  if (vc12 >= 0.0) then
    u12 = U_c(i,j)
    v12 = V_c(i,j)
    p12 = P_c(i,j+1)
  else
    u12 = U_c(i,j+1)
    v12 = V_c(i,j+1)
    p12 = P_c(i,j)
  end if
             
             
  Mass = v12
  Imp_x = vc12*u12 - NU*(U_c(i,j+1)-U_c(i,j))/dy
  Imp_y = vc12*v12 + p12 - NU*(V_c(i,j+1)-V_c(i,j))/dy

  !коррекция на необходимой границе
  if(j == 0) then
    Mass = vc12
  end if


             
  !i,j
  Res_p(i,j) = Res_p(i,j) + Mass/dy
  Res_u(i,j) = Res_u(i,j) + Imp_x/dy
  Res_v(i,j) = Res_v(i,j) + Imp_y/dy
             
  !i,j+1
  Res_p(i,j+1) = Res_p(i,j+1) - Mass/dy
  Res_u(i,j+1) = Res_u(i,j+1) - Imp_x/dy
  Res_v(i,j+1) = Res_v(i,j+1) - Imp_y/dy
             
end do
end do
         
do i = 0,NI
do j = 0,NJ
  U_c(i,j) = U_c(i,j) - dt*Res_u(i,j)
  V_c(i,j) = V_c(i,j) - dt*Res_v(i,j)
  P_c(i,j) = P_c(i,j) - dt/A*Res_p(i,j)
end do

Res_p(i,0) = 0.0
Res_u(i,0) = 0.0
Res_v(i,0) = 0.0
Res_p(i,NJ) = 0.0
Res_u(i,NJ) = 0.0
Res_v(i,NJ) = 0.0

end do
         
iter_num = iter_num + 1
         
write(1,*) iter_num,maxval(abs(Res_p))/A, & 
maxval(abs(Res_u)),maxval(abs(Res_v))
         
end do
    
close(1)
       
END  SUBROUTINE

!************************************************************************************************
                    
