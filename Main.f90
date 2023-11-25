Program Nv_ageev
Implicit none
 
INTEGER, PARAMETER:: IO = 12 ! input-output unit
INTEGER I,J,NI,NJ, NITER
double precision :: L,H,U0,MU,NU,R0,P0
double precision :: dx,dy,CFL,Uref,EPS
double precision :: k
double precision,ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
double precision,ALLOCATABLE :: X_Cell(:,:),Y_Cell(:,:)
double precision,ALLOCATABLE :: P_c(:,:),U_c(:,:),V_c(:,:)
double precision,dimension (:), allocatable :: Sc_fr,Rex, delta
 
write(*,*) 'Read input file'
open(IO,FILE='Input.txt')
read(IO,*) L
read(IO,*) H
read(IO,*) NI
read(IO,*) NJ
read(IO,*) NITER
read(IO,*) CFL
read(IO,*) Uref
read(IO,*) EPS
 
read(IO,*) U0
read(IO,*) MU
read(IO,*) R0
read(IO,*) P0
read(IO,*) k
CLOSE(IO)

allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
allocate(X_Cell(0:NI,0:NJ)) ! cell centers X-coordinates
allocate(Y_Cell(0:NI,0:NJ)) ! cell centers Y-coordinates

!------ Cell-centered variables
allocate(U_c(0:NI,0:NJ))   ! Velocity U
allocate(V_c(0:NI,0:NJ))   ! Velocity V
allocate(P_c(0:NI,0:NJ))   ! Pressure
allocate(Sc_fr(NI),Rex(NI),delta(NI))

dx=L/(NI-1)
dy=H/(NJ-1)

!------ Coordinate of nodes
DO I=1,NI
DO J=1,NJ
  X_Node(I,J)=(I-1)*dx
  Y_Node(I,J)=(J-1)*dy
END DO
END DO

!------ Coordinate of cell centers
X_Cell(0,1:NJ)=-dx/2
Y_Cell(0,1:NJ)=Y_Node(1,1:NJ)+dy/2
X_Cell(1:NI,0)=X_Node(1:NI,1)+dx/2
Y_Cell(1:NI,0)=-dy/2
DO I=1,NI
DO J=1,NJ
  X_Cell(I,J)=X_Node(I,J)+dx/2
  Y_Cell(I,J)=Y_Node(I,J)+dy/2
END DO
END DO

!----------------- Parameters ------------------------

NU=MU/R0

DO I = 1, NI
	Rex(I) = U0 * X_Node(I,1)/NU
END DO

write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
write(*,*)'ReL= ', U0*L/NU
!        pause    

!----------------- Initial fields -----------------------------

DO I=0,NI
DO J=0,NJ
  U_c(I,J)=U0
  V_c(I,J)=1.0e-5
  P_c(I,J)=0.0
ENDDO
ENDDO

!---------------- Solve Navier-Stokes equations ---------------------
 
write(*,*) 'Solve Navier-Stokes equations' 
Call NavierStokes(NI,NJ,dx,dy,U_c,V_c,P_c,CFL,Uref,NITER,EPS,U0,NU, &
k,X_Cell)

call cf_function(Sc_fr,NI,NJ,dy,U_c,Nu,R0,U0,MU)
 !----------------- Output data ------------------------------
 
write(*,*) 'Output data cell (Navier-Stokes)' 
Open(IO,FILE='Results_Navier_vduv.plt')
Call OutputFields_Cell(IO,NI,NJ,X_Node,Y_Node,U_c,V_c,P_c)
Close(IO)

Open(IO,FILE='CF_Navier_vduv.plt')
call Output_CF(IO,NI,NJ,X_Cell,Sc_fr,Rex)  
close(IO)     

do I = 1, NI - 1
j = 0
do while (U_c(i,j) < 0.99 * U0)
j = j + 1
end do
delta(i) = Y_Cell(i,j)

end do

Open(IO,FILE='profiles_Navier_vduv.plt')
call Output_profiles_vduv(IO,NI,NJ,Y_Cell,U_c,V_c,delta,U0) 
close(IO)

Open(IO,FILE='profiles_Navier_vduv_sec.plt')
call Output_profiles_vduv_sec(IO,NI,NJ,Y_Cell,U_c,V_c,delta,U0)
close(IO)

END PROGRAM
   
