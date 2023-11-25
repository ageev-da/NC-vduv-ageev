SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
IMPLICIT NONE

INTEGER NI,NJ,IO
double precision,DIMENSION(NI,NJ):: X,Y
double precision,DIMENSION(0:NI,0:NJ)::U,V,P
       
Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
Write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
Write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
Write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)

END SUBROUTINE 


subroutine Output_CF(IO,NI,NJ,X_Cell,Sc_fr,Rex)
implicit none
       
integer :: IO,NI,NJ, NI_C, i
double precision,dimension(NI) :: Sc_fr,Rex, Cf_th
double precision, dimension(0:NI,0:NJ) :: X_Cell
       
NI_C = NI - 1

write(IO,*) 'Variables = "Re<sub>x</sub>", "Cf","Blasius"'
Write(IO,*) 'ZONE I=',NI-2,', J=',1, ',DATAPACKING=BLOCK'
write(IO,'(100E25.16)') Rex(2:NI_C)
write(IO,'(100E25.16)') Sc_fr(2:NI_C)

Cf_th(1) = 0
do i = 2, NI_C
	Cf_th(i) = 0.664/sqrt(Rex(i))
end do

write(IO,'(100E25.16)') Cf_th(2:NI_C)
       
end subroutine


subroutine Output_profiles_vduv(IO,NI,NJ,Y_Cell,U_c,V_c,delta,U0)
implicit none

integer :: IO,NI,NJ, NI_C, NI_S, NI_N, j
double precision :: U0
double precision,dimension(NI) :: delta
double precision, dimension(0:NI,0:NJ) :: U_c,V_c
double precision, dimension(0:NI,0:NJ) :: Y_Cell

NI_N = 44
NI_S = 46
NI_C = NI - 2

write(IO,*) 'Variables = "Velocity x", "Y", "Velocity y", "phi_nac", "nu_nac","phi_ser", "nu_ser", "phi_con", "nu_con"'
Write(IO,*) 'ZONE I=',NJ,', J=',1, ',DATAPACKING=BLOCK'

write(IO,'(100E25.16)') U_c(NI_C,0:NJ-1)
write(IO,'(100E25.16)') Y_Cell(NI_C,0:NJ-1)
write(IO,'(100E25.16)') V_c(NI_C,0:NJ-1)

write(IO,'(100E25.16)') U_c(NI_N,0:NJ-1)/U0
write(IO,'(100E25.16)') Y_Cell(NI_N,0:NJ-1)/delta(NI_N)

write(IO,'(100E25.16)') U_c(NI_S,0:NJ-1)/U0
write(IO,'(100E25.16)') Y_Cell(NI_S,0:NJ-1)/delta(NI_S)

write(IO,'(100E25.16)') U_c(NI_C,0:NJ-1)/U0
write(IO,'(100E25.16)') Y_Cell(NI_C,0:NJ-1)/delta(NI_C)



end subroutine

subroutine Output_profiles_vduv_sec(IO,NI,NJ,Y_Cell,U_c,V_c,delta,U0)
implicit none

integer :: IO,NI,NJ, NI_C, NI_S, NI_N, j
double precision :: U0
double precision,dimension(NI) :: delta
double precision, dimension(0:NI,0:NJ) :: U_c,V_c
double precision, dimension(0:NI,0:NJ) :: Y_Cell

NI_N = 10
NI_S = NI/2
NI_C = NI - 2

write(IO,*) 'Variables = "Velocity x (nac)", "Y", "Velocity y (nac)", "Velocity x (ser)", "Velocity y (ser)", & 
"Velocity x (con)", "Velocity y (con)"'
Write(IO,*) 'ZONE I=',NJ,', J=',1, ',DATAPACKING=BLOCK'

write(IO,'(100E25.16)') U_c(NI_N,0:NJ-1)
write(IO,'(100E25.16)') Y_Cell(NI_C,0:NJ-1)
write(IO,'(100E25.16)') V_c(NI_N,0:NJ-1)
write(IO,'(100E25.16)') U_c(NI_S,0:NJ-1)
write(IO,'(100E25.16)') V_c(NI_S,0:NJ-1)
write(IO,'(100E25.16)') U_c(NI_C,0:NJ-1)
write(IO,'(100E25.16)') V_c(NI_C,0:NJ-1)
end subroutine
