SUBROUTINE cf_function(Sc_fr,NI,NJ,dy,U_c,Nu,R0,U0,MU)
      
implicit none
integer :: NI, NJ, i, j
double precision :: dy, Nu, R0, U0, MU
double precision,dimension(NI) :: Sc_fr
double precision, dimension(0:NI,0:NJ) :: U_c

do i = 1, NI - 1      
  Sc_fr(i) = -(MU*((U_c(i,0) - U_c(i,1))/dy)) / &
(0.5*R0*U0**2)

end do
END SUBROUTINE
