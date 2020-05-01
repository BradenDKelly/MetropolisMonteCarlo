function Adjust(move_accept, attempt, naccept, L)
# This code adjusts the maximum translation such that the probability of a successful move equals move_accept (fraction)
# Written by Braden Kelly June 19, 2016 / modified for Julia April 2020
# Based on the code by Frenkel and Smith (Fortran 77/90)

#integer, intent(in)             :: attempt, naccept
#integer                         :: attempp, naccepp  ! counts number of attempts and accepts from previous call
#!real, intent(inout)             :: dr ! new adjusted displacement
#real, intent(in)                :: L, move_accept
#real                            :: ratio, dr_ratio, dr_old    ! ratio and place holder for old displacement
#save naccepp, attempp

 if(attempt == 0 ) #! .or. attempp .ge. attempt
     naccepp = naccept
     attempp = attempt
 else
     ratio =  real(naccept - naccepp)/real(attempt-attempp)
     dr_old = dr
     dr = dr*ratio/move_accept
     dr_ratio = dr/dr_old
     #! Place limits on changes
     if(dr_ratio > 1.5) dr = dr_old*1.5 end
     if(dr_ratio < 0.5) dr = dr_old*0.5 end
     if(dr > L/2) dr = L/2 end
     #!write(*,*) dr, dr_old, ratio, attempt - attempp, naccept - naccepp
     naccepp = naccept
     attempp = attempt
 end
 # If you are seeing this, then perhaps it is too late.
end #translation
