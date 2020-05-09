function Adjust!(saved::Moves, L)
# This code adjusts the maximum translation such that the probability of a successful move equals move_accept (fraction)
# Written by Braden Kelly June 19, 2016 / modified for Julia April 2020
# Based on the code by Frenkel and Smith (Fortran 77/90)

#integer, intent(in)             :: attempt, naccept
#integer                         :: attempp, naccepp  ! counts number of attempts and accepts from previous call
#!real, intent(inout)             :: dr ! new adjusted displacement
#real, intent(in)                :: L, move_accept
#real                            :: ratio, dr_ratio, dr_old    ! ratio and place holder for old displacement
#save naccepp, attempp


#=
---------------------------
    naccepp is number of moves accepted as of previous call
    naccept is number of moves accepted as of this call
    attempp is number of attempts as of previous call
    attemt is number of attempts as of this call

    ratio is a numerical derivative

=#
 if(saved.attempp == 0 ) #! .or. attempp .ge. attempt
     saved.naccepp = saved.naccept
     saved.attempp = saved.attempt
 else
     ratio =  real(saved.naccept - saved.naccepp)/real(saved.attempt-saved.attempp)
     dr_old = saved.dr_max
     saved.dr_max = saved.dr_max*ratio/saved.set_value #move_accept
     dr_ratio = saved.dr_max/dr_old
     #! Place limits on changes
     if(dr_ratio > 1.5) saved.dr_max = dr_old*1.5 end
     if(dr_ratio < 0.5) saved.dr_max = dr_old*0.5 end
     if(saved.dr_max > L/2) saved.dr_max = L/2 end
     #!write(*,*) dr, dr_old, ratio, attempt - attempp, naccept - naccepp
     saved.naccepp = saved.naccept
     saved.attempp = saved.attempt
 end
 # If you are seeing this, then perhaps it is too late.
 return saved
end #translation
