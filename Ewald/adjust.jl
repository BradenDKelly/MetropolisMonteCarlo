function Adjust!(saved::Moves, L)
# This code adjusts the maximum translation such that the probability of a successful move equals move_accept (fraction)
# Written by Braden Kelly June 19, 2016 / modified for Julia April 2020
# Based on the code by Frenkel and Smith (Fortran 77/90)
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
     dr_old = saved.d_max
     saved.d_max = saved.d_max*ratio/saved.set_value #move_accept
     dr_ratio = saved.d_max/dr_old
     #! Place limits on changes
     if(dr_ratio > 1.5) saved.d_max = dr_old*1.5 end
     if(dr_ratio < 0.5) saved.d_max = dr_old*0.5 end
     if(saved.d_max > L/2) saved.d_max = L/2 end
     #!write(*,*) dr, dr_old, ratio, attempt - attempp, naccept - naccepp
     saved.naccepp = saved.naccept
     saved.attempp = saved.attempt
 end
 # If you are seeing this, then perhaps it is too late.
 return saved
end #translation

function Adjust_rot!(saved::Moves, L)
# This code adjusts the maximum rotation such that the probability of a successful move equals move_accept (fraction)
# Written by Braden Kelly June 19, 2016 / modified for Julia April 2020
# Based on the code by Frenkel and Smith (Fortran 77/90)
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
     dϕ_old = saved.d_max
     saved.d_max = saved.d_max*ratio/saved.set_value #move_accept
     dϕ_ratio = saved.d_max/dϕ_old
     #! Place limits on changes
     if(dϕ_ratio > 1.5) saved.d_max = dϕ_old*1.5 end
     if(dϕ_ratio < 0.5) saved.d_max = dϕ_old*0.5 end
     if(saved.d_max > L/2) saved.d_max = L/2 end
     #!write(*,*) dr, dr_old, ratio, attempt - attempp, naccept - naccepp
     saved.naccepp = saved.naccept
     saved.attempp = saved.attempt
 end
 # If you are seeing this, then perhaps it is too late.
 return saved
end #translation
