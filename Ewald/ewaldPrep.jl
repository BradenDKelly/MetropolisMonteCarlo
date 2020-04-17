struct EWALD{I}
    kappa::Float64
    nk::I
    k_sq_max::I
    NKVECS::I
end

#                           --------------------------
##########################!!!!!!                 !!!!!!#########################
#                                Prepare EWALDS                               #
##########################!!!!!!                !!!!!!##########################
#                           --------------------------
function PrepareEwaldVariables(ewald::EWALD, boxSize::SVector{3,T} where T)
    kappa = ewald.kappa
    nk = ewald.nk
    k_sq_max = ewald.k_sq_max
    box = min(boxSize...)
    b = 1.0 / 4.0 / kappa / kappa / box / box # kappa already divided by box, this undoes that...
    twopi = 2 * pi
    twopi_sq = twopi^2
    NKVECS = 0

    for kx = 0:nk
        for ky = -nk:nk
            for kz = -nk:nk
                k_sq = kx^2 + ky^2 + kz^2
                if ( ( k_sq <= k_sq_max ) && ( k_sq != 0 ) ) # Test to ensure within range
    			    NKVECS += 1
                end # End test to ensure within range
            end
        end
    end
    ewald = EWALD(ewald.kappa, ewald.nk, ewald.k_sq_max, NKVECS)
    kxyz = Vector{SVector{3,Int32}}(undef,NKVECS)
    cfac = Vector{Float64}(undef,NKVECS)
    #cfac=0; kxx = 0; kyy=0;kzz=0
    NKVECS = 0
    for kx = 0:nk
        for ky = -nk:nk
            for kz = -nk:nk

                k_sq = kx^2 + ky^2 + kz^2

                if ( ( k_sq <= k_sq_max ) && ( k_sq != 0 ) ) # Test to ensure within range
    			    NKVECS += 1
                    kxyz[NKVECS] = SVector{3,Int32}(kx,ky,kz)
                    kr_sq        = twopi_sq * float( k_sq )           # k**2 in real units
                    cfac[NKVECS] = twopi * exp( -b * kr_sq ) / kr_sq / box# Stored expression for later use
    			    if kx > 0 cfac[NKVECS] = cfac[NKVECS] * 2.0 end

                end # End test to ensure within range

            end # kz
        end # ky
    end # kx
    return cfac, kxyz, ewald
end # function
