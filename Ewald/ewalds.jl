using SpecialFunctions
using Setfield
using StaticArrays

#include("boundaries.jl")


struct EWALD{I}
    kappa::Float64
    nk::I
    k_sq_max::I
    NKVECS::I
    factor::Real
end

"Vector between two coordinate values, accounting for mirror image seperation"
@inline function vector1D(c1::Float64, c2::Float64, box_size::Float64)
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) : (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) : (c2 - c1 + box_size)
    end
end

#                           --------------------------
##########################!!!!!!                 !!!!!!#########################
#                                Prepare EWALDS                               #
##########################!!!!!!                !!!!!!##########################
#                           --------------------------
function PrepareEwaldVariables(ewald::EWALD, boxSize::Real where T)
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

"
struct EWALD{I}
    kappa::Float64
    nk::I
    k_sq_max::I
    NKVECS::I
end
"

"""Secondary Preparation for EWALD Summation"""
#function SetupKVecs(nk, kappa, boxSize)
function SetupKVecs(ewald::EWALD, boxSize::Real where T)

    kappa = ewald.kappa
    nk = ewald.nk
    k_sq_max = ewald.k_sq_max


    #k_sq_max = nk^2 + 2
    println("k_sq_max: ", k_sq_max)
    #kfac = zeros(SVector{k_sq_max + 2, Float64})
    kfacTemp = Vector{Float64}(undef,k_sq_max)
    b = 1.0 / 4.0 / kappa / kappa / boxSize / boxSize  # boxSize may be issue
    #b = 1.0 / 4.0 / kappa / kappa

    @inbounds for kx = 0:(nk)
        for ky = 0:(nk)
            for kz = 0:(nk)
                k_sq = kx^2 + ky^2 + kz^2
                if k_sq <= k_sq_max && k_sq > 0
                    kr_sq = (2*pi)^2 * float(k_sq)
                    kfacTemp[k_sq] = 2 * pi * exp( -b * kr_sq) / kr_sq / boxSize
                    #println(i," ", j," ", k, " ", kfacTemp[k_sq])
                end
            end
        end
    end
    kfac = [kfacTemp[i] for i=1:length(kfacTemp) ]

    return kfac, ewald
end


"Real Ewald contribution"
#function EwaldReal(diff,coord1, coord2,q1,q2, L, rcut2, kappa)
function EwaldReal(qq_r::Vector, qq_q::Vector, kappa::Real, box::Real,
                    thisMol_thisAtom::Vector, chosenOne::Int, r_cut::Real)
    ####
    #
    #    Some prep stuff
    #
    #############
    start_i = thisMol_thisAtom[chosenOne][1]
    end_i   = thisMol_thisAtom[chosenOne][2]
    #r_cut = 10.0
    r_cut_sq = r_cut ^ 2
    pot = 0.0
    rᵢⱼ =rij= SVector{3}(0.0,0.0,0.0)

    """Loop over all atoms in molecule A"""
    for i = start_i:end_i
        rᵢ = qq_r[i]
        """Loop over all atoms in molecule B"""
        for (j,rⱼ) in enumerate(qq_r)

            if j in start_i:end_i continue end # keep mol i from self interacting
            @inbounds for k=1:3
                #println(i," ",j, " ", k, " ", rᵢ[k], " ", rⱼ[k] )
                 rij  = @set rij[k] = vector1D(rᵢ[k], rⱼ[k], box)
            end
            rᵢⱼ² = rij[1]*rij[1] + rij[2]*rij[2] + rij[3]*rij[3]

            if rᵢⱼ² < r_cut_sq
                rᵢⱼ = sqrt( rᵢⱼ² )
                pot += qq_q[i] * qq_q[j] * erfc( kappa * rᵢⱼ ) / rᵢⱼ
                #println(pot, " ", rᵢⱼ, " ",kappa)
            else
                pot += 0.0
            end # potential cutoff
        end # loop over atom in molecule b
    end # loop over atoms in molecule a
    return pot
end

#function RecipLong(boxSize, n, kappa, kfac, nk, k_sq_max, r, qq_q)
""" Recipricol space Ewald energy """
function RecipLong(system::Requirements, ewald::EWALD, r::Vector, qq_q::Vector, kfac::Vector)
    energy = 0.0
    L = system.box
    twopi = 2.0 * π
    nk = ewald.nk
    n = length(r)
    k_sq_max = ewald.k_sq_max

    eikx = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)
    eikz = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)

    @inbounds for j = 1:n
        # calculate kx,ky,kz =1,2 explicitly
        eikx[j, 0] = 1.0 + 0.0 * im
        eiky[j, 0] = 1.0 + 0.0 * im
        eikz[j, 0] = 1.0 + 0.0 * im

        eikx[j, 1] =
            cos(twopi * (r[j][1]) / L) + sin(twopi * (r[j][1]) / L) * im
        eiky[j, 1] =
            cos(twopi * (r[j][2]) / L) + sin(twopi * (r[j][2]) / L) * im
        eikz[j, 1] =
            cos(twopi * (r[j][3]) / L) + sin(twopi * (r[j][3]) / L) * im

        #eikx[j,-1] = conj(eikx[j,1])
        eiky[j, -1] = conj(eiky[j, 1])
        eikz[j, -1] = conj(eikz[j, 1])
    end
    # calculate remaining positive kx, ky, kz by recurrence
    for k = 2:(nk)
        @inbounds for j = 1:n
            eikx[j, k] = eikx[j, k-1] * eikx[j, 1]
            eiky[j, k] = eiky[j, k-1] * eiky[j, 1]
            eikz[j, k] = eikz[j, k-1] * eikz[j, 1]

            #eikx[j,-k] = conj( eikx[j,k] )
            eiky[j, -k] = conj(eiky[j, k])
            eikz[j, -k] = conj(eikz[j, k])
        end
    end

    term = 0.0 + 0.0 * im
    #fill!(term,0.0+0.0*im) #term = zeros(ComplexF64,n)
    # fun party time - lets bring out the loops
    @inbounds for kx = 0:nk
        if kx == 0
            factor = 1.0
        else
            factor = 2.0
        end

        for ky = -nk:nk
            for kz = -nk:nk
                k_sq = kx * kx + ky * ky + kz * kz

                if k_sq <= k_sq_max && k_sq > 0
                    term = 0.0 + 0.0 * im
                    @inbounds for l = 1:n
                        term +=
                            qq_q[l] * eikx[l, kx] * eiky[l, ky] * eikz[l, kz]
                    end
                    energy += factor * kfac[k_sq] * real(conj(term) * term)  # / L
                end # if
            end # kz for loop
        end # ky for loop
    end # kx for loop
    #energy *= qq #/ L
    return energy
end # function

function EwaldSelf(ewald::EWALD, qq_q::Vector)
    kappa = ewald.kappa
    factor = ewald.factor
    return self = - kappa * sum( qq_q.^2 ) / sqrt(π) * factor
end

function TinfoilBoundary(system::Requirements, ewald::EWALD, qq_q::Vector, qq_r::Vector)
    vol = system.box ^ 3
    return 2 * π / 3 / vol * dot(qq_q .* qq_r,qq_q .* qq_r )
end
