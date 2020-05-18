using SpecialFunctions
using Setfield
using StaticArrays

#include("boundaries.jl")

# TODO make another struct that holds qq_q, qq_r, sumQExpOld, sumQExpNew since
# they are mutable, and the rest are not.
mutable struct EWALD{I}
    kappa::Float64
    nk::I
    k_sq_max::I
    NKVECS::I
    kxyz::Vector{SVector{3,Int32}}
    cfac::Vector{Float64}
    sumQExpOld::Vector{ComplexF64}
    sumQExpNew::Vector{ComplexF64}
    factor::Float64
end
#=
struct EWALD{I}
    kappa::Float64
    nk::I
    k_sq_max::I
    cfac::Vector
    factor::Real
end
=#
"Vector between two coordinate values, accounting for mirror image seperation"
@inline function vector1D(c1::Float64, c2::Float64, box_size::Float64)
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) :
               (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) :
               (c2 - c1 + box_size)
    end
end

#                           --------------------------
##########################!!!!!!                 !!!!!!#########################
#                                Prepare EWALDS                               #
##########################!!!!!!                !!!!!!##########################
#                           --------------------------
function PrepareEwaldVariables(ewald::EWALD, boxSize::Real where {T})
    kappa = ewald.kappa
    nk = ewald.nk
    k_sq_max = ewald.k_sq_max
    @assert k_sq_max == 27
    fact = ewald.factor
    box = min(boxSize...)
    b = 1.0 / 4.0 / kappa / kappa / box / box # kappa already divided by box, this undoes that...
    twopi = 2.0 * pi
    twopi_sq = twopi^2
    NKVECS = 0

    for kx = 0:nk
        for ky = -nk:nk
            for kz = -nk:nk
                k_sq = kx^2 + ky^2 + kz^2
                if ((k_sq < k_sq_max) && (k_sq != 0)) # Test to ensure within range
                    NKVECS += 1
                end # End test to ensure within range
            end
        end
    end
    kxyz = Vector{SVector{3,Int32}}(undef, NKVECS)
    cfac = Vector{Float64}(undef, NKVECS)
    #cfac=0; kxx = 0; kyy=0;kzz=0
    NKVECS = 0
    for kx = 0:nk
        for ky = -nk:nk
            for kz = -nk:nk
                k_sq = kx^2 + ky^2 + kz^2

                if ((k_sq < k_sq_max) && (k_sq != 0)) # Test to ensure within range
                    NKVECS += 1
                    kxyz[NKVECS] = SVector{3,Int32}(kx, ky, kz)
                    kr_sq = twopi_sq * float(k_sq)           # k**2 in real units
                    cfac[NKVECS] = twopi * exp(-b * kr_sq) / kr_sq / box# Stored expression for later use
                    if kx > 0
                        cfac[NKVECS] = cfac[NKVECS] * 2.0
                    end

                end # End test to ensure within range

            end # kz
        end # ky
    end # kx
    #ewald = EWALD(ewald.kappa, ewald.nk, ewald.k_sq_max, cfac, fact)
    ewald = EWALD(
        ewald.kappa,
        ewald.nk,
        ewald.k_sq_max,
        NKVECS,
        kxyz,
        cfac,
        zeros(ComplexF64, NKVECS),    # dummy values
        zeros(ComplexF64, NKVECS),
        fact,
    )
    return ewald #cfac, kxyz, ewald
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
function SetupKVecs(ewald::EWALD, boxSize::Real where {T})

    kappa = ewald.kappa
    nk = ewald.nk
    k_sq_max = ewald.k_sq_max
    fact = ewald.factor


    #k_sq_max = nk^2 + 2
    println("k_sq_max: ", k_sq_max)
    #kfac = zeros(SVector{k_sq_max + 2, Float64})
    kfacTemp = Vector{Float64}(undef, k_sq_max)
    b = 1.0 / 4.0 / kappa / kappa / boxSize / boxSize  # boxSize may be issue
    #b = 1.0 / 4.0 / kappa / kappa

    @inbounds for kx = 0:(nk)
        for ky = 0:(nk)
            for kz = 0:(nk)
                k_sq = kx^2 + ky^2 + kz^2
                if k_sq <= k_sq_max && k_sq > 0
                    kr_sq = (2 * pi)^2 * float(k_sq)
                    kfacTemp[k_sq] = 2 * pi * exp(-b * kr_sq) / kr_sq / boxSize
                    #println(i," ", j," ", k, " ", kfacTemp[k_sq])
                end
            end
        end
    end
    kfac = [kfacTemp[i] for i = 1:length(kfacTemp)]
    ewald = EWALD(ewald.kappa, ewald.nk, ewald.k_sq_max, kfac, fact)

    return kfac, ewald
end

#=
"Real Ewald contribution - without molecular cutoff"
#function EwaldReal(diff,coord1, coord2,q1,q2, L, rcut2, kappa)
@fastmath function EwaldReal(
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    kappa::Real,
    box::Real,
    thisMol_thisAtom::Vector,
    chosenOne::Int,
    r_cut::Real,
)
    ####
    #
    #    Some prep stuff
    #
    #############
    start_i = thisMol_thisAtom[chosenOne][1]
    end_i = thisMol_thisAtom[chosenOne][2]
    #r_cut = 10.0
    r_cut_sq = r_cut^2
    pot = 0.0
    rij = SVector{3}(0.0, 0.0, 0.0)
    rᵢ = similar(rij)
    rⱼ = similar(rij)



    """Loop over all atoms in molecule A"""
    @inbounds for i = start_i:end_i
        rᵢ = qq_r[i]
        """Loop over all atoms in molecule B"""
        @inbounds for (j, rⱼ) in enumerate(qq_r)
            if j in start_i:end_i
                continue
            end # keep mol i from self interacting
            @inbounds for k = 1:3

                rij = @set rij[k] = vector1D(rᵢ[k], rⱼ[k], box)
            end
            rᵢⱼ² = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

            if rᵢⱼ² < r_cut_sq
                rᵢⱼ = sqrt(rᵢⱼ²)
                pot += qq_q[i] * qq_q[j] * erfc(kappa * rᵢⱼ) / rᵢⱼ

            else
                pot += 0.0
            end # potential cutoff
        end # loop over atom in molecule b
    end # loop over atoms in molecule a
    return pot
end
=#
"Real Ewald contribution with molecular cutoff"
#function EwaldReal(diff,coord1, coord2,q1,q2, L, rcut2, kappa)
function EwaldReal(
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    kappa::Real,
    box::Float64,
    thisMol_thisAtom::Vector{SVector{2,Int64}},
    chosenOne::Int64,
    system::Requirements,
)
    ####
    #
    #    Some prep stuff
    #
    #############
    ri = system.rm[chosenOne]
    start_a = thisMol_thisAtom[chosenOne][1]
    end_a = thisMol_thisAtom[chosenOne][2]
    rMol = system.rm
    #r_cut = 10.0
    r_cut = system.r_cut
    #@assert r_cut == 10.0
    diameter = 0 #r_cut * 0.25  + 5.0#2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq = r_cut^2
    #@assert r_cut_sq == 100.0
    pot = float(0.0)
    rab = SVector{3,Float64}(0.0, 0.0, 0.0)
    rij = SVector{3,Float64}(0.0, 0.0, 0.0)
    #ri = similar(rab)
    rj = SVector{3,Float64}(0.0, 0.0, 0.0)
    ra = SVector{3,Float64}(0.0, 0.0, 0.0)
    rb = SVector{3,Float64}(0.0, 0.0, 0.0)

    overlap = false
    ovr = 1.0

    @inbounds for (j, rj) in enumerate(rMol) # cycle through all molecules
        if j == chosenOne # skip self interactions
            continue
        end

        @inbounds for k = 1:3
            rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij2 = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

        if rij2 < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
            @inbounds for a = start_a:end_a
                ra = qq_r[a]
                """Loop over all atoms in molecule B"""
                #@inbounds for (b, rb) in enumerate(qq_r)
                start_b = thisMol_thisAtom[j][1]
                end_b = thisMol_thisAtom[j][2]
                @inbounds for b = start_b:end_b
                    rb = qq_r[b]
                    #if b in start_a:end_a
                    #    continue
                    #end # keep mol i from self interacting
                    @inbounds for k = 1:3
                        rab = @set rab[k] = vector1D(ra[k], rb[k], box)
                    end
                    rab2 = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    if (rab2 < ovr) && (qq_q[a] * qq_q[b] < 0)
                        return 0.0, true

                    elseif rab2 < (r_cut_sq + 100)
                        # this uses only molecular cutoff, hence the +100
                        # TODO remove this if statement
                        rab_mag = sqrt(rab2)
                        pot +=
                            qq_q[a] * qq_q[b] * erfc(kappa * rab_mag) / rab_mag
                    else
                        pot += 0.0
                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end # if statement for molecular cutoff
    end # loop over molecules
    return pot, overlap
end

##########################################

function EwaldReal(chosenOne::Int64,
                    moa::StructArray,
                    soa::StructArray,
                    ewald::EWALD,
                    r_cut::Float64,
                    box::Float64
    )
    ####
    #
    #    Some prep stuff
    #
    #############
    ri = moa.COM[chosenOne]
    start_a = moa.firstAtom[chosenOne]
    end_a = moa.lastAtom[chosenOne]
    rMol = moa.COM
    #r_cut = 10.0
    r_cut = r_cut
    #@assert r_cut == 10.0
    diameter = 0 #r_cut * 0.25  + 5.0#2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq = r_cut^2
    #@assert r_cut_sq == 100.0
    pot = float(0.0)
    kappa = ewald.kappa
    rab = SVector{3,Float64}(0.0, 0.0, 0.0)
    rij = SVector{3,Float64}(0.0, 0.0, 0.0)
    #ri = similar(rab)
    rj = SVector{3,Float64}(0.0, 0.0, 0.0)
    ra = SVector{3,Float64}(0.0, 0.0, 0.0)
    rb = SVector{3,Float64}(0.0, 0.0, 0.0)

    overlap = false
    ovr = 0.5

    @inbounds for (j, rj) in enumerate(rMol) # cycle through all molecules
        if j == chosenOne # skip self interactions
            continue
        end

        @inbounds for k = 1:3
            rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij2 = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

        if rij2 < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
            for a = start_a:end_a
                ra = soa.coords[a]
                """Loop over all atoms in molecule B"""
                #@inbounds for (b, rb) in enumerate(qq_r)
                start_b = moa.firstAtom[j]
                end_b = moa.lastAtom[j]
                for b = start_b:end_b
                    rb = soa.coords[b]
                    #if b in start_a:end_a
                    #    continue
                    #end # keep mol i from self interacting
                    @inbounds for k = 1:3
                        rab = @set rab[k] = vector1D(ra[k], rb[k], box)
                    end
                    rab2 = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    if (rab2 < ovr) && (soa.charge[a] * soa.charge[b] < 0)
                        return 0.0, true

                    elseif rab2 < (r_cut_sq + 100)
                        # this uses only molecular cutoff, hence the +100
                        # TODO remove this if statement
                        rab_mag = sqrt(rab2)
                        pot +=
                            soa.charge[a] * soa.charge[b] * erfc(kappa * rab_mag) / rab_mag
                    else
                        pot += 0.0
                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end # if statement for molecular cutoff
    end # loop over molecules
    return pot, overlap
end


#function RecipLong(boxSize, n, kappa, kfac, nk, k_sq_max, r, qq_q)
""" Recipricol space Ewald energy """
function RecipLong(
    system::Requirements,
    ewald::EWALD,
    r::Vector,
    qq_q::Vector{Float64},
    kfac::Vector,
)
    energy = 0.0
    L = system.box
    twopi = 2.0 * π
    nk = ewald.nk
    n = length(r)
    k_sq_max = ewald.k_sq_max

    eikx = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)
    eikz = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)

    @inbounds @simd for j = 1:n
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
        @inbounds @simd for j = 1:n
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
    # TODO use other code which is single loop and @simd it
    # TODO only evaluate moved charges - huge speedup will occur from this
    # ~ 50x increase in speed for a 300 spce water system. Bigger speedups within
    # larger systems.
    exit()
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
                    @inbounds @simd for l = 1:n
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

""" Recipricol space Ewald energy parallel and molecular version"""
function RecipLong(
    system::Requirements,
    ewald::EWALD,
    r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
)
    #=
    Intended to calculate the recipricol for the molecule that was moved, not then
    entire system.
    =#
    #cfac = ewald.cfac
    kxyz = ewald.kxyz[:]
    energy = 0.0
    L = system.box
    twopi = 2.0 * π
    nk = ewald.nk
    n = length(r)
    k_sq_max = ewald.k_sq_max

    eikx = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)
    eikz = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)

    for j = 1:n
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

        eiky[j, -1] = conj(eiky[j, 1])
        eikz[j, -1] = conj(eikz[j, 1])
    end
    # calculate remaining positive kx, ky, kz by recurrence
    for k = 2:(nk)
        for j = 1:n
            eikx[j, k] = eikx[j, k-1] * eikx[j, 1]
            eiky[j, k] = eiky[j, k-1] * eiky[j, 1]
            eikz[j, k] = eikz[j, k-1] * eikz[j, 1]

            #eikx[j,-k] = conj( eikx[j,k] )
            eiky[j, -k] = conj(eiky[j, k])
            eikz[j, -k] = conj(eikz[j, k])
        end
    end
    # TODO This can be sped up. It is generating too many temp Arrays
    # It should be faster if devectorized, but that makes even more temp arrays
    term = 0.0 + 0.0 * im
    for i = 1:length(ewald.cfac)
        term = 0.0 + 0.0 * im
        for l = 1:n
            term +=
                qq_q[l] *
                eikx[l, kxyz[i][1]] *
                eiky[l, kxyz[i][2]] *
                eikz[l, kxyz[i][3]]#* eiky[l,kxyz[i][2]] * eikz[l,kxyz[i][3]]
        end
        #term = sum( qq_q[:] .* eikx[:,kxyz[i][1]] .* eiky[:,kxyz[i][2]] .* eikz[:,kxyz[i][3]] )
        energy += ewald.cfac[i] * real(conj(term) * term)
        ewald.sumQExpNew[i] = term
        ewald.sumQExpOld[i] = term
    end
    return energy, ewald
end # function


""" Recipricol space Ewald energy Move"""
function RecipMove(
    system::Requirements,
    ewalds::EWALD,
    r_old::Vector,
    r_new::Vector,
    qq_q::Vector,
)
    #=
    Intended to calculate the recipricol for the molecule that was moved, not then
    entire system.

    This is twice as long as it should be.
    It calculates the fourier contibution for both old and new positions.
    TODO store old contributions
    =#
    #cfac = ewald.cfac
    kxyz = ewalds.kxyz[:]
    energy = float(0.0)
    L = system.box
    twopi = 2.0 * π
    nk = ewalds.nk
    n = length(r_old)
    @assert n == 3
    k_sq_max = ewalds.k_sq_max
    @assert k_sq_max == 27
    @assert nk == 5

    eikx_n = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky_n = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)

    eikz_n = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #similar(eiky_n) #  OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)
    eikx_o = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #similar(eikx_n) #  OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky_o = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #similar(eiky_n) # OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)
    eikz_o = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #similar(eiky_n) # OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)

    for j = 1:n
        # calculate kx,ky,kz =1,2 explicitly
        eikx_n[j, 0] = 1.0 + 0.0 * im
        eiky_n[j, 0] = 1.0 + 0.0 * im
        eikz_n[j, 0] = 1.0 + 0.0 * im
        eikx_o[j, 0] = 1.0 + 0.0 * im
        eiky_o[j, 0] = 1.0 + 0.0 * im
        eikz_o[j, 0] = 1.0 + 0.0 * im

        eikx_n[j, 1] =
            cos(twopi * (r_new[j][1]) / L) + sin(twopi * (r_new[j][1]) / L) * im
        eiky_n[j, 1] =
            cos(twopi * (r_new[j][2]) / L) + sin(twopi * (r_new[j][2]) / L) * im
        eikz_n[j, 1] =
            cos(twopi * (r_new[j][3]) / L) + sin(twopi * (r_new[j][3]) / L) * im
        eikx_o[j, 1] =
            cos(twopi * (r_old[j][1]) / L) + sin(twopi * (r_old[j][1]) / L) * im
        eiky_o[j, 1] =
            cos(twopi * (r_old[j][2]) / L) + sin(twopi * (r_old[j][2]) / L) * im
        eikz_o[j, 1] =
            cos(twopi * (r_old[j][3]) / L) + sin(twopi * (r_old[j][3]) / L) * im

        eiky_n[j, -1] = conj(eiky_n[j, 1])
        eikz_n[j, -1] = conj(eikz_n[j, 1])
        eiky_o[j, -1] = conj(eiky_o[j, 1])
        eikz_o[j, -1] = conj(eikz_o[j, 1])
    end
    # calculate remaining positive kx, ky, kz by recurrence
    for k = 2:(nk)
        for j = 1:n
            eikx_n[j, k] = eikx_n[j, k-1] * eikx_n[j, 1]
            eiky_n[j, k] = eiky_n[j, k-1] * eiky_n[j, 1]
            eikz_n[j, k] = eikz_n[j, k-1] * eikz_n[j, 1]
            eikx_o[j, k] = eikx_o[j, k-1] * eikx_o[j, 1]
            eiky_o[j, k] = eiky_o[j, k-1] * eiky_o[j, 1]
            eikz_o[j, k] = eikz_o[j, k-1] * eikz_o[j, 1]

            #eikx[j,-k] = conj( eikx[j,k] )
            eiky_n[j, -k] = conj(eiky_n[j, k])
            eikz_n[j, -k] = conj(eikz_n[j, k])
            eiky_o[j, -k] = conj(eiky_o[j, k])
            eikz_o[j, -k] = conj(eikz_o[j, k])
        end
    end
    # TODO This can be sped up. It is generating too many temp Arrays
    # It should be faster if devectorized, but that makes even more temp arrays
    # TODO can probably vectorize better if I swap i and l loops.

    # Calculate difference between the twopi

    for i = 1:length(ewalds.cfac)
        for l = 1:n
            ewalds.sumQExpNew[i] +=
                qq_q[l] * (
                    eikx_n[l, kxyz[i][1]] *
                    eiky_n[l, kxyz[i][2]] *
                    eikz_n[l, kxyz[i][3]] -
                    eikx_o[l, kxyz[i][1]] *
                    eiky_o[l, kxyz[i][2]] *
                    eikz_o[l, kxyz[i][3]]
                )
        end


        energy +=
            ewalds.cfac[i] * (
                real(conj(ewalds.sumQExpNew[i]) * ewalds.sumQExpNew[i]) -
                real(conj(ewalds.sumQExpOld[i]) * ewalds.sumQExpOld[i])
            )
    end

    #
    return energy * ewalds.factor, ewalds
end # function Recipricol Move

function EwaldSelf(ewald::EWALD, qq_q::Vector)
    kappa = ewald.kappa
    factor = ewald.factor
    return self = -kappa * sum(qq_q .^ 2) / sqrt(π) * factor
end

"""Calculates the surface term to remove tinfoil boundary assumption"""
function TinfoilBoundary(
    system::Requirements,
    ewald::EWALD,
    qq_q::Vector,
    qq_r::Vector,
)
    vol = system.box^3
    return 2 * π / 3 / vol * dot(qq_q .* qq_r, qq_q .* qq_r)
end

"""Calculates real, self and tinfoil contributions to Coulombic energy for
a single molecule interacting with the rest of the system"""
function EwaldShort(
    i::Int64,
    system::Requirements,
    ewald::EWALD,
    box::Float64,
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    tinfoil = false,
)
    partial_e = 0.0
    partial_v = 0.0
    overlap = false
    # Calculate new ewald REAL energy
    realEwald, overlap = EwaldReal(
        qq_r,
        qq_q,
        ewald.kappa,
        box,
        system.thisMol_theseAtoms,
        i,
        system,
    )
    realEwald *= ewald.factor
    partial_e += realEwald
    partial_v += (realEwald / 3)

    # Calculate Self interaction energy
    #selfEnergy = EwaldSelf(ewald, qq_q)
    #partial_e += selfEnergy
    #partial_v += (selfEnergy / 3)

    #Calculate tinfoil Boundaries
    #=
    if tinfoil
        tfb = TinfoilBoundary(system, ewald, qq_q, qq_r) * ewald.factor
        partial_e += tfb
        partial_e += (tfb / 3)
    end
    =#
    return partial_e, partial_v, overlap
end

"""Calculates real, self and tinfoil contributions to Coulombic energy for
a single molecule interacting with the rest of the system"""
function EwaldShort(
    i::Int64,
    moa::StructArray,
    soa::StructArray,
    sim_props::Properties2,
    ewald::EWALD,
    box::Float64,
)
    partial_e = 0.0
    partial_v = 0.0
    overlap = false
    # Calculate new ewald REAL energy
    realEwald, overlap = EwaldReal(i, moa, soa, ewald,sim_props.qq_rcut, box)
    realEwald *= ewald.factor
    partial_e += realEwald
    partial_v += (realEwald / 3)

    return partial_e, partial_v, overlap
end
