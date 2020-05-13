#####
#
#            Hodgepodge of stuff
#
##########
using StaticArrays

abstract type ForceField end
abstract type Gromacs <: ForceField end
abstract type GAFF <: ForceField end
abstract type AngleType <: ForceField end
abstract type DihedralType <: ForceField end

function potential_lrc(ρ, r_cut)
    """Calculates long-range correction for Lennard-Jones potential per atom."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * ((8.0 / 9.0) * sr3^3 - (8.0 / 3.0) * sr3) * ρ
end

function pressure_lrc(ρ, r_cut)
    """Calculates long-range correction for Lennard-Jones pressure."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * ((32.0 / 9.0) * sr3^3 - (16.0 / 3.0) * sr3) * ρ^2
end

function pressure_delta(ρ, r_cut)
    """Calculates correction for Lennard-Jones pressure due to discontinuity in the potential at r_cut."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * (8.0 / 3.0) * (sr3^3 - sr3) * ρ^2
end

mutable struct Properties
    energy::Float64
    virial::Float64
    coulomb::Float64
    recip::Float64
    recipOld::Float64
    old_e::Float64
    old_v::Float64
end

"""Struct for tracking and optimizing translation moves"""
mutable struct Moves
    naccepp::Int           # previous number of attempts
    naccept::Int           # current number of attempts
    attempp::Int
    attempt::Int
    set_value::Float64     # desired success rate
    dr_max::Float64
end


"""Struct for FF parameters (currently LJ and EXP6) """
struct Tables #{T<:Vector} #<: ForceField
    # passed two 1D arrays, convert them both to 2D matrices and
    # apply geometric and arithmetic mixing rules
    ϵᵢⱼ::Array{Float64,2}
    σᵢⱼ::Array{Float64,2}
    function Tables(a::Vector{T}, b::Vector{T}) where {T}
        e = [sqrt(a[i] * a[j]) for i = 1:length(a), j = 1:length(a)]
        s = [(b[i] + b[j]) / 2 for i = 1:length(b), j = 1:length(b)]
        new(e, s)
    end
end

"""
function Tables(a::Vector{T}, b::Vector{T}) where T
    e = [sqrt(a[i]*a[j]) for i=1:length(a), j=1:length(a)]
    s = [(b[i] + b[j]) / 2 for i=1:length(b), j=1:length(b)]
    return Tables(e,s)
end
"""

"""structure for energy calculation properties """

mutable struct Requirements
    rm::Vector{SVector{3,Float64}}
    ra::Vector{SVector{3,Float64}}
    nMols::Int
    nAtoms::Int
    nCharges::Int
    thisMol_theseAtoms::Vector{SVector{2,Int64}}
    molNames::Vector     # strings
    molTypes::Vector     # integers
    atomNames::Vector    # strings
    atomTypes::Vector    # integers
    #ϵ::Array{Float64,2}
    #σ::Array{Float64,2}
    table::Tables
    box::Float64
    r_cut::Float64
end

""" structure for generic simulation properties"""
mutable struct Properties2
    temperature::Float64
    ρ::Float64
    pressure::Float64
    dr_max::Float64
    dϕ_max::Float64
    move_accept::Float64
    numTranAccepted::Int
    totalStepsTaken::Int
    quat::Vector{SVector{4,Float64}}
end


function random_translate_vector(dr_max::Float64, old::SVector, box::Float64)

    # In: dr_max, old
    # Out: SVector{3,Float64}

    zeta = rand(Float64, 3)    # Three uniform random numbers in range (0,1)
    zeta = zeta .- 0.5         #! Now in range (-1,+1)
    return PBC(old + zeta .* dr_max, box) #! Move to new position

end #random_translate_vector

"""Evaluates if we keep a translation move or not"""
function Metropolis(delta)
    if delta < 0.0
        return true
    elseif exp(-delta) > rand()
        return true
    else
        return false
    end
end
""" Calculates pressure including the tail correction"""
function Pressure(vir, ρ, T, vol, r_cut)
    return ρ * T + vir.virial / vol + pressure_lrc(ρ, r_cut)
end

"""Calculates pressure without tail correction"""
function Pressure(vir, ρ, T, vol)
    return ρ * T + vir.virial / vol
end

#test_LJ()

""" returns max and min values in an array of SVectors"""
function maxmin(array::Vector)
    maxVal = maximum(array[1])
    minVal = minimum(array[1])

    for svector in array
        if maximum(svector) > maxVal
            maxVal = maximum(svector)
        end
        if minimum(svector) < minVal
            minVal = minimum(svector)
        end
    end

    return minVal, maxVal
end

"""Returns the center of mass of an array of atoms and array of masses"""
function COM(atoms::Vector, masses)
    totalMass = sum(masses)
    numerator = sum(atoms .* masses)
    #println("From inside COM: ", numerator ./ totalMass)
    return numerator ./ totalMass
end
#COM([[1,2,3],[2,3,4],[0,1,2]],[1,1,100])

"""Returns the Fortran equivalent of the same name"""
function MATMUL(ai, db)
    # In: ai:3x3 SMatrix
    #     db: 3x1 vector
    # Out: SVector{3}
    return SVector(dot(db, ai[:, 1]), dot(db, ai[:, 2]), dot(db, ai[:, 3]))
end
