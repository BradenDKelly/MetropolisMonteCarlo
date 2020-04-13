
using Distributions
using Random
using Base.Threads
using LinearAlgebra: norm, normalize, dot, ×
using BenchmarkTools
using Setfield
using StaticArrays
using StructArrays
using SpecialFunctions
using OffsetArrays
using Printf

import JSON

Random.seed!(11234)

################################################################################
# TODO (BDK) create JSON input file with starting parameters:
# TODO (BDK) add proper sampling
# TODO (BDK) Add some timing (~ 30 times faster than numpy)
# until then, manually enter at top
################################################################################
temperature = 1.0 #0.8772  # 1.2996
ρ = 0.75
nAtoms = 256
#ϵ = 1.0
#σ = 1.0
r_cut = 2.5  # box / 2
nSteps = 1000
nblock = 100
outputInterval = 100
initialConfiguration = "crystal"  # place atoms in a crystal structure

# Units
const na = 6.02214129e23      # mol^-1
const R = 8.3144621e-3        # kJ mol^-1 K^-1
const kb = R / na             # kJ K^-1
const h = 0.399031271         # kJ mol^-1 ps
const hbar = 0.0635077993     # kJ mol^-1 ps
const c = 299792.458          # nm ps^-1
const e = 1.6021765653e-19     # C
ϵ₀ = 8.854187817e-12 # C²/J/m
ϵ₀ *= 1e-9  # C²/(J ⋅ nm)
ϵ₀ *= 1000   # C²/(kJ ⋅ nm)
const qq_convert = 138.935458   # kJ mol^-1 nm e^-2 # 1 / (4πϵ₀)
factor2 =  e^2 / ϵ₀ / 4 / pi / kb
conv = e^2 / ϵ₀ * na / 4 / pi

kb1 = 1.3806488e-23 # J/K
ϵ₀1 = 8.854187817e-12 # C²/J/m
ϵ₀1 *= 1e-10
e1 = 1.602176565e-19
factor = e1^2 / ϵ₀1 / 4 / pi / kb1

# Set default values, check keys and typecheck values
defaults = Dict("nblock"=> 10, "nstep"=> 1000, "temperature"=> 1.0, "r_cut"=> 2.5, "dr_max"=> 0.15,
            "natoms"=> 256, "initConfig"=> "crystal", "rho"=>0.75, "ϵ"=> 1.0,
             "σ"=>1.0)

""" Generates .pdb files """
function PrintPDB(r, box, step=1, filename="pdbOutput")

    open(filename * "_" * string(step) * ".pdb", "w") do file

        line = @sprintf("%-7s %7.3f %7.3f %7.3f %30s \n","CRYST1", 10.0 * box,
               10.0 * box, 10.0 * box, "90.00  90.00  90.00 P 1           1")
        write(file,line)

        for (i,atom) in enumerate(r)

            atomName = "Ar"
            molName = "Ar"
            #atomName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].atomnm
            #molName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].resnm

            line = @sprintf("%-6s %4d %3s %4s %5d %3s %7.3f %7.3f %7.3f %5.2f %5.2f \n",
            "ATOM",i, atomName, molName, i, " ", 10.0 * r[i][1],
            10.0 * r[i][2],10.0 * r[i][3], 1.00, 0.00   )
            write(file,line)
        end
    end
    println("finished making pdb")
end

function InitCubicGrid(n::Int,rho::Real)

#!------------------------------------------------------------------------
#! Created by Braden Kelly
#!------------------------------------------------------------------------
#! Creates an initial configuration
#!------------------------------------------------------------------------
#! input:  n       number of particles
#!         rho     density
#! output: coords  coordinates of n particles
#!------------------------------------------------------------------------

        #! Calculate box length (L)
        L = (n/rho)^(1.0/3.0)

        #! Calculate the lowest perfect cube that will contain all of the particles
        nCube = 2

        while (nCube^3 < n)
           nCube = nCube+1
        end
        coords = Array{Float64,2}(undef, 3, n)
        # initial position of particle 1
        posit=[0, 0, 0]

        #!begin assigning particle positions
        for i=1:n
           coords[:,i] = (posit+[0.01,0.01,0.01])*(L/nCube)
           #! Advancing the index (posit)
                posit[1]=posit[1]+1
                if (posit[1] == nCube)
                   posit[1] = 0
                   posit[2] = posit[2]+1
                   if (posit[2] == nCube)
                      posit[2] = 0
                      posit[3] = posit[3] + 1
                   end
                end
        end
    return [SVector{3}(coords[1,i],coords[2,i],coords[3,i]) for i = 1:size(coords,2)]
end #function

function potential_lrc( ρ, r_cut )
    """Calculates long-range correction for Lennard-Jones potential per atom."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * ( (8.0/9.0)  * sr3^3  - (8.0/3.0)  * sr3 ) * ρ
end

function pressure_lrc( ρ, r_cut )
    """Calculates long-range correction for Lennard-Jones pressure."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut ^ 3
    return π * ( (32.0/9.0) * sr3 ^ 3  - (16.0/3.0) * sr3 ) * ρ ^ 2
end

function pressure_delta( ρ, r_cut )
    """Calculates correction for Lennard-Jones pressure due to discontinuity in the potential at r_cut."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * (8.0/3.0) * ( sr3 ^3  - sr3 ) * ρ ^ 2
end

mutable struct Properties
    energy::Float64
    virial::Float64
    old_e::Float64
    old_v::Float64
end

"""structure for energy calculation properties """
mutable struct Requirements
    r::Vector{SVector{3,Float64}}
    ϵ::Vector{Float64}
    σ::Vector{Float64}
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
end
function PBC(in::SVector{3,Float64}, box::Real)
    x, y, z = in[1], in[2], in[3]
    if x > box x -= box end
    if x < 0 x += box end
    if y > box y -= box end
    if y < 0 y += box end
    if z > box z -= box end
    if z < 0 z += box end

    return SVector{3,Float64}(x, y, z)
end

function random_translate_vector( dr_max::Float64, old::SVector, box::Float64 )

# In: dr_max, old
# Out: SVector{3,Float64}

  zeta = rand(Float64,3)    # Three uniform random numbers in range (0,1)
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
    return ρ * T + vir.virial / vol + pressure_lrc( ρ, r_cut )
end

"""Calculates pressure without tail correction"""
function Pressure(vir, ρ, T, vol)
    return ρ * T + vir.virial / vol
end

"Vector between two coordinate values, accounting for mirror image seperation"
@inline function vector1D(c1::Float64, c2::Float64, box_size::Float64)
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) : (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) : (c2 - c1 + box_size)
    end
end

function LJ_ΔU(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    r = system.r
    ϵ = system.ϵ
    σ = system.σ
    box = system.box
    diff = SVector{3}(0.0,0.0,0.0)
    rcut_sq = system.r_cut^2
    pot, vir = 0.0, 0.0

    for (j, atom) in enumerate(r)

        if j == i continue end

        @inbounds for k=1:3
             diff = @set diff[k] = vector1D(r[i][k], atom[k], box)
        end

        rij_sq = diff[1]*diff[1] + diff[2]*diff[2] + diff[3]*diff[3]

        if  rij_sq > rcut_sq
            pot += 0.0
            vir += 0.0
        else
            # possibly needed in the case of random starting config to avoid
            # numerical overflow from overlapped particles. Not needed
            # if starting from lattice.
            #if rij_sq < 0.6^2
            #    rij_sq = 0.6^2
            #end

            sr2 = σ[j]^2 / rij_sq
            sr6 = sr2 ^ 3
            sr12 = sr6 ^ 2
            pot += ϵ[j] * (sr12 - sr6)  # LJ pair potentials (cut but not shifted)
            vir += ϵ[j] * (2*sr12 - sr6)  # LJ pair virials
        end

    end

    return pot * 4.0, vir * 24.0 / 3.0
end

""" Calculates total potential energy of the system. Double Counts. """
function potential(system, tot)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0

    for i = 1:length(system.r)
        ener, vir = LJ_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    return tot

end

"""unit test of lj"""
function test_LJ()
    # this tests the LJ potential and indirectly the mirror image Seperation
    #
    # the point is to calculate the energy of 1 particle interacting with
    # 2 other particles at known positions.
    # Part 1: tests the case where all particles are within cutoff
    # Part 2: tests the case where 1 particle is outside cutoff
    # Part 2: tests mirror image seperation since particle outside cutoff
    #        has its mirror image inside the cutoff.

    #### Part 1
    box_test = 5.0
    r_cut = box_test / 2
    r = [SVector{3, Float64}(0,0,0),SVector{3, Float64}(0,0,2),SVector{3, Float64}(0,1.5,0)]
    test = Requirements(r, ones(3), ones(3), box_test, r_cut)
    enn, virr = LJ_ΔU(1,test)
    println("Testing 3 particles (0,0,0), (0,0,2), (0,1.5,0) where rcut is ",r_cut," None are outside of rcut")
    println("Calculated: ", enn, " manual: ", LennardJones(2.0)+LennardJones(1.5) )

    #### Part 2:
    r = [SVector{3, Float64}(0,0,0),SVector{3, Float64}(0,0,4),SVector{3, Float64}(0,1.5,0)]
    test = Requirements(r, ones(3), ones(3), box_test, r_cut)
    enn, virr = LJ_ΔU(1,test)
    lj1 = LennardJones(4.0)+LennardJones(1.5)
    lj2 = LennardJones(1.0)+LennardJones(1.5)
    println("Testing 3 particles (0,0,0), (0,0,4), (0,1.5,0) where rcut is ",r_cut,"  One  is outside of rcut")
    println("Calculated: ", enn, " incorrect manual: ", lj1 )
    println("Calculated: ", enn, " correct manual: ", lj2 )
    if abs(lj2 - enn) < 0.001
        println("Energy is correct in unit test.")
        println("Mirror Image Seperation passed this simple test.")
    end

end

function LennardJones(rij)
    return 4 * 1 * ( (1/rij)^12 - (1/rij)^6 )
end
test_LJ()

""" returns max and min values in an array of SVectors"""
function maxmin(array::Vector)
    maxVal = maximum(array[1])
    minVal = minimum(array[1])

    for svector in array
        if maximum(svector) > maxVal maxVal = maximum(svector) end
        if minimum(svector) < minVal minVal = minimum(svector) end
    end

    return minVal, maxVal
end

################################################################################
#
#                        Start of Configuration
#
################################################################################
box = (nAtoms / ρ ) ^ (1 / 3)
dr_max = box / 30   # at 256 particles, ρ=0.75, T=1.0 this is 48% acceptance

if lowercase(initialConfiguration) == "crystal"
    r = InitCubicGrid(nAtoms,ρ)
else
    r = [SVector{3,Float64}(rand(), rand(), rand()) .* box for i = 1:nAtoms]
end

PrintPDB(r, box, 1, "pdbOutput")
ϵ = ones(nAtoms)
σ = ones(nAtoms)

total     = Properties(0.0, 0.0, 0.0, 0.0)
system    = Requirements(r, ϵ, σ, box, r_cut)
total     = potential(system, Properties(0.0,0.0, 0.0, 0.0))
averages  = Properties(total.energy, total.virial,total.energy, total.virial) # initialize struct with averages
totProps  = Properties2(temperature, ρ, Pressure(total, ρ, temperature, box^3),
                            dr_max, 0.0, 0.3, 0, 0)

println("TEST ", total.energy, " ", totProps.pressure)

#""" Main loop of simulation. Sweep over all atoms, then Steps, then Blocks."""
for blk = 1:nblock
    for step =1:nSteps
        for i = 1:nAtoms

            partial_old_e, partial_old_v = LJ_ΔU(i, system)
            rold = deepcopy(system.r[i])
            rnew = random_translate_vector(totProps.dr_max, system.r[i], box)
            system.r[i] = rnew
            partial_new_e, partial_new_v = LJ_ΔU(i, system)

            delta = partial_new_e - partial_old_e

            if Metropolis(delta/temperature)
                total.energy += delta
                total.virial += (partial_new_v - partial_old_v)
                totProps.numTranAccepted += 1
                ne = averages.old_e + delta
                nv = averages.old_v + partial_new_v - partial_old_v
                averages.energy += ne
                averages.virial += nv
                averages.old_e = ne
                averages.old_v = nv
            else
                system.r[i] = rold
                averages.energy += averages.old_e
                averages.virial += averages.old_v
            end

            # for troubleshooting checks that particles are in box
            minV, maxV = maxmin(system.r)
            #println(minV, maxV)
            if minV < 0.0
                println("Shit, particle is less than 0")
            end
            if maxV > box
                println("Shit, particle is outside box")
            end

            totProps.totalStepsTaken += 1

        end # i to nAtoms


    end # step to nSteps
    #PrintPDB(system.r, box, blk, "pdbOutput")
    # Hella ugly output
    # TODO (BDK) modify to formatted output
    println(blk, " ", averages.energy / totProps.totalStepsTaken / nAtoms,
         " ", totProps.numTranAccepted / totProps.totalStepsTaken,
         "   ", ρ*temperature + averages.virial / box^3 / totProps.totalStepsTaken, #  Pressure(total, ρ, temperature, box^3),
         " Pcut: ", pressure_lrc( ρ, system.r_cut ), " Ecorr: ",
         potential_lrc( ρ, r_cut ), "p delta: ", pressure_delta(ρ,system.r_cut ),
         " A & T p_c: " , ρ*temperature + averages.virial / box^3 /
         totProps.totalStepsTaken + pressure_delta(ρ,system.r_cut ) )
end # blk to nblock

