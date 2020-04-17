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
using Dates

import JSON

include("tests.jl")
include("auxillary.jl")
include("initialConfigurations.jl")
include("boundaries.jl")
include("energy.jl")
include("quaternions.jl")
include("volumeChange.jl")

Random.seed!(11234)
start=Dates.now()
println(Dates.now())

################################################################################
# TODO (BDK) create JSON input file with starting parameters:
# TODO (BDK) add proper sampling
# TODO (BDK) Add some timing (~ 30 times faster than numpy)
# TODO (BDK) write more unit tests
# TODO (BDK) extend to mixtures
# TODO (BDK) add ewalds (routines written, need to incorporate)
# TODO (BDK) make tables for parameters
# until then, manually enter at top
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# Things to organize for EWALD 
# 1) prefactor: need k-vectors - this needs prep (code written already)
# 2) need array for charges, their coords, atomID, molID, k-fac max, alpha
# 3) need real, recipricol, self and tinfoil boundary routines (code already written)
# 4) test vs. NIST SPC/E database:
#   https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10a-cutoff
# 5) DOUBLE CHECK ALL FREAKING UNITS!!!
################################################################################
temperature = 0.6 #0.8772  # 1.2996
ρ = 0.30533
nMol = 256 * 4
nAtoms = nMol * 3
#ϵ = 1.0
#σ = 1.0
r_cut = 2.5  # box / 2
nSteps = 2000
nblock = 100
outputInterval = 100
initialConfiguration = "cnf"  # place atoms in a crystal structure
dϕ_max = 0.05
dr_max = 0.05

# Set default values, check keys and typecheck values
defaults = Dict("nblock"=> 10, "nstep"=> 1000, "temperature"=> 1.0, "r_cut"=> 2.5, "dr_max"=> 0.15,
            "natoms"=> 256, "initConfig"=> "crystal", "rho"=>0.75, "ϵ"=> 1.0,
             "σ"=>1.0)

################################################################################
#
#                        Start of Configuration
#
################################################################################
box = (nMol / ρ ) ^ (1 / 3)
 #box / 30   # at 256 particles, ρ=0.75, T=1.0 this is 48% acceptance

# Generate molecular COM coordinates
if lowercase(initialConfiguration) == "crystal"
    rm = InitCubicGrid(nMol,ρ)
elseif occursin(lowercase(initialConfiguration),"cnf") #"cnf"  lowercase(initialConfiguration)
    rm, quat, box = ReadCNF("cnf_input.inp")
    nMol = length(rm)
    nAtoms = nMol * 3
    ρ = nMol / (box ^ 3)
    println(" Initial configuration is from file: cnf_input.inp")
    println("boxsize is: ", box)
    println("density is: ", ρ)
    sleep(3)

    """
    A&T use a box centered at [0,0,0] whereas we use a box from 0->box in all
    dimensions. This next part shift all coordinates by the magnitude of the
    smallest coordinate so all coordinates are not between 0->box
    """
    xl=yl=zl=0.0
    for (i,mol) in enumerate(rm)
        global xl,yl,zl
        if i == 1
            xl = mol[1]
            yl = mol[2]
            zl = mol[3]
        else
            if mol[1] < xl xl = mol[1] end
            if mol[2] < yl yl = mol[2] end
            if mol[3] < zl zl = mol[3] end
        end
    end
    xl = abs(xl)
    yl = abs(yl)
    zl = abs(zl)
    for (i,mol) in enumerate(rm)
        rm[i] = mol .+ [xl,yl,zl]
    end

else
    rm = [SVector{3,Float64}(rand(), rand(), rand()) .* box for i = 1:nMol]
end

# Generate atomic coordinates
ra = [SVector{3,Float64}(  0.0, 0.0, 0.0) for i = 1:nAtoms]
ra = []
finish = 0
thisMol_thisAtom = []
initQuaternions = []

for (i,com) in enumerate(rm)
    global finish
    start = finish + 1
    finish = start + 2
    if occursin(lowercase(initialConfiguration), "cnf")
        ei = quat[i]
    else
        ei = random_quaternion()
    end
    push!(initQuaternions,ei)
    ai = q_to_a( ei ) # Rotation matrix for i
    for a = 1:at_per_mol # Loop over all atoms
       #di(:,a) = MATMUL ( db(:,a), ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai
       push!(ra,com + SVector( MATMUL(ai , db[:,a]) ))
    end # End loop over all atoms
    push!(thisMol_thisAtom, SVector(start, finish))
end

PrintPDB(ra, box, 0, "pdbOutput")
ϵ = σ = ones(nAtoms)

total     = Properties(0.0, 0.0, 0.0, 0.0)
system    = Requirements(rm, ra, thisMol_thisAtom, ϵ, σ, box, r_cut)
totProps  = Properties2(temperature, ρ, Pressure(total, ρ, temperature, box^3),
                            dr_max, dϕ_max, 0.3, 0, 0, initQuaternions)

"""If not using a premade initial configuration, do Energy Minimization"""

if lowercase(initialConfiguration)  == "crystal"
    totProps.quat = @time EnergyMinimize(system,db, totProps.quat)
    initQuaternions = totProps.quat
end

for (i,com) in enumerate(system.rm)

    start = system.thisMol_theseAtoms[i][1]
    finish = system.thisMol_theseAtoms[i][2]
    ai = q_to_a( totProps.quat[i] ) # Rotation matrix for i
    j = 0
    for a = start:finish # Loop over all atoms
        j += 1
       system.ra[a] = com + SVector( MATMUL(ai , db[:,j]) )
    end # End loop over all atoms
end

total     = potential(system, Properties(0.0,0.0, 0.0, 0.0))
averages  = Properties(total.energy, total.virial,total.energy, total.virial) # initialize struct with averages
totProps  = Properties2(temperature, ρ, Pressure(total, ρ, temperature, box^3),
                            dr_max, dϕ_max, 0.3, 0, 0, initQuaternions)

println("TEST ", total.energy, " ", totProps.pressure)
#""" Main loop of simulation. Sweep over all atoms, then Steps, then Blocks."""

for blk = 1:nblock
    for step =1:nSteps
        for i = 1:length(system.rm)

            partial_old_e, partial_old_v = LJ_poly_ΔU(i, system)
            rm_old = deepcopy(system.rm[i])
            ra_old = deepcopy(system.ra[
                        system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]]
                                )
            rnew = random_translate_vector(totProps.dr_max, system.rm[i], box)
            ei = random_rotate_quaternion(totProps.dϕ_max, totProps.quat[i]) #quaternion()
            ai = q_to_a( ei ) # Rotation matrix for i
            ra_new = []

            for a = 1:at_per_mol # Loop over all atoms
               push!(ra_new,rnew + SVector( MATMUL(ai , db[:,a]) ))
            end # End loop over all atoms

            system.rm[i] = rnew

            system.ra[
                    system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]
                    ] = ra_new

            partial_new_e, partial_new_v = LJ_poly_ΔU(i, system)

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
                totProps.quat[i] = ei
            else
                system.rm[i] = rm_old
                system.ra[
                        system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]
                        ] = ra_old
                averages.energy += averages.old_e
                averages.virial += averages.old_v
            end

            # for troubleshooting checks that particles are in box
            minV, maxV = maxmin(system.rm)
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
    total2     = potential(system, Properties(0.0,0.0, 0.0, 0.0))
    if abs(total2.energy - total.energy) > 0.001
        print("SHITTTTTT, things aren't adding up")
    end
    PrintPDB(system.ra, box, blk, "pdbOutput")
    # Hella ugly output
    # TODO (BDK) modify to formatted output
    println(blk, " Energy: ", averages.energy / totProps.totalStepsTaken / nMol,
         " Ratio: ", totProps.numTranAccepted / totProps.totalStepsTaken,
         " Pressure: ", ρ*temperature + averages.virial / box^3 / totProps.totalStepsTaken, #  Pressure(total, ρ, temperature, box^3),
         " Pcut: ", pressure_lrc( ρ, system.r_cut ), " Ecorr: ",
         potential_lrc( ρ, r_cut ), "p delta: ", pressure_delta(ρ,system.r_cut ),
         " A & T p_c: " , ρ*temperature + averages.virial / box^3 /
         totProps.totalStepsTaken + pressure_delta(ρ,system.r_cut ),
         " instant energy: ", total.energy / nMol,
         " instant pressure: ", Pressure(total, ρ, temperature, box^3) +
         pressure_delta(ρ,system.r_cut ))
         println("box: ", box, "  density: ", ρ )
end # blk to nblock

finish = Dates.now()
difference =  finish - start


println("start: ", start)
println("finish: ", finish)
#println("difference: Seconds: ", difference)
println("Total runtime was: ", Dates.canonicalize(Dates.CompoundPeriod(difference) ) )
#Dates.format(DateTime("2017-10-01T01:02:03"), "H:M:S.s"
#Dates.canonicalize(Dates.CompoundPeriod(t2-t1))
"""
println("difference: Seconds: ", parse(Float64,split(difference)[1]) / 1000)
println("difference: Minutes: ", parse(Float64,split(difference)[1]) / 1000 / 60)
println("difference: Hours  : ", parse(Float64,split(difference)[1]) / 1000 / 3600)
"""
