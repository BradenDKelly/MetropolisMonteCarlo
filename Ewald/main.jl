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
using Profile


include("tests.jl")
include("auxillary.jl")
include("initialConfigurations.jl")
include("boundaries.jl")
include("energy.jl")
include("quaternions.jl")
include("volumeChange.jl")
include("ewalds.jl")
include("constants.jl")
include("banners.jl")
include("adjust.jl")


function PrintLine(s::String, n::Int64) # handy tool for outputting lines
    println(repeat(s, n)) # prints n copies of whatever s is.
end

Logo()
Random.seed!(11234)
start = Dates.now()
println(Dates.now())

################################################################################
# TODO (BDK) create JSON input file with starting parameters:
# TODO (BDK) add proper sampling
# TODO (BDK) Add some timing (~ 30 times faster than numpy)
# TODO (BDK) write more unit tests
# TODO (BDK) NPT volume change (implement NPT in general)
# TODO (BDK) Implement REMC
# TODO (BDK) Make benchmarks - print out configurations and associated properties
################################################################################
temperature = 298.15 #0.6 #0.8772  # 1.2996
ρ = 0.033101144   #0.015047707 #0.003633451 #0.00375000533 0.015047712

nMol = 300 #256 #100 * 10
nAtoms = nMol * 3
r_cut = 10.0 #2.5  # box / 2 Angstrom
nSteps = 10
nblock = 200
outputInterval = 100
initialConfiguration = "crystal"  # place atoms in a crystal structure
dϕ_max = 0.05
dr_max = 0.05
coulombStyle = "ewald"

# Set default values, check keys and typecheck values
defaults = Dict(
    "nblock" => 10,
    "nstep" => 1000,
    "temperature" => 1.0,
    "r_cut" => 2.5,
    "dr_max" => 0.15,
    "natoms" => 256,
    "initConfig" => "crystal",
    "rho" => 0.75,
    "ϵ" => 1.0,
    "σ" => 1.0,
)


################################################################################
#
#                        Start of Configuration
#
################################################################################
box = (nMol / ρ)^(1 / 3)
dr_max = 0.316555789 * 10.0 * 0.1
println("BoxSize: ", box)
#box / 30   # at 256 particles, ρ=0.75, T=1.0 this is 48% acceptance
warnings = []

# Shift COM of body-fixed reference molecules to [0,0,0]
push!(
    warnings,
    "Overwriting triatomic glass to spce body-fixed. line 184 tests.jl",
)
a = []
push!(a, SVector(db[:, 1]...))
push!(a, SVector(db[:, 2]...))
push!(a, SVector(db[:, 3]...))
com_db = COM(a, [15.9994, 1.008, 1.008]) # COM takes Vector{SVector}
db = db .- com_db

# Generate molecular COM coordinates
if lowercase(initialConfiguration) == "crystal"
    push!(
        warnings,
        "hardcoded ϵ and σ for crystal to spc/e . Line 98, main.jl.",
    )
    # make LJ table of values.
    σ_O = 0.316555789 * 10.0 # Å
    σ_H = 0.0 # nm
    ϵ_O = 78.1974311 # K   (ϵ/kᵦ)
    ϵ_H = 0.0 # K   (ϵ/kᵦ)
    rm = InitCubicGrid(nMol, ρ)
    molNames = ["Wat" for i = 1:nMol]
    molType = [1 for i = 1:nMol]
    atomName = [
        if (i - 1) % 3 == 0
            "O"
        else
            "H"
        end for i = 1:nAtoms
    ]
    atomType = [
        if (i - 1) % 3 == 0
            1
        else
            2
        end for i = 1:nAtoms
    ]
    qq_q = [
        if (i - 1) % 3 == 0
            -0.42380 * 2
        else
            0.42380
        end for i = 1:nAtoms
    ]

elseif occursin(lowercase(initialConfiguration), "cnf") #"cnf"  lowercase(initialConfiguration)
    rm, quat, box = ReadCNF("cnf_input.inp")
    nMol = length(rm)
    nAtoms = nMol * 3
    ρ = nMol / (box^3)
    println(" Initial configuration is from file: cnf_input.inp")
    println("boxsize is: ", box)
    println("density is: ", ρ)

    """
    A&T use a box centered at [0,0,0] whereas we use a box from 0->box in all
    dimensions. This next part shift all coordinates by the magnitude of the
    smallest coordinate so all coordinates are not between 0->box
    """
    xl = yl = zl = 0.0
    for (i, mol) in enumerate(rm)
        global xl, yl, zl
        if i == 1
            xl = mol[1]
            yl = mol[2]
            zl = mol[3]
        else
            if mol[1] < xl
                xl = mol[1]
            end
            if mol[2] < yl
                yl = mol[2]
            end
            if mol[3] < zl
                zl = mol[3]
            end
        end
    end
    xl = abs(xl)
    yl = abs(yl)
    zl = abs(zl)
    for (i, mol) in enumerate(rm)
        rm[i] = mol .+ [xl, yl, zl]
    end

elseif occursin(lowercase(initialConfiguration), "nist")
    # THis file stores a single configuration of 100-750 SPC/E water molecules.
    # This is for testing only. NIST has results for energy contributions.
    # load the configuration, get the atom types, charges, coordinates.
    # Return to here and test.
    filename = "spce_sample_config_periodic1.txt"
    filename = "coord750.txt"
    qq_r, qq_q, rm, ra, atomTracker, box, atomName, atomType =
        ReadNIST(filename)

    # make LJ table of values.
    σ_O = 0.316555789 * 10.0 # convert nm to Å
    σ_H = 0.0 # nm
    ϵ_O = 78.1974311 # K   (ϵ/kᵦ)
    ϵ_H = 0.0 # K   (ϵ/kᵦ)

    xl = yl = zl = 0.0
    for (i, mol) in enumerate(rm)
        global xl, yl, zl
        if i == 1
            xl = mol[1]
            yl = mol[2]
            zl = mol[3]
        else
            if mol[1] < xl
                xl = mol[1]
            end
            if mol[2] < yl
                yl = mol[2]
            end
            if mol[3] < zl
                zl = mol[3]
            end
        end
    end
    xl = abs(xl)
    yl = abs(yl)
    zl = abs(zl)
    for (i, mol) in enumerate(rm)
        rm[i] = mol .+ [xl, yl, zl]
    end
    for (i, part) in enumerate(ra)
        ra[i] = part .+ [xl, yl, zl]
        qq_r[i] = part .+ [xl, yl, zl]
    end
else
    rm = [SVector{3,Float64}(rand(), rand(), rand()) .* box for i = 1:nMol]
end

###########################
#
#                       Test electrostatics
#
########################################

ewald = EWALD(
    5.6 / box,
    5,
    27,
    1,
    [SVector{3,Int32}(i, i, i) for i = 1:3],
    [0.0, 0.0],     # dummy values
    zeros(ComplexF64, 2),     # dummy values
    zeros(ComplexF64, 2),
    factor,
)  # kappa, nk, k_sq_max, NKVECS
#ewald = EWALD(5.6 / box, 5, 27, [0.0,0.0],factor)
println("Check initialization of EWALD on line 162. 5.6 / box?")
ewald = PrepareEwaldVariables(ewald, box) # better one # cfac, kxyz,
#kfacs, ewald = SetupKVecs(ewald, box)
println("Set up initial Ewald k-vectors")

if occursin(lowercase(initialConfiguration), "cnf")
    ϵ = σ = ones(1, 1)
else
    ϵ = [ϵ_O, ϵ_H]
    σ = [σ_O, σ_H]
    molNames = ["Wat" for i = 1:length(rm)]
    molTypes = [1 for i = 1:length(rm)]

end

if occursin(lowercase(initialConfiguration), "cnf") ||
   occursin(lowercase(initialConfiguration), "crystal")
    ra = []
end
thisMol_thisAtom = []
initQuaternions = []
finish = 0

for (i, com) in enumerate(rm)
    global finish
    start = finish + 1
    finish = start + 2

    if occursin(lowercase(initialConfiguration), "cnf")
        ei = quat[i]
    else
        ei = random_quaternion()
    end
    push!(initQuaternions, ei)
    if occursin(lowercase(initialConfiguration), "cnf") ||
       occursin(lowercase(initialConfiguration), "crystal")
        ai = q_to_a(ei) # Rotation matrix for i
        for a = 1:at_per_mol # Loop over all atoms
            # di(:,a) = MATMUL ( db(:,a), ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai
            push!(ra, com + SVector(MATMUL(ai, db[:, a])))
        end # End loop over all atoms
    end
    push!(thisMol_thisAtom, SVector(start, finish))
end
ra = [SVector(ra[i]...) for i = 1:length(ra)]
#qq_r::Vector{SVector{3,Float64}}
qq_r = ra  # this assumes charges are atom centered

# TODO add qq_q and qq_r to system::Requirements

# check that simulation box is charge neutral
@assert sum(qq_q) == 0.0

PrintPDB(ra, box, 0, "pdbOutput")
total = Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
system = Requirements(
    rm,
    ra,
    length(rm), #nMols,
    length(ra), #nAtoms,
    length(qq_r), #nCharges,
    thisMol_thisAtom,
    molNames,
    molTypes,
    atomName,
    atomType,
    Tables(ϵ, σ),
    box,
    r_cut,
)

# set up struct with general system properties like T and move acceptance
totProps = Properties2(
    temperature,
    ρ,
    Pressure(total, ρ, temperature, box^3),
    dr_max,
    dϕ_max,
    0.3,
    0,
    0,
    initQuaternions,
)

trans_moves = Moves(0,0,0,0,0.5, dr_max)

println("TEst rcut, press_corr ", system.r_cut, "...........", factor)
println(ener_corr(system, 2, [system.nMols, system.nAtoms]))
println(press_corr(system, 2, [system.nMols, system.nAtoms]))
#sleep(20)
"""If not using a premade initial configuration, do Energy Minimization"""
#=
if lowercase(initialConfiguration) == "crystal"
    totProps.quat = @time EnergyMinimize(system, db, totProps.quat,qq_q, qq_r, ewald)
    initQuaternions = totProps.quat
else
    totProps.quat = initQuaternions
end
=#
#=
for (i, com) in enumerate(system.rm)
    start = system.thisMol_theseAtoms[i][1]
    finish = system.thisMol_theseAtoms[i][2]
    ai = q_to_a(totProps.quat[i]) # Rotation matrix for i
    j = 0
    for a = start:finish # Loop over all atoms
        j += 1
        system.ra[a] = com + SVector(MATMUL(ai, db[:, j]))
    end # End loop over all atoms
end
=#
LJ, reall, recipEnergy = 0.0, 0.0, 0.0
if coulombStyle == "bare"
    total, LJ, reall = potential(
        system,
        Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        ewald,
        qq_q,
        qq_r,
        "bare"
    )
else
    total, LJ, reall, recipEnergy = potential(
        system,
        Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        ewald,
        qq_q,
        qq_r,
        #kfacs,
    )
end
initial = total.energy
averages = Properties(
    total.energy,
    total.virial,
    0.0,
    0.0,
    0.0,
    total.energy,
    total.virial,
) # initialize struct with averages
totProps = Properties2(
    temperature,
    ρ,
    Pressure(total, ρ, temperature, box^3),
    dr_max,
    dϕ_max,
    0.3,
    0,
    0,
    initQuaternions,
)

println("TEST ", total.energy, " ", totProps.pressure)
#""" Main loop of simulation. Sweep over all atoms, then Steps, then Blocks."""

if occursin(lowercase(initialConfiguration), "nist")
    error("NIST can only do starting configuration, Stopping now.")
end

PrintOutput(
    system,
    totProps,
    atomType,
    atomName,
    qq_r,
    qq_q,
    box,
    1,
    "xyz_quat",
)
#=
@code_warntype RecipMove(
   system,
   ewald,
   system.ra[system.thisMol_theseAtoms[1][1]:system.thisMol_theseAtoms[1][2]],
   system.ra[system.thisMol_theseAtoms[1][1]:system.thisMol_theseAtoms[1][2]].* 1.01,
   qq_q[system.thisMol_theseAtoms[1][1]:system.thisMol_theseAtoms[1][2]]
   )
=#
#=
@code_warntype EwaldReal(
        qq_r,
        qq_q,
        ewald.kappa,
        box,
        system.thisMol_theseAtoms,
        2,
        system.r_cut,
    )
=#
#@code_warntype EwaldReal(qq_r, qq_q, ewald.kappa, box, thisMol_thisAtom, 2, system)
#total.recipOld = RecipLong(system, ewald, qq_r, qq_q) * ewald.factor
#total.recipOld = RecipLong(system, ewald, qq_r, qq_q, ewald.cfac) * ewald.factor
ovr_count = 0

testing = false

#@code_warntype CoulombReal(qq_r, qq_q, box, 3, system)
function Loop(system, totProps, ovr_count, box, temperature, total, trans_moves,
                qq_q, qq_r, ewald, averages, ρ, atomName, atomType, coulombStyle)
    for blk = 1:nblock
        #global ovr_count, trans_moves
        @time for step = 1:nSteps
            for i = 1:length(system.rm)
                partial_old_e, partial_old_v = LJ_poly_ΔU(i, system)
                LLJ1 = partial_old_e
                # calculates all short range ewald energies
                if coulombStyle == "bare"
                    partial_ewald_e,  overlap1 =
                        CoulombReal(qq_r, qq_q, box, i, system)
                else

                    partial_ewald_e, partial_ewald_v, overlap1 =
                        EwaldShort(i, system, ewald, box, qq_r, qq_q, false)
                        partial_old_v += partial_ewald_v #+ total.recipOld / 3.0
                end
                reall1 = partial_ewald_e * ewald.factor
                partial_old_e += partial_ewald_e * ewald.factor #+ total.recipOld


                #=
                first, LJ1, real1, recipEnergy1 = potential(
                    system,
                    Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    ewald,
                    qq_q,
                    qq_r,
                    #kfacs,
                )
                =#
                #################################################
                #
                #      Move and Rotate a particle
                #
                #################################################

                rm_old = deepcopy(system.rm[i])
                ra_old =
                    deepcopy(system.ra[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]])

                ##recipEnergy = RecipLong(system, ewald, ra_old, qq_q[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]], ewald.cfac) * ewald.factor
                ##partial_old_e += total.recipOld #  recipEnergy
                ##partial_old_v += total.recipOld / 3.0 #recipEnergy / 3


                # move particle
                rnew = random_translate_vector(totProps.dr_max, system.rm[i], box)
                system.rm[i] = rnew
                # rotate molecule and update atom positions
                ei = random_rotate_quaternion(totProps.dϕ_max, totProps.quat[i]) #quaternion()
                ai = q_to_a(ei) # Rotation matrix for i
                ra_new = []
                for a = 1:at_per_mol # Loop over all atoms
                    push!(ra_new, SVector(rnew + SVector(MATMUL(ai, db[:, a]))))
                end # End loop over all atoms
                ra_new = [SVector(item...) for item in ra_new]

                # Update atom coords
                system.ra[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]] =
                    ra_new
                # Update charge coords
                qq_r[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]] =
                    ra_new

                # Calculate new LJ energy
                partial_new_e, partial_new_v = LJ_poly_ΔU(i, system)
                LLJ2 = partial_new_e
                # calculate new real contribution to ewalds
                if coulombStyle == "bare"
                    partial_ewald_e,  overlap2 =
                    CoulombReal(qq_r, qq_q, box, i, system)
                else
                    partial_ewald_e, partial_ewald_v, overlap2 =
                        EwaldShort(i, system, ewald, box, qq_r, qq_q, false)
                        partial_new_v += partial_ewald_v
                end
                reall2 = partial_ewald_e * ewald.factor
                partial_new_e += partial_ewald_e * ewald.factor


                # Calculate new recipricol energy
                # note: this is a total energy, not the energy of i with the system
                ##recipEnergy =
                #RecipLong(system, ewald, qq_r, qq_q, kfacs) * ewald.factor
                ##    RecipLong(system, ewald, qq_r, qq_q) * ewald.factor
                #recipEnergy = RecipLong(system, ewald, ra_new, qq_q[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]], ewald.cfac) * ewald.factor
                ##partial_new_e += recipEnergy
                ##partial_new_v += recipEnergy / 3
                if overlap1 || overlap2
                    overlap = true
                else
                    overlap = false
                end

                if overlap == false && coulombStyle != "bare"
                    deltaRecip =
                        RecipMove(
                            system,
                            ewald,
                            ra_old,
                            ra_new,
                            qq_q[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]],
                        ) * ewald.factor
                else
                    deltaRecip = 0.0
                end
                #println("outer: ", deltaRecip)
                # Calculate difference in old and new system energy
                delta = (partial_new_e) - (partial_old_e) + deltaRecip
                #=
                second, LJ2, real2, recipEnergy2 = potential(
                    system,
                    Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    ewald,
                    qq_q,
                    qq_r,
                    #kfacs,
                )
                # #=
                println("total diff: ", second.energy - first.energy)
                println("single diff: ", delta)
                println("diff LJ full: ", LJ2 - LJ1)
                println("diff real full: ", real2 - real1)
                println("diff recip full: ", recipEnergy2 - recipEnergy1)
                println("diff recip single: ", deltaRecip)
                println("diff real single: ", reall2 - reall1)
                println("diff LJ single: ", LLJ2 - LLJ1)
                sleep(1)
                # =#
                delta = second.energy - first.energy
                =#
                #sleep(1)
                #(partial_new_e + recipEnergy) - (partial_old_e)
                #println(delta , "      ", deltaRecip)
                # Test if we accept the move

                if overlap
                    ovr_count += 1
                end
                if Metropolis(delta / temperature) && overlap == false# make sure units work
                    total.energy += delta
                    total.virial += (partial_new_v - partial_old_v) + deltaRecip / 3 # + recipEnergy / 3
                    #total.recipOld = recipEnergy
                    #total.recip = recipEnergy
                    totProps.numTranAccepted += 1
                    trans_moves.naccept += 1
                    ne = averages.old_e + delta
                    nv =
                        averages.old_v + partial_new_v - partial_old_v +
                        deltaRecip / 3#+ recipEnergy / 3
                    averages.energy += ne
                    averages.virial += nv
                    averages.old_e = ne
                    averages.old_v = nv
                    #averages.recipOld = recipEnergy
                    #averages.recip = recipEnergy
                    totProps.quat[i] = ei
                else
                    system.rm[i] = rm_old
                    system.ra[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]] =
                        ra_old
                    qq_r[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]] =
                        ra_old
                    averages.energy += averages.old_e
                    averages.virial += averages.old_v
                    averages.recip = averages.recipOld
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
                trans_moves.attempt += 1

            end # i to nAtoms
            trans_moves.dr_max = totProps.dr_max
            trans_moves = Adjust!(trans_moves,box)
            totProps.dr_max = trans_moves.dr_max


        end # step to nSteps
        println(trans_moves.naccept / trans_moves.attempt, "   ", trans_moves.dr_max)
        #total2 = potential(system, Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        #if abs(total2.energy - total.energy) > 0.001
        #    println("SHITTTTTT, things aren't adding up")
        #end
        PrintPDB(system.ra, box, blk, "pdbOutput_atoms")
        #PrintPDB(qq_r, box, blk, "pdbOutput_qq")
        # Hella ugly output
        # TODO (BDK) modify to formatted output
        println(
            blk,
            " Energy: ",
            averages.energy / totProps.totalStepsTaken / nMol,
            " Ratio: ",
            totProps.numTranAccepted / totProps.totalStepsTaken,
            " Pressure: ",
            ρ * temperature + averages.virial / box^3 / totProps.totalStepsTaken, #  Pressure(total, ρ, temperature, box^3),
            " Pcut: ",
            pressure_lrc(ρ, system.r_cut),
            " Ecorr: ",
            potential_lrc(ρ, r_cut),
            "p delta: ",
            pressure_delta(ρ, system.r_cut),
            " A & T p_c: ",
            ρ * temperature +
            averages.virial / box^3 / totProps.totalStepsTaken +
            pressure_delta(ρ, system.r_cut),
            " instant energy: ",
            total.energy / nMol,
            " instant pressure: ",
            Pressure(total, ρ, temperature, box^3) +
            pressure_delta(ρ, system.r_cut),
            " overlap count: ", ovr_count
        )
        println("box: ", box, "  density: ", ρ)

        PrintOutput(
            system,
            totProps,
            atomType,
            atomName,
            qq_r,
            qq_q,
            box,
            blk,
            "xyz_quat",
        )
    end # blk to nblock
end

Loop(system, totProps, ovr_count, box, temperature, total, trans_moves,
                qq_q, qq_r, ewald, averages, ρ, atomName, atomType, coulombStyle)

PrintOutput(
    system,
    totProps,
    atomType,
    atomName,
    qq_r,
    qq_q,
    box,
    1,
    "xyz_quat_final",
)

finish = Dates.now()
difference = finish - start


println("start: ", start)
println("finish: ", finish)
#println("difference: Seconds: ", difference)
println(
    "Total runtime was: ",
    Dates.canonicalize(Dates.CompoundPeriod(difference)),
)

Completion()
