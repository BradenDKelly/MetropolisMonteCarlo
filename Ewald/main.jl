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
include("setup.jl")
include("structs.jl")


function PrintLine(s::String, n::Int64) # handy tool for outputting lines
    println(repeat(s, n)) # prints n copies of whatever s is.
end

Logo()
#Random.seed!(11234)
start = Dates.now()
println(Dates.now())

################################################################################
# TODO (BDK) Add in setup.jl and structs.jl
# TODO (BDK) Make initial configuration generation w/ Topology struct
# TODO (BDK) Implement pressure /w ewalds and wolf summation
# TODO (BDK) Make functions in place so no temp arrays are generated
# TODO (BDK) make input density in kg/m3 and convert
# TODO (BDK) read in checkpoint file/ make restart file
# TODO (BDK) Separate translation and rotation moves (done ✔)
# TODO (BDK) Tidy up code
# TODO (BDK) create JSON input file with starting parameters:
# TODO (BDK) make fully generalized setup
# TODO (BDK) add proper sampling
# TODO (BDK) Add some timing (~ 30 times faster than numpy)
# TODO (BDK) write more unit tests
# TODO (BDK) NPT volume change (implement NPT in general)
# TODO (BDK) Implement REMC
# TODO (BDK) Make benchmarks - print out configurations and associated properties
# TODO (BDK) Add 'smart mc' where certain particles are moved more
# TODO (BDK) Add configurational bias
# TODO (BDK) Add alchemical free energy perturbation
# TODO (BDK) Make simulation structs for object-oriented scripting
################################################################################
temperature = 298.15 #0.6 #0.8772  # 1.2996
ρ = 0.033101144   #0.015047707 #0.003633451 #0.00375000533 0.015047712

nMol = 1000 #256 #100 * 10
nAtoms = nMol * 3
r_cut = 10.0 #2.5  # box / 2 Angstrom
nSteps = 100
nblock = 100
outputInterval = 100
initialConfiguration = "crystal"  # place atoms in a crystal structure
dr_max = 0.15
dϕ_max = 0.05
coulombStyle = "ewald"
Wolf = false

if Wolf
    println("This simulation uses Wolf Summation")
else
    println("This simulation uses Ewald Summation")
end
# Set default values, check keys and typeche ck values
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
    "mass" => [15.998, 1.008]
)

probability_of_move = Dict("translation" => 0.5, "rotation" => 0.5)
sum_value = 0.0
for (key, value) in probability_of_move
    global sum_value
    sum_value += value
    probability_of_move[key] = sum_value
end
# ensure probabilities sum to 1.0
for (key, value) in probability_of_move
    probability_of_move[key] /= sum_value
end
println(probability_of_move)



################################################################################
#
#                        Start of Configuration
#
################################################################################
box = (nMol / ρ)^(1 / 3)
dr_max = 0.316555789
println("BoxSize: ", box)
#box / 30   # at 256 particles, ρ=0.75, T=1.0 this is 48% acceptance
warnings = []

# Shift COM of body-fixed reference molecules to [0,0,0]
push!(
    warnings,
    "Overwriting triatomic glass to spce body-fixed. line 120 tests.jl",
)
#=
a = []
push!(a, SVector(db[:, 1]...))
push!(a, SVector(db[:, 2]...))
push!(a, SVector(db[:, 3]...))
com_db = COM(a, [15.9994, 1.008, 1.008]) # COM takes Vector{SVector}
db = db .- com_db
=#
# Generate molecular COM coordinates
if lowercase(initialConfiguration) == "crystal"
    push!(
        warnings,
        "hardcoded ϵ and σ for crystal to spc/e . Line 98, main.jl.",
    )
    #=
    what the deal is:
    ----------------------
        Generate center-of-mass coords as per a crystal
        Assign random quaternion to each, and add in atoms

    Outcome:
    ---------------------
    Have a starting structure and all parameters initialized
    =#

    moleculeList = []
    bodyFixed = []

    top_file = "water.top"
    specieList = ["tip3p.pdb"]  # "mea.pdb",
    @time systemTop = ReadTopFile(top_file) # returns struct FFParameters # located in Setup.jl

    for i in eachindex(specieList)
        push!( moleculeList,ReadPDB(specieList[i])  )  # get molecule info from pdb. Returns struct Topology # located in Setup.jl
        temp = BodyFixed(moleculeList[i], systemTop)
        tempCOM = Center_of_Mass(temp.r, temp.mass)
        Shift_COM_to_Zero!(temp.r, tempCOM )

        temp = BodyFixed(temp.r,temp.mass,temp.atype)
        push!( bodyFixed,temp) #BodyFixed(moleculeList[i], systemTop) )  # get atom coords for quaternions # located in quaternions.jl
    end
    db = [coords.r for coords in bodyFixed]

topology, initQuaternions = Initialize(systemTop, moleculeList, bodyFixed,box)

#println(topology)


###########################################
#      Main Data Structures !!
##########################################
soa, moa = MakeAtomArrays(systemTop,topology,initQuaternions, "kmc")  # located in setup.jl
#intraFF, vdwTable, qqTable, nonbonded_matrix, scaled_pairs = MakeTables(systemTop,atomsPDB) # located in Setup.jl
@time PrintPDB(soa, moa, topology.box, 333, "kmc_output")

intraFF, vdwTable, nonbonded_matrix, scaled_pairs = MakeTables(systemTop,topology) #  qqTable, located in Setup.jl
num_atom_types = length(systemTop.atomTypes)
vdwTable.ϵᵢⱼ /= R # convert to K from kJ/mol
vdwTable.σᵢⱼ *= 10.0 # convert to Å from nm

# this stucture holds information on the number of atoms, molecules, atom types, molecule types, charges
numbers = Numbers(length(soa), length(moa), num_atom_types, length(systemTop.molParams), count(!iszero, soa.charge) )

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
if Wolf
    alpha = 5.6
else
    alpha = 5.6
end
ewald = EWALD(
    alpha / box,
    #alpha,
    5,
    27,
    1,
    [SVector{3,Int32}(i, i, i) for i = 1:3],
    [0.0, 0.0],     # dummy values
    zeros(ComplexF64, 2),     # dummy values
    zeros(ComplexF64, 2),
    factor,
)  # kappa, nk, k_sq_max, NKVECS

ewald = PrepareEwaldVariables(ewald, box) # better one # cfac, kxyz,
#kfacs, ewald = SetupKVecs(ewald, box)
println("Set up initial Ewald k-vectors")

if occursin(lowercase(initialConfiguration), "cnf")
    ϵ = σ = ones(1, 1)
else
    #=
    ϵ = [ϵ_O, ϵ_H]
    σ = [σ_O, σ_H]
    molNames = ["Wat" for i = 1:length(rm)]
    molTypes = [1 for i = 1:length(rm)]
    =#
end

if occursin(lowercase(initialConfiguration), "cnf") ||
   occursin(lowercase(initialConfiguration), "crystal")
    ra = []
end

try
    thisMol_thisAtom = []
    initQuaternions = []
    finish = 0
    for (i, com) in enumerate(rm)
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
catch e
    println("Using soa and moa rather than cnf or nist")
end


# TODO add qq_q and qq_r to system::Requirements

# check that simulation box is charge neutral
@assert isapprox(sum(soa.charge[:]),0.0,atol=0.00001)

#PrintPDB(ra, box, 0, "pdbOutput_molecular")
total = Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

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
    r_cut,
    r_cut,
    box
)

trans_moves = Moves(0, 0, 0, 0, 0.5, dr_max)
rot_moves = Moves(0, 0, 0, 0, 0.5, dϕ_max)

println("TEst rcut, press_corr ", totProps.LJ_rcut, "...........", factor)
# TODO new press_cor and ener_cor
#println(ener_corr(system, 2, [system.nMols, system.nAtoms]))
#println(press_corr(system, 2, [system.nMols, system.nAtoms]))

LJ, reall, recipEnergy = 0.0, 0.0, 0.0
if coulombStyle == "bare"
    total, LJ, reall = potential(
        system,
        Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        ewald,
        qq_q,
        qq_r,
        "bare",
    )
elseif Wolf
    total = potential(
                    moa,
                    soa,
                    Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),#triggers wolf summations using double strings
                    ewald,
                    vdwTable,
                    totProps,
                    )

else
    total  = potential(moa,
                        soa,
                        Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                        ewald,
                        vdwTable,
                        totProps,
                        "ewald"
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
    r_cut,
    r_cut,
    box
)

println("TEST ", total.energy, " ", totProps.pressure)
#""" Main loop of simulation. Sweep over all atoms, then Steps, then Blocks."""

if occursin(lowercase(initialConfiguration), "nist")
    error("NIST can only do starting configuration, Stopping now.")
end

ovr_count = 0

testing = false

#@code_warntype CoulombReal(qq_r, qq_q, box, 3, system)
#@code_warntype LJ_poly_ΔU(2, moa, soa, vdwTable, r_cut, box)
#@code_warntype EwaldReal(2, moa, soa, ewald,totProps.qq_rcut, box)


#################################################
"""soa and moa loop for MMC"""
function Loop(
    soa,
    moa,
    systemTop,
    vdwTable,
    numbers,
    totProps,
    ovr_count,
    box,
    temperature,
    total,
    trans_moves,
    rot_moves,
    #qq_q,
    #qq_r,
    ewald,
    averages,
    ρ,
    #atomName,
    #atomType,
    coulombStyle,
    db
)
    @assert totProps.LJ_rcut < box / 2
    @assert ρ > 0.0
    @assert box > 0.0

    for blk = 1:nblock
        #global ovr_count, trans_moves
        for step = 1:nSteps
            for i = 1:numbers.molecules
                partial_old_e, partial_old_v = LJ_poly_ΔU(i, moa, soa, vdwTable, r_cut, box)
                LLJ1 = partial_old_e
                # calculates all short range ewald energies
                if coulombStyle == "bare"
                    partial_ewald_e, overlap1 =
                        CoulombReal(qq_r, qq_q, box, i, system)
                    reall1 = partial_ewald_e * ewald.factor
                    partial_old_e += partial_ewald_e * ewald.factor #+ total.recipOld
                else

                    partial_ewald_e, partial_ewald_v, overlap1 =
                    EwaldShort( i,moa, soa, totProps, ewald, box)      # qq_factor already included
                    partial_old_v += partial_ewald_v #+ total.recipOld / 3.0
                    reall1 = partial_ewald_e
                    partial_old_e += partial_ewald_e
                end

                #################################################
                #
                #      Move or Rotate a particle
                #
                #################################################

                rm_old = deepcopy(moa.COM[i])
                ra_old = deepcopy(soa.coords[moa[i].firstAtom:moa[i].lastAtom])
                chose_move = rand()
                at_per_mol = length(systemTop.molParams[moa[i].molType].atoms)
                #println(" $i  $at_per_mol  $chose_move")
                if chose_move < probability_of_move["translation"]
                    # move particle
                    trans_moves.attempt += 1
                    rnew = random_translate_vector(
                        totProps.dr_max,
                        moa.COM[i],
                        box,
                    )
                    moa.COM[i] = rnew
                    ei = totProps.quat[i]
                    ai = q_to_a(ei)
                elseif chose_move <= probability_of_move["rotation"]
                    rot_moves.attempt += 1
                    rnew = moa.COM[i]
                    # rotate molecule and update atom positions
                    ei = random_rotate_quaternion(
                        totProps.dϕ_max,
                        totProps.quat[i],
                    ) #quaternion()
                    ai = q_to_a(ei) # Rotation matrix for i
                else
                    println("No move selected, exiting main.jl ~ line 490")
                    exit()
                end
                ra_new = []
                #if i == 1 @assert(at_per_mol == 11) end
                for a = 1:at_per_mol # Loop over all atoms
                    # di(:,a) = MATMUL ( db(:,a), ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai
                    push!(ra_new, moa[i].COM + SVector(MATMUL(ai, db[soa[i].molType][a])))
                end # End loop over all atoms
                ra_new = [SVector(item...) for item in ra_new]

                # Update atom coords
                soa.coords[moa[i].firstAtom:moa[i].lastAtom] = ra_new
                # Update charge coords
                #qq_r[moa[i].firstAtom:moa[i].lastAtom] =  ra_new

                # Calculate new LJ energy
                partial_new_e, partial_new_v = LJ_poly_ΔU(i, moa, soa, vdwTable, r_cut, box)
                LLJ2 = partial_new_e
                # calculate new real contribution to ewalds
                if coulombStyle == "bare"
                    partial_ewald_e, overlap2 =
                        CoulombReal(qq_r, qq_q, box, i, system)
                    reall2 = partial_ewald_e * ewald.factor
                    partial_new_e += partial_ewald_e * ewald.factor
                else
                    partial_ewald_e, partial_ewald_v, overlap2 =
                        EwaldShort( i,moa, soa, totProps, ewald, box)
                    partial_new_v += partial_ewald_v
                    reall2 = partial_ewald_e
                    partial_new_e += partial_ewald_e
                end


                if overlap1 || overlap2
                    overlap = true
                else
                    overlap = false
                end

                if overlap == false && coulombStyle != "bare" && Wolf != true
                    deltaRecip, ewald = RecipMove(
                        box,
                        ewald,
                        ra_old,
                        ra_new,
                        soa.charge[moa.firstAtom[i]:moa.lastAtom[i]]
                    )
                else
                    deltaRecip = 0.0
                end

                # Calculate difference in old and new system energy
                delta = (partial_new_e) - (partial_old_e) + deltaRecip

                if overlap
                    ovr_count += 1
                end
                if Metropolis(delta / temperature) && overlap == false# make sure units work
                    total.energy += delta
                    total.virial +=
                        (partial_new_v - partial_old_v) + deltaRecip / 3 # + recipEnergy / 3
                    #total.recipOld = recipEnergy
                    #total.recip = recipEnergy
                    totProps.numTranAccepted += 1
                    if chose_move < probability_of_move["translation"]
                        trans_moves.naccept += 1
                    elseif chose_move <= probability_of_move["rotation"]
                        rot_moves.naccept += 1
                    end
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
                    ewald.sumQExpOld = [item for item in ewald.sumQExpNew]
                else
                    moa.COM[i] = rm_old
                    soa.coords[moa[i].firstAtom:moa[i].lastAtom] =  ra_old
                    averages.energy += averages.old_e
                    averages.virial += averages.old_v
                    averages.recip = averages.recipOld
                    ewald.sumQExpNew = [item for item in ewald.sumQExpOld]
                end

                # for troubleshooting checks that particles are in box
                minV, maxV = maxmin(moa.COM)
                #println(minV, maxV)
                if minV < 0.0
                    println("Shit, particle is less than 0")
                end
                if maxV > box
                    println("Shit, particle is outside box")
                end

                totProps.totalStepsTaken += 1


            end # i to nAtoms
            trans_moves.d_max = totProps.dr_max
            trans_moves = Adjust!(trans_moves, box)
            totProps.dr_max = trans_moves.d_max

            rot_moves.d_max = totProps.dϕ_max
            rot_moves = Adjust_rot!(rot_moves, box)
            totProps.dϕ_max = rot_moves.d_max

            #@assert qq_r == system.ra


        end # step to nSteps
        #println(trans_moves.naccept / trans_moves.attempt, "   ", trans_moves.d_max)
        #total2 = potential(system, Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        #if abs(total2.energy - total.energy) > 0.001
        #    println("SHITTTTTT, things aren't adding up")
        #end
        #PrintPDB(system.ra, box, blk, "pdbOutput_molecular_rcut")
        #PrintPDB(qq_r, box, blk, "pdbOutput_qq")
        # Hella ugly output
        # TODO (BDK) modify to formatted output
        PrintPDB(soa,moa, topology.box, blk, "final")
        line = @sprintf(
            "Block: %4d, Energy: %8.2f, Ratio trans: %4.2f, dr_max: %4.2f, Ratio rot: %4.2f, dϕ_max: %4.2f, instant energy: %8.2f, overlap count: %4d, pressure: %8.2f",
            blk,
            averages.energy / totProps.totalStepsTaken / numbers.molecules,
            trans_moves.naccept / trans_moves.attempt,
            trans_moves.d_max,
            rot_moves.naccept / rot_moves.attempt,
            rot_moves.d_max,
            total.energy / numbers.molecules,
            ovr_count,
            4.60453 + total.virial / box / box / box
        )
        println(line)

        #println("box: ", box, "  density: ", ρ)
        #=
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
        =#
    end # blk to nblock
end
#=
Loop(
    system,
    totProps,
    ovr_count,
    box,
    temperature,
    total,
    trans_moves,
    rot_moves,
    qq_q,
    qq_r,
    ewald,
    averages,
    ρ,
    atomName,
    atomType,
    coulombStyle,
)
=#
Loop(
    soa,
    moa,
    systemTop,
    vdwTable,
    numbers,
    totProps,
    ovr_count,
    box,
    temperature,
    total,
    trans_moves,
    rot_moves,
    #qq_q,
    #qq_r,
    ewald,
    averages,
    ρ,
    #atomName,
    #atomType,
    coulombStyle,
    db
)
PrintPDB(soa,moa, topology.box, 100, "final")
#=
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
=#
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

#=
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
=#
