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
# TODO (BDK) extend to mixtures
# TODO (BDK) add ewalds (routines written, need to incorporate)
# TODO (BDK) make tables for parameters
# TODO (BDK) need arrays for ϵ, σ
# TODO (BDK) NPT volume change (implement NPT in general)
# TODO (BDK) Implement REMC
# TODO (BDK) Make benchmarks - print out configurations and associated properties
# until then, manually enter at top
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# Things to organize for EWALD
# 1) prefactor: need k-vectors - this needs prep (code written already)
# 2) need array for charges, their coords, atomID, molID, k-fac max, alpha
# 3) need real, recipricol, self and tinfoil boundary routines (code already written)
# 4) need to read in NIST database format for testing
# 5) test vs. NIST SPC/E database:
#   https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10a-cutoff
# 6) DOUBLE CHECK ALL FREAKING UNITS!!!
################################################################################
temperature = 298.15 #0.6 #0.8772  # 1.2996
ρ = 0.00375000533
nMol = 100 * 3
nAtoms = nMol * 3
r_cut = 10 #2.5  # box / 2
nSteps = 20
nblock = 100
outputInterval = 100
initialConfiguration = "crystal"  # place atoms in a crystal structure
dϕ_max = 0.05
dr_max = 0.05

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
dr_max = box * 0.01
println("BoxSize: ", box)
#box / 30   # at 256 particles, ρ=0.75, T=1.0 this is 48% acceptance
global warnings = []

# Shift COM of body-fixed reference molecules to [0,0,0]
push!(warnings, "Overwriting triatomic glass to spce body-fixed. line 184 tests.jl")
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
            -0.8436
        else
            0.4218
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
ewald = EWALD(5.6 / box, 5, 27, 1, [SVector{3,Int32}(i, i, i) for i=1:3], [0.0,0.0],factor)  # kappa, nk, k_sq_max, NKVECS
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
            push!(ra, SVector( [com + SVector(MATMUL(ai, db[:, a]) )]... )
        end # End loop over all atoms
    end
    push!(thisMol_thisAtom, SVector(start, finish))
end
qq_r = ra  # this assumes charges are atom centered

# TODO add qq_q and qq_r to system::Requirements

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

println("TEst rcut, press_corr ", system.r_cut, "...........", factor)
println(ener_corr(system, 2, [system.nMols, system.nAtoms]))
println(press_corr(system, 2, [system.nMols, system.nAtoms]))
#sleep(20)
"""If not using a premade initial configuration, do Energy Minimization"""

if lowercase(initialConfiguration) == "crystal"
    totProps.quat = @time EnergyMinimize(system, db, totProps.quat,qq_q, qq_r, ewald)
    initQuaternions = totProps.quat
else
    totProps.quat = initQuaternions
end

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

total = potential(
    system,
    Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    ewald,
    qq_q,
    qq_r,
    #kfacs,
)
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

@code_warntype EwaldReal(
        qq_r,
        qq_q,
        ewald.kappa,
        box,
        system.thisMol_theseAtoms,
        2,
        system.r_cut,
    )

#total.recipOld = RecipLong(system, ewald, qq_r, qq_q) * ewald.factor
#total.recipOld = RecipLong(system, ewald, qq_r, qq_q, ewald.cfac) * ewald.factor
for blk = 1:nblock
    @time for step = 1:nSteps
        for i = 1:length(system.rm)
            partial_old_e, partial_old_v = LJ_poly_ΔU(i, system)

            # calculates all short range ewald energies
            partial_ewald_e, partial_ewald_v = EwaldShort(
                i,
                system,
                ewald,
                temperature,
                box,
                qq_r,
                qq_q,
                false,
            )
            partial_old_e += partial_ewald_e #+ total.recipOld
            partial_old_v += partial_ewald_v #+ total.recipOld / 3.0

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



            rnew = random_translate_vector(totProps.dr_max, system.rm[i], box)
            ei = random_rotate_quaternion(totProps.dϕ_max, totProps.quat[i]) #quaternion()
            ai = q_to_a(ei) # Rotation matrix for i
            ra_new = []

            for a = 1:at_per_mol # Loop over all atoms
                push!(ra_new, rnew + SVector(MATMUL(ai, db[:, a])))
            end # End loop over all atoms

            system.rm[i] = rnew

            # Update atom coords
            system.ra[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]] =
                ra_new
            # Update charge coords
            qq_r[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]] =
                ra_new

            # Calculate new LJ energy
            partial_new_e, partial_new_v = LJ_poly_ΔU(i, system)

            partial_ewald_e, partial_ewald_v = EwaldShort(
                i,
                system,
                ewald,
                temperature,
                box,
                qq_r,
                qq_q,
                false,
            )
            partial_new_e += partial_ewald_e
            partial_new_v += partial_ewald_v

            # Calculate new recipricol energy
            # note: this is a total energy, not the energy of i with the system
            ##recipEnergy =
                #RecipLong(system, ewald, qq_r, qq_q, kfacs) * ewald.factor
            ##    RecipLong(system, ewald, qq_r, qq_q) * ewald.factor
            #recipEnergy = RecipLong(system, ewald, ra_new, qq_q[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]], ewald.cfac) * ewald.factor
            ##partial_new_e += recipEnergy
            ##partial_new_v += recipEnergy / 3

             deltaRecip = RecipMove(
                system,
                ewald,
                ra_old,
                ra_new,
                qq_q[system.thisMol_theseAtoms[i][1]:system.thisMol_theseAtoms[i][2]]
                )
            #println("outer: ", deltaRecip)
            # Calculate difference in old and new system energy
            delta = (partial_new_e ) - (partial_old_e) + deltaRecip
                #(partial_new_e + recipEnergy) - (partial_old_e)
            #println(delta)
            # Test if we accept the move
            #println(delta, "    ", delta / temperature)
            if Metropolis(delta / temperature) # make sure units work
                total.energy += delta
                total.virial += (partial_new_v  - partial_old_v  ) + deltaRecip/3 # + recipEnergy / 3
                #total.recipOld = recipEnergy
                #total.recip = recipEnergy
                totProps.numTranAccepted += 1
                ne = averages.old_e + delta
                nv =
                    averages.old_v + partial_new_v  - partial_old_v + deltaRecip / 3#+ recipEnergy / 3
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

        end # i to nAtoms


    end # step to nSteps
    #total2 = potential(system, Properties(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    #if abs(total2.energy - total.energy) > 0.001
    #    println("SHITTTTTT, things aren't adding up")
    #end
    PrintPDB(system.ra, box, blk, "pdbOutput_atoms")
    PrintPDB(qq_r, box, blk, "pdbOutput_qq")
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
