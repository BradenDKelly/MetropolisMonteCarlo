##########
#
#                         Energy Routines
#
###############
using Setfield
using Distributions
using BenchmarkTools

include("auxillary.jl")
include("ewalds.jl")
""" Calculates LJ potential between particle 'i' and the other N-1 particles.
 Cut & Shifted Potential"""
function LJ_poly_ΔU_shifted(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # Uses σ = ϵ = 1 units
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    # Cutoff distance and force-shift parameters (all private) chosen as per the reference:
    # S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)
    r_cut = 2.612 # in sigma=1 units, where r_cut = 1.2616 nm, sigma = 0.483 nm
    sr_cut = 1.0 / r_cut
    sr_cut6 = sr_cut^6
    sr_cut12 = sr_cut6^2
    lambda1 = 4.0 * (7.0 * sr_cut6 - 13.0 * sr_cut12)
    lambda2 = -24.0 * (sr_cut6 - 2.0 * sr_cut12) * sr_cut
    sr2_ovr = 1.77

    diameter = 1.327441  #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq = r_cut^2                   # Potential cutoff squared in sigma=1 units

    """
     ri   = deepcopy(system.rm[i])
     rMol = deepcopy(system.rm)
     rMol = deleteat!(rMol,i)

     startAtom = deepcopy(system.thisMol_theseAtoms[i][1])
     endAtom   = deepcopy(system.thisMol_theseAtoms[i][2])

     sF = deepcopy(system.thisMol_theseAtoms)
     sF = deleteat!(sF,i)

     ra = deepcopy(system.ra[startAtom:endAtom])
     rb = deepcopy(system.ra)
     rb = deleteat!(rb,startAtom:endAtom)
     #println(length(system.rm))
    """

    ri = system.rm[i]
    rMol = system.rm

    startAtom = system.thisMol_theseAtoms[i][1]
    endAtom = system.thisMol_theseAtoms[i][2]

    sF = system.thisMol_theseAtoms

    ra = system.ra[startAtom:endAtom]
    rb = system.ra

    ϵ = system.ϵ
    σ = system.σ
    box = system.box
    rcut_sq = system.r_cut^2

    rij = rab = SVector{3}(0.0, 0.0, 0.0)

    pot, vir = 0.0, 0.0

    for (j, rj) in enumerate(rMol)
        if j == i
            continue
        end

        @inbounds for k = 1:3
            rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij_sq = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

        if rij_sq < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
            for a = 1:3
                """Loop over all atoms in molecule B"""
                for b = sF[j][1]:sF[j][2]
                    if (sF[j][2] - sF[j][1]) != 2
                        println("INdex wrong for molecule B")
                    end

                    @inbounds for k = 1:3
                        rab = @set rab[k] = vector1D(ra[a][k], rb[b][k], box)
                    end

                    rab_sq = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    sr2 = 1.0 / rab_sq             # (sigma/rab)**2
                    ovr = sr2 > sr2_ovr
                    if ovr
                        println("OVERLAP, OVERLAP!!!!!")
                    end
                    if rab_sq < r_cut_sq
                        #println(i, " ", j)
                        sr2 = 1.0 / rab_sq
                        rmag = sqrt(rab_sq)
                        sr6 = sr2^3
                        sr12 = sr6^2
                        pot += 4 * (sr12 - sr6) + lambda1 + lambda2 * rmag
                        virab = 24.0 * (2.0 * sr12 - sr6) - lambda2 * rmag
                        fab = rab * virab * sr2
                        vir += dot(rij, fab)

                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end
    end

    return pot, vir / 3.0
end # lj_poly_ΔU

""" LJ energy of particle i with all other particles. Real units of σ and ϵ """
function LJ_poly_ΔU(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    ri = system.rm[i]
    rMol = system.rm

    startAtom = system.thisMol_theseAtoms[i][1]
    endAtom = system.thisMol_theseAtoms[i][2]
    list_type = system.atomTypes
    moli_type = system.atomTypes[startAtom:endAtom]

    sF = system.thisMol_theseAtoms

    ra = system.ra[startAtom:endAtom]
    rb = system.ra

    table = system.table
    box = system.box
    r_cut = system.r_cut
    diameter = 0 #r_cut * 0.25 + 5 #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq = r_cut^2                   # Potential cutoff squared in sigma=1 units

    rᵢⱼ = rab = SVector{3}(0.0, 0.0, 0.0)

    pot, vir = 0.0, 0.0

    # cycle through all molecules
    @inbounds for (j, rj) in enumerate(rMol)
        if j == i
            continue # skip if molecule j is same as i
        end

        """Molecular mirror image separation"""
        @inbounds for k = 1:3
            rᵢⱼ = @set rᵢⱼ[k] = vector1D(ri[k], rj[k], box)
        end

        rᵢⱼ² = rᵢⱼ[1] * rᵢⱼ[1] + rᵢⱼ[2] * rᵢⱼ[2] + rᵢⱼ[3] * rᵢⱼ[3]

        if rᵢⱼ² < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
            @inbounds for a = 1:3

                """Loop over all atoms in molecule B"""
                @inbounds for b = sF[j][1]:sF[j][2]

                    """Atomic mirror image separation"""
                    @inbounds for k = 1:3
                        rab = @set rab[k] = vector1D(ra[a][k], rb[b][k], box)
                    end

                    rab² = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    ϵᵢⱼ = table.ϵᵢⱼ[moli_type[a], list_type[b]]
                    if rab² < (r_cut_sq + 100) && ϵᵢⱼ > 0.001
                        # this uses only molecular cutoff

                        σᵢⱼ = table.σᵢⱼ[moli_type[a], list_type[b]]

                        σ² = σᵢⱼ^2 / rab²             # (sigma/rab)**2
                        σ⁶ = σ²^3
                        σ¹² = σ⁶^2
                        pot += ϵᵢⱼ * (σ¹² - σ⁶)
                        virab = ϵᵢⱼ * (2.0 * σ¹² - σ⁶)
                        fab = rab * virab * σ²
                        vir += dot(rᵢⱼ, fab)

                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end
    end

    return pot * 4, vir * 24 / 3.0
end # lj_poly_ΔU
################################################################################
"""polyatomic LJ with mos and soa"""
function LJ_poly_ΔU(i, moa::StructArray, soa::StructArray,
                            vdwTable, r_cut, box)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    ri = moa.COM[i]   # molecule i
    rMol = moa.COM    # all molecules

    startAtom = moa.firstAtom[i] # first atom in molecule i
    endAtom = moa.lastAtom[i]    # lastatom in molecule i
    list_type = soa.atype
    moli_type = soa.atype[startAtom:endAtom]

    sFs, sFe = moa.firstAtom, moa.lastAtom

    ra = soa.coords[startAtom:endAtom]
    rb = soa.coords

    table = vdwTable
    box = box
    r_cut = r_cut
    diameter = 0 #r_cut * 0.25 + 5 #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq = r_cut^2                   # Potential cutoff squared in sigma=1 units

    rᵢⱼ = rab = SVector{3}(0.0, 0.0, 0.0)

    pot, vir = 0.0, 0.0

    # cycle through all molecules
     for (j, rj) in enumerate(rMol)
        if j == i
            continue # skip if molecule j is same as i
        end

        """Molecular mirror image separation"""
        @inbounds for k = 1:3
            rᵢⱼ = @set rᵢⱼ[k] = vector1D(ri[k], rj[k], box)
        end

        rᵢⱼ² = rᵢⱼ[1] * rᵢⱼ[1] + rᵢⱼ[2] * rᵢⱼ[2] + rᵢⱼ[3] * rᵢⱼ[3]

        if rᵢⱼ² < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
             for a = 1:(endAtom-startAtom+1)

                """Loop over all atoms in molecule B"""
                for b = sFs[j]:sFe[j]

                    """Atomic mirror image separation"""
                    @inbounds for k = 1:3
                        rab = @set rab[k] = vector1D(ra[a][k], rb[b][k], box)
                    end

                    rab² = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    ϵᵢⱼ = table.ϵᵢⱼ[moli_type[a], list_type[b]]
                    if rab² < (r_cut_sq + 100) && ϵᵢⱼ > 0.001
                        # this uses only molecular cutoff

                        σᵢⱼ = table.σᵢⱼ[moli_type[a], list_type[b]]

                        σ² = σᵢⱼ^2 / rab²             # (sigma/rab)**2
                        σ⁶ = σ²^3
                        σ¹² = σ⁶^2
                        pot += ϵᵢⱼ * (σ¹² - σ⁶)
                        virab = ϵᵢⱼ * (2.0 * σ¹² - σ⁶)
                        fab = rab * virab * σ²
                        vir += dot(rᵢⱼ, fab)

                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end
    end

    return pot * 4, vir * 24 / 3.0
end # lj_poly_ΔU


""" Calculates LJ potential between particle 'i' and the other N-1 particles"""
function LJ_ΔU(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    r = system.r
    ϵ = system.ϵ
    σ = system.σ
    box = system.box
    diff = SVector{3}(0.0, 0.0, 0.0)
    rcut_sq = system.r_cut^2
    pot, vir = 0.0, 0.0

    for (j, atom) in enumerate(r)
        if j == i
            continue
        end

        @inbounds for k = 1:3
            diff = @set diff[k] = vector1D(r[i][k], atom[k], box)
        end

        rij_sq = diff[1] * diff[1] + diff[2] * diff[2] + diff[3] * diff[3]

        if rij_sq > rcut_sq
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
            sr6 = sr2^3
            sr12 = sr6^2
            pot += ϵ[j] * (sr12 - sr6)  # LJ pair potentials (cut but not shifted)
            vir += ϵ[j] * (2 * sr12 - sr6)  # LJ pair virials
        end

    end

    return pot * 4.0, vir * 24.0 / 3.0
end

""" Calculates total LJ potential energy of the system. Double Counts. """
function potential(system, tot)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0

    #########
    #
    #     Calculate LJ
    #
    ##################

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2

    return tot

end

""" Calculates total LJ && Electrostatic potential energy of the system.
 Double Counts. """
function potential(
    system,
    tot,
    ewalds::EWALD,
    qq_q::Vector,
    qq_r::Vector,
    #kfacs::Vector,
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = system.r_cut
    box = system.box
    thisMol_thisAtom = system.thisMol_theseAtoms
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real, recip = 0.0, 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real
    #ener = EwaldReal(qq_r, qq_q, ewald.kappa, box, thisMol_thisAtom, 2, system)
    #=
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    kappa::Real,
    box::Float64,
    thisMol_thisAtom::Vector,
    chosenOne::Int64,
    system::Requirements
    =#
    totReal = 0.0
    for i = 1:length(rm)
        ener, overlap =
            EwaldReal(qq_r, qq_q, ewalds.kappa, box, thisMol_thisAtom, i, system)
        totReal += ener
        if overlap
            exit()
        end
    end
    tot.energy += totReal / 2 * ewalds.factor # divide by 2 to account for double counting
    tot.coulomb += totReal / 2 * ewalds.factor
    reall = totReal / 2 * ewalds.factor

    println("Total real Ewald is: ", totReal / 2 * ewalds.factor)

    # Recipricol
    #@time recipEnergy = RecipLong(system, ewald, qq_r, qq_q) * ewald.factor #RecipLong(system, ewald, qq_r, qq_q, kfacs)
    recipEnergy, ewalds = RecipLong(system, ewalds, qq_r, qq_q)  #RecipLong(system, ewald, qq_r, qq_q, kfacs)
    recipEnergy *= ewalds.factor
    #@btime RecipLong(system, ewald, qq_r, qq_q) * ewald.factor
    #println("first: ", recipEnergy)

    println("Total recipricol Ewald is: ", recipEnergy)

    tot.energy += recipEnergy #* ewald.factor
    tot.coulomb += recipEnergy #* ewald.factor
    tot.recipOld = recipEnergy #* ewald.factor
    tot.recip = recipEnergy #* ewald.factor
    # Self

    println("Self energy: ", EwaldSelf(ewalds, qq_q))
    tot.energy += EwaldSelf(ewalds, qq_q)
    # tinfoil boundary
    #@time tinfoil = TinfoilBoundary(system, ewald, qq_q, qq_r)
    #println("tinfoil boundaries: ", tinfoil * ewald.factor)
    #tot.energy += tinfoil * ewald.factor
    #tot.coulomb += tinfoil * ewald.factor
    return tot, LJ, reall, recipEnergy, ewalds

end

"""Rotate each molecule X times, save the lowest energy rotation. Loop Through
All Molecules. Can do this loop Y times. """
function EnergyMinimize(thisSystem::Requirements, db, quatVec::Vector, qq_q, qq_r, ewald::EWALD)
    sys = deepcopy(thisSystem)
    qV = deepcopy(quatVec)
    nLoops = 5
    nMols = length(sys.rm)
    nTrials = 3
    ra = [SVector(0.0, 0.0, 0.0) for i = 1:3]
    println("EnergyMinimize in energy.jl, line 360, has some hardcoding")
    @inbounds for i = 1:nLoops
        for j = 1:nMols
            com = sys.rm[j]
            lowE, lowV = LJ_poly_ΔU(j, sys)
            partial_ewald_e, partial_ewald_v = EwaldShort(
                j, system,ewald, temperature, box, qq_r, qq_q, false )
            lowE += partial_ewald_e
            lowV += partial_ewald_v
            savedQuat = qV[j]
            for k = 1:nTrials
                ei = random_rotate_quaternion(0.05, savedQuat)
                ai = q_to_a(ei)
                for a = 1:at_per_mol # Loop over all atoms
                    ra[a] = com + MATMUL(ai, db[:, a])
                end # End loop over all atoms
                sys.ra[sys.thisMol_theseAtoms[j][1]:sys.thisMol_theseAtoms[j][2]] =
                    ra
                qq_r[sys.thisMol_theseAtoms[j][1]:sys.thisMol_theseAtoms[j][2]] =
                        ra
                newE, newV = LJ_poly_ΔU(j, sys)
                partial_ewald_e, partial_ewald_v = EwaldShort(
                    j, system,ewald, temperature, box, qq_r, qq_q, false )
                newE += partial_ewald_e
                newV += partial_ewald_v
                if newE < lowE
                    savedQuat = ei
                end
            end
            qV[j] = savedQuat
        end
        println("loop: ", i)
    end

    return deepcopy(qV)

end

function LennardJones(rij)
    return 4 * 1 * ((1 / rij)^12 - (1 / rij)^6)
end

# ================================================================================================================================
function press_corr(system::Requirements, num_atom_types = 2, b = [100, 200])
    # ================================================================================================================================

    # Adds lennard-jones or Buckingham (EXP-6) tail correction to pressure
    # ======================================================================================================
    #                 Correction to Each atom type
    # ======================================================================================================

    corp = 0.0
    corpi = 0.0
    Rc = system.r_cut
    vol = system.box^3
    for i = 1:num_atom_types # cycle through all of the atoms, this counts as being atom 1, cycle "k" is atom 2

        for j = 1:num_atom_types #for each atom, we must cycle through all of the molecules it interacts with

            #    			write(*,*) " coru pre-error error : ", i,coru, cori(i)
            epsij = system.table.ϵᵢⱼ[i, j]  #EPS( i, j )
            sigij = system.table.σᵢⱼ[i, j]  #SIG( i, j )

            # ============================================================================
            #					Lennard-Jones Force-Field
            # ============================================================================

            sig3 = (sigij * sigij * sigij)
            sigor3 = sig3 / (Rc * Rc * Rc)
            sigor9 = sigor3 * sigor3 * sigor3

            corp =
                corp +
                b[i] * b[j] * epsij * sig3 * ((2.0 / 3.0) * sigor9 - sigor3)
            #corpi(i) = corpi(i) + b(i)*b(j) * epsij * sig3 * ( (2.0/3.0) * sigor9 - sigor3 ) # individual molecule corrections (Not used, it is here to trick you)

        end

    end


    # ============================================================================
    #       Add on constants kept out of loop for speed purposes
    # ============================================================================

    #if(FF_Flag .eq. 1) then # use Lennard-Jones force-field

    corp = 16 * π / (3.0 * vol * vol) * corp
    #corpi = 16*π/(3.0*vol*vol) * corpi
    return corp
    #elseif(FF_Flag .eq. 2) then # use EXP-6 force-field

    #    corp  = 4.0 * 3.141592654 / (3.0 * vol * vol) * corp
    #    corpi = 4.0 * 3.141592654 / (3.0 * vol * vol) * corpi

    #end
end #Subroutine press_corr

# ================================================================================================================================
function ener_corr(system::Requirements, num_atom_types = 2, b = [100, 200])
    # ================================================================================================================================

    # -Adds lennard-jones or EXP-6 tail correction to potential energy.
    # -Works for polyatomic molecules. Calculates total potential for a pure
    # monatomic molecule as well as individual species tail corrections for
    # mixtures.

    # ============================================================
    #         Written by Braden Kelly a.k.a. Zarathustra sometime
    #                  prior to your reading this.
    # ============================================================

    coru = 0.0
    Rc = system.r_cut
    vol = system.box^3
    #println(system.table.ϵᵢⱼ)
    for i = 1:num_atom_types # cycle through all of the atoms, this counts as being atom 1, cycle "k" is atom 2

        for j = 1:num_atom_types #for each atom, we must cycle through all of the molecules it interacts with

            #    			write(*,*) " coru pre-error error : ", i,coru, cori(i)
            epsij = system.table.ϵᵢⱼ[i, j]  #EPS( i, j )
            sigij = system.table.σᵢⱼ[i, j]  #SIG( i, j )

            # ============================================================================
            #					Lennard-Jones Force-Field
            # ============================================================================


            sig3 = (sigij * sigij * sigij)
            sigor3 = sig3 / (Rc * Rc * Rc)
            sigor9 = sigor3 * sigor3 * sigor3

            coru += b[i] * b[j] * epsij * sig3 * ((1.0 / 3.0) * sigor9 - sigor3)

        end

    end

    coru = 8.0 * π / (3.0 * vol) * coru
    print("ener_corr value: ", coru, "vol: ", vol)
    return coru

end #Subroutine ener_corr

"Real Ewald contribution with molecular cutoff"
#function EwaldReal(diff,coord1, coord2,q1,q2, L, rcut2, kappa)
function CoulombReal(
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    box::Float64,
    chosenOne::Int64,
    system::Requirements
)
    ####
    #
    #    Some prep stuff
    #    --------------------------------
    #    system: a struct with arrays for molecular and atomic coords and other stuff
    #    chosenOne: molecule that was moved
    #    qq_q: array of atomic partial charges
    #    qq_r: array of atomic coordinates. Each index is a Vector with x, y, z
    #
    #############
    ri = system.rm[chosenOne]                       # x, y, z coordinates for COM of molecule "chosenOne"
    thisMol_thisAtom = system.thisMol_theseAtoms    # length number of molecules
    start_a = thisMol_thisAtom[chosenOne][1]        # this atom starts molecule "chosenOne"
    end_a = thisMol_thisAtom[chosenOne][2]          # this atom ends molecule "chosenOne"
    rMol = system.rm                                # array of molecule COM coords
    r_cut = system.r_cut                            # site-site cutoff

    diameter = r_cut * 0.25  + 5.0        # made up diameter of molecule
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff
    rm_cut_box_sq = rm_cut_box^2          # Molecular cutoff squared
    r_cut_sq = r_cut^2                    # atomic cutoff squared

    # quick double check on atomic cutoffs
    @assert r_cut == 10.0
    #@assert r_cut_sq == 100.0

    pot = 0.0
    # Allocate for speed
    rab = SVector{3, Float64}(0.0, 0.0, 0.0)
    rij = SVector{3, Float64}(0.0, 0.0, 0.0)
    rj = SVector{3, Float64}(0.0, 0.0, 0.0)
    ra = SVector{3, Float64}(0.0, 0.0, 0.0)
    rb = SVector{3, Float64}(0.0, 0.0, 0.0)

    overlap = false
    ovr = 1.0

    for (j, rj) in enumerate(rMol) # cycle through all molecules
        if j == chosenOne # skip self interactions
            continue
        end

        #Molecular mirror image separation
        for k = 1:3
            rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij2 = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

        #Molecular cutoff check
        if rij2 < rm_cut_box_sq

            #Loop over all atoms in molecule A
            for a = start_a:end_a
                ra = qq_r[a]

                start_b = thisMol_thisAtom[j][1]  # first atom in molecule j
                end_b = thisMol_thisAtom[j][2]    # last atom in molecule j

                # Loop over all atoms in molecule B
                for b=start_b:end_b
                    rb = qq_r[b]       # coordinates of charge b

                    # Atomic mirror image separation
                    for k = 1:3
                        rab = @set rab[k] = vector1D(ra[k], rb[k], box)
                    end
                    rab2 = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    # Overlap Check
                    if (rab2 < ovr) && (qq_q[a] * qq_q[b] < 0)
                        println("got one")
                        return 0.0, true
                    # Atomic cutoff check
                    elseif rab2 < r_cut_sq
                        #println(a, "   ", b)
                        #rab_mag =
                        pot += qq_q[a] * qq_q[b] / sqrt(rab2)
                    else
                        pot += 0.0
                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end # if statement for molecular cutoff
    end # loop over molecules
    return pot , overlap
end

""" Calculates total LJ && Bare Coulomb potential energy of the system.
 Double Counts. """
function potential(
    system,
    tot,
    ewald::EWALD,
    qq_q::Vector,
    qq_r::Vector,
    string,
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = system.r_cut
    #@assert r_cut == 10.0
    box = system.box
    thisMol_thisAtom = system.thisMol_theseAtoms
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real

    totReal = 0.0
    for i = 1:length(rm)
        ener, overlap =
            CoulombReal(qq_r, qq_q, box, i, system)
        totReal += ener
        if overlap
            exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb += totReal
    reall = totReal
    println("Coulomb contribution is: ", reall)


    return tot, LJ, reall

end

""" Calculates total LJ && Wolf Summation potential energy of the system.
 Double Counts. """
function potential(
    system,
    tot,
    ewald::EWALD,
    qq_q::Vector,
    qq_r::Vector,
    string,
    string2 #triggers wolf summations using double strings
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = system.r_cut
    #@assert r_cut == 10.0
    box = system.box
    thisMol_thisAtom = system.thisMol_theseAtoms
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real
    totReal = 0.0
    for i = 1:length(rm)
        ener, overlap =
            EwaldReal(qq_r, qq_q,ewald.kappa, box,thisMol_thisAtom, i, system)
        totReal += ener
        if overlap
            exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb += totReal
    reall = totReal
    println("Coulomb contribution is: ", reall)

    prefactor = 0.0
    for i=1:length(qq_q)
        for j=1:length(qq_q)
            prefactor += qq_q[i]*qq_q[j] * erfc(ewald.kappa*r_cut) / r_cut
        end
    end
    prefactor *= -1
    prefactor2 = (erfc(ewald.kappa*r_cut) / 2 / r_cut +ewald.kappa / sqrt(π) ) *
                    dot(qq_q,qq_q)
    tot.energy += (prefactor - prefactor2)*ewald.factor
    tot.coulomb += (prefactor - prefactor2)*ewald.factor

    println("total energy is: ", tot.energy)
    println("total coulomb is: ", tot.coulomb)
    println("total LJ energy is: ", tot.energy - tot.coulomb)


    return tot

end

"""Calculate total potential energy using moa and soa for WolfSummation"""
function potential(
    moa::StructArray,
    soa::StructArray,
    tot::Properties,
    ewald::EWALD,
    vdwTable::Tables,
    sim_props::Properties2#triggers wolf summations using double strings
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = sim_props.LJ_rcut
    #@assert r_cut == 10.0
    box = sim_props.box
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(moa.COM)
        ener, vir = LJ_poly_ΔU(i, moa, soa, vdwTable, r_cut, box)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real
    totReal = 0.0
    for i = 1:length(moa.COM)
        ener, overlap =
            #EwaldReal(qq_r, qq_q,ewald.kappa, box,thisMol_thisAtom, i, system)
             EwaldReal(i, moa, soa, ewald,totProps.qq_rcut, box)
        totReal += ener
        if overlap
            println("overlap after EwaldReal")
            #exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb += totReal
    reall = totReal
    println("Coulomb contribution is: ", reall)

    prefactor = 0.0
    for i=1:length(soa.charge)
        for j=1:length(soa.charge)
            prefactor += soa.charge[i]*soa.charge[j] * erfc(ewald.kappa*r_cut) / r_cut
        end
    end
    prefactor *= -1
    prefactor2 = (erfc(ewald.kappa*r_cut) / 2 / r_cut +ewald.kappa / sqrt(π) ) *
                    dot(soa.charge,soa.charge)
    tot.energy += (prefactor - prefactor2)*ewald.factor
    tot.coulomb += (prefactor - prefactor2)*ewald.factor

    println("total energy is: ", tot.energy)
    println("total coulomb is: ", tot.coulomb)
    println("total LJ energy is: ", tot.energy - tot.coulomb)


    return tot

end

"""Calculate total potential energy using moa and soa for Ewald Summation"""
function potential(
    moa::StructArray,
    soa::StructArray,
    tot::Properties,
    ewalds::EWALD,
    vdwTable::Tables,
    sim_props::Properties2,
    coulomb_style::String #triggers wolf summations using double strings
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = sim_props.LJ_rcut
    #@assert r_cut == 10.0
    box = sim_props.box
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(moa.COM)
        ener, vir = LJ_poly_ΔU(i, moa, soa, vdwTable, r_cut, box)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real
    totReal = 0.0
    for i = 1:length(moa.COM)
        ener, overlap =
            #EwaldReal(qq_r, qq_q,ewald.kappa, box,thisMol_thisAtom, i, system)
             EwaldReal(i, moa, soa, ewald,totProps.qq_rcut, box)
        totReal += ener
        if overlap
            println("overlap after EwaldReal")
            #exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb += totReal
    tot.virial += totReal / 3.0
    reall = totReal
    println("Coulomb contribution is: ", reall)

    recipEnergy, ewalds = RecipLong(ewalds, soa.coords, soa.charge, box)  #RecipLong(system, ewald, qq_r, qq_q, kfacs)
    recipEnergy *= ewalds.factor


    println("Total recipricol Ewald is: ", recipEnergy)
    tot.energy += recipEnergy   # divide by 2 to account for double counting
    tot.coulomb += recipEnergy
    tot.virial += recipEnergy / 3.0

    selfEnergy = EwaldSelf(ewalds, soa.charge)
    println("Self energy: ", selfEnergy)
    tot.energy += selfEnergy
    tot.coulomb += selfEnergy
    tot.virial += selfEnergy / 3.0



    println("total energy is: ", tot.energy)
    println("total coulomb is: ", tot.coulomb)
    println("total LJ energy is: ", tot.energy - tot.coulomb)


    return tot

end
