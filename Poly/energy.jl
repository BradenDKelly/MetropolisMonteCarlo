#########
#
#                         Energy Routines
#
###############
""" Calculates LJ potential between particle 'i' and the other N-1 particles"""
function LJ_poly_ΔU(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    # Cutoff distance and force-shift parameters (all private) chosen as per the reference:
    # S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)
    r_cut   = 2.612 # in sigma=1 units, where r_cut = 1.2616 nm, sigma = 0.483 nm
    sr_cut  = 1.0/r_cut; sr_cut6 = sr_cut^6; sr_cut12 = sr_cut6^2
    lambda1 = 4.0*(7.0*sr_cut6-13.0*sr_cut12)
    lambda2 = -24.0*(sr_cut6-2.0*sr_cut12)*sr_cut
    sr2_ovr = 1.77

    diameter      = 1.327441  #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box    = ( r_cut + diameter )       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq      = r_cut^2                   # Potential cutoff squared in sigma=1 units

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

    ri   = system.rm[i]
    rMol = system.rm

    startAtom = system.thisMol_theseAtoms[i][1]
    endAtom   = system.thisMol_theseAtoms[i][2]

    sF = system.thisMol_theseAtoms

    ra = system.ra[startAtom:endAtom]
    rb = system.ra

    ϵ = system.ϵ
    σ = system.σ
    box = system.box
    rcut_sq = system.r_cut^2

    rij = rab = SVector{3}(0.0,0.0,0.0)

    pot, vir = 0.0, 0.0

    for (j,rj) in enumerate(rMol)

        if j == i continue end

        @inbounds for k=1:3
             rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij_sq = rij[1]*rij[1] + rij[2]*rij[2] + rij[3]*rij[3]

        if rij_sq < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
            for a = 1:3
                """Loop over all atoms in molecule B"""
                for b = sF[j][1]:sF[j][2]

                    if (sF[j][2]-sF[j][1]) != 2
                        println("INdex wrong for molecule B")
                    end

                    @inbounds for k=1:3
                         rab = @set rab[k] = vector1D(ra[a][k], rb[b][k], box)
                    end

                    rab_sq = rab[1]*rab[1] + rab[2]*rab[2] + rab[3]*rab[3]

                    sr2      = 1.0 / rab_sq             # (sigma/rab)**2
                    ovr      = sr2 > sr2_ovr
                    if ovr println("OVERLAP, OVERLAP!!!!!") end
                    if rab_sq < r_cut_sq
                        #println(i, " ", j)
                        sr2  = 1.0 / rab_sq
                        rmag = sqrt(rab_sq)
                        sr6  = sr2 ^ 3
                        sr12 = sr6 ^ 2
                        pot  += 4 * (sr12-sr6) + lambda1 + lambda2*rmag
                        virab = 24.0 * (2.0*sr12-sr6 ) - lambda2*rmag
                        fab  = rab * virab * sr2
                        vir  += dot(rij,fab)

                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end
    end

    return pot , vir / 3.0
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

""" Calculates total LJ potential energy of the system. Double Counts. """
function potential(system, tot)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    return tot

end

"""Rotate each molecule X times, save the lowest energy rotation. Loop Through
All Molecules. Can do this loop Y times. """
function EnergyMinimize(thisSystem::Requirements,db, quatVec::Vector)
    sys = deepcopy(thisSystem)
    qV = deepcopy(quatVec)
    nLoops = 500
    nMols = length(sys.rm)
    nTrials = 15
    ra = [SVector(0.,0.,0.) for i=1:3]
    @inbounds for i=1:nLoops
        for j=1:nMols
            com = sys.rm[j]
            lowE, lowV = LJ_poly_ΔU(j,sys)
            savedQuat = qV[j]
            for k=1:nTrials
                ei = random_rotate_quaternion(0.05,savedQuat)
                ai = q_to_a( ei )
                for a = 1:at_per_mol # Loop over all atoms
                   ra[a] = com +  MATMUL(ai , db[:,a])
                end # End loop over all atoms
                sys.ra[sys.thisMol_theseAtoms[j][1]:sys.thisMol_theseAtoms[j][2]] = ra
                newE, newV = LJ_poly_ΔU(j,sys)
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
    return 4 * 1 * ( (1/rij)^12 - (1/rij)^6 )
end
