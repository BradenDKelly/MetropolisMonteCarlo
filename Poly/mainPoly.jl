
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
# TODO (BDK) write unit test to make sure atoms and molecules COM match
# until then, manually enter at top
################################################################################
temperature = 1.0 #0.8772  # 1.2996
ρ = 0.32655
nMol = 256
nAtoms = nMol * 3
#ϵ = 1.0
#σ = 1.0
r_cut = 2.5  # box / 2
nSteps = 1000
nblock = 100
outputInterval = 100
initialConfiguration = "crystal"  # place atoms in a crystal structure
dϕ_max = 0.05
dr_max = 0.05

################################################################################
#
# ! Bond vectors in body-fixed frame (na and db are public)
# ! Isosceles triangle, 3 sites, with unit bond length and bond angle alpha, which we set to 75 degrees here

alpha = 75.0 * π / 180.0
alpha2 = alpha / 2.0
at_per_mol = 3
"""
 REAL, DIMENSION(3,na), PARAMETER, PUBLIC :: db = RESHAPE ( [ &
      & -SIN(alpha2), 0.0,    -COS(alpha2)/3.0, &
      &  0.0,         0.0, 2.0*COS(alpha2)/3.0, &
      &  SIN(alpha2), 0.0,    -COS(alpha2)/3.0 ], [3,na] )
"""
db = reshape([-sin(alpha2), 0.0 , -cos(alpha2)/3.0 ,
              0.0,          0.0 , 2*cos(alpha2)/3.0,
              sin(alpha2),  0.0 , -cos(alpha2)/3.0],3,at_per_mol)

println(db[:,1])

"""Returns the Fortran equivalent of the same name"""
function MATMUL(ai,db)
    # In: ai:3x3 SMatrix
    #     db: 3x1 vector
    # Out: SVector{3}
    return SVector(dot(db,ai[:,1]),dot(db,ai[:,2]),dot(db,ai[:,3]) )
end

function test_quaternion(db)
    println("Here is the fortran matrix for db")
    println("-0.608761489       0.00000000      0.608761489")
    println("0.00000000       0.00000000       0.00000000" )
    println("-0.264451116      0.528902233     -0.264451116 ")
    println("Here is the Julia matrix for db")
    println(db[1,:])
    println(db[2,:])
    println(db[3,:])
    println("Testing Quaternion MATMUL")
    println("Testing MATMUL(mat_a,mat_b) as per fortran")
    println("Calculated in Fortran as MATMUL(db(:,1),db)")
    answer = [0.440524936,-0.139868781, -0.300656140 ]

    #di(:,a) = MATMUL ( db(:,a), ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai
    println("Method 1: transpose(db) * db[:,1]              ", transpose(db) * db[:,1])
    println("Method 2: db * db[:,1]                         ", db * db[:,1])
    println("Method 3: manual dot prod db[:,1] with db[:,j] ", MATMUL(db,db[:,1]))
    println("the answer as per Fortran is:                  ", answer)
end
#test_quaternion(db)

################################################################################

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

function q_to_a( q )

  # INTENT(out) DIMENSION(3,3) :: a ! Returns a 3x3 rotation matrix calculated from
  # INTENT(in) :: q ! supplied quaternion
  # The rows of the rotation matrix correspond to unit vectors of the molecule in the space-fixed frame
  # The third row  a(3,:) is "the" axis of the molecule, for uniaxial molecules
  # Use a to convert space-fixed to body-fixed axes thus: db = matmul(a,ds)
  # Use transpose of a to convert body-fixed to space-fixed axes thus: ds = matmul(db,a)
  # The supplied quaternion should be normalized and we check for this
  tol = 1.e-6
  norm = dot(q,q) # Quaternion squared length
  if  abs( norm - 1.0 ) > tol
     print( "quaternion normalization error ", norm, " ", tol)
     exit()
 end
  # Write out row by row, for clarity
  #a[1,:] = [ q[0]^2+q[1]^2-q[2]^2-q[3]^2,   2*(q[1]*q[2]+q[0]*q[3]),       2*(q[1]*q[3]-q[0]*q[2])     ] # 1st row
  #a[2,:] = [     2*(q[1]*q[2]-q[0]*q[3]),   q[0]^2-q[1]^2+q[2]^2-q[3]^2,   2*(q[2]*q[3]+q[0]*q[1])     ] # 2nd row
  #a[3,:] = [     2*(q[1]*q[3]+q[0]*q[2]),       2*(q[2]*q[3]-q[0]*q[1]),   q[0]^2-q[1]^2-q[2]^2+q[3]^2 ] # 3rd row

  #a = SMatrix{3,3}(q[1]^2+q[2]^2-q[3]^2-q[4]^2, 2*(q[2]*q[3]-q[1]*q[3]), 2*(q[2]*q[3]+q[1]*q[3]),
  #             2*(q[2]*q[3]+q[1]*q[3]),q[1]^2-q[2]^2+q[3]^2-q[3]^2, 2*(q[3]*q[3]-q[1]*q[2]),
  #             2*(q[2]*q[3]-q[1]*q[3]), 2*(q[3]*q[3]+q[1]*q[2]), q[1]^2-q[2]^2-q[3]^2+q[3]^2)

 a = SMatrix{3,3}([ q[1]^2+q[2]^2-q[3]^2-q[4]^2   2*(q[2]*q[3]+q[1]*q[4])       2*(q[2]*q[4]-q[1]*q[3])   ; # 1st row
               2*(q[2]*q[3]-q[1]*q[4])          q[1]^2-q[2]^2+q[3]^2-q[4]^2   2*(q[2]*q[4]+q[1]*q[2])     ; # 2nd row
               2*(q[2]*q[4]+q[1]*q[3])          2*(q[3]*q[4]-q[1]*q[2])   q[1]^2-q[2]^2-q[3]^2+q[4]^2 ]) # 3rd row
  return a
end #q_to_a

function random_vector()

  #REAL, DIMENSION(3) :: e ! Returns a uniformly sampled unit vector

  # The vector is chosen uniformly within the cube surrounding the unit sphere
  # Vectors lying outside the unit sphere are rejected
  # Having found a vector within the unit sphere, it is normalized
  #! Essentially the same routine will work in 2d, or for quaternions in 4d

  #REAL :: norm
 norm = 0.0
  while true# Loop until within unit sphere
     e = rand(Float64,3) # 3 random numbers uniformly sampled in range (0,1)
     e    = 2.0 .* e .- 1.0     # Now in range (-1,+1) i.e. within containing cube
     norm = dot( e, e )      # Square modulus
     if norm < 1.0  break end  # Within unit sphere
  end # End loop until within unit sphere

  e = e ./ sqrt( norm ) # Normalize
  return e
end # random_vector_1

function quatmul( a, b )

    #REAL, DIMENSION(0:3)             :: c    ! Returns quaternion product of
    #REAL, DIMENSION(0:3), INTENT(in) :: a, b ! two supplied quaternions

    #c(0) = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
    #c(1) = a(1)*b(0) + a(0)*b(1) - a(3)*b(2) + a(2)*b(3)
    #c(2) = a(2)*b(0) + a(3)*b(1) + a(0)*b(2) - a(1)*b(3)
    #c(3) = a(3)*b(0) - a(2)*b(1) + a(1)*b(2) + a(0)*b(3)

    c0 = a[1]*b[1] - a[2]*b[2] - a[3]*b[3] - a[4]*b[4]
    c1 = a[2]*b[1] + a[1]*b[2] - a[4]*b[3] + a[3]*b[4]
    c2 = a[3]*b[1] + a[4]*b[2] + a[1]*b[3] - a[2]*b[4]
    c3 = a[4]*b[1] - a[3]*b[2] + a[2]*b[3] + a[1]*b[4]
    return SVector{4, Float64}(c0, c1, c2, c3)
end #quatmul

function rotate_quaternion( angle, axis, old )
  #IMPLICIT NONE
  #REAL, DIMENSION(0:3)             :: e     ! Returns a quaternion rotated by
  #REAL,                 INTENT(in) :: angle ! specified rotation angle (in radians) about
  #REAL, DIMENSION(3),   INTENT(in) :: axis  ! specified rotation axis relative to
  #REAL, DIMENSION(0:3), INTENT(in) :: old   ! old quaternion

  #! Note that the axis vector should be normalized and we test for this
  #! In general, the old quaternion need not be normalized, and the same goes for the result
  #! although in our applications we only ever use unit quaternions (to represent orientations)

  #REAL                 :: norm
  #REAL, DIMENSION(0:3) :: rot
  tol = 1.e-6
  norm = dot( axis, axis ) #! Axis squared length
  if abs( norm - 1.0 ) > tol
     print("axis normalization error", norm," ", tol)
     exit()
 end
 rot = zeros(4)
  #! Standard formula for rotation quaternion, using half angles
  rot[1]   = cos(0.5*angle)
  rot[2:4] = sin(0.5*angle).* axis

  e = quatmul( rot, old ) # Apply rotation to old quaternion

  return e
end # rotate_quaternion

function random_quaternion()
  #IMPLICIT NONE
  #REAL, DIMENSION(0:3) :: e ! Returns a uniformly sampled unit quaternion

  #REAL, DIMENSION(2) :: ζ
  #REAL               :: norm1, norm2, f
  ζ = zeros(2)
  norm1 = norm2 = 0.0
  while true #! Loop until within unit disk
     ζ = rand(Float64,2) #! Two uniform random numbers between 0 and 1
     ζ = 2.0 .* ζ  .- 1.0     #! Now each between -1 and 1
     norm1 = dot(ζ,ζ)       #! Squared magnitude
     if ( norm1 < 1.0 ) break end     #! Test for within unit disk
 end #! End loop until within unit disk

  e0 = ζ[1]
  e1 = ζ[2]

  while true# ! Loop until within unit disk
     ζ  = rand(Float64,2) #! Two uniform random numbers between 0 and 1
     ζ  = 2.0 .* ζ  .- 1.0     #! Now each between -1 and 1
     norm2 = dot(ζ,ζ)       #! Squared magnitude
     if ( norm2 < 1.0 ) break end   #! Test for within unit disk
 end #! End loop until within unit disk

  f = sqrt( ( 1.0 - norm1 ) / norm2 )
  e2 = ζ[1]*f
  e3 = ζ[2]*f
  
 return SVector{4, Float64}(e0, e1, e2, e3)
end # random_quaternion

function random_rotate_quaternion( angle_max, old )

  #REAL, DIMENSION(0:3)             :: e         ! Returns a unit quaternion rotated by a
  #REAL,                 INTENT(in) :: angle_max ! maximum angle (in radians) relative to
  #REAL, DIMENSION(0:3), INTENT(in) :: old       ! the old quaternion

  # Note that the reference quaternion should be normalized and we test for this

  #REAL, DIMENSION(3) :: axis
  #REAL               :: zeta, angle, norm
  tol = 1.e-6
  norm = dot(old, old ) # Old squared length
  if ( abs( norm - 1.0 ) > tol )
     print( "old normalization error", norm, tol)
     exit()
 end

  axis = random_vector()               # Choose random unit vector
  zeta = rand()           # Random number between 0 and 1
  angle = ( 2.0*zeta - 1.0 ) * angle_max # Uniform random angle in desired range

  e = rotate_quaternion( angle, axis, old ) # General rotation function

  return e
end # random_rotate_quaternion


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
    rm::Vector{SVector{3,Float64}}
    ra::Vector{SVector{3,Float64}}
    thisMol_theseAtoms::Vector{SVector{2,Int64}}
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
    quat::Vector{SVector{4,Float64}}
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

    diameter      = 1.0 #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box    = ( r_cut + diameter )       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq      = r_cut^2                   # Potential cutoff squared in sigma=1 units

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
    
    ϵ = system.ϵ
    σ = system.σ
    box = system.box
    rcut_sq = system.r_cut^2
        
    rij = rab = SVector{3}(0.0,0.0,0.0)

    pot, vir = 0.0, 0.0

    for (j,rj) in enumerate(rMol)

        @inbounds for k=1:3
             rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij_sq = rij[1]*rij[1] + rij[2]*rij[2] + rij[3]*rij[3]

        if rij_sq < rm_cut_box_sq

            for a = 1:3
                for b = sF[j][1]:sF[j][2]

                    @inbounds for k=1:3
                         rab = @set rab[k] = vector1D(ra[a][k], rb[b][k], box)
                    end

                    rab_sq = rab[1]*rab[1] + rab[2]*rab[2] + rab[3]*rab[3]

                    if rab_sq < r_cut_sq
                        if rab_sq < 0.9^2 rab_sq = 0.9^2 end
                        sr2  = 1.0 / rab_sq
                        rmag = sqrt(rab_sq)
                        sr6  = sr2 ^ 3
                        sr12 = sr6 ^ 2
                        pot  += (sr12-sr6) + lambda1 + lambda2*rmag 
                        virab = (2.0*sr12-sr6 ) - lambda2*rmag 
                        fab  = rab * virab * sr2 
                        vir  += dot( rij,fab)

                    end # potential cutoff
                end # loop over atom in molecule b 
            end # loop over atoms in molecule a
                

        end

    end

    return pot * 4.0, vir * 24.0 / 3.0
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

""" Calculates total potential energy of the system. Double Counts. """
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

"""unit test of monatomic lj"""
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

    test = Requirements(r,ones(3), ones(3), box_test, r_cut)
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
#test_LJ()

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
box = (nMol / ρ ) ^ (1 / 3)
 #box / 30   # at 256 particles, ρ=0.75, T=1.0 this is 48% acceptance

# Generate molecular COM coordinates
if lowercase(initialConfiguration) == "crystal"
    rm = InitCubicGrid(nMol,ρ)
else
    rm = [SVector{3,Float64}(rand(), rand(), rand()) .* box for i = 1:nMol]
end

# Generate atomic coordinates
ra = [SVector{3,Float64}(  0.0, 0.0, 0.0) for i = 1:nAtoms]
ra = []
finish = 0
thisMol_thisAtom = []
initQuaternions = []

for com in rm
    global finish
    start = finish + 1
    finish = start + 2
    ei = random_quaternion()
    push!(initQuaternions,ei)
    ai = q_to_a( ei ) # Rotation matrix for i
    for a = 1:at_per_mol # Loop over all atoms
       #di(:,a) = MATMUL ( db(:,a), ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai
       push!(ra,com + SVector( MATMUL(ai , db[:,a]) ))
    end # End loop over all atoms
    push!(thisMol_thisAtom, SVector(start, finish))
end


PrintPDB(ra, box, 0, "pdbOutput")
ϵ = ones(nAtoms)
σ = ones(nAtoms)

total     = Properties(0.0, 0.0, 0.0, 0.0)
system    = Requirements(rm, ra, thisMol_thisAtom, ϵ, σ, box, r_cut)
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
end # blk to nblock