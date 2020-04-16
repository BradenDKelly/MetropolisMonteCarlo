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

Random.seed!(11234)
start=Dates.now()
println(Dates.now())

################################################################################
# TODO (BDK) create JSON input file with starting parameters:
# TODO (BDK) add proper sampling
# TODO (BDK) Add some timing (~ 30 times faster than numpy)
# TODO (BDK) write unit test to make sure atoms and molecules COM match
# TODO (BDK)
# until then, manually enter at top
################################################################################
temperature = 0.6 #0.8772  # 1.2996
ρ = 0.1832655
nMol = 256
nAtoms = nMol * 3
#ϵ = 1.0
#σ = 1.0
r_cut = 2.5  # box / 2
nSteps = 20000
nblock = 10
outputInterval = 100
initialConfiguration = "cnf"  # place atoms in a crystal structure
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

"""Returns the center of mass of an array of atoms and array of masses"""
function COM(atoms,masses)
    totalMass = sum(masses)
    numerator = sum(atoms .* masses)
    println("From inside COM: ", numerator ./ totalMass)
    return numerator ./ totalMass
end
#COM([[1,2,3],[2,3,4],[0,1,2]],[1,1,100])


function test_two_LJ_triangles()

 """ two dummy molecules each of 3 atoms, interact via LJ"""
    alpha = 75.0 * π / 180.0
    alpha2 = alpha / 2.0
    at_per_mol = 3
    boxSize = 1000

    db = reshape([-sin(alpha2), 0.0 , -cos(alpha2)/3.0 ,
                  0.0,          0.0 , 2*cos(alpha2)/3.0,
                  sin(alpha2),  0.0 , -cos(alpha2)/3.0],3,at_per_mol)

    a = []
    push!(a,SVector(db[:,1]...) )
    push!(a,SVector(db[:,2]...) )
    push!(a,SVector(db[:,3]...) )
    b_b = deepcopy(db) .+ [0,0,2]
    b = []
    push!(b,SVector(b_b[:,1]...) )
    push!(b,SVector(b_b[:,2]...) )
    push!(b,SVector(b_b[:,3]...) )

    mass = [1.,1.,1.]
    r=[]

    push!(r,SVector{3}(COM(a,mass)...))
    push!(r,SVector{3}(COM(b,mass)...))

    ra = []
    push!(ra,a[1])
    push!(ra,a[2])
    push!(ra,a[3])
    push!(ra,b[1])
    push!(ra,b[2])
    push!(ra,b[3])

    system_test    = Requirements(r, ra, [[1,3],[4,6]], ones(length(ra)),
                         ones(length(ra)), boxSize, boxSize/2)
    potential = 0
    potential += LennardJones(sqrt(sum( (ra[1] - ra[4]) .* (ra[1] - ra[4]) ) ) )  #1 -> 4
    potential += LennardJones(sqrt(sum( (ra[1] - ra[5]) .* (ra[1] - ra[5]) ) ) ) #1 -> 5
    potential += LennardJones(sqrt(sum( (ra[1] - ra[6]) .* (ra[1] - ra[6]) ) ) ) #1 -> 6
    potential += LennardJones(sqrt(sum( (ra[2] - ra[4]) .* (ra[2] - ra[4]) ) ) ) #2 -> 4
    potential += LennardJones(sqrt(sum( (ra[2] - ra[5]) .* (ra[2] - ra[5]) ) ) ) #2 -> 5
    potential += LennardJones(sqrt(sum( (ra[2] - ra[6]) .* (ra[2] - ra[6]) ) ) ) #2 -> 6
    potential += LennardJones(sqrt(sum( (ra[3] - ra[4]) .* (ra[3] - ra[4]) ) ) ) #3 -> 4
    potential += LennardJones(sqrt(sum( (ra[3] - ra[5]) .* (ra[3] - ra[5]) ) ) ) #3 -> 5
    potential += LennardJones(sqrt(sum( (ra[3] - ra[6]) .* (ra[3] - ra[6]) ) ) ) #3 -> 6
    calculated, virial = LJ_poly_ΔU(1,system_test)
    println("Testing two triangles made of LJ atoms")
    println("Presenting a: ")
    println(a)
    println("Presenting b: ")
    println(b)
    println("The answer is: ", potential)
    println("The test is  : ", calculated)
    println("1-4 ", sum( (ra[1] - ra[4]) .* (ra[1] - ra[4]) ))
    println("1-5 ", sum( (ra[1] - ra[5]) .* (ra[1] - ra[5]) ))
    println("1-6 ", sum( (ra[1] - ra[6]) .* (ra[1] - ra[6]) ))
    println("2-4 ", sum( (ra[2] - ra[4]) .* (ra[2] - ra[4]) ))
    println("2-5 ", sum( (ra[2] - ra[5]) .* (ra[2] - ra[5]) ))
    println("2-6 ", sum( (ra[2] - ra[6]) .* (ra[2] - ra[6]) ))
    println("3-4 ", sum( (ra[3] - ra[4]) .* (ra[3] - ra[4]) ))
    println("3-5 ", sum( (ra[3] - ra[5]) .* (ra[3] - ra[5]) ))
    println("3-6 ", sum( (ra[3] - ra[6]) .* (ra[3] - ra[6]) ))

    println("Warning, LJ may be cut and shifted - look inside LJ_poly_ΔU")

    if abs(potential - calculated) < 0.0001
        println("LJ_poly_ΔU Passed the test")
    else
        println("LJ_poly_ΔU Fails the test")
    end

end
#test_two_LJ_triangles()

"""Unit test for the center-of-mass calculation"""
function test_COM()
    answer = [1., 2., 3.]
    coords = [SVector([1, 2,3]...),SVector([2,3,4]...),SVector([0,1,2]...)]
    masses = SVector([1,1,1]...)
    result = COM(coords, masses)
    println("Testing COM calculation")
    println("Testing center of mass calculation")
    println("sending [[1,2,3],[2,3,4],[0,1,2]] as coords")
    println("sending [1,1,1] as masses")
    println("Answer is:      ", answer)
    println("Calculation is: ", result)
    if abs(sum( (1. ./ answer) .* result) - 3) < 0.0001
        println("Test for COM calculation Passed")
    else
        println("Test for COM calculation Failed")
    end
end
#test_COM()

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
        mol = 1
        atomName = "O"
        f = 1.0  # scaling factor for coordinates
        for (i,atom) in enumerate(r)
            if i==2
                atomName = "H"
            end

            molName = "Sol"
            #atomName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].atomnm
            #molName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].resnm

            line = @sprintf("%-6s %4d %3s %4s %5d %3s %7.3f %7.3f %7.3f %5.2f %5.2f \n",
            "ATOM",i, atomName, molName, mol, " ", f*r[i][1],
            f * r[i][2],f * r[i][3], 1.00, 0.00   )
            write(file,line)
            if i % 3 == 0
                mol += 1
                atomName = "O"
            else
                atomName = "H"
            end
        end
    end
    println("finished making pdb")
end

function ReadCNF(input="cnf_input.inp")
    r = []
    e = []
    box=9.42953251
    open(input) do file
        for line in eachline(file)

            if length(split(line)) == 1 # && typeof(split(line)) == Float64
                box = 9.42953251 #parse(Float64, split(line)) # this should be on the 2nd Line
                println("hardcoded box at line 257 in ReadCNF")
            end

            if length(split(line)) > 4

                lin = split(line)

                push!(r, [ parse(Float64, strip(lin[1])),
                    parse(Float64, strip(lin[2])),
                    parse(Float64, strip(lin[3])) ] )

                push!(e, [ parse(Float64, strip(lin[4])),
                    parse(Float64, strip(lin[5])),
                    parse(Float64, strip(lin[6])),
                    parse(Float64, strip(lin[7])) ] )
            end
        end
    end
    return r, e, box
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

 #a = SMatrix{3,3}([ q[1]^2+q[2]^2-q[3]^2-q[4]^2   2*(q[2]*q[3]+q[1]*q[4])       2*(q[2]*q[4]-q[1]*q[3])   ; # 1st row
#               2*(q[2]*q[3]-q[1]*q[4])          q[1]^2-q[2]^2+q[3]^2-q[4]^2   2*(q[2]*q[4]+q[1]*q[2])     ; # 2nd row
#               2*(q[2]*q[4]+q[1]*q[3])          2*(q[3]*q[4]-q[1]*q[2])   q[1]^2-q[2]^2-q[3]^2+q[4]^2 ]) # 3rd row
  return [ q[1]^2+q[2]^2-q[3]^2-q[4]^2   2*(q[2]*q[3]+q[1]*q[4])       2*(q[2]*q[4]-q[1]*q[3])   ; # 1st row
                2*(q[2]*q[3]-q[1]*q[4])          q[1]^2-q[2]^2+q[3]^2-q[4]^2   2*(q[2]*q[4]+q[1]*q[2])     ; # 2nd row
                2*(q[2]*q[4]+q[1]*q[3])          2*(q[3]*q[4]-q[1]*q[2])   q[1]^2-q[2]^2-q[3]^2+q[4]^2 ] #a
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
    return [c0,c1,c2,c3] #SVector{4, Float64}(c0, c1, c2, c3)
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

 return [e0,e1,e2,e3] #SVector{4, Float64}(e0, e1, e2, e3)
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

function MIS(vec::SVector,box::Real)
    x, y, z = vec[1], vec[2], vec[3]
    
end

# Generate molecular COM coordinates
if lowercase(initialConfiguration) == "crystal"
    rm = InitCubicGrid(nMol,ρ)
elseif occursin(lowercase(initialConfiguration),"cnf") #"cnf"  lowercase(initialConfiguration)
    rm, quat, box = ReadCNF("cnf_input.inp")
    ρ = length(rm) / (box ^ 3)    
    println(" Initial configuration is from file: cnf_input.inp")


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
ϵ = ones(nAtoms)
σ = ones(nAtoms)

total     = Properties(0.0, 0.0, 0.0, 0.0)
system    = Requirements(rm, ra, thisMol_thisAtom, ϵ, σ, box, r_cut)

"""Rotate each molecule X times, save the lowest energy rotation. Loop Through
All Molecules. Can do this loop Y times. """
function EnergyMinimize(thisSystem::Requirements,db, quatVec::Vector)
    sys = deepcopy(thisSystem)
    qV = deepcopy(quatVec)
    nLoops = 100
    nMols = length(sys.rm)
    nTrials = 10
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

    return qV

end
totProps  = Properties2(temperature, ρ, Pressure(total, ρ, temperature, box^3),
                            dr_max, dϕ_max, 0.3, 0, 0, initQuaternions)

# totProps.quat = @time EnergyMinimize(system,db, totProps.quat)

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

"""
function VolumeChange()
    #Subroutine MC_vol(vol_attempt,vol_accept,vmax,L)
    #!coords,atom_coords,,
    # !use MT19937 ! Uses the Mersenne Twister random number generator grnd()
    #use random_generators
    #use Global_Variables
    #use energy_routines

    #implicit none

    #integer, intent(inout)                   :: vol_accept, vol_attempt
    #real, intent(inout)                      :: L
    #! real, dimension(num_mol)			     :: Xb,Xnew,Yb,Ynew,Zb,Znew
    #real, dimension(3,num_mol)		         :: coords_new
    #real, intent(in)                         :: vmax
    #!integer, intent(in)                      :: Em
    #real                                     :: vol_old,vol_new
    #!real                                     :: lnvol_new
    #real                                     :: L_new,L_old,energy_old,vir_old
    #real                                     :: test,f,energy_new,vir_new,rho_new,coru_new
    #!integer                                  :: Nb
    #!real                                     :: BoxSize
    #integer	                                 :: i,j, start,endd,begin,finish

    # ! cori is a vector containing the energy tail correction for each molecule type
    #! TypeAtom is a matrix containing the types of each atom in each molecule
    #! maxlenthmol is the number of atoms in the largest molecule
    #! atoms_per_molecule a vector containing the number of atoms in each type of molecule
    #! i.e. if CO2 is molecule type #1 then the first entry in atoms_per_molecule is 3 (C+O+O = 3 atoms)
    #! real, dimension(num_atom_types)        :: cori
    #! integer, dimension(maxlengthmol,nm)     :: TypeAtom
    #! integer	                             :: maxlengthmol
    #! integer, dimension(num_mol_types)                  :: atoms_per_molecule
    #! integer, intent(in)									:: FF_Flag ! 1 = Lennard-Jones, 2 = EXP-6
    #! ALPHA is the matrix containing the cross parameters of the EXP-6 ALPHA parameters
    #!real, dimension (num_atom_types,num_atom_types), intent(in)							:: ALPHA
    #real	                        :: enn, virn
    #real, dimension(3,num_mol)                               :: change_in_coords
    #real, dimension(3,num_atoms)                             ::atom_XYZ !,atom_YT,atom_ZT
    #!==============================================================================================

      vol_attempt = vol_attempt + 1

      L_old = L
      vir_old = vir
      energy_old = energy

      vol_old = L_old^3

    #!  lnvol_new = LOG(vol_old) + (rranf()-0.5)*LOG(vmax) #! if using lnV modify acceptance to (num_mol + 1)
      vol_new = vol_old + (rand()-0.5)*vmax #!EXP(lnvol_new)
      L_new = vol_new^(1.0/3.0)

      f = L_new/L_old #! Scaling factor

      coords_new(:,:) = f*coords(:,:)

     change_in_coords(1,:) = coords_new(1,:) - coords(1,:)
     change_in_coords(2,:) = coords_new(2,:) - coords(2,:)
     change_in_coords(3,:) = coords_new(3,:) - coords(3,:)

      for i = 1:num_mol

        start = ThisMol_startAtom_endAtom(i,1) #! from the array of atoms, this is the atom that starts molecule 1
        finish = ThisMol_startAtom_endAtom(i,2)

         for j = start:finish

             atom_XYZ(:,j) = atom_coords(:,j) + change_in_coords(:,i)

         end
     end
      #!Need to insert IF Statement for SHIFT
      start = 0;finish=0;endd=0;begin=0
      energy_new = 0.0
      vir_new = 0.0
      enn = 0.0
      virn = 0.0
     #!=======================================================================================================
     #!				Calculate total energy of new configuration
     #!=======================================================================================================

       for i = 1:num_mol-1

        start = ThisMol_startAtom_endAtom(i,1) #! from the array of atoms, this is the atom that starts molecule 1
        finish = ThisMol_startAtom_endAtom(i,2) #! from the array of atoms, this is the atom that ends molecule 1

        for j = i+1:num_mol

    #!        if(i == j) cycle ! That's right, skidaddle if you can't be different from i!

            begin = ThisMol_startAtom_endAtom(j,1) #! from the array of atoms, this is the atom that starts molecule 2
            endd = ThisMol_startAtom_endAtom(j,2) #! from the array of atoms, this is the atom that ends molecule 2

            call ener_single(finish-start+1,atom_XYZ(:,start:finish),coords_new(:,i),endd-begin+1,
                atom_XYZ(:,begin:endd),coords_new(:,j),L_new,enn,virn,
    	        atom_type_list(start:finish),atom_type_list(begin:endd)) #!energy of new coord

            energy_new = energy_new + enn
            vir_new = vir_new + virn

        end
     end


     #!===============================================================================================
     rho_new = num_mol/vol_new
     #!===============================================================================================
     #!				Add tail correction to potential energy
     #!===============================================================================================

      if( cut_type == 'tail_corr')
         call ener_corr(rho_new,vol_new,coru_new)
         energy_new = energy_new + coru_new
      end


     #!=====================================================================================================
     #! 					Acceptance criteria
     #!=====================================================================================================
      test = exp(-beta*(P*(vol_new-vol_old)-(num_mol*log(vol_new/vol_old)/beta)
           + (energy_new-energy_old) ) )

      if (rand() < test)
          #accept volume move
          vol_accept = vol_accept + 1

          if( cut_type == 'tail_corr')
             energy_new = energy_new - coru_new
          elseif( cut_type(1:9) == 'cut_shift') then
              energy = energy_new
          end
          vir = vir_new
          #Rescale box
     #!     f = L_new/L_old ! Scaling factor
          coords(:,:) = f*coords(:,:)
          atom_coords(:,:) = atom_XYZ(:,:)
          L = L_new
      end

end #Subroutine MC_vol
"""
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
println("difference: Seconds: ", difference)

"""
println("difference: Seconds: ", parse(Float64,split(difference)[1]) / 1000)
println("difference: Minutes: ", parse(Float64,split(difference)[1]) / 1000 / 60)
println("difference: Hours  : ", parse(Float64,split(difference)[1]) / 1000 / 3600)
"""
