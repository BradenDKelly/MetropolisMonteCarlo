# Unit tests for certain operations:
#
#  1) test_two_LJ_triangles tests the LJ interaction of 2 triatomic Molecules
#
#  2) test_COM  tests that the center-of-mass is calculated accurately
#
#  3) test_quaternion tests that quaternions are calculated correctly
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

#= Overwriting to spc/e water =#
alpha = 75.0 * π / 180.0
alpha2 = alpha / 2.0
at_per_mol = 3
db = reshape([0.816, 0.0 , 0.577 ,
              0.0,  0.0 , 0.0,
              1.633,  0.0 , 0.0],3,at_per_mol)
#=
db = reshape([x1, y1 , z1 ,
              x2, y2 , z2 ,
              x3, y3 , z3 ],3,at_per_mol)
=#
