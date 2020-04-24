#####
#
#                   Quaternion Routines
#
#####

function q_to_a(q)

    # INTENT(out) DIMENSION(3,3) :: a ! Returns a 3x3 rotation matrix calculated from
    # INTENT(in) :: q ! supplied quaternion
    # The rows of the rotation matrix correspond to unit vectors of the molecule in the space-fixed frame
    # The third row  a(3,:) is "the" axis of the molecule, for uniaxial molecules
    # Use a to convert space-fixed to body-fixed axes thus: db = matmul(a,ds)
    # Use transpose of a to convert body-fixed to space-fixed axes thus: ds = matmul(db,a)
    # The supplied quaternion should be normalized and we check for this
    tol = 1.e-6
    norm = dot(q, q) # Quaternion squared length
    if abs(norm - 1.0) > tol
        print("quaternion normalization error ", norm, " ", tol)
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
    return [
        q[1]^2 + q[2]^2 - q[3]^2 - q[4]^2 2 * (q[2] * q[3] + q[1] * q[4]) 2 * (
            q[2] * q[4] - q[1] * q[3]
        ) # 1st row
        2 * (q[2] * q[3] - q[1] * q[4]) q[1]^2 - q[2]^2 + q[3]^2 - q[4]^2 2 * (
            q[2] * q[4] + q[1] * q[2]
        ) # 2nd row
        2 * (q[2] * q[4] + q[1] * q[3]) 2 * (q[3] * q[4] - q[1] * q[2]) q[1]^2 -
                                                                        q[2]^2 -
                                                                        q[3]^2 +
                                                                        q[4]^2
    ] #a
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
        e = rand(Float64, 3) # 3 random numbers uniformly sampled in range (0,1)
        e = 2.0 .* e .- 1.0     # Now in range (-1,+1) i.e. within containing cube
        norm = dot(e, e)      # Square modulus
        if norm < 1.0
            break
        end  # Within unit sphere
    end # End loop until within unit sphere

    e = e ./ sqrt(norm) # Normalize
    return e
end # random_vector_1

function quatmul(a, b)

    #REAL, DIMENSION(0:3)             :: c    ! Returns quaternion product of
    #REAL, DIMENSION(0:3), INTENT(in) :: a, b ! two supplied quaternions

    #c(0) = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
    #c(1) = a(1)*b(0) + a(0)*b(1) - a(3)*b(2) + a(2)*b(3)
    #c(2) = a(2)*b(0) + a(3)*b(1) + a(0)*b(2) - a(1)*b(3)
    #c(3) = a(3)*b(0) - a(2)*b(1) + a(1)*b(2) + a(0)*b(3)

    c0 = a[1] * b[1] - a[2] * b[2] - a[3] * b[3] - a[4] * b[4]
    c1 = a[2] * b[1] + a[1] * b[2] - a[4] * b[3] + a[3] * b[4]
    c2 = a[3] * b[1] + a[4] * b[2] + a[1] * b[3] - a[2] * b[4]
    c3 = a[4] * b[1] - a[3] * b[2] + a[2] * b[3] + a[1] * b[4]
    return [c0, c1, c2, c3] #SVector{4, Float64}(c0, c1, c2, c3)
end #quatmul

function rotate_quaternion(angle, axis, old)
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
    norm = dot(axis, axis) #! Axis squared length
    if abs(norm - 1.0) > tol
        print("axis normalization error", norm, " ", tol)
        exit()
    end
    rot = zeros(4)
    #! Standard formula for rotation quaternion, using half angles
    rot[1] = cos(0.5 * angle)
    rot[2:4] = sin(0.5 * angle) .* axis

    e = quatmul(rot, old) # Apply rotation to old quaternion

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
        ζ = rand(Float64, 2) #! Two uniform random numbers between 0 and 1
        ζ = 2.0 .* ζ .- 1.0     #! Now each between -1 and 1
        norm1 = dot(ζ, ζ)       #! Squared magnitude
        if (norm1 < 1.0)
            break
        end     #! Test for within unit disk
    end #! End loop until within unit disk

    e0 = ζ[1]
    e1 = ζ[2]

    while true# ! Loop until within unit disk
        ζ = rand(Float64, 2) #! Two uniform random numbers between 0 and 1
        ζ = 2.0 .* ζ .- 1.0     #! Now each between -1 and 1
        norm2 = dot(ζ, ζ)       #! Squared magnitude
        if (norm2 < 1.0)
            break
        end   #! Test for within unit disk
    end #! End loop until within unit disk

    f = sqrt((1.0 - norm1) / norm2)
    e2 = ζ[1] * f
    e3 = ζ[2] * f

    return [e0, e1, e2, e3] #SVector{4, Float64}(e0, e1, e2, e3)
end # random_quaternion

function random_rotate_quaternion(angle_max, old)

    #REAL, DIMENSION(0:3)             :: e         ! Returns a unit quaternion rotated by a
    #REAL,                 INTENT(in) :: angle_max ! maximum angle (in radians) relative to
    #REAL, DIMENSION(0:3), INTENT(in) :: old       ! the old quaternion

    # Note that the reference quaternion should be normalized and we test for this

    #REAL, DIMENSION(3) :: axis
    #REAL               :: zeta, angle, norm
    tol = 1.e-6
    norm = dot(old, old) # Old squared length
    if (abs(norm - 1.0) > tol)
        print("old normalization error", norm, tol)
        exit()
    end

    axis = random_vector()               # Choose random unit vector
    zeta = rand()           # Random number between 0 and 1
    angle = (2.0 * zeta - 1.0) * angle_max # Uniform random angle in desired range

    e = rotate_quaternion(angle, axis, old) # General rotation function

    return e
end # random_rotate_quaternion
