#######
#
#           periodic Boundaries, Mirror Image, etc...
#
#############

"Vector between two coordinate values, accounting for mirror image seperation"
@inline function vector1D(c1::Float64, c2::Float64, box_size::Float64)
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) : (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) : (c2 - c1 + box_size)
    end
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
