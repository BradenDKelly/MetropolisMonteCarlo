
using BenchmarkTools
using StaticArrays
using OffsetArrays
using Profile

function test(a,b)
    Num=353
    num2 = 30
    energy = 0.0
    kxyz = a
    qq_q = b
    # do not read into this, it is dummy indices, not a pattern to be played with

    eikx = OffsetArray{Complex{Float64}}(undef, 1:num2,0:Num)
    eikx = OffsetArray{Complex{Float64}}(undef, 1:num2,-Num:Num)
    eikx = OffsetArray{Complex{Float64}}(undef, 1:num2,-Num:Num)
    # initialize manually
    for j=1:num2
        eikx[j, 0] = 1.0 + 0.0 * im
        eikx[j, 1] =
                cos(2 * π  * ( rand() / j ) ) + sin(2 * π * rand() / j ) * im
    end
    # Calculate via recursion
    for k = 2:(Num)
        @inbounds for j = 1:num2
            eikx[j, k] = eikx[j, k-1] * eikx[j, 1]
        end
    end
    energy = Recip(kxyz,eikx,qq_q, Num, num2)
    #=
    @inbounds for i=1:Num
        term = 0.0 + 0.0*im
        kx, ky, kz = kxyz[i][1], kxyz[i][2], kxyz[i][3]
        @inbounds for l = 1:num2
                term +=  qq_q[l] * ( eikx[l,kx] * eikx[l,ky]* eikx[l,kz]
                - 2*(eikx[l,kx] * eikx[l,ky]* eikx[l,kz]) )
                #term +=  eikx[l,kxyz[i,1]] * eikx[l,kxyz[i,2]]* eikx[l,kxyz[i,3]]
                #term +=  eikx[l,i] * eikx[l,i]* eikx[l,i]
        end
        energy += real(conj(term) * term)
    end
    =#
    return energy
end

function Recip(kxyz, eikx, qq_q, Num, num2)
    energy = 0.0
    @inbounds for i=1:Num
        term = 0.0 + 0.0*im
        kx, ky, kz = kxyz[i][1], kxyz[i][2], kxyz[i][3]
        @inbounds for l = 1:num2
                term +=  qq_q[l] * ( eikx[l,kx] * eikx[l,ky]* eikx[l,kz]
                - 2*(eikx[l,kx] * eikx[l,ky]* eikx[l,kz]) )
                #term +=  eikx[l,kxyz[i,1]] * eikx[l,kxyz[i,2]]* eikx[l,kxyz[i,3]]
                #term +=  eikx[l,i] * eikx[l,i]* eikx[l,i]
        end
        energy += real(conj(term) * term)
    end
    return energy
end


function test2(a,b)
    Num=353
    num2 = 30
    energy = 0.0
    kxyz = a
    qq_q = b
    # do not read into this, it is dummy indices, not a pattern to be played with

    eikx = OffsetArray{Complex{Float64}}(undef, 0:Num, 1:num2)
    # initialize manually
    for j=1:num2
        eikx[0, j] = 1.0 + 0.0 * im
        eikx[1, j] =
                cos(2 * π  * ( rand() / j ) ) + sin(2 * π * rand() / j ) * im
    end
    # Calculate via recursion
    for k = 2:(Num)
        @inbounds for j = 1:num2
            eikx[k, j] = eikx[k-1, j] * eikx[1, j]
        end
    end

    @inbounds for i=1:Num
        term = 0.0 + 0.0*im
        kx, ky, kz = kxyz[i][1], kxyz[i][2], kxyz[i][3]
        @inbounds for l = 1:num2
                term +=  qq_q[l] * ( eikx[kx,l] * eikx[ky,l]* eikx[kz,l]
                - 2*(eikx[kx,l] * eikx[ky,l]* eikx[kz,l]) )
                #term +=  eikx[l,kxyz[i,1]] * eikx[l,kxyz[i,2]]* eikx[l,kxyz[i,3]]
                #term +=  eikx[l,i] * eikx[l,i]* eikx[l,i]
        end
        energy += real(conj(term) * term)
    end
    return energy
end


Num=353
num2 = 30

kxyz = Vector{SVector{3,Int32}}(undef, Num)
#kxyz=[]
for i=1:Num
    kxyz[i] = SVector(i,i,i)
    #push!(kxyz,[i,i,i])
end
#kxyz = Tuple(kxyz)
qq_q = rand(num2)
test(kxyz,qq_q)
@btime test(kxyz,qq_q)
@btime test2(kxyz,qq_q)


nStep = 300
#@btime test(kxyz,qq_q)
#@time for step = 1:nStep
#    test(kxyz,qq_q)
#end
