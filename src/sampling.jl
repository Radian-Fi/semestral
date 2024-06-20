using Makie

mutable struct random_index_struct
    # see http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/

    N::Int
    ind::Int  # index of low discrepancy sequence
    acoeff::Vector{Float64}

    function random_index_struct(n::Int)
        i = 0
        ϕ = (MathConstants.φ, 1.324717957244746, 1.2207440846057596)[n]
        a = ϕ.^(-(1:n))
        return new(n, i, a)
    end
end

s = nothing

function reset_point_index()
    global s = nothing    
end

"""
    get_new_point_index(resolution, mask)

Implemetation of pseudo-random sampling.
Adapted from the initial implementation of streamlines.
"""
function get_new_point_index(resolution, mask)
    # see http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    global s
    if isnothing(s)
        s = random_index_struct(length(resolution))
    end
    while true
        c = CartesianIndex(ntuple(s.N) do i
            j = ceil(Int, ((0.5 + s.acoeff[i]*s.ind) % 1)*resolution[i])
            clamp(j, 1, size(mask, i))
        end)
        s.ind += 1
        if mask[c]; return c; end
    end
end

mutable struct fpss_struct
    M::Array{Int64,3}

    function fpss_struct(resolution::Vec{3, Int64})
        m = ones(Int64, resolution...)
        for i in 1:(maximum(resolution)-1)÷2
            m[1+i:end-i, 1+i:end-i, 1+i:end-i] .+= ones(Int64, (resolution .- 2*i)...)
        end
        return new(m)
    end
end

fs = nothing

function reset_point_index_fpss()
    global fs = nothing
end

"""
    get_new_point_index_fpss(resolution, mask)

Approximation of the farthest point sampling (using Menhettan distance).
It should provide a slightly better sampling than the pseudo-random approach,
but it is also much slower (both asymptotically and in practice).
"""
function get_new_point_index_fpss(resolution, mask)
    global fs
    if isnothing(fs)
        fs = fpss_struct(resolution)
    end

    # update "Voronoi diagram" matrix M with new occupied cells (approximately)
    fs.M .*= mask
    indices = vcat(CartesianIndex.(-1:2:1, 0, 0), CartesianIndex.(0, -1:2:1, 0), CartesianIndex.(0, 0, -1:2:1))
    for iter in 1:(maximum(resolution)÷2)-1
        for i in CartesianIndices((2:resolution[1]-1, 2:resolution[2]-1, 2:resolution[3]-1))
            fs.M[i] = min(fs.M[i], minimum(fs.M[i .+ indices]) + 1)
        end
    end

    # find farthest point
    c = argmax(fs.M)

    return c
end