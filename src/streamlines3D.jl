using MakieCore
using MakieCore: automatic
using GLMakie
import GLMakie: plot, plot!, convert_arguments

include("cone.jl")      # include an implementation of a cone mesh/marker
include("sampling.jl")  # include an implementation of pseudo-random sampling
"""
This file contains a specialized implementation of 3D streamline visualization as a recipe for the GLMakie library.
This implementation is based on the original streamlines implementation by Moritz Schauer 
and borrows parts of the code from the library.
"""



"""
    streamlines3d(f::function, xinterval, yinterval; color = norm, kwargs...)

f must either accept `f(::Point)` or `f(x::Number, y::Number, z::Number)`.
f must return a Point3.

Example:
```julia
v(x::Point3{T}) where T = Point3f(x[3], 4*x[2], x[1])
streamlines3d(v, -4..4, -4..4, -4..4)
```
Mesh{3, Float32, Triangle}
## Implementation
See the function `streamlines3d_impl` for implementation details.
"""
@recipe(StreamLines3D, f, limits) do scene
    attr = Attributes(
        stepsize = 0.01,
        gridsize = (32, 32, 32),
        maxsteps = 500,
        color = norm,

        arrow_size = automatic,
        # arrow_head = automatic,  # replaced by a cone
        density = 1.0,
        quality = 16,

        linewidth = theme(scene, :linewidth),
        linecap = theme(scene, :linecap),
        joinstyle = theme(scene, :joinstyle),
        miter_limit = theme(scene, :miter_limit),
        linestyle = nothing,
    )
    MakieCore.colormap_attributes!(attr, theme(scene, :colormap))
    MakieCore.generic_plot_attributes!(attr)
    return attr
end

function convert_arguments(::Type{<: StreamLines3D}, f::Function, xrange, yrange, zrange)
    xmin, xmax = extrema(xrange)
    ymin, ymax = extrema(yrange)
    zmin, zmax = extrema(zrange)
    mini = Vec3(xmin, ymin, zmin)
    maxi = Vec3(xmax, ymax, zmax)
    return (f, Rect(mini, maxi .- mini))
end

"""
    streamlines3d_impl(CallType, f, limits::Rect{N, T}, resolutionND, stepsize)

Code adapted from an example implementation by Moritz Schauer (@mschauer)
from https://github.com/MakieOrg/Makie.jl/issues/355#issuecomment-504449775

Background: The algorithm puts an arrow somewhere and extends the
streamline in both directions from there. Then, it chooses a new
position (from the remaining ones), repeating the the exercise until the
streamline gets blocked, from which on a new starting point, the process
repeats.

So, ideally, the new starting points for streamlines are not too close to
current streamlines.

Links:

[Quasirandom sequences](http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/)
"""
function streamlines3d_impl(CallType, f, limits::Rect{N, T}, resolutionND, stepsize, maxsteps=500, dens=1.0, color_func = norm) where {N, T}
    # pase parameters and variables
    resolution = to_ndim(Vec{N, Int}, resolutionND, last(resolutionND))
    mask = trues(resolution...)  # unvisited squares
    arrow_pos = Point{N, Float32}[]
    arrow_dir = Vec{N, Float32}[]
    line_points = Point{N, Float32}[]
    _cfunc = x-> to_color(color_func(x))
    ColorType = typeof(_cfunc(Point{N,Float32}(0.0)))
    line_colors = ColorType[]
    colors = ColorType[]
    dt = Point{N, Float32}(stepsize)
    mini, maxi = minimum(limits), maximum(limits)
    r = ntuple(N) do i
        LinRange(mini[i], maxi[i], resolution[i] + 1)
    end
    #apply_f(x0) = type(f) <: Observable ? f(x0)[] : f(x0)
    apply_f(x0, P) = P <: Point ? f(x0) : f(x0...)

    reset_point_index()
    n_points = 0  # count visited squares
    while n_points < prod(resolution)*min(one(dens), dens)  # fill up to 100*dens% of markersize
        c = get_new_point_index(resolution, mask)

        x0 = Point{N}(ntuple(N) do i
            first(r[i]) + (c[i] - 0.5) * step(r[i])
        end)
        point = apply_f(x0, CallType)
        if typeof(point) <: Observable
            point = point[]
        end
        #println(point, " ", typeof(point))
        #println(apply_f, " ", typeof(apply_f))
        if !(point isa Point3)
            error("Function passed to streamlines3d must return Point3")
        end
        pnorm = norm(point)
        color = _cfunc(point)
        push!(arrow_pos, x0)
        push!(arrow_dir, point ./ pnorm)
        push!(colors, color)
        mask[c] = false
        n_points += 1
        for d in (-1, 1)
            n_linepoints = 1
            x = x0
            ccur = c
            push!(line_points, Point{N, Float32}(NaN), x)
            push!(line_colors, color, color)
            while x in limits && n_linepoints < maxsteps
                point = apply_f(x, CallType)
                if typeof(point) <: Observable
                    point = point[]
                end
                pnorm = norm(point)
                x = x .+ d .* dt .* point ./ pnorm
                if !(x in limits)
                    break
                end
                # WHAT? Why does point behave different from tuple in this
                # broadcast
                idx = CartesianIndex(searchsortedlast.(r, Tuple(x)))
                if idx != ccur
                    if !mask[idx]
                        break
                    end
                    mask[idx] = false
                    n_points += 1
                    ccur = idx
                end
                push!(line_points, x)
                push!(line_colors, _cfunc(point))
                n_linepoints += 1
            end
        end
    end

    return (
        arrow_pos,
        arrow_dir,
        line_points,
        colors,
        line_colors,
    )
end

function plot!(p::StreamLines3D)
    data = lift(p, p.f, p.limits, p.gridsize, p.stepsize, p.maxsteps, p.density, p.color) do f, limits, resolution, stepsize, maxsteps, density, color_func
        P = if applicable(f, Point3f(0))
            Point
        else
            Number
        end
        streamlines3d_impl(P, f, limits, resolution, stepsize, maxsteps, density, color_func)
    end
    colormap_args = MakieCore.colormap_attributes(p)
    generic_plot_attributes = MakieCore.generic_plot_attributes(p)
    println("streamlines working 8")
    lines!(
        p,
        lift(x->x[3], p, data),
        color = lift(last, p, data),
        linestyle = p.linestyle,
        linecap = p.linecap,
        joinstyle = p.joinstyle,
        miter_limit = p.miter_limit,
        linewidth = p.linewidth;
        colormap_args...,
        generic_plot_attributes...
    )
    println("streamlines working 9")
    rotations = map(x -> x[2], data)
    println("streamlines working 10")
    arrow_size = map(p, p.arrow_size) do arrow_size
        if arrow_size === automatic
            return 0.2 * minimum(p.limits[].widths) / minimum(p.gridsize[])
        else
            return arrow_size
        end
    end
    println("streamlines working 11")
    meshscatter!(
        p,
        lift(first, p, data);
        markersize=arrow_size, rotation=rotations,
        color=lift(x -> x[4], p, data),
        marker = lift(q -> Cone(q), p, p.quality),
        colormap_args...,
        generic_plot_attributes...
    )
end