# import Pkg; Pkg.add("Revise"); Pkg.add("GLMakie"); Pkg.add("Makie"); Pkg.add("MakieCore"); Pkg.add("LinearAlgebra"); Pkg.add("Interpolations"); Pkg.add("DelimitedFiles"); Pkg.add("ZipFile"); Pkg.add("Statistics")
using GLMakie

using LinearAlgebra
using Interpolations
using DelimitedFiles
using ZipFile
using Statistics
using Base.Threads
using NativeFileDialog

export run_demo;
export run_main;

include("streamlines3D.jl")

function read_dataset(filename::String, interval::Int, status_bar::Observable{Float32})
    #println("Loading in...")
    m = 512
    n = m ÷ interval
    field = zeros(Point3f, (n, n, n))
    xmin = Inf
    xmax = -Inf
    ymin = Inf
    ymax = -Inf
    zmin = Inf
    zmax = -Inf
    bhpoints = Vector{Point3f}()
    bhcenter = Point3f((0, 0, 0))
    bhr = 0.0f0
    if last(filename, 4) == ".zip"
        r = ZipFile.Reader(filename)
        for f in r.files
            # file_matrix = readdlm(f, ',', Float64)
            for k in 1:m
                for j in 1:m
                    for i in 1:m
                        line = split(readline(f), ",")
                        x, y, z, vx, vy, vz = parse.(Float64, line)
                        if any(isnan, (x, y, z, vx, vy, vz))
                            # println("Blackhole found at: ($x, $y, $z)")
                            xmin = min(isnan(x) ? Inf : x, xmin)
                            xmax = max(isnan(x) ? -Inf : x, xmax)
                            ymin = min(isnan(y) ? Inf : y, ymin)
                            ymax = max(isnan(y) ? -Inf : y, ymax)
                            zmin = min(isnan(z) ? Inf : z, zmin)
                            zmax = max(isnan(z) ? -Inf : z, zmax)
                            xi = ceil(Int, i / interval)
                            yi = ceil(Int, j / interval)
                            zi = ceil(Int, k / interval)
                            field[xi, yi, zi] = convert(Point3f, (vx, vy, vz))
                            push!(bhpoints, convert(Point3f, (x, y, z)))
                        elseif i % interval == 0 && j % interval == 0 && k % interval == 0
                            xmin = min(isnan(x) ? Inf : x, xmin)
                            xmax = max(isnan(x) ? -Inf : x, xmax)
                            ymin = min(isnan(y) ? Inf : y, ymin)
                            ymax = max(isnan(y) ? -Inf : y, ymax)
                            zmin = min(isnan(z) ? Inf : z, zmin)
                            zmax = max(isnan(z) ? -Inf : z, zmax)
                            xi, yi, zi = convert.(Int, div.((i, j, k), interval))
                            if field[xi, yi, zi][1] != NaN
                                field[xi, yi, zi] = convert(Point3f, (vx, vy, vz))
                            end
                        end
                        #println("status")
                        #println(Float32(Float32(i + j*m + k*m*m) / Float32(m + m*m + m*m*m)))
                        status_bar[] = Float32(Float32(i + j * m + k * m * m) / Float32(m + m * m + m * m * m))
                        #notify(status_bar)
                        #println("Loading...\t$(round(100*status_bar[]))%")
                    end
                end
            end
        end
        close(r)
    else
        # file_matrix = readdlm(filename, ',', Float64)
        open(filename) do f
            for k in 1:m
                for j in 1:m
                    for i in 1:m
                        line = split(readline(f), ",")
                        x, y, z, vx, vy, vz = parse.(Float64, line)
                        if any(isnan, (x, y, z, vx, vy, vz))
                            # println("Blackhole found at: ($x, $y, $z)")
                            xmin = min(isnan(x) ? Inf : x, xmin)
                            xmax = max(isnan(x) ? -Inf : x, xmax)
                            ymin = min(isnan(y) ? Inf : y, ymin)
                            ymax = max(isnan(y) ? -Inf : y, ymax)
                            zmin = min(isnan(z) ? Inf : z, zmin)
                            zmax = max(isnan(z) ? -Inf : z, zmax)
                            xi = ceil(Int, i / interval)
                            yi = ceil(Int, j / interval)
                            zi = ceil(Int, k / interval)
                            field[xi, yi, zi] = convert(Point3f, (vx, vy, vz))
                            push!(bhpoints, convert(Point3f, (x, y, z)))
                        elseif i % interval == 0 && j % interval == 0 && k % interval == 0
                            xmin = min(isnan(x) ? Inf : x, xmin)
                            xmax = max(isnan(x) ? -Inf : x, xmax)
                            ymin = min(isnan(y) ? Inf : y, ymin)
                            ymax = max(isnan(y) ? -Inf : y, ymax)
                            zmin = min(isnan(z) ? Inf : z, zmin)
                            zmax = max(isnan(z) ? -Inf : z, zmax)
                            xi, yi, zi = convert.(Int, div.((i, j, k), interval))
                            if field[xi, yi, zi][1] != NaN
                                field[xi, yi, zi] = convert(Point3f, (vx, vy, vz))
                            end
                        end
                        #println("status")
                        #println(Float32(Float32(i + j*m + k*m*m) / Float32(m + m*m + m*m*m)))
                        status_bar[] = Float32(Float32(i + j * m + k * m * m) / Float32(m + m * m + m * m * m))
                        #notify(status_bar)
                        #println("Loading...\t$(round(100*status_bar[]))%")
                    end
                end
            end
        end
    end
    # println(field[1, 1, 1])
    # println(field[n, n, n])

    xinterval = (xmax - xmin) / (n - 1)
    yinterval = (ymax - ymin) / (n - 1)
    zinterval = (zmax - zmin) / (n - 1)

    bound = (xmin:xinterval:xmax, ymin:yinterval:ymax, zmin:zinterval:zmax)

    if length(bhpoints) > 1
        bhcenter = mean(bhpoints)
        bhr = maximum(map(x -> norm(bhcenter - x), bhpoints))
    elseif length(bhpoints) == 1
        bhcenter = bhpoints[1]
        bhr = xinterval / 2
    end

    return bound, field, bhcenter, bhr
end;

function divergence(field, scale)
    l, m, n = size(field)
    volumedata = zeros(l, m, n)
    for k in 1:n
        for j in 1:m
            for i in 1:l
                dvx, dvy, dvz = 0, 0, 0
                if i == 1
                    #dvx = (field[i+1, j, k][1] - field[i, j, k][1]) / scale
                elseif i == l
                    #dvx = (field[i, j, k][1] - field[i-1, j, k][1]) / scale
                else
                    dvx = (field[i+1, j, k][1] - field[i-1, j, k][1]) / (2 * scale)
                end
                if j == 1
                    #dvy = (field[i, j+1, k][2] - field[i, j, k][2]) / scale
                elseif j == m
                    #dvy = (field[i, j, k][2] - field[i, j-1, k][2]) / scale
                else
                    dvy = (field[i, j+1, k][2] - field[i, j-1, k][2]) / (2 * scale)
                end
                if k == 1
                    #dvz = (field[i, j, k+1][3] - field[i, j, k][3]) / scale
                elseif k == n
                    #dvz = (field[i, j, k][3] - field[i, j, k-1][3]) / scale
                else
                    dvz = (field[i, j, k+1][3] - field[i, j, k-1][3]) / (2 * scale)
                end
                volumedata[i, j, k] = dvx + dvy + dvz
            end
        end
    end
    return volumedata
end;

"""
  link_cameras_lscene(f; step=0.01)

Solution to link cameras in subplots (different axes) as suggested by user `Thomberger` taken from: 
https://discourse.julialang.org/t/makie-linking-cameras-in-multiple-3d-plots/83309/6
Slightly edited by Adrian Filcík.
"""
function link_cameras_lscene(f; step=0.01)
    scenes = f.content[findall(x -> typeof(x) == LScene, f.content)]
    cameras = [x.scene.camera_controls for x in scenes]

    for (i, camera1) in enumerate(cameras)
        on(camera1.eyeposition) do x
            for j in collect(1:length(cameras))[1:end.!=i]
                if sum(abs.(x - cameras[j].eyeposition[])) > step
                    cameras[j].lookat[] = camera1.lookat[]
                    cameras[j].eyeposition[] = camera1.eyeposition[]
                    cameras[j].upvector[] = camera1.upvector[]
                    # cameras[j].zoom_mult[]      = cameras[i].zoom_mult[]
                    # cameras[j].fov              = cameras[i].fov
                    # cameras[j].near             = cameras[i].near
                    # cameras[j].far              = cameras[i].far
                    # cameras[j].bounding_sphere  = cameras[i].bounding_sphere

                    update_cam!(scenes[j].scene, cameras[j])
                end
            end
        end
    end
    return f
end;

function update_visualization_parameters(bound, field)
    println("begin update")
    itp = interpolate(bound, field, Gridded(Linear()))
    etp = extrapolate(itp, fill(NaN, 3))
    flow_field_fun(x, y, z) = convert(Point3f, etp(x, y, z))
    println("flow_field update done")

    offx = (bound[1][end] - bound[1][begin]) / 8
    boundx = bound[1][begin] - offx .. bound[1][end] + offx

    offy = (bound[2][end] - bound[2][begin]) / 8
    boundy = bound[2][begin] - offy .. bound[2][end] + offy

    offz = (bound[3][end] - bound[3][begin]) / 8
    boundz = bound[3][begin] - offz .. bound[3][end] + offz
    println("bounds updated")

    arrow_size = convert(Float64, 2 * bound[1].step)
    println("arrow_size updated")

    vnorm = norm.(field)
    minmaxnorm = extrema(vnorm .|> x -> (isnan(x) ? Inf : x))
    println("minmaxnorm updated")

    limits = (minmaxnorm[1], maximum(vnorm .|> x -> (!isfinite(x) ? minmaxnorm[1] : x)))
    println("limits updated")

    return flow_field_fun, boundx, boundy, boundz, arrow_size, vnorm, minmaxnorm, limits
end;

function visualize_field(ax, bound, field, bhcenter, bhr, limits, alphas=[0.5, 1.0, 0.5], volume_threshold=2, volume_gamma=2)
    inside = 0
    println("executed in ", inside += 1)
    flow_field(x, y, z) = @lift begin
        itp = interpolate($bound, $field, Gridded(Linear()))
        etp = extrapolate(itp, fill(NaN, 3))
        convert(Point3f, etp(x, y, z))
    end
    println(typeof(flow_field))

    println("executed in ", inside += 1)
    boundx = @lift begin
        offx = ($bound[1][end] - $bound[1][begin]) / 8
        $bound[1][begin] - offx .. $bound[1][end] + offx
    end
    boundy = @lift begin
        offy = ($bound[2][end] - $bound[2][begin]) / 8
        $bound[2][begin] - offy .. $bound[2][end] + offy
    end
    boundz = @lift begin
        offz = ($bound[3][end] - $bound[3][begin]) / 8
        $bound[3][begin] - offz .. $bound[3][end] + offz
    end

    println("executed in ", inside += 1)
    #println(flow_field)
    #println(boundx)
    #println(boundy)
    #println(boundz)
    println(flow_field, boundx, boundy, boundz, alphas[1])

    arrow_size = @lift(convert(Float64, 2 * $bound[1].step))

    println(arrow_size)
    streamlines3d!(ax, flow_field, boundx, boundy, boundz, alpha=alphas[1],
        color=norm, colormap=:viridis, gridsize=(16, 16, 16), arrow_size=arrow_size, linewidth=2,
        transparency=true)

    """
    if bhr > 0
        mesh!(ax, Sphere(bhcenter, bhr), alpha=alphas[2], color=:black)
    end
    """

    println("executed in ", inside += 1)
    mesh!(ax, @lift(Sphere($bhcenter, $bhr)), alpha=alphas[2], color=:black, visible=@lift($bhr > 0), transparency=true)

    println("executed in ", inside += 1)
    volumecolormap = @lift begin
        colormap = to_colormap(:coolwarm)
        for (i, v) in enumerate(abs.(-1:2/(length(colormap)-1):1))
            color = colormap[i]
            colormap[i] = RGBAf(color.r, color.g, color.b, color.alpha * v^$volume_gamma)
        end
        colormap
    end

    println("executed in ", inside += 1)
    vnorm = lift(x -> norm.(x), field)
    # minnorm, maxnorm = extrema(vnorm .|> x -> (isnan(x) ? Inf : x))
    println("executed in ", inside += 1)
    minmaxnorm = lift(y -> extrema(y .|> x -> (isnan(x) ? Inf : x)), vnorm)
    # fieldnorm = 2 .* (vnorm .- minnorm)
    # fieldnorm = lift(x -> clamp.(x .* vnorm, 0, 1), volume_factor)
    println("executed in ", inside += 1)
    fieldnorm = lift((x, y, z) -> Float64.(ifelse.(y .>= x, z[2], y) ./ z[2]), volume_threshold, vnorm, minmaxnorm)
    # fieldnorm = lift(x -> clamp.(vnorm ./ x, 0, 1), volume_threshold)
    # clamp!(fieldnorm, 0, 1)
    println("executed in ", inside += 1)
    field_sign = @lift(sign.(divergence($field, convert(Float64, $bound[1].step))))
    println("executed in ", inside += 1)
    volumedata = lift((x, y) -> ((y .* (1 .- x) .+ 1) ./ 2), fieldnorm, field_sign)
    println("executed in ", inside += 1)
    volume!(ax, @lift($bound[1]), @lift($bound[2]), @lift($bound[3]), volumedata, algorithm=:mip, alpha=alphas[3], colormap=volumecolormap, colorrange=(0, 1), transparency=true)

    println("executed in ", inside += 1)
    # limits = @lift begin 
    #     maxnorm = maximum($vnorm .|> x -> (!isfinite(x) ? $minmaxnorm[1] : x))
    #     ($minmaxnorm[1], maxnorm)
    # end
    println(typeof(minmaxnorm), typeof(vnorm))
    #GLMakie.Observables.@map!(limits, ((&minmaxnorm)[1], maximum(&vnorm .|> x -> (!isfinite(x) ? (&minmaxnorm)[1] : x))))
    #limits[] = (minmaxnorm[][1], maximum(vnorm[] .|> x -> (!isfinite(x) ? minmaxnorm[][1] : x)))

    println("executed in ", inside += 1)
    println(limits)

    return ax, vnorm, minmaxnorm
end;

function main()
    filename1 = "/media/radian-fi/KINGSTON/VIZ/Blackhole/el3_512_512_512.csv"
    filename2 = "/media/radian-fi/KINGSTON/VIZ/Blackhole/mag3_512_512_512.csv"
    subsample = 4

    # elbound, elfield, elbhcenter, elbhr = read_dataset(filename1, subsample)
    # magbound, magfield, magbhcenter, magbhr = read_dataset(filename2, subsample)

    # filenames = Observable([filename1, filename2])

    GLMakie.activate!()

    fig = Figure(size=(600, 800), fontsize=14)

    # status_bar = Vector{Observable{Float32}}(fill(Observable(0.0f0), 2))
    status_bar = [Observable(0.0f0) for _ in 1:2]
    # status_bar_string = [Observable("Loading... 0%") for _ in 1:2]
    el_status_string = Observable("Loading... 0%")
    mag_status_string = Observable("Loading... 0%")
    status_bar_visible = [Observable(false) for _ in 1:2]

    r, c = 1, 1
    Label(fig[r, c, Top()], "Datasets:", color=:dodgerblue, halign=:left, font=:bold, padding=(0, 0, 5, 5))
    folder = Observable(joinpath(pwd(), "data"))
    btn_dir = Button(fig[r, c], label=" Change folder... ", halign=:left, tellwidth=false)
    on(btn_dir.clicks) do n
        @async begin
            folder[] = fetch(@spawn pick_folder())
        end
    end

    el_options = Observable([("no valid datasets in given folder", "")])
    mag_options = Observable([("no valid datasets in given folder", "")])

    on(folder) do n
        files = filter!(contains(r"(mag|el)\d+(_\d+){3}\.(zip|csv)"), readdir(folder[]))
        data_paths = joinpath.(folder[], files)
        data_names = [file[1:end-4] for file in files]
        data_extensions = [file[end-2:end] for file in files]
        el_dict = Dict{String,String}()
        mag_dict = Dict{String,String}()

        for i in 1:length(files)
            if data_names[i][1:2] == "el" && (data_extensions[i] == "csv" || !haskey(el_dict, data_names[i]))
                el_dict[data_names[i]] = data_paths[i]
            elseif data_names[i][1:3] == "mag" && (data_extensions[i] == "csv" || !haskey(mag_dict, data_names[i]))
                mag_dict[data_names[i]] = data_paths[i]
            end
        end

        new_el_options = sort([(k, v) for (k, v) in el_dict], by=x -> parse(Int, match(r"\d+", x[1][3:end]).match))
        if isempty(el_dict)
            el_options[] = [("no valid datasets in given folder", "")]
        elseif el_options[] != new_el_options
            el_options[] = new_el_options
        end
        new_mag_options = sort([(k, v) for (k, v) in mag_dict], by=x -> parse(Int, match(r"\d+", x[1][4:end]).match))
        if isempty(mag_dict)
            mag_options[] = [("no valid datasets in given folder", "")]
        elseif mag_options[] != new_mag_options
            mag_options[] = new_mag_options
        end
    end

    btn_width = btn_dir.layoutobservables.computedbbox[].widths[1]

    Label(fig[r, c], "Folder: ", halign=:left, font=:bold, padding=(btn_width + 14, 0, 5, 5), tellwidth=false)
    Label(fig[r, c], folder, halign=:left, padding=(btn_width + 70, 0, 5, 5), tellwidth=false)

    r += 1
    el_label = Label(fig[r, c], "Electric field: ", halign=:left, font=:bold, padding=(0, 0, 5, 5), tellwidth=false)
    el_label_width = lift(x -> x.widths[1], el_label.layoutobservables.computedbbox)
    el_parent_width = lift(x -> x.widths[1], el_label.layoutobservables.suggestedbbox)
    el_menu_width = lift((x, y) -> (x - y - 7), el_parent_width, el_label_width)

    el_menu = Menu(fig[r, c],
        options=el_options,
        width=el_menu_width[],
        halign=:right,
        tellwidth=false)
    on(el_menu_width) do n
        if el_menu.width[] ≉ el_menu_width[]
            el_menu.width = el_menu_width[]
        end
    end

    r += 1
    el_status_bar = Label(fig[r, c], el_status_string, halign=:left, font=:bold, padding=(0, 0, 5, 5), tellwidth=false, visible=status_bar_visible[1])
    el_dot_update = 3
    # el_status = GLMakie.Observables.throttle(0.1, status_bar[1])
    el_status = status_bar[1]
    el_last_status = 0.0f0
    on(el_status) do el_percent
        if !isapprox(el_percent, el_last_status; atol=0.01)
            println("el. status update: $el_percent")
            el_dot_update = el_dot_update % 3
            el_dot_update += 1
            # status_bar_string[1][] = el_percent == 1.0f0 ? "Dataset loaded." : "Loading$('.'^(el_dot_update))$(' '^(3-el_dot_update))\t$(round(Int, el_percent*100))%"
            # notify(status_bar_string[1])
            # el_status_string[] = el_percent == 1.0f0 ? "Dataset loaded." : "Loading$('.'^(el_dot_update))$(' '^(3-el_dot_update))\t$(round(Int, el_percent*100))%"
            el_status_bar.text[] = el_percent == 1.0f0 ? "Dataset loaded." : "Loading$('.'^(el_dot_update))$(' '^(3-el_dot_update))\t$(round(Int, el_percent*100))%"
            el_last_status = el_percent
        elseif el_percent == 1.0f0
            # status_bar_string[1][] = "Dataset loaded."
            # notify(status_bar_string[1])
            # el_status_string[] = "Dataset loaded."
            el_status_bar.text[] = "Dataset loaded."
        end
        if !status_bar_visible[1][]
            status_bar_visible[1][] = true
        end
        #sleep(1)
        #yield()
    end
    """
    on(status_bar_string[1]) do s
        display(fig)
    end
    """

    r += 1
    mag_label = Label(fig[r, c], "Magnetic field: ", halign=:left, font=:bold, padding=(0, 0, 5, 5), tellwidth=false)
    mag_label_width = lift(x -> x.widths[1], mag_label.layoutobservables.computedbbox)
    mag_parent_width = lift(x -> x.widths[1], mag_label.layoutobservables.suggestedbbox)
    mag_menu_width = lift((x, y) -> (x - y - 7), mag_parent_width, mag_label_width)

    mag_menu = Menu(fig[r, c],
        options=mag_options,
        width=mag_menu_width[],
        halign=:right,
        tellwidth=false)
    on(mag_menu_width) do n
        if mag_menu.width[] ≉ mag_menu_width[]
            mag_menu.width = mag_menu_width[]
        end
    end

    r += 1
    mag_status_bar = Label(fig[r, c], mag_status_string, halign=:left, font=:bold, padding=(0, 0, 5, 5), tellwidth=false, visible=status_bar_visible[2])
    mag_dot_update = 3
    # mag_status = GLMakie.Observables.throttle(0.1, status_bar[2])
    mag_status = status_bar[2]
    mag_last_status = 0.0f0
    on(mag_status) do mag_percent
        if !isapprox(mag_percent, mag_last_status; atol=0.01)
            println("mag. status update: $mag_percent")
            mag_dot_update = mag_dot_update % 3
            mag_dot_update += 1
            # status_bar_string[2][] = mag_percent == 1.0f0 ? "Dataset loaded." : "Loading$('.'^(mag_dot_update))$(' '^(3-mag_dot_update))\t$(round(Int, mag_percent*100))%"
            # notify(status_bar_string[2])
            # mag_status_string[] = mag_percent == 1.0f0 ? "Dataset loaded." : "Loading$('.'^(mag_dot_update))$(' '^(3-mag_dot_update))\t$(round(Int, mag_percent*100))%"
            mag_status_bar.text[] = mag_percent == 1.0f0 ? "Dataset loaded." : "Loading$('.'^(mag_dot_update))$(' '^(3-mag_dot_update))\t$(round(Int, mag_percent*100))%"
            mag_last_status = mag_percent
        elseif mag_percent == 1.0f0
            # status_bar_string[2][] = "Dataset loaded."
            # notify(status_bar_string[2])
            # mag_status_string[] = "Dataset loaded."
            mag_status_bar.text[] = "Dataset loaded."
        end
        if !status_bar_visible[2][]
            status_bar_visible[2][] = true
        end
        #sleep(1)
        #yield()
    end
    """
    on(status_bar_string[2]) do s
        display(GLMakie.Screen(), fig)
    end
    """

    r += 1
    Label(fig[r, c, Top()], "Streamlines visualization parameters:", color=:dodgerblue, halign=:left, font=:bold, padding=(0, 0, 5, 5))
    sg_streamlines = SliderGrid(
        fig[r, c],
        (label="Streamlines opacity", range=0:0.01:1, format="{:.2f}", startvalue=0.5),
        (label="Streamlines threshold", range=0:0.01:1, format="{:.2f}", startvalue=1.0),
        halign=:left
    )
    r += 1
    Label(fig[r, c, Top()], "Zero-point visualization parameters:", color=:dodgerblue, halign=:left, font=:bold, padding=(0, 0, 5, 5))
    sg_zeropoint = SliderGrid(
        fig[r, c],
        (label="Zero-point opacity", range=0:0.01:1, format="{:.2f}", startvalue=0.5),
        (label="Zero-point threshold", range=0:0.1:20, format="{:.1f}", startvalue=0.1),
        (label="Zero-point border smoothing", range=0:0.01:1, format="{:.2f}", startvalue=0.33),
        halign=:left
    )
    r += 1
    Label(fig[r, c, Top()], "Balck hole visualization parameters:", color=:dodgerblue, halign=:left, font=:bold, padding=(0, 0, 5, 5))
    sg_bh = SliderGrid(
        fig[r, c],
        (label="Black hole opacity", range=0:0.01:1, format="{:.2f}", startvalue=1.0),
        halign=:left
    )

    alphas = [sg_streamlines.sliders[1].value, sg_bh.sliders[1].value, sg_zeropoint.sliders[1].value]

    volume_threshold = sg_zeropoint.sliders[2].value

    volume_gamma = lift(x -> (1 - x) / x, sg_zeropoint.sliders[3].value)
    volumecolormap = @lift begin
        colormap = to_colormap(:coolwarm)
        for (i, v) in enumerate(abs.(-1:2/(length(colormap)-1):1))
            color = colormap[i]
            colormap[i] = RGBAf(color.r, color.g, color.b, color.alpha * v^$volume_gamma)
        end
        colormap
    end

    fig_labels = ["Electric field", "Magnetic field"]

    vis_fig = Figure(size=(1200, 800), fontsize=26)
    axs = [LScene(vis_fig[1, i]) for i = 1:2]

    vis_fig[0, :] = Label(vis_fig, "Black Hole Visualization")

    # bound = Vector{Observable{Tuple{StepRangeLen{Number,Number,Number,Int},StepRangeLen{Number,Number,Number,Int},StepRangeLen{Number,Number,Number,Int}}}}(undef, 2)
    # bound = Vector{Observable{Tuple{StepRangeLen{Number,Number,Number,Int},StepRangeLen{Number,Number,Number,Int},StepRangeLen{Number,Number,Number,Int}}}}(fill(Observable((-4.:4., -4.:4., -4.:4.)), 2))
    bound = [Observable((-4.0:4.0, -4.0:4.0, -4.0:4.0)) for _ in 1:2]
    # field = Vector{Observable{Array{Point{3,Float32},3}}}(undef, 2)
    # field = Vector{Observable{Array{Point{3,Float32},3}}}(fill(Observable(reshape([Point3f(0) for _ in 1:9*9*9], (9, 9, 9))), 2))
    field = [Observable(reshape([Point3f(0) for _ in 1:9*9*9], (9, 9, 9))) for _ in 1:2]
    # bhcenter = Vector{Observable{Point3f}}(undef, 2)
    # bhcenter = Vector{Observable{Point3f}}(fill(Observable(Point3f(0)), 2))
    bhcenter = [Observable(Point3f(0)) for _ in 1:2]
    # bhr = Vector{Observable{Float32}}(undef, 2)
    # bhr = Vector{Observable{Float32}}(fill(Observable(0.0f0), 2))
    bhr = [Observable(0.0f0) for _ in 1:2]
    # limits = Vector{Observable{Tuple{Float32,Float32}}}(undef, 2)
    # limits = Vector{Observable{Tuple{Float32,Float32}}}(fill(Observable((0.0f0, 1.0f0)), 2))
    limits = [Observable((0.0f0, 1.0f0)) for _ in 1:2]
    # colorbar_visible = Vector{Observable{Bool}}(fill(Observable(false), 2))
    fun_itp = interpolate(bound[1][], field[1][], Gridded(Linear()))
    fun_etp = extrapolate(fun_itp, fill(NaN, 3))
    fun(x, y, z) = convert(Point3f, fun_etp(x, y, z))
    flow_field = [Observable{Function}((x, y, z) -> fun(x, y, z)) for _ in 1:2]
    boundx = [Observable(-5.0 .. 5.0) for _ in 1:2]
    boundy = [Observable(-5.0 .. 5.0) for _ in 1:2]
    boundz = [Observable(-5.0 .. 5.0) for _ in 1:2]
    arrow_size = [Observable(2.0) for _ in 1:2]
    colorbar_visible = [Observable(false) for _ in 1:2]
    vnorm = [Observable(norm.(field[i][])) for i in 1:2]
    minmaxnorm = [Observable(extrema(vnorm[i][] .|> x -> (isnan(x) ? Inf : x))) for i in 1:2]
    field_sign = [Observable(sign.(divergence(field[i][], convert(Float64, bound[i][][1].step)))) for i in 1:2]
    volumedata = [Observable(norm.(field[i][]) .+ (1 / 2)) for i in 1:2]

    """
    @threads for i in 1:2
        @inbounds bound[i], field[i], bhcenter[i], bhr[i] = read_dataset(filenames[i], subsample)
    end
    """

    """
    on(bhr[1]) do n
        println("vis_figure update bhr[1]")
        sleep(1)  # hopefully this is enough time
        display(vis_fig)
    end

    on(bhr[2]) do n
        println("vis_figure update bhr[2]")
        sleep(1)  # hopefully this is enough time
        display(vis_fig)
    end

    onany(bound[1], bound[2], field[1], field[2], bhcenter[1], bhcenter[2], bhr[1], bhr[2]) do el_bound, mag_bound, el_field, mag_field, el_bhcenter, mag_bhcenter, el_bhr, mag_bhrb
        println("vis_figure update")
        sleep(1)  # hopefully this is enough time
        display(vis_fig)
    end
    """

    onany(vnorm[1], minmaxnorm[1]) do el_vnorm, el_minmaxnorm
        nothing
    end

    #status_bar[1][] = Float32(Float32(1 + 1*512 + 1*512*512) / Float32(512 + 512*512 + 512*512*512))

    el_valid = Observable(false)
    # el_flow_field = nothing
    on(el_menu.selection) do filename
        if !isnothing(filename) && filename != ""
            @spawn begin
                el_bound, el_field, el_bhcenter, el_bhr = fetch(@spawn read_dataset(filename, subsample, status_bar[1]))
                #println(typeof(bound[1][]) == typeof(el_bound))
                #println(typeof(field[1][]) == typeof(el_field))
                #println(typeof(bhcenter[1][]) == typeof(el_bhcenter))
                #println(typeof(bhr[1][]) == typeof(el_bhr))
                bound[1].val = el_bound
                field[1].val = el_field
                bhcenter[1].val = el_bhcenter
                bhr[1][] = el_bhr
                println("Loading el. field done: ", filename)
                #notify(bound[1])
                #notify(field[1])
                #notify(bhcenter[1])
                #notify(bhr[1])
                el_flow_field, el_boundx, el_boundy, el_boundz, el_arrow_size, el_vnorm, el_minmaxnorm, el_limits = update_visualization_parameters(el_bound, el_field)
                #println(typeof(el_flow_field), " ", typeof(el_boundx), " ", typeof(el_boundy), " ", typeof(el_boundz), " ", typeof(el_arrow_size))
                #println(typeof(el_flow_field) == typeof(flow_field[1][]))
                println(typeof(el_flow_field))
                println(typeof((x, y, z) -> el_flow_field(x, y, z)))
                println((x, y, z) -> el_flow_field(x, y, z))
                println(typeof(flow_field[1][]))
                println(typeof(flow_field[1]))
                println(flow_field[1][])
                #println(typeof(el_boundx) == typeof(boundx[1][]))
                #println(typeof(el_arrow_size) == typeof(arrow_size[1][]))
                #flow_field[1].val = (x, y, z) -> el_flow_field(x, y, z)
                boundx[1].val = el_boundx
                boundy[1].val = el_boundy
                boundz[1].val = el_boundz
                println("el. streamlines parameters loaded")
                #arrow_size[1][] = el_arrow_size
                arrow_size[1].val = el_arrow_size
                limits[1][] = el_limits
                println("el. vis. parameters updated")
                vnorm[1].val = el_vnorm
                minmaxnorm[1].val = el_minmaxnorm
                field_sign[1][] = sign.(divergence(el_field, convert(Float64, el_bound[1].step)))
                println("El. field renderer notified.")
                el_valid[] = true
            end
        end
    end

    # mag_valid = Observable(false)
    on(mag_menu.selection) do filename
        if !isnothing(filename) && filename != ""
            @spawn begin
                mag_bound, mag_field, mag_bhcenter, mag_bhr = fetch(@spawn read_dataset(filename, subsample, status_bar[2]))
                bound[2].val = mag_bound
                field[2].val = mag_field
                bhcenter[2].val = mag_bhcenter
                bhr[2].val = mag_bhr
                println("Loading mag. field done: ", filename)
                notify(bound[2])
                notify(field[2])
                notify(bhcenter[2])
                notify(bhr[2])
                println("Mag. field renderer notified.")
                # mag_valid[] = true
            end
        end
    end

    # while !(el_valid[] && mag_valid[]) end

    """
    on(limits[1]) do limit
        if limit[1] == limit[2]
            limits[1][] = (0.0f0, 1.0f0)
            colorbar_visible[1][] = false
        end
    end

    on(limits[2]) do limit
        if limit[1] == limit[2]
            limits[2][] = (0.0f0, 1.0f0)
            colorbar_visible[2][] = false
        end
    end
    """

    on(limits[1]) do n
        println("el. limits updated")
    end

    on(limits[2]) do n
        println("mag. limits updated")
    end

    for i in 1:2
        on(limits[i]) do limit
            println("notified")
            if limit[1] == limit[2]
                limits[i][] = (0.0f0, 1.0f0)
                colorbar_visible[i][] = false
            else
                colorbar_visible[i][] = true
            end
        end

        println("executed 2")
        Label(vis_fig[1, i, Top()], fig_labels[i], padding=(0, 0, 5, 0))
        println("executed 3")
        # ax, vnorm[i], minmaxnorm[i] = visualize_field(axs[i], bound[i], field[i], bhcenter[i], bhr[i], limits[i], alphas, volume_threshold, volume_gamma)
        streamlines3d!(axs[i], flow_field[i], boundx[i], boundy[i], boundz[i], alpha=alphas[1], color=norm, colormap=:viridis, gridsize=(16, 16, 16), arrow_size=arrow_size[i], linewidth=2, transparency=true)
        mesh!(axs[i], @lift(Sphere($(bhcenter[i]), $(bhr[i]))), alpha=alphas[2], color=:black, visible=@lift($(bhr[i]) > 0), transparency=true)
        volume!(axs[i], lift(x -> x[1], bound[i]), lift(x -> x[2], bound[i]), lift(x -> x[3], bound[i]), volumedata[i], algorithm=:mip, alpha=alphas[3], colormap=volumecolormap, colorrange=(0, 1), transparency=true)
        #GLMakie.Observables.@map!(limits[i], ((&minmaxnorm[i])[1], maximum(&vnorm[i] .|> x -> (!isfinite(x) ? (&minmaxnorm[i])[1] : x))))
        #notify(limits[i])
        #println(limits[i][])
        println("executed 4")

        Colorbar(vis_fig[1, i][2, 1], limits=limits[i], colormap=:viridis, vertical=false, flipaxis=false,
            bottomspinevisible=colorbar_visible[i], labelvisible=colorbar_visible[i], leftspinevisible=colorbar_visible[i],
            rightspinevisible=colorbar_visible[i], ticklabelsvisible=colorbar_visible[i], ticksvisible=colorbar_visible[i],
            topspinevisible=colorbar_visible[i], label="Vector magnitude")
    end

    println("executed 5")
    vis_fig = link_cameras_lscene(vis_fig)
    println("executed 6")
    display(GLMakie.Screen(), vis_fig)
    println("executed 7")
    display(GLMakie.Screen(), fig)

    onany(vnorm[1], minmaxnorm[1], field_sign[1], volume_threshold) do el_vnorm, el_minmaxnorm, el_field_sign, el_volume_threshold
        println("begin volume data update")
        el_fieldnorm = Float64.(ifelse.(el_vnorm .>= el_volume_threshold, el_minmaxnorm[2], el_vnorm) ./ el_minmaxnorm[2])
        volumedata[1][] = (el_field_sign .* (1 .- el_fieldnorm) .+ 1) ./ 2
        println("el. volumedata updated")
        # display(vis_fig)
    end

    on(bhr[1]) do el_bhr
        println("el. bhr updated: ", el_bhr)
    end

    on(el_valid) do n
        println("rendering new vis_figure")
        sleep(1.0)
        # display(vis_fig)
        println("new vis_figure rendered")
    end

    println("end of main")
end;

main();