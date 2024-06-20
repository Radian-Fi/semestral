# import Pkg; Pkg.add("Revise"); Pkg.add("GLMakie"); Pkg.add("LinearAlgebra"); Pkg.add("Interpolations"); Pkg.add("DelimitedFiles"); Pkg.add("ZipFile"); 

using GLMakie

using LinearAlgebra
using Interpolations
using DelimitedFiles
using ZipFile
using Statistics

export run_demo;
export run_main;

include("streamlines3D.jl")

function run_demo(n = 64)
  GLMakie.activate!()

  n -= 1;

  ps = [Point3f(x, y, z) for x in -n:2:n for y in -n:2:n for z in -n:2:n]
  ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)
  lengths = norm.(ns)
  arrows(
      ps, ns, fxaa=true, # turn on anti-aliasing
      color=lengths,
      linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.4),
      align = :center, axis=(type=Axis3,)
  )
end;

"""
function flow_field_interpolation(n, points, x, y, z)
  x = clamp(x, 1, n)
  y = clamp(y, 1, n)
  z = clamp(z, 1, n)

  x_idx_l = clamp(floor(Int64, x), 1, n)
  x_idx_u = clamp(ceil(Int64, x), 1, n)
  y_idx_l = clamp(floor(Int64, y), 1, n)
  y_idx_u = clamp(ceil(Int64, y), 1, n)
  z_idx_l = clamp(floor(Int64, z), 1, n)
  z_idx_u = clamp(ceil(Int64, z), 1, n)

  x_diff = x - x_idx_l
  y_diff = y - y_idx_l
  z_diff = z - z_idx_l

   = points[x_idx_l]
end;
"""

function read_dataset(filename::String, interval::Int)
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
  bhr = 0
  if last(filename, 4) == ".zip"
    r = ZipFile.Reader(filename)
    for f in r.files
      # file_matrix = readdlm(f, ',', Float64)
      for k in 1:m
        for j in 1:m
          for i in 1:m
            line = split(readline(f), ",")
            x, y, z, vx, vy, vz = parse.(Float64, line)
            if vx == NaN
              xmin = min(isnan(x) ? Inf : x, xmin)
              xmax = max(isnan(x) ? -Inf : x, xmax)
              xi = round(Int, x / interval) + 1
              yi = round(Int, y / interval) + 1
              zi = round(Int, z / interval) + 1
              field[xi, yi, zi] = convert(Point3f, (vx, vy, vz))
            elseif x % interval == 0 && y % interval == 0 && z % interval == 0
              xmin = min(isnan(x) ? Inf : x, xmin)
              xmax = max(isnan(x) ? -Inf : x, xmax)
              xi, yi, zi = convert.(Int, div.((x, y, z), interval)) .+ 1
              if field[xi, yi, zi][1] != NaN
                field[xi, yi, zi] = convert(Point3f, (vx, vy, vz))
              end
            end
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
          end
        end
      end
    end
  end
  println(field[1, 1, 1])
  println(field[n, n, n])

  xinterval = (xmax - xmin) / (n - 1)
  yinterval = (ymax - ymin) / (n - 1)
  zinterval = (zmax - zmin) / (n - 1)

  bound = (xmin:xinterval:xmax, ymin:yinterval:ymax, zmin:zinterval:zmax)

  if length(bhpoints) > 1
    bhcenter = mean(bhpoints)
    bhr = maximum(map(x->norm(bhcenter-x), bhpoints))
  elseif length(bhpoints) == 1
    bhcenter = bhpoints[1]
    bhr = xinterval / 2
  end

  return bound, field, bhcenter, bhr
end;

function streamplot_in_3d()
  GLMakie.activate!()

  ps = [Point3f(x, y, z) for x=-3:1:3 for y=-3:1:3 for z=-3:1:3]
  ns = map(p -> 0.1 * rand() * Vec3f(p[2], p[3], p[1]), ps)
  ns = reshape(ns, 7, 7, 7)
  # lengths = norm.(ns)
  # flowField(x, y, z) = Point3f(-y + x * (-1 + x^2 + y^2)^2, x + y * (-1 + x^2 + y^2)^2, z + x * (y - z^2))
  itp = interpolate((-3:1:3, -3:1:3, -3:1:3), ns, Gridded(Linear()))
  etp = extrapolate(itp, fill(NaN, 3))
  flowField(x, y, z) = convert(Point3f, etp(x, y, z))

  fig = Figure(size=(1200, 800), fontsize=26)
  axs = [LScene(fig[1, i]) for i=1:2]

  fig[0, :] = Label(fig, "Black Hole Visualization")

  Label(fig[1, 1, Top()], "Electric field", padding = (0, 0, 5, 0))
  streamlines3d!(axs[1], flowField, -4 .. 4, -4 .. 4, -4 .. 4, alpha=0.5,
      color=norm, colormap=:viridis, gridsize=(16, 16, 16), arrow_size=0.1, linewidth=2)
  
  Label(fig[1, 2, Top()], "Magnetic field", padding = (0, 0, 5, 0))
  streamlines3d!(axs[2], flowField, -4 .. 4, -4 .. 4, -4 .. 4, alpha=0.5,
      color=norm, colormap=:viridis, gridsize=(16, 16, 16), arrow_size=0.1, linewidth=2)
  fig = link_cameras_lscene(fig)
  fig
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

function visualize_main(elbound, elfield, elbhcenter, elbhr, magbound, magfield, magbhcenter, magbhr)
  GLMakie.activate!()

  # ps = [Point3f(x, y, z) for x=-3:1:3 for y=-3:1:3 for z=-3:1:3]
  # ns = map(p -> 0.1 * rand() * Vec3f(p[2], p[3], p[1]), ps)
  # ns = reshape(ns, 7, 7, 7)
  # lengths = norm.(ns)
  # flowField(x, y, z) = Point3f(-y + x * (-1 + x^2 + y^2)^2, x + y * (-1 + x^2 + y^2)^2, z + x * (y - z^2))
  # itp = interpolate((-3:1:3, -3:1:3, -3:1:3), ns, Gridded(Linear()))
  elitp = interpolate(elbound, elfield, Gridded(Linear()))
  eletp = extrapolate(elitp, fill(NaN, 3))
  el_flow_field(x, y, z) = convert(Point3f, eletp(x, y, z))

  magitp = interpolate(magbound, magfield, Gridded(Linear()))
  magetp = extrapolate(magitp, fill(NaN, 3))
  mag_flow_field(x, y, z) = convert(Point3f, magetp(x, y, z))

  fig = Figure(size=(1200, 800), fontsize=26)
  axs = [LScene(fig[1, i]) for i=1:2]

  fig[0, :] = Label(fig, "Black Hole Visualization")

  Label(fig[1, 1, Top()], "Electric field", padding = (0, 0, 5, 0))

  eloffx = (elbound[1][end] - elbound[1][begin]) / 8
  elboundx = elbound[1][begin]-eloffx .. elbound[1][end]+eloffx
  eloffy = (elbound[2][end] - elbound[2][begin]) / 8
  elboundy = elbound[2][begin]-eloffy .. elbound[2][end]+eloffy
  eloffz = (elbound[3][end] - elbound[3][begin]) / 8
  elboundz = elbound[3][begin]-eloffz .. elbound[3][end]+eloffz

  println((elboundx, elboundy, elboundz))

  streamplot!(axs[1], el_flow_field, elboundx, elboundy, elboundz, alpha=0.5,
  color=norm, colormap=:viridis, gridsize=(16, 16, 16),  arrow_size=convert(Float64, 2*elbound[1].step), linewidth=2,
  transparency=true)

  if elbhr > 0
    mesh!(axs[1], Sphere(elbhcenter, elbhr), alpha=1.0, color=:black)
  end

  volumecolormap = to_colormap(:coolwarm)
  for (i, v) in enumerate(-1:2/(length(volumecolormap)-1):1)
    color = volumecolormap[i]
    volumecolormap[i] = RGBAf(color.r, color.g, color.b, color.alpha * v^2)
  end

  elvnorm = norm.(elfield)
  elminnorm, elmaxnorm = extrema(elvnorm .|> x->(isnan(x) ? Inf : x))
  elnorm = 2 .* (elvnorm .- elminnorm)
  clamp!(elnorm, 0, 1)
  elvolumedata = (sign.(divergence(elfield, convert(Float64, elbound[1].step))) .* (1 .- elnorm) .+1) ./ 2
  volume!(axs[1], elbound[1], elbound[2], elbound[3], elvolumedata, algorithm=:mip, alpha=0.5, colormap=volumecolormap, colorrange=(0, 1), transparency=true)

  elmaxnorm = maximum(elvnorm .|> x->(!isfinite(x) ? elminnorm : x))
  println((elminnorm, elmaxnorm))
  Colorbar(fig[1, 1][2, 1], limits = (elminnorm, elmaxnorm), colormap = :viridis, vertical = false, flipaxis=false, label="Vector magnitude")
  

  Label(fig[1, 2, Top()], "Magnetic field", padding = (0, 0, 5, 0))

  magoffx = (magbound[1][end] - magbound[1][begin]) / 8
  magboundx = magbound[1][begin]-magoffx .. magbound[1][end]+magoffx
  magoffy = (magbound[2][end] - magbound[2][begin]) / 8
  magboundy = magbound[2][begin]-magoffy .. magbound[2][end]+magoffy
  magoffz = (magbound[3][end] - magbound[3][begin]) / 8
  magboundz = magbound[3][begin]-magoffz .. magbound[3][end]+magoffz

  streamplot!(axs[2], mag_flow_field, magboundx, magboundy, magboundz, alpha=0.5,
  color=norm, colormap=:viridis, gridsize=(16, 16, 16), arrow_size=convert(Float64, 2*magbound[1].step), linewidth=2,
  transparency=true)

  if magbhr > 0
    mesh!(axs[2], Sphere(magbhcenter, magbhr), alpha=1.0, color=:black)
  end

  normcolormap = to_colormap(:binary)
  for (i, v) in enumerate(0:1/(length(normcolormap)-1):1)
    color = normcolormap[i]
    normcolormap[i] = RGBAf(color.r, color.g, color.b, color.alpha * v)
  end

  magvnorm = norm.(magfield)
  magminnorm, magmaxnorm = extrema(magvnorm .|> x->(isnan(x) ? Inf : x))
  # magnorm = (magvnorm .- magminnorm) ./ (magmaxnorm - magminnorm)
  magnorm = 2 .* (magvnorm .- magminnorm)
  clamp!(magnorm, 0, 1)
  #magvolumedata = divergence(magfield, convert(Float64, magbound[1].step)) .* (1 .- magnorm)
  #magmindiv, magmaxdiv = extrema(magvolumedata)
  #magvolumedata .= clamp.(Float32.((magvolumedata) ./ (2 * min(abs(magmaxdiv), abs(magmindiv))) .+ 0.5), 0, 1)
  magvolumedata = (sign.(divergence(magfield, convert(Float64, magbound[1].step))) .* (1 .- magnorm) .+1) ./ 2
  volume!(axs[2], magbound[1], magbound[2], magbound[3], magvolumedata, algorithm=:mip, alpha=0.5, colormap=volumecolormap, colorrange=(0, 1), transparency=true)
  #svolume!(axs[2], magbound[1], magbound[2], magbound[3], (1 .- magnorm), algorithm=:mip, alpha=1.0, colormap=normcolormap)

  magmaxnorm = maximum(magvnorm .|> x->(!isfinite(x) ? magminnorm : x))
  println((magminnorm, magmaxnorm))
  Colorbar(fig[1, 2][2, 1], limits = (magminnorm, magmaxnorm), colormap = :viridis, vertical = false, flipaxis=false, label="Vector magnitude")

  magmindiv, magmaxdiv = extrema(magvolumedata)
  println((magmindiv, magmaxdiv))
  #println(magvolumedata[ceil(Int64, (0.723-magbound[1][begin])/convert(Float64, magbound[1].step)), ceil(Int64, (-1.117-magbound[2][begin])/convert(Float64, magbound[2].step)), ceil(Int64, (-0.97-magbound[3][begin])/convert(Float64, magbound[3].step))])
  #println(magvnorm[ceil(Int64, (0.723-magbound[1][begin])/convert(Float64, magbound[1].step)), ceil(Int64, (-1.117-magbound[2][begin])/convert(Float64, magbound[2].step)), ceil(Int64, (-0.97-magbound[3][begin])/convert(Float64, magbound[3].step))])
  fig = link_cameras_lscene(fig)
  fig
end;

"""
  link_cameras_lscene(f; step=0.01)

Solution to link cameras in subplots (different axes) as suggested by user `Thomberger` taken from: 
https://discourse.julialang.org/t/makie-linking-cameras-in-multiple-3d-plots/83309/6
Slightly edited by Adrian Filcík.
"""
function link_cameras_lscene(f; step=0.01)
  scenes = f.content[findall(x -> typeof(x) == LScene,f.content)]
  cameras = [x.scene.camera_controls for x in scenes]

  for (i, camera1) in enumerate(cameras)
      on(camera1.eyeposition) do x
          for j in collect(1:length(cameras))[1:end .!= i]
              if sum(abs.(x - cameras[j].eyeposition[])) > step
                  cameras[j].lookat[]         = camera1.lookat[]
                  cameras[j].eyeposition[]    = camera1.eyeposition[]
                  cameras[j].upvector[]       = camera1.upvector[]
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

function run_main(filename1="/media/radian-fi/KINGSTON/VIZ/Blackhole/el3_512_512_512.csv", filename2="/media/radian-fi/KINGSTON/VIZ/Blackhole/mag3_512_512_512.csv", subsample=4)
  elbound, elfield, elbhcenter, elbhr = read_dataset(filename1, subsample)
  magbound, magfield, magbhcenter, magbhr = read_dataset(filename2, subsample)
  println(elbound)
  # streamplot_in_3d()
  visualize_main(elbound, elfield, elbhcenter, elbhr, magbound, magfield, magbhcenter, magbhr)
end

#elbound, elfield, elbhcenter, elbhr = read_dataset("/media/radian-fi/KINGSTON/VIZ/Blackhole/el3_512_512_512.csv", 4)
#magbound, magfield, magbhcenter, magbhr = read_dataset("/media/radian-fi/KINGSTON/VIZ/Blackhole/mag3_512_512_512.csv", 4)
#println(elbound)
# arrows_and_streamplot_in_3d()
#visualize_main(elbound, elfield, elbhcenter, elbhr, magbound, magfield, magbhcenter, magbhr)

streamplot_in_3d()
