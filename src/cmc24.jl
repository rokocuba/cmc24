# Computational Modeling Challenge 2024
# CMC24 solution evaluation script
# Author: Hrvoje Abraham, hrvoje.abraham@avl.com

using FileIO
using ImageIO
using Measures
using Plots;
gr();
using UUIDs
using JLD2

temple_string =
#1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0
"O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O
 O  .  .  .  .  O  .  .  .  .  .  .  .  .  O  .  .  .  .  O
 O  .  .  .  .  .  .  .  O  .  .  O  .  .  .  .  .  .  .  O
 O  .  .  .  .  .  O  .  .  .  .  .  .  O  .  .  .  .  .  O
 O  .  .  O  .  .  .  .  .  O  O  .  .  .  .  .  O  .  .  O
 O  O  .  .  .  .  .  .  .  O  O  .  .  .  .  .  .  .  O  O
 O  .  .  .  O  .  .  .  .  .  .  .  .  .  .  O  .  .  .  O
 O  .  .  .  .  .  .  O  .  .  .  .  O  .  .  .  .  .  .  O
 O  .  .  .  .  .  .  .  .  O  O  .  .  .  .  .  .  .  .  O
 O  .  O  .  .  O  O  .  .  .  .  .  .  O  O  .  .  O  .  O
 O  .  O  .  .  O  O  .  .  .  .  .  .  O  O  .  .  O  .  O
 O  .  .  .  .  .  .  .  .  O  O  .  .  .  .  .  .  .  .  O
 O  .  .  .  .  .  .  O  .  .  .  .  O  .  .  .  .  .  .  O
 O  .  .  .  O  .  .  .  .  .  .  .  .  .  .  O  .  .  .  O
 O  O  .  .  .  .  .  .  .  O  O  .  .  .  .  .  .  .  O  O
 O  .  .  O  .  .  .  .  .  O  O  .  .  .  .  .  O  .  .  O
 O  .  .  .  .  .  O  .  .  .  .  .  .  O  .  .  .  .  .  O
 O  .  .  .  .  .  .  .  O  .  .  O  .  .  .  .  .  .  .  O
 O  .  .  .  .  O  .  .  .  .  .  .  .  .  O  .  .  .  .  O
 O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O  O"

block_size = 1
mirror_length = 0.5
light_halfwidth = 1
ε = 1e-12

"""Float64 infinity"""
∞ = Inf

"""
	⋅(v, w)

Dot product of two 2D vectors.
"""
function ⋅(v, w)
   return v[1] * w[1] + v[2] * w[2]
end

"""
	×(v, w)

Cross product of two 2D vectors.
"""
function ×(v, w)
   return v[1] * w[2] - v[2] * w[1]
end

"""Last 12 digits of hexadecimal format of the input integer."""
function hex12(x::Integer)
   return last(string(x, base=16), 12)
end

"""
	Convert integer into string with digit grouping.

	E.g. 1234567 => "1,234,567"
"""
function commas(num::Integer)
   str = string(num)
   return replace(str, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

function load_temple(temple_string, block_size)
   # println(stderr, " " * temple_string)
   rows = split(replace(strip(temple_string), " " => ""), '\n')
   temple_shape = length(rows[1]), length(rows)

   temple = Set()
   for (j, row) ∈ enumerate(rows)
      for (i, c) ∈ enumerate(row)
         if c == 'O'
            x = (i - 1) * block_size
            y = temple_shape[2] - j * block_size

            v1 = [x, y]
            v2 = [x + block_size, y]
            v3 = [x + block_size, y + block_size]
            v4 = [x, y + block_size]

            block = (
               v1=v1,  # bottom left corner
               v2=v2,
               v3=v3,  # up right corner
               v4=v4, s1=(v1, block_size, 0),
               s2=(v2, block_size, π / 2),
               s3=(v3, block_size, π),
               s4=(v4, block_size, 3π / 2),
            )

            push!(temple, block)
         end
      end
   end
   # display(temple)

   println(stderr, "The temple of size $temple_shape is loaded.")

   return (
      blocks=temple,
      shape=temple_shape,
      size=block_size .* temple_shape,
   )
end

function load_solution(cmc24_solution, mirror_length)
   if size(cmc24_solution) ≠ (9, 3)
      println(stderr, "ERROR! The solution isn't 9x3 size matrix.")
      finalize()
   end

   if !(eltype(cmc24_solution) <: Number)
      println(stderr, "ERROR! The solution contains non-numerical inputs.")
      finalize()
   end

   try
      cmc24_solution = float(cmc24_solution)
   catch
      println(stderr, "ERROR! The solution can't be converted to double precision floating point format.")
      finalize()
   end

   # preprocess the lamp
   α = cmc24_solution[1, 3]
   lamp = (
      v=cmc24_solution[1, 1:2],
      α=α,
      e=[cos(α), sin(α)],
   )

   # preprocess the mirrors
   mirrors = []
   for m ∈ 1:8
      α = cmc24_solution[m+1, 3]

      v = cmc24_solution[m+1, 1:2]
      e = [cos(α), sin(α)]
      n = [-sin(α), cos(α)]  # normal

      mirror = (
         v1=v,
         v2=v + mirror_length * e,
         s=(v, mirror_length, α),
         α=α,
         e=e,
         n=n,
      )

      push!(mirrors, mirror)
   end

   println(stderr, "The solution is loaded.")

   return (lamp, mirrors)
end

function load_solution_partial(cmc24_solution, mirror_length)
   cmc24_solution = float(cmc24_solution)

   # preprocess the lamp
   α = cmc24_solution[1, 3]
   lamp = (
      v=cmc24_solution[1, 1:2],
      α=α,
      e=[cos(α), sin(α)],
   )

   # preprocess the mirrors
   mirrors = []
   for m ∈ 2:size(cmc24_solution, 1)
      α = cmc24_solution[m, 3]

      v = cmc24_solution[m, 1:2]
      e = [cos(α), sin(α)]
      n = [-sin(α), cos(α)]  # normal

      mirror = (
         v1=v,
         v2=v + mirror_length * e,
         α=α,
         e=e,
         n=n,
      )

      push!(mirrors, mirror)
   end

   println(stderr, "The solution is loaded.")

   return (lamp, mirrors)
end

function point_in_block(temple, point)
   for cell ∈ temple.blocks
      # if the point is within bottom-left and top-right vertex
      if all(cell.v1 .≤ point .≤ cell.v3)
         return true
      end
   end

   return false
end

function point_sector(temple, point)
   sx = floor(Int, 3 * point[1] / temple.size[1])
   sy = floor(Int, 3 * point[2] / temple.size[2])

   return 3 * sy + sx + 1
end

function ray_ray_intersection(ray1, ray2)
   (p, r) = ray1
   (q, s) = ray2

   rs = r × s
   qpr = (q - p) × r

   # CASE 1 - rays are collinear and maybe overlap
   if (rs == 0) && (qpr == 0)
      t0 = (q - p) ⋅ r / (r ⋅ r)
      t1 = (q + s - p) ⋅ r / (r ⋅ r)
      return (1, t0, t1)
   end

   # CASE 2 - rays are parallel so they don't intersect
   if (rs == 0) && (qpr ≠ 0)
      return (2, 0, 0)
   end

   # CASE 3 - rays intersect
   qps = (q - p) × s
   t = qps / rs
   u = qpr / rs
   if (rs ≠ 0) && (t ≥ 0) && (u ≥ 0)
      return (3, t, u)
   end

   # CASE 4 - rays don't intersect
   return (4, 0, 0)
end

function ray_segment_intersection(ray, segment)
   (p, r) = ray
   (q, l, β) = segment

   s = l * [cos(β), sin(β)]

   (case, t, u) = ray_ray_intersection((p, r), (q, s))

   # CASE 1 - No intersection
   if case == 1 && t < 0 && u < 0
      return (1, 0, 0)
   end
   if case == 2
      return (1, 0, 0)
   end
   if case == 3 && (t ≤ 0 || u < 0 || u > 1)
      return (1, 0, 0)
   end
   if case == 4
      return (1, 0, 0)
   end

   # CASE 2 - Ray and segment are collinear and they intersect
   if case == 1
      if t > 0 && u ≥ 0
         return (2, min(t, u), 0)
      end
      if t ≥ 0
         return (2, t, 0)
      end
      if u ≥ 0
         return (2, 0, 0)
      end
   end

   # CASE 3 - Ray and segment intersect in ordinary way
   return (3, t, u)
end

function segment_segment_intersection(segment1, segment2)
   (p, la, α) = segment1
   (q, lb, β) = segment2

   r = la * [cos(α), sin(α)]
   s = lb * [cos(β), sin(β)]

   (case, t, u) = ray_ray_intersection((p, r), (q, s))

   if case == 1 && r ⋅ s > 0 && t ≤ u && (0 ≤ t ≤ 1 || 0 ≤ u ≤ 1)
      return true
   end

   if case == 1 && r ⋅ s < 0 && t ≥ u && (0 ≤ t ≤ 1 || 0 ≤ u ≤ 1)
      return true
   end

   if case == 2
      return false
   end

   if case == 3 && (0 ≤ t ≤ 1 && 0 ≤ u ≤ 1)
      return true
   end

   if case == 4
      return false
   end

   return false
end

function segment_block_intersection(segment, block)
   return any((
      segment_segment_intersection(segment, block.s1),
      segment_segment_intersection(segment, block.s2),
      segment_segment_intersection(segment, block.s3),
      segment_segment_intersection(segment, block.s4),
   ))
end

function temple_segment_intersection(temple, segment)
   return any(segment_block_intersection(segment, block) for block ∈ temple.blocks)
end

function temple_ray_intersection(temple, ray)
   t_min = ∞
   for block ∈ temple.blocks
      for segment ∈ [block.s1, block.s2, block.s3, block.s4]
         (case, t, u) = ray_segment_intersection(ray, segment)
         if (case == 2 || case == 3) && (t < t_min) && (t > ε)
            t_min = t
         end
      end
   end

   return t_min
end

function check_solution(temple, lamp, mirrors)
   # check the lamp is within the temple
   if !all([0, 0] .≤ lamp.v .≤ temple.size)
      println(stderr, "ERROR! The lamp isn't placed within temple limits which is of size $(temple.size).")
      finalize(temple, lamp, mirrors)
   end

   # check mirrors' ends are within the temple
   if !all(all([0, 0] .≤ mirror.v1 .≤ temple.size) for mirror ∈ mirrors)
      println(stderr, "ERROR! Some mirror isn't placed within temple of size $(temple.size).")
      finalize(temple, lamp, mirrors)
   end

   if !all(all([0, 0] .≤ mirror.v2 .≤ temple.size) for mirror ∈ mirrors)
      println(stderr, "ERROR! Some mirror isn't placed within temple of size $(temple.size).")
      finalize(temple, lamp, mirrors)
   end

   # check the lamp isn't in some building block
   if point_in_block(temple, lamp.v)
      println(stderr, "ERROR! Lamp is placed in a building block.")
      finalize(temple, lamp, mirrors)
   end

   # check some mirror end isn't in some building block
   for (m, mirror) ∈ enumerate(mirrors)
      if point_in_block(temple, mirror.v1) || point_in_block(temple, mirror.v2)
         println(stderr, "ERROR! Mirror $m has one of its ends inside a building block.")
         finalize(temple, lamp, mirrors)
      end
   end

   # check some mirror doesn't overlap with some building block
   for (m, mirror) ∈ enumerate(mirrors)
      if temple_segment_intersection(temple, (mirror.v1, mirror_length, mirror.α))
         println(stderr, "ERROR! Mirror $m intersects with a building block.")
         finalize(temple, lamp, mirrors)
      end
   end

   # check if some mirrors intersect
   for (m1, mirror1) ∈ enumerate(mirrors[1:end-1]), (m2, mirror2) ∈ enumerate(mirrors[m1+1:end])
      if segment_segment_intersection((mirror1.v1, mirror_length, mirror1.α), (mirror2.v1, mirror_length, mirror2.α))
         println(stderr, "ERROR! Mirrors $m1 & $m2 intersect.")
         finalize(temple, lamp, mirrors)
      end
   end

   println(stderr, "The solution geometry is correct.")
end

function raytrace(temple, lamp, mirrors)
   local hit_mirror

   path = (
      points=[lamp.v],
      directions=[],
   )

   ray = (
      v=lamp.v,
      e=lamp.e,
   )

   hit_mirrors = []
   while true
      # check if ray can hit some mirror
      t_mirror = ∞
      for (m, mirror) ∈ enumerate(mirrors)
         (case, t, u) = ray_segment_intersection(ray, (mirror.v1, mirror_length, mirror.α))
         if ((case == 2) || (case == 3)) && (t < t_mirror) && (t > ε)
            t_mirror = t
            hit_mirror = mirror
            push!(hit_mirrors, m)
         end
      end

      # check where ray would hit the temple
      t_temple = temple_ray_intersection(temple, ray)

      # closest hit point
      t = min(t_mirror, t_temple)
      hitting_point = ray.v + t * ray.e
      push!(path.directions, ray.e)
      push!(path.points, hitting_point)

      # ray hit a mirror, calculate new direction
      if t_mirror < t_temple
         ray = (
            v=hitting_point,
            e=ray.e - 2 * (ray.e ⋅ hit_mirror.n) * hit_mirror.n,
         )
         continue
      end

      # ray hit the temple
      break
   end

   return path
end

function correlate_color1(value::Float64)

   # Interpolate between red (1, 0, 0) and green (0, 1, 0)
   r = 1.0 - value  # Red decreases from 1 to 0
   g = value        # Green increases from 0 to 1
   b = 0.0         # Blue is constant at 0

   # Return the RGB color
   return RGBA(r, g, b, 1)
end

function correlate_color2(value::Float64)

   if value <= 0.25
      # Transition from Red to Orange
      r = 1
      g = (value / 0.25)
      b = 0
   elseif value <= 0.5
      # Transition from Orange to Yellow
      r = 1 - (value - 0.25) / 0.25
      g = 1
      b = 0
   elseif value <= 0.75
      # Transition from Yellow to Light Green
      r = 0
      g = 1
      b = (value - 0.5) / 0.25
   else
      # Transition from Light Green to Green
      r = 0
      g = 1
      b = 1 - (value - 0.75) / 0.25
   end
   return RGBA(r, g, b, 1)
end

function map_to_red_green(value::Float64)
   # Ensure x is clamped between -1 and 1
   value = clamp(x, -1, 1)

   # Calculate red and green components based on the value of x
   red = (1 + value) / 2  # Red increases from 0 to 1 as x goes from -1 to 1
   green = (1 - value) / 2  # Green decreases from 1 to 0 as x goes from -1 to 1

   return RGBA(red, green, 0, 1)  # Return the RGB color, with blue set to 0
end

function batlow_palette()
   return [
      RGB(0.042, 0.167, 0.478),  # Dark Blue
      RGB(0.031, 0.385, 0.675),  # Blue
      RGB(0.263, 0.635, 0.878),  # Light Blue
      RGB(0.407, 0.866, 0.886),  # Cyan
      RGB(0.745, 0.886, 0.675),  # Light Green
      RGB(0.922, 0.933, 0.471),  # Yellow Green
      RGB(0.992, 0.682, 0.314),  # Orange
      RGB(0.973, 0.330, 0.109),  # Red
      RGB(0.784, 0.063, 0.189),  # Dark Red
      RGB(0.356, 0.019, 0.152),   # Dark Purple
   ]
end

# Step 4: Function to interpolate between colors
function interpolate_color(c1::RGB, c2::RGB, t::Float64)
   return RGB(
      (1 - t) * red(c1) + t * red(c2),
      (1 - t) * green(c1) + t * green(c2),
      (1 - t) * blue(c1) + t * blue(c2),
   )
end

# Step 5: Create a continuous Batlow color function
function batlow_color(value::Float64)
   # Ensure the value is within the range [0, 1]
   if value < 0 || value > 1
      throw(ArgumentError("Input value must be in the range [0, 1]."))
   end
   value = 1 - value
   # Get the Batlow color palette
   palette = batlow_palette()
   n = length(palette)

   # Calculate the index and fractional part
   index = Int(floor(value * (n - 1)))  # Find the lower index
   t = (value * (n - 1)) - index         # Fractional part

   # Interpolate between the two colors
   if index < n - 1
      return interpolate_color(palette[index+1], palette[index+2], t)
   else
      return palette[end]  # Return the last color if at the end
   end
end

function grid_plot(temple, heatgrid_matrix, resgrid, n, fi_ε)
   plot_scale = 150
   plot_size = plot_scale .* temple.shape

   #=
   if minimum(heatgrid_matrix[:, 3]) < 0
      clr = batlow_color
   else
      clr = map_to_red_green
   end
   =#
   clr = batlow_color

   println("DONE")

   plot(
      size=plot_size,
      xlims=(0, temple.size[1]),
      ylims=(0, temple.size[2]),
      background_color=RGBA(0.95, 0.95, 0.95, 1),
      label=false,
      showaxis=false,
      grid=false,
      legend=false,
      aspect_ratio=1,
      bottom_margin=-20mm,
      right_margin=-10mm,
      top_margin=-10mm,
      left_margin=-20mm,
   )



   function circleShape(x, y, r, n)
      θ = LinRange(0, 2π, n + 1)
      return Shape(x .+ r * cos.(θ), y .+ r * sin.(θ))
   end

   for block ∈ temple.blocks
      (x, y) = block.v1
      plot!(
         Shape(
            x .+ [0, block_size, block_size, 0],
            y .+ [0, 0, block_size, block_size]),
         color=RGBA(0.50, 0.48, 0.47, 1),
         linecolor=RGBA(0, 0, 0, 0),
         linewidth=0,
      )
   end

   dims = size(heatgrid_matrix)
   println("STEPS: $(dims[1]*dims[2])")
   c = 0
   for i in 1:dims[1]
      for j in 1:dims[2]
         x = i * resgrid + 1
         y = j * resgrid + 1
         plot!(
            circleShape(x, y, resgrid / 2, 10),
            color=clr(heatgrid_matrix[i, j][1]),
            linecolor=RGBA(0, 0, 0, 0),
            linewidth=0,
         )
      end
      c += 1
      print("$(c)|")
   end
   println()


   savefig("heatgrid_plot$(resgrid)_$(fi_ε)_$(n)diff.png")
   return


end

function cmc24_plot(temple; lamp=nothing, mirrors=nothing, path=nothing)
   plot_scale = 150
   plot_size = plot_scale .* temple.shape

   plot(
      size=plot_size,
      xlims=(0, temple.size[1]),
      ylims=(0, temple.size[2]),
      background_color=RGBA(0.9, 0.87, 0.7, 1),
      label=false,
      showaxis=false,
      grid=false,
      legend=false,
      aspect_ratio=1,
      bottom_margin=-20mm,
      right_margin=-10mm,
      top_margin=-10mm,
      left_margin=-20mm,
   )

   function circleShape(x, y, r, n)
      θ = LinRange(0, 2π, n + 1)
      return Shape(x .+ r * cos.(θ), y .+ r * sin.(θ))
   end

   # plot the lightened area
   if path ≠ nothing
      # circle parts of the lightened area
      for p ∈ path.points
         plot!(
            circleShape(p[1], p[2], light_halfwidth, 1000),
            color=RGBA(1, 0.7, 0.6, 1),
            linecolor=RGBA(0, 0, 0, 0),
            linewidth=0,
         )
      end

      # rectangle parts of the lightened area
      for (p1, p2, e) ∈ zip(path.points, path.points[2:end], path.directions)
         (x1, y1) = p1
         (x2, y2) = p2
         (nx, ny) = (e[2], -e[1])

         xs = [x1 - nx, x2 - nx, x2 + nx, x1 + nx]
         ys = [y1 - ny, y2 - ny, y2 + ny, y1 + ny]

         plot!(
            Shape(xs, ys),
            color=RGBA(1, 0.7, 0.6, 1),
            linecolor=RGBA(0, 0, 0, 0),
            linewidth=0,
         )
      end
   end

   # plot the mirrors
   if mirrors ≠ nothing
      for mirror ∈ mirrors
         (x1, y1) = mirror.v1
         (x2, y2) = mirror.v2
         (nx, ny) = 0.05 * mirror.n

         xs = [x1 - nx, x2 - nx, x2 + nx, x1 + nx]
         ys = [y1 - ny, y2 - ny, y2 + ny, y1 + ny]

         plot!(
            Shape(xs, ys),
            color=RGBA(0, 0, 1, 1),
            linecolor=RGBA(0, 0, 0, 0),
            linewidth=0,
         )
      end
   end

   # plot the ray
   if path ≠ nothing
      plot!(
         first.(path.points),
         last.(path.points),
         linecolor=RGBA(1, 0, 0, 1),
         linewidth=0.04 * plot_scale,
      )
   end

   # plot the lamp
   if lamp ≠ nothing
      plot!(
         circleShape(lamp.v[1], lamp.v[2], 0.2, 6),
         color=RGBA(0.9, 0, 1, 1),
         linecolor=RGBA(1, 1, 1, 1),
         linewidth=5,
      )
   end

   # plot the building blocks
   for block ∈ temple.blocks
      (x, y) = block.v1
      plot!(
         Shape(
            x .+ [0, block_size, block_size, 0],
            y .+ [0, 0, block_size, block_size]),
         color=RGBA(0.50, 0.48, 0.47, 1),
         linecolor=RGBA(0, 0, 0, 0),
         linewidth=0,
      )
   end

   solution_hash = hex12(hash([temple, lamp, mirrors, path]))
   uuid = hex12(UUIDs.uuid4().value)
   filename = "cmc24_solution_" * solution_hash * "_" * uuid * ".png"
   savefig(filename)

   return filename
end

function evaluate(temple, path)
   fplot1 = cmc24_plot(temple)
   fplot2 = cmc24_plot(temple, path=path)

   img1 = FileIO.load(fplot1)
   img2 = FileIO.load(fplot2)

   # count the total number of the plot pixels
   total = length(img1)

   # count the number of vacant pixels recognized by being bright
   vacant = sum(p.r > 0.7 for p ∈ img1)

   # count the number of pixels changed due to the light ray
   score = sum(p1 ≠ p2 for (p1, p2) ∈ zip(img1, img2))

   return total, vacant, score
end

function finalize(temple=nothing, lamp=nothing, mirrors=nothing, path=nothing)
   if temple ≠ nothing
      cmc24_plot(temple, lamp=lamp, mirrors=mirrors, path=path)
   end

   print(0)  # stdout print of a fallback result in a case of early exit

   exit()
end

function print_parameters(lamp, mirrors)
   print(lamp.v[1], " ", lamp.v[2], " ", lamp.α)
   for mirror in mirrors
      println(";")
      print(mirror.v1[1], " ", mirror.v1[2], " ", mirror.α)
   end
   println()
end

function create_pars(lamp, mirrors)
   println(mirrors)
   println(lamp)
   cmc24_solution = zeros(Float64, 9, 3)
   cmc24_solution[1, 1] = lamp.v[1]
   cmc24_solution[1, 2] = lamp.v[2]
   cmc24_solution[1, 3] = lamp.α

   for i in 2:9  # Loop over rows
      cmc24_solution[i, 1] = mirrors[i-1].v1[1]
      cmc24_solution[i, 2] = mirrors[i-1].v1[2]
      cmc24_solution[i, 3] = mirrors[i-1].α
   end

   return cmc24_solution
end

function main()
   #=
   cmc24_solution = [
   	15.45 18.91 4.84;
   	17.8442373553 1.5260603781 -2.43434867321;
   	11.9 16.5 0.95;
   	15.2 17.6 2.45;
   	13.8 12.0 0.92;
   	1.6 6.2 2.53;
   	2.2 14.7 0.7;
   	8.5 14.2 2.325;
   	8.7 3.05 2.525
   ]=#
   #=cmc24_solution = [
   			   15.44 18.740000000000002 4.838999999999999;
   			   17.44483078674054 1.0547062779840373 0.3388146928204122;
   			   10.075717649968665 13.92520787619475 3.1618146928204123;
   			   2.3394450081080054 1.152258526200875 -0.2813706143591739;
   			   3.351580742743956 14.933603017834102 3.818629385640827;
   			   18.777266219955518 11.874165569569342 1.2204440784612416;
   			   11.156337101419227 16.97389816157055 3.6654440784612423]=#
   #cmc24_solution = [5.0 5.0 0.26; 10.5 6.5 0.9; 10.9 16.5 0.95; 14.2 17.6 2.45; 13.8 12.0 0.92; 0.6000000000000001 6.2 2.53; 2.2 14.7 0.7; 7.5 14.2 2.325; 7.699999999999999 3.05 2.525]
   # cmc24_solution = [5.0 5.0 0.26; 11.5 6.5 0.89; 11.9 16.5 0.96; 15.2 17.6 2.45; 13.8 12.0 0.92; 1.6 6.2 2.53; 2.2 14.7 0.7; 8.5 14.2 2.325; 8.7 3.05 2.538]
   # cmc24_solution = [5.0 5.0 0.26; 11.5 6.5 0.89; 11.9 16.5 0.96; 15.2 17.6 2.45; 13.8 12.0 0.92; 1.6 6.2 2.53; 2.2 14.7 0.7; 8.5 14.2 2.325; 8.7 3.05 2.538]
   # load the temple

   #=cmc24_solution = [
   			   15.44 18.740000000000002 4.838999999999999;
   			   17.94144862808972 2.0992674118612387 -2.41734867321
   			]=#
   cmc24_solution = [
      18.9899966959408 2.0456944803577 3.1501342250355;
      10.0120712164148 1.9689771571414 3.1415652715598;
      2.4063041114551 2.1136692357806 -0.7280893154205;
      1.0226335567006 13.5984354661066 0.388350749844;
      18.0172580316316 12.8567267177671 5.931103847365;
      17.5776224809363 5.1281172810842 0.8746528405693;
      3.9837798002119613 17.20617704051384 0.1522652160650959;
      8.674412896651 7.3683022170594 -0.6903934384405;
      12.998514834419021 16.45356313240268 2.530702623740865
   ]
   #lamp_mirrors = JLD2.load("914-328.jld2", "lamp_mirrors")
   #=cmc24_solution = [
      8.2 1.05 1.955;
      1.511282553966332 18.777019514155974 4.0025;
      17.999382047667247 14.246360007657092 2.0184073464102075;
      12.00162535314205 2.030116236592267 -0.2040926535897929;
      11.975290046657713 13.282072734711987 2.307739240961644;
      2.25474834444133 14.487338958039668 0.7936465873718507;
      4.693509218519093 1.0000435609535732 3.0543010339512917;
      7.33225422573161 10.083618226045012 0.49020838036149944;
      17.98606137702609 7.201722591937579 2.2309073464102074
   ]=#

   #mirrors = JLD2.load("BEST_SCORE_15870_0000001_0.1.jld2", "mirrors")
   temple = load_temple(temple_string, block_size)

   # load the solution
   #lamp, mirrors = (lamp_mirrors[1], lamp_mirrors[2:end])
   lamp, mirrors = load_solution_partial(cmc24_solution, mirror_length)
   print_parameters(lamp, mirrors)
   check_solution(temple, lamp, mirrors)
   # compute the ray path
   path = raytrace(temple, lamp, mirrors)
   println(path)
   # evaluate the solution
   total, vacant, score = evaluate(temple, path)
   println(stderr, "Base plot has $(commas(vacant)) vacant of total $(commas(total)) pixels.")
   println(stderr, "Your CMC24 score is $(commas(score)) / $(commas(vacant)) = $(100. * score / vacant) %.")
   println(score)

   # create the presentation plot
   cmc24_plot(temple, lamp=lamp, mirrors=mirrors, path=path)
end
t = @timed begin
   main()
end

println("Execution time: ", t.time)
println("Memory allocated: ", t.bytes * 1e-6, " MB")
#julia --track-allocation=user cmc24edit.jl
