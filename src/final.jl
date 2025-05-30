using JLD2

#DEFINITIONS

#=
lamp = (
		v = [x, y],
		α = α,
		e = [cos(α), sin(α)],
	)
=#

#=
mirror = (
			v1 = v,
			v2 = v + mirror_length * e,
			α = α,
			e = e,
			n = n,
		)
=#

#=
path = (
		points = [lamp.v, A, B,...],
		directions = [lamp.e, rayA.e, rayB.e,...],
=#

#=
	ray = (
		v = v,
		e = e,
	)
=#

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
mirror_length_half = 0.25
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

#LOADERS

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
      # finalize()
   end

   if !(eltype(cmc24_solution) <: Number)
      println(stderr, "ERROR! The solution contains non-numerical inputs.")
      # finalize()
   end

   try
      cmc24_solution = float(cmc24_solution)
   catch
      println(stderr, "ERROR! The solution can't be converted to double precision floating point format.")
      # finalize()
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
         v2=v + mirror_length * e, α=α,
         e=e,
         n=n,
      )

      push!(mirrors, mirror)
   end

   println(stderr, "The solution is loaded.")

   return (lamp, mirrors)
end

#PRINTERS/TRANSFORMERS

function print_mirror(mirror)
   println("MIRROR: (x1 y1||x2 y2)=$(mirror.v1[1]) $(mirror.v1[2])||$(mirror.v2[1]) $(mirror.v2[2]), α=$(mirror.α)")
end

function print_mirrors(mirrors)
   if length(mirrors) == 0
      return
   end
   i = 1
   println("MIRRORS-----------")
   for mirror in mirrors
      println("MIRROR ($(i)):")
      print_mirror(mirror)
      i += 1
   end
   println("------------------\n")
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

#CONDITION CHECKING/RAY CALCULATING

function is_within_one_unit_of_segment(A, B, P)
   # Unpack tuples
   (x1, y1) = A
   (x2, y2) = B
   (x0, y0) = P

   # Vector AB
   ABx = x2 - x1
   ABy = y2 - y1

   # Vector AP
   APx = x0 - x1
   APy = y0 - y1

   # Calculate the squared length of AB
   AB_squared = squared_distance(A, B)

   if AB_squared == 0
      # A and B are the same point, check distance to A (or B)
      return squared_distance(A, P) <= 1.0
   end

   # Project point P onto the line extending through A and B
   t = (APx * ABx + APy * ABy) / AB_squared

   # Check the closest point on the segment
   if t < 0
      closest_x, closest_y = x1, y1  # Closest to A
   elseif t > 1
      closest_x, closest_y = x2, y2  # Closest to B
   else
      closest_x = x1 + t * ABx
      closest_y = y1 + t * ABy
   end

   # Check if the distance to the closest point is within 1 unit
   return squared_distance((closest_x, closest_y), P) <= light_halfwidth
end

function is_within_epsilon_of_segment(A, B, P)
   # Unpack tuples
   (x1, y1) = A
   (x2, y2) = B
   (x0, y0) = P

   # Vector AB
   ABx = x2 - x1
   ABy = y2 - y1

   # Vector AP
   APx = x0 - x1
   APy = y0 - y1

   # Calculate the squared length of AB
   AB_squared = squared_distance(A, B)

   if AB_squared == 0
      # A and B are the same point, check distance to A (or B)
      return squared_distance(A, P) <= 1.0
   end

   # Project point P onto the line extending through A and B
   t = (APx * ABx + APy * ABy) / AB_squared

   # Check the closest point on the segment
   if t < 0
      closest_x, closest_y = x1, y1  # Closest to A
   elseif t > 1
      closest_x, closest_y = x2, y2  # Closest to B
   else
      closest_x = x1 + t * ABx
      closest_y = y1 + t * ABy
   end

   # Check if the distance to the closest point is within 1 unit
   return squared_distance((closest_x, closest_y), P) <= 1e-4
end

function point_illuminated(path, point)
   points = path.points
   return any(is_within_one_unit_of_segment(A, B, point) for (A, B) in [(points[i], points[i+1]) for i in 1:length(points)-1])
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

function mirrors_oversect(mirrors)
   # check if some mirrors intersect
   for (m1, mirror1) ∈ enumerate(mirrors[1:end-1]), (m2, mirror2) ∈ enumerate(mirrors[m1+1:end])
      if segment_segment_intersection((mirror1.v1, mirror_length, mirror1.α), (mirror2.v1, mirror_length, mirror2.α))
         println(stderr, "ERROR! Mirrors $m1 & $m2 intersect.")
         return true
         # finalize(temple, lamp, mirrors)
      end
   end
   return false
end

function segment_mirrors_oversect(mirrors, segment)
   return any(segment_segment_intersection(segment, (mirror.v1, mirror_length, mirror.α)) for mirror in mirrors)
end

function check_solution(temple, lamp, mirrors)
   # check the lamp is within the temple
   if !all([0, 0] .≤ lamp.v .≤ temple.size)
      println(stderr, "ERROR! The lamp isn't placed within temple limits which is of size $(temple.size).")
      # finalize(temple, lamp, mirrors)
   end

   # check mirrors' ends are within the temple
   if !all(all([0, 0] .≤ mirror.v1 .≤ temple.size) for mirror ∈ mirrors)
      println(stderr, "ERROR! Some mirror isn't placed within temple of size $(temple.size).")
      # finalize(temple, lamp, mirrors)
   end

   if !all(all([0, 0] .≤ mirror.v2 .≤ temple.size) for mirror ∈ mirrors)
      println(stderr, "ERROR! Some mirror isn't placed within temple of size $(temple.size).")
      # finalize(temple, lamp, mirrors)
   end

   # check the lamp isn't in some building block
   if point_in_block(temple, lamp.v)
      println(stderr, "ERROR! Lamp is placed in a building block.")
      # finalize(temple, lamp, mirrors)
   end

   # check some mirror end isn't in some building block
   for (m, mirror) ∈ enumerate(mirrors)
      if point_in_block(temple, mirror.v1) || point_in_block(temple, mirror.v2)
         println(stderr, "ERROR! Mirror $m has one of its ends inside a building block.")
         # finalize(temple, lamp, mirrors)
      end
   end

   # check some mirror doesn't overlap with some building block
   for (m, mirror) ∈ enumerate(mirrors)
      if temple_segment_intersection(temple, (mirror.v1, mirror_length, mirror.α))
         println(stderr, "ERROR! Mirror $m intersects with a building block.")
         # finalize(temple, lamp, mirrors)
      end
   end

   # check if some mirrors intersect
   for (m1, mirror1) ∈ enumerate(mirrors[1:end-1]), (m2, mirror2) ∈ enumerate(mirrors[m1+1:end])
      if segment_segment_intersection((mirror1.v1, mirror_length, mirror1.α), (mirror2.v1, mirror_length, mirror2.α))
         println(stderr, "ERROR! Mirrors $m1 & $m2 intersect.")
         # finalize(temple, lamp, mirrors)
      end
   end

   # println(stderr, "The solution geometry is correct.")
end

function mirror_visible(temple, mirror, P)

   for k in -0.3:0.2:0.3
      M = mirror.v1 + (k * mirror_length) * mirror.e
      d = segment_length(P, M)
      e = (M .- P) ./ d
      rayP = (
         v=P,
         e=e,
      )
      tP = temple_ray_intersection(temple, rayP)

      rayM = (
         v=P,
         e=-e,
      )
      tM = temple_ray_intersection(temple, rayM)

      if (tP + tM) > d
         return true
      end
   end

   return false
end

#SIMPLE CALCULATING FUNCTIONS

function line_angle(x1, y1, x2, y2)
   return atan(y2 - y1, x2 - x1)
end

function segment_length(A, B)
   (x1, y1) = A
   (x2, y2) = B
   return sqrt((x2 - x1)^2 + (y2 - y1)^2)
end

function squared_distance(p1, p2)
   return (p2[1] - p1[1])^2 + (p2[2] - p1[2])^2
end

function grid_floatpoint(point, resgrid)
   x, y = point
   return [round((x - 1) / resgrid) * resgrid + 1, round((y - 1) / resgrid) * resgrid + 1]
end

function grid_indexpoint(point, resgrid, dim_M=calculate_dim_M(resgrid))
   pointx, pointy = point
   (gridpointx, gridpointy) = (Int(round((pointx - 1) / resgrid)), Int(round((pointy - 1) / resgrid)))

   gridpointx = gridpointx == 0 ? 1 : (gridpointx == dim_M + 1 ? dim_M : gridpointx)
   gridpointy = gridpointy == 0 ? 1 : (gridpointy == dim_M + 1 ? dim_M : gridpointy)

   return [gridpointx, gridpointy]
end

function floatpoint_gridpoint(gridpoint, resgrid)
   gridpointx, gridpointy = gridpoint
   return [gridpointx * resgrid + 1, gridpointy * resgrid + 1]
end

function calculate_dim_M(resgrid)
   return Int(round(18 / resgrid) - 1)
end

function distance_from_segment(A, B, P)
   # Unpack tuples
   (x1, y1) = A
   (x2, y2) = B
   (x0, y0) = P

   # Vector AB
   ABx = x2 - x1
   ABy = y2 - y1

   # Vector AP
   APx = x0 - x1
   APy = y0 - y1

   # Calculate the squared length of AB
   AB_squared = squared_distance(A, B)

   if AB_squared == 0
      # A and B are the same point, check distance to A (or B)
      return 0
   end

   # Project point P onto the line extending through A and B
   t = (APx * ABx + APy * ABy) / AB_squared

   # Check the closest point on the segment
   if t < 0
      closest_x, closest_y = x1, y1  # Closest to A
   elseif t > 1
      closest_x, closest_y = x2, y2  # Closest to B
   else
      closest_x = x1 + t * ABx
      closest_y = y1 + t * ABy
   end

   # Check if the distance to the closest point is within 1 unit
   return segment_length((closest_x, closest_y), P)
end

function score_w_lin(score, mirror_count, gridpoint_v)
   return score * (1 + gridpoint_v / (mirror_count * 2 + 1))
end

function score_w_lin_rnd(score, mirror_count, gridpoint_v)
   return score * ((1 + gridpoint_v / (mirror_count * 2 + 1)) + rand() / 10)
end

function score_w(score, mirror_count, gridpoint_v)
   return score
end
#LAMP OPERATIONS

function rotate_lamp(lamp, α)
   lamp_rotated = (
      v=lamp.v,
      α=lamp.α + α,
      e=[cos(lamp.α + α), sin(lamp.α + α)],
   )
   return lamp_rotated
end



#MIRROR OPERATIONS

function rotate_mirror(mirror, α)
   # Step 1: Calculate the midpoint
   mid_x = (mirror.v1[1] + mirror.v2[1]) / 2
   mid_y = (mirror.v1[2] + mirror.v2[2]) / 2

   # Step 2: Rotate (x1, y1) around the midpoint
   x1_rot = mid_x + (mirror.v1[1] - mid_x) * cos(α) - (mirror.v1[2] - mid_y) * sin(α)
   y1_rot = mid_y + (mirror.v1[1] - mid_x) * sin(α) + (mirror.v1[2] - mid_y) * cos(α)

   # Step 3: Rotate (x2, y2) around the midpoint
   x2_rot = mid_x + (mirror.v2[1] - mid_x) * cos(α) - (mirror.v2[2] - mid_y) * sin(α)
   y2_rot = mid_y + (mirror.v2[1] - mid_x) * sin(α) + (mirror.v2[2] - mid_y) * cos(α)

   mirror_rotated = (
      v1=[x1_rot, y1_rot],
      v2=[x2_rot, y2_rot],
      α=mirror.α + α,
      e=[cos(mirror.α + α), sin(mirror.α + α)],
      n=[-sin(mirror.α + α), cos(mirror.α + α)],
   )
   return mirror_rotated
end

function translate_mirror(mirror, len)
   mirror_translated = (
      v1=mirror.v1 + len * mirror.e,
      v2=mirror.v2 + len * mirror.e,
      α=mirror.α,
      e=mirror.e,
      n=mirror.n,
   )
   return mirror_translated
end

function mirrors_different(mirror1, mirror2)
   mirror1_s = (mirror1.v1 .+ mirror1.v2) / 2
   mirror2_s = (mirror2.v1 .+ mirror2.v2) / 2
   return squared_distance(mirror1_s, mirror2_s) > mirror_length_half || abs(mirror1.α - mirror2.α) > 0.1
end

function adjust_mirror_pars(temple, mirrors, path, point, α)
   for k in -mirror_length+0.001:0.002:0
      if !temple_segment_intersection(temple, (point + k .* [cos(α), sin(α)], mirror_length, α)) &&
         !segment_mirrors_oversect(mirrors, (point + k .* [cos(α), sin(α)], mirror_length, α)) &&
         !any(segment_segment_intersection((point + k .* [cos(α), sin(α)], mirror_length, α),
            (A, segment_length(A, B), line_angle(A[1], A[2], B[1], B[2]))) for (A, B) in zip(path.points[1:end-2], path.points[2:end-1]))
         return (true, point + k .* [cos(α), sin(α)])
      end
   end
   return (false, [0, 0])
end

function adjust_mirror_pars_nopath(temple, mirrors, point, α)
   for k in -mirror_length+0.001:0.002:0
      if !temple_segment_intersection(temple, (point + k .* [cos(α), sin(α)], mirror_length, α)) &&
         !segment_mirrors_oversect(mirrors, (point + k .* [cos(α), sin(α)], mirror_length, α))
         return (true, point + k .* [cos(α), sin(α)])
      end
   end
   return (false, [0, 0])
end

function place_mirror_near_gridpoint(temple, mirrors, path, gridpointx, gridpointy, resgrid, α, A, B, α_P=line_angle(A[1], A[2], B[1], B[2]), d=segment_length(A, B))
   (pointx, pointy) = floatpoint_gridpoint([gridpointx, gridpointy], resgrid)
   mirror = nothing
   fi = (α + α_P) / 2
   #println("fi=$fi, α=$α")
   for i in -1:2:1
      ray = (
         v=[pointx, pointy],
         e=i .* [cos(α), sin(α)],
      )
      (case, t, u) = ray_segment_intersection(ray, (A, d, α_P))
      #println("RAY: $ray\n CASE:$case\n T:$t")
      if case == 1
         continue

      elseif case == 3

         floatpoint = ray.v + t .* ray.e
         (placed, point) = adjust_mirror_pars(temple, mirrors, path, floatpoint, fi)
         if placed
            mirror = (
               v1=point,
               v2=point + mirror_length * [cos(fi), sin(fi)],
               α=fi,
               e=[cos(fi), sin(fi)],
               n=[-sin(fi), cos(fi)],
            )
         end
      end
   end

   return mirror
end

function mirrors_hitcounts(path, mirrors)
   mirror_count = length(mirrors)
   hitcounts = [0 for k in 1:mirror_count]
   for point in path.points
      for i in 1:mirror_count
         if is_within_epsilon_of_segment(mirrors[i].v1, mirrors[i].v2, point)
            hitcounts[i] += 1
         end
      end
   end
   return hitcounts
end

function hitpoint_index(path, hitpoints)
   index = zeros(Int, length(hitpoints))
   for i in 1:length(path.points)
      for j in 1:length(hitpoints)
         if !isnothing(hitpoints[j]) && segment_length(hitpoints[j], path.points[i]) < 1e-4
            index[j] = i
            break
         end
      end
   end
   return index
end

function moveable_mirrors_hitpoints(path, mirrors, hitcounts)
   mirror_count = length(mirrors)
   hitpoints = []

   for i in 1:mirror_count
      if hitcounts[i] == 1
         for point in path.points
            if is_within_epsilon_of_segment(mirrors[i].v1, mirrors[i].v2, point)
               push!(hitpoints, point)
            end
         end
      else
         push!(hitpoints, nothing)
      end
   end
   return hitpoints
end

function mirrors_three_points(temple, mirrors, hitpoints)
   mirror_count = length(mirrors)
   mirrors_three_points = []
   for i in 1:mirror_count

      if isnothing(hitpoints[i])
         push!(mirrors_three_points, nothing)
         continue
      end

      three_points = [[0.0, 0.0], hitpoints[i], [0.0, 0.0]]

      for j in -1:2:1
         ray = (
            v=hitpoints[i],
            e=j .* mirrors[i].e,
         )

         t_mirror = ∞
         for mirror ∈ mirrors
            (case, t, u) = ray_segment_intersection(ray, (mirror.v1, mirror_length, mirror.α))
            if ((case == 2) || (case == 3)) && (t < t_mirror) && (t > ε)
               t_mirror = t
            end
         end
         t_temple = temple_ray_intersection(temple, ray)
         three_points[2+j] = hitpoints[i] + (min(mirror_length, t_mirror, t_temple)) .* ray.e
      end

      push!(mirrors_three_points, three_points)
   end

   return mirrors_three_points
end

#PATH OPERATIONS

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

#points ... [weighted_score, dirs, [gridpointx, gridpointy]]
function find_grid_points_of_segment(matrix, A, B, resgrid, n, mirror_count)

   dim_M = calculate_dim_M(resgrid)
   #println("M_2($(dim_M ))")

   #matrix = [[0, ([0.0, 0] for k in 1:n)...] for i in 1:dim_M, j in 1:dim_M] #REMOVE THIS IN FINAL
   #for k in eachindex(matrix)
   #   println(k, "\n")
   #end

   points = [[0, [0.0, 0.0], [0, 0]] for _ in 1:n]

   #println(typeof(points))
   #println(points[1][3][1])

   (x1, y1) = A
   (pointx, pointy) = grid_floatpoint([x1, y1], resgrid)
   #println("x:$(pointx), y:$(pointy)")
   (gridpointx, gridpointy) = grid_indexpoint([pointx, pointy], resgrid, dim_M)

   #println("i:$(gridpointx), j:$(gridpointy)")


   (x2, y2) = B
   (endpointx, endpointy) = grid_floatpoint([x2, y2], resgrid)
   #("x:$(endpointx), y:$(endpointy)")
   (endgridpointx, endgridpointy) = grid_indexpoint([endpointx, endpointy], resgrid, dim_M)

   #println("i:$(endgridpointx), j:$(endgridpointy)")

   angle_AB = line_angle(x1, y1, x2, y2)

   direction = collect(B .- A) ./ segment_length(A, B)
   #println("DIRECTION_fi: ", direction)

   direction .= (x -> x < 0 ? -1 : 1).(direction)
   direction = Int.(direction)
   #println("DIRECTION: ", direction)

   #c = 0
   while gridpointx != endgridpointx && gridpointy != endgridpointy
      (pointx, pointy) = floatpoint_gridpoint([gridpointx, gridpointy], resgrid)
      for dirs in matrix[gridpointx, gridpointy][2:end]
         if abs(rem(dirs[1] - angle_AB, pi)) < 0.1
            continue
         end
         weighted_score = dirs[2] + segment_length(A, [pointx, pointy])
         #println("\n\n $points \n\n")
         (same_point_found, req_sort) = direct_upgrade(gridpointx, gridpointy, dirs, points, n, resgrid, weighted_score)
         if req_sort
            #println("UNSORTED: $points \n\n")
            points = sort(points, by=x -> -x[1])
            #println("SORTED: $points \n\n")
         end
         if same_point_found
            continue
         end
         for i in 1:n
            if points[i][1] < weighted_score
               insert!(points, i, [weighted_score, dirs, [gridpointx, gridpointy]])
               pop!(points)
               #println("$gridpoints,\n $gridpointx,\n $gridpointy")
               #println(matrix[gridpointx, gridpointy])
               break
            end
         end
      end
      #c += 1

      d1 = distance_from_segment(A, B, [pointx + resgrid * direction[1], pointy])
      d2 = distance_from_segment(A, B, [pointx, pointy + resgrid * direction[2]])
      d1 < d2 ? gridpointx += direction[1] : gridpointy += direction[2]
      #println("STEP: $c (GRIDX,GRIDY) = ($(gridpointx),$(gridpointy)) || ($(d1),$(d2))")
   end
   return (points)
end

function return_mirror_from_gridpoints(temple, lamp, mirrors, matrix, A, path, resgrid, reseval, n, mirror_count, fi, d, calculator)

   points = find_grid_points_of_segment(matrix, A, path.points[end], resgrid, n, mirror_count)

   score_points_best = 0
   mirror_best = nothing

   c = 1
   for point in points
      print("$c|")
      mirror_try = place_mirror_near_gridpoint(temple, mirrors, path, point[3][1], point[3][2], resgrid, point[2][1], A, path.points[end], fi, d)
      #println("POINTS:$(floatpoint_gridpoint(point[3][1],point[3][2],resgrid)), ANGLE:$(point[2][1]), C:$c")
      if !isnothing(mirror_try)
         #println("!")
         #println(mirror_try)
         mirrors_try = vcat(mirrors, mirror_try)
         path_try = raytrace(temple, lamp, mirrors_try)

         score_try = evaluate_path(temple, path_try, reseval)

         if mirror_count < 7
            (i, j) = point[3]
            score_try = calculator(score_try, mirror_count, matrix[i, j][1])

         end
         if score_try > score_points_best
            print("!")
            score_points_best = score_try
            mirror_best = mirror_try
         end
      end

      c += 1
   end

   return (mirror_best, score_points_best)
end

function direct_upgrade(gridpointx, gridpointy, dirs, points, n, resgrid, weighted_score)
   for k in 1:n
      if (abs(points[k][3][1] - gridpointx) < (0.09 / resgrid) || abs(points[k][3][2] - gridpointy) < (0.09 / resgrid)) && abs(rem((points[k][2][1] - dirs[1]), pi)) < 0.001
         if weighted_score > points[k][1]
            points[k] = [weighted_score, dirs, [gridpointx, gridpointy]]
            return (true, true)
         end
         return (true, false)
      end
   end
   return (false, false)
end

function next_mirror_pars(temple, path, mirrors)
   A = path.points[end-1]
   dir = path.directions[end]
   rayA = (
      v=A,
      e=dir,
   )
   B = A + temple_ray_intersection(temple, rayA) .* dir
   fi = line_angle(A[1], A[2], B[1], B[2])
   d = segment_length(A, B)
   mirror_count = length(mirrors)

   return (A, dir, rayA, B, fi, d, mirror_count)
end

function next_mirror_pars_arr(temple, path, mirrors)
   hitcounts = mirrors_hitcounts(path, mirrors)
   hitpoints = moveable_mirrors_hitpoints(path, mirrors, hitcounts)
   mirrors_three_points_array = mirrors_three_points(temple, mirrors, hitpoints)
   return (hitcounts, hitpoints, mirrors_three_points_array)
end

function next_mirror(matrix, resgrid, n, temple, lamp, mirrors, resstep, reseval, resrot, calculator)
   path = raytrace(temple, lamp, mirrors)
   if length(path.points) < 2
      println("		PATH IS EMPTY")
      return
   end
   #score_initial = evaluate_path(temple, path, reseval)

   (A, dir, _, B, fi, d, mirror_count) = next_mirror_pars(temple, path, mirrors)
   (hitcounts, hitpoints, mirrors_three_points_array) = next_mirror_pars_arr(temple, path, mirrors)
   (score_mirrors_moveable_best, score_mirrors_unmoveable_best, mirrors_moveable_best, mirror_unmoveable_best) =
      check_reflections(temple, matrix, resgrid, lamp, mirrors, hitcounts, mirrors_three_points_array, A, B, dir, d, fi, reseval, resstep, resrot, calculator)

   (mirror_points_best, score_points_best) = return_mirror_from_gridpoints(temple, lamp, mirrors, matrix, A, path, resgrid, reseval, n, mirror_count, fi, d, calculator)
   println()
   println("		WEIGHTED SCORES: MOVABLE: $score_mirrors_moveable_best, UNMOVEABLE: $score_mirrors_unmoveable_best, POINTS: $score_points_best")
   return pick_solution(mirrors, score_mirrors_moveable_best, score_mirrors_unmoveable_best, score_points_best, mirrors_moveable_best, mirror_unmoveable_best, mirror_points_best)
end

function check_reflections(temple, matrix, resgrid, lamp, mirrors, hitcounts, mirrors_three_points_array, A, B, dir, d, fi, reseval, resstep, resrot, calculator)
   mirror_count = length(mirrors)
   B_step = B - resstep .* dir
   score_mirrors_moveable_best = 0
   score_mirrors_unmoveable_best = 0
   mirrors_moveable_best = nothing
   mirror_unmoveable_best = nothing

   steps = d / resstep
   println("		CHECKING MIRROR REFLECTIONS:,\n		A:$A\n		B:$B")
   c = 1
   while ⋅(dir, B_step - A - 2 .* dir) > 0
      print("$c|")
      for k in 1:mirror_count
         #print("$k%")
         if !mirror_visible(temple, mirrors[k], B_step)
            continue
         end
         if hitcounts[k] == 1
            #println("MOVEABLE!")
            (score_try, mirror_try, mirror_i_try) = mirror_hit_pars_moveable(temple, lamp, mirrors, k, mirrors_three_points_array[k], B_step, fi, reseval, resrot)
            if isnothing(mirror_try)
               continue
            end

            if mirror_count < 7 && !isnothing(mirror_try) && !isnothing(mirror_i_try)
               (i, j) = grid_indexpoint(grid_floatpoint(B_step, resgrid), resgrid)
               score_try = calculator(score_try, mirror_count, matrix[i, j][1])
            end
            if score_try > score_mirrors_moveable_best && !isnothing(mirror_try) && !isnothing(mirror_i_try)
               mirrors_moveable_best = vcat(mirrors[1:k-1], mirror_i_try, mirrors[k+1:end], mirror_try)
               score_mirrors_moveable_best = score_try
               print("#")
            end
         else
            #println("UNMOVEABLE!")
            (score_try, mirror_try) = mirror_hit_pars_unmoveable(temple, lamp, mirrors, k, B_step, fi, reseval, resrot)
            if isnothing(mirror_try)
               continue
            end
            if mirror_count < 7 && !isnothing(mirror_try)
               (i, j) = grid_indexpoint(grid_floatpoint(B_step, resgrid), resgrid)
               score_try = calculator(score_try, mirror_count, matrix[i, j][1])
            end
            if score_try > score_mirrors_unmoveable_best && !isnothing(mirror_try)
               mirror_unmoveable_best = mirror_try
               score_mirrors_unmoveable_best = score_try
               print("%")
            end
         end
      end

      B_step = B_step + (-resstep * ((40 + c) / (40 + steps))) .* dir

      c += 1
   end

   return (score_mirrors_moveable_best, score_mirrors_unmoveable_best, mirrors_moveable_best, mirror_unmoveable_best)
end



function pick_solution(mirrors, score_mirrors_moveable_best, score_mirrors_unmoveable_best, score_points_best, mirrors_moveable_best, mirror_unmoveable_best, mirror_points_best)
   if score_mirrors_moveable_best >= score_mirrors_unmoveable_best && score_mirrors_moveable_best >= score_points_best

      return [mirrors_moveable_best, score_mirrors_unmoveable_best > score_points_best ? vcat(mirrors, mirror_unmoveable_best) : vcat(mirrors, mirror_points_best)]
   elseif score_mirrors_unmoveable_best >= score_mirrors_moveable_best && score_mirrors_unmoveable_best >= score_points_best

      return [vcat(mirrors, mirror_unmoveable_best), score_mirrors_moveable_best > score_points_best ? mirrors_moveable_best : vcat(mirrors, mirror_points_best)]
   elseif score_points_best != 0

      return [vcat(mirrors, mirror_points_best), score_mirrors_moveable_best > score_mirrors_unmoveable_best ? mirrors_moveable_best : vcat(mirrors, mirror_unmoveable_best)]
   end

   println("\n\n\nERROR, DID NOT FIND A MIRROR\n\n\n")
   return
end

#M1, M2 and M3 (three_points) are, in the order, the "leftmost", the single point on a path a mirror is hitting and the "rightmost" point a mirror can be placed on
#i is the index of a checked mirror
function mirror_hit_pars_moveable(temple, lamp, mirrors, i, three_points, P, α_P, reseval, resrot)
   mirror_i_best = mirrors[i]
   mirror_best = nothing
   (M1, M2, M3) = three_points
   α_M = line_angle(M1[1], M1[2], M3[1], M3[2])
   d_M = segment_length(M1, M3)
   path = raytrace(temple, lamp, mirrors)
   score_best = evaluate_path(temple, path, reseval)
   base_points_num = length(path.points)
   resrot_initial = resrot

   fi_1 = line_angle(P[1], P[2], M1[1], M1[2])
   fi_2 = line_angle(P[1], P[2], M3[1], M3[2])
   if fi_1 > fi_2
      fi_help = fi_1
      fi_1 = fi_2
      fi_2 = fi_help
   end

   α = fi_1
   c = 1
   while α < fi_2
      rayP = (
         v=P,
         e=[cos(α), sin(α)],
      )
      (case, tM, u) = ray_segment_intersection(rayP, (M1, d_M, α_M))
      tT = temple_ray_intersection(temple, rayP)

      if tT > tM && case == 3
         α_try = (α_P + α) / 2
         (can_place, p_place) = adjust_mirror_pars(temple, mirrors, path, P, α_try)
         if can_place

            mirror_try = (
               v1=p_place,
               v2=p_place + mirror_length * [cos(α_try), sin(α_try)],
               α=α_try,
               e=[cos(α_try), sin(α_try)],
               n=[-sin(α_try), cos(α_try)],
            )
            hit_point = P + tM .* rayP.e
            α_i = line_angle(hit_point[1], hit_point[2], M2[1], M2[2])
            mirror_i_try = (
               v1=hit_point - 0.0001 .* [cos(α_i), sin(α_i)],
               v2=hit_point + (mirror_length - 0.0001) * [cos(α_i), sin(α_i)],
               α=α_i,
               e=[cos(α_i), sin(α_i)],
               n=[-sin(α_i), cos(α_i)],
            )
            if temple_segment_intersection(temple, (mirror_i_try.v1, mirror_length, mirror_i_try.α))
               α_i = α_i + pi
               mirror_i_try = (
                  v1=[M2[1], M2[2]] - 0.0001 .* [cos(α_i), sin(α_i)],
                  v2=[M2[1], M2[2]] + (mirror_length - 0.0001) * [cos(α_i), sin(α_i)],
                  α=α_i,
                  e=[cos(α_i), sin(α_i)],
                  n=[-sin(α_i), cos(α_i)],
               )
               if temple_segment_intersection(temple, (mirror_i_try.v1, mirror_length, mirror_i_try.α))
                  α += resrot
                  continue
               end
            end
            path_try = raytrace(temple, lamp, vcat(mirrors[1:i-1], mirror_i_try, mirrors[i+1:end], mirror_try))
            score_try = evaluate_path(temple, path_try, reseval)
            if score_try > score_best
               mirror_i_best = mirror_i_try
               mirror_best = mirror_try
               score_best = score_try
            end
            resrot = resrot_initial / (min((2.0^(length(path_try.points) - base_points_num - 2)), 16))
         end
      end
      α += resrot
      c += 1
      if c > 20000
         break
      end
   end
   return (score_best, mirror_best, mirror_i_best)
end

function mirror_hit_pars_unmoveable(temple, lamp, mirrors, i, P, α_P, reseval, resrot)
   mirror_best = nothing

   fi_1 = line_angle(P[1], P[2], mirrors[i].v1[1], mirrors[i].v1[2])
   fi_2 = line_angle(P[1], P[2], mirrors[i].v2[1], mirrors[i].v2[2])
   path = raytrace(temple, lamp, mirrors)
   score_best = evaluate_path(temple, path, reseval)
   base_points_num = length(path.points)
   resrot_initial = resrot

   if fi_1 > fi_2
      fi_help = fi_1
      fi_1 = fi_2
      fi_2 = fi_help
   end

   α = fi_1
   c = 1
   while α < fi_2
      rayP = (
         v=P,
         e=[cos(α), sin(α)],
      )
      (case, tM, u) = ray_segment_intersection(rayP, (mirrors[i].v1, mirror_length, mirrors[i].α))
      tT = temple_ray_intersection(temple, rayP)

      if tT > tM && case == 3
         α_try = (α_P + α) / 2
         (can_place, p_place) = adjust_mirror_pars(temple, mirrors, path, P, α_try)

         if can_place
            mirror_try = (
               v1=p_place,
               v2=p_place + mirror_length * [cos(α_try), sin(α_try)],
               α=α_try,
               e=[cos(α_try), sin(α_try)],
               n=[-sin(α_try), cos(α_try)],
            )
            path_try = raytrace(temple, lamp, vcat(mirrors, mirror_try))
            score_try = evaluate_path(temple, path_try, reseval)
            if score_try > score_best
               mirror_best = mirror_try
               score_best = score_try
            end
            resrot = resrot_initial / (min((2.0^(length(path_try.points) - base_points_num - 2)), 16))
            #println("			resrot:$resrot")
         end
      end
      α += resrot
      #println("			resrot:$resrot")
      c += 1
      if c > 20000
         break
      end
   end
   return (score_best, mirror_best)
end

#EVALUATOR

function evaluate(temple, path, res)
   score = 0
   st_coo = 1 + min(0.3, res / 3)
   p = (st_coo, st_coo)
   # Double for loop to increment p.x in the inner loop and p.y in the outer loop
   for y ∈ p[2]:res:(temple.size[2]-1)
      for x ∈ p[1]:res:(temple.size[1]-1)
         if point_in_block(temple, (x, y))
            continue
         end
         if point_illuminated(path, (x, y))
            score += 1
         end
      end
   end

   return score
end

function evaluate_path(temple, path, res)
   score = 0
   points = path.points
   len = length(points)
   s = 0.95
   c = 0

   for i in 1:1:len-1
      seg = [points[i], points[i+1]]
      fi = line_angle(seg[1][1], seg[1][2], seg[2][1], seg[2][2])
      e = [cos(fi), sin(fi)]
      n = [-sin(fi), cos(fi)]
      seg[1] -= e
      seg[2] += e
      steps_l = floor((segment_length(seg[1], seg[2])) / res)
      steps_w = floor(max(2, 1 / res))

      gen_values = collect(s * k for k in -1:1/steps_w:1)

      starting_points = [(A, B) for (A, B) in [seg[1] + n .* value for value in gen_values]]

      for j in 0:1:steps_l
         for point in starting_points
            point_curr = (j * res) .* e .+ point
            c += 1

            if point_in_block(temple, point_curr) || any(is_within_one_unit_of_segment(A, B, point_curr) for (A, B) in [(points[s], points[s+1]) for s in 1:i-1])
               continue
            end

            if is_within_one_unit_of_segment(points[i], points[i+1], point_curr)
               score += 1
            end

         end
      end
   end

   return score
end

#TUNERS

function finetune_α1d(temple, lamp, mirrors, mirror_index, score_best, res, fi_ε, fi)
   if mirror_index > length(mirrors)
      println("MIRROR_INDEX IN FINETUNE OUT OF BOUNDS")
      return (false, nothing, 0)
   end
   improv = false
   mirror_try = mirrors[mirror_index]
   mirror_best = mirrors[mirror_index]
   steps = floor(2 * fi / fi_ε)

   for k in 0:steps
      mirror_try = rotate_mirror(mirror_try, (-1)^k * k * fi_ε)
      if temple_segment_intersection(temple, (mirror_try.v1, mirror_length, mirror_try.α))
         continue
      end
      if mirror_index == length(mirrors)
         mirrors_try = vcat(
            mirrors[1:mirror_index-1],
            mirror_try,
         )
      else
         mirrors_try = vcat(
            mirrors[1:mirror_index-1],
            mirror_try,
            mirrors[mirror_index+1:end],
         )
      end
      path = raytrace(temple, lamp, mirrors_try)
      score = evaluate_path(temple, path, res)

      if score > score_best
         improv = true
         score_best = score
         mirror_best = mirror_try
      end
   end
   mirrors[mirror_index] = mirror_best
   return (improv, mirrors, score_best)
end

function finetune_lamp_α1d(temple, lamp, mirrors, score_best, res, fi_ε, fi)
   improv = false
   lamp_best = lamp
   lamp_try = lamp
   steps = floor(2 * fi / fi_ε)

   for k in 0:steps
      lamp_try = rotate_lamp(lamp_try, (-1)^k * k * fi_ε)
      path = raytrace(temple, lamp_try, mirrors)
      score = evaluate_path(temple, path, res)

      if score > score_best
         improv = true
         score_best = score
         lamp_best = lamp_try
      end
   end
   return (improv, lamp_best, score_best)
end

function last_optimization(temple, lamp, mirrors, tolerance, reseval)
   path = raytrace(temple, lamp, mirrors)
   hitcounts = mirrors_hitcounts(path, mirrors)
   hitpoints = moveable_mirrors_hitpoints(path, mirrors, hitcounts)
   index = hitpoint_index(path, hitpoints)
   #println(index)
   println("	HITCOUNTS AND HITPOINTS\n", hitcounts, "\n", hitpoints)
   len = length(hitcounts)
   for i in 1:len-2
      print("	MIRROR OPT: $i")
      if hitcounts[i] == 1 && hitcounts[i+1] == 1 && index[i+1] - index[i] == 1
         if hitcounts[i+2] != 1 || index[i+2] - index[i+1] != 1
            println(", 0 tolerance")
            mirrors = optimize_last_2(temple, lamp, mirrors, i, i + 1, path, hitpoints[i:i+1], 0, reseval)
         else
            println(", with tolerance")
            mirrors = optimize_last_2(temple, lamp, mirrors, i, i + 1, path, hitpoints[i:i+1], tolerance, reseval)
         end
      else
         println(", not able to optimize")
      end
      path = raytrace(temple, lamp, mirrors)

   end

   return mirrors
end

function optimize_last_2(temple, lamp, mirrors, i, j, path, hitpoints, tolerance, reseval)

   score_best = evaluate_path(temple, path, reseval)
   #println("!!!!!!!!!!!!	$score_best")
   p_len = length(path.points)
   if p_len == 0
      println("		ADJUSTMENT UNSUCCESSFUL")
      return mirrors
   end

   mirror_directions = [[0.0, 0.0], [0.0, 0.0]] #dir1,dir2,dir3
   for k in 2:p_len
      if segment_length(path.points[k], hitpoints[1]) < 1e-4
         mirror_directions[1] = -path.directions[k-1]
      end
      if segment_length(path.points[k], hitpoints[2]) < 1e-4
         mirror_directions[2] = path.directions[k]
      end
   end
   #println(mirror_directions)
   #println(hitpoints[1], "\n", hitpoints[2])

   An1 = [0.0, 0.0]
   An2 = [0.0, 0.0]

   k = 0.01
   while true
      point_dir = hitpoints[1] - k .* mirror_directions[1]
      if point_in_block(temple, point_dir)
         An1 = point_dir + 0.01 .* mirror_directions[1]
         break
      end
      k += 0.01
      if k > 30
         println("		ADJUSTMENT UNSUCCESSFUL")
         return mirrors
      end
   end
   k = 0.01

   while true
      point_dir = hitpoints[2] - k .* mirror_directions[2]
      if point_in_block(temple, point_dir)
         An2 = point_dir + (tolerance + 0.01) .* mirror_directions[2]
         break
      end
      k += 0.01
      if k > 30
         println("		ADJUSTMENT UNSUCCESSFUL")
         return mirrors
      end
   end
   #println(An1, "\n", An2)

   px1 = An1
   px2 = An2
   k = 0.05


   while ⋅(px1 - hitpoints[1] + 0.5 .* mirror_directions[1], mirror_directions[1]) < 0

      px2 = An2
      step = 1
      print("		START")
      while ⋅(px2 - hitpoints[2] - 0.5 .* mirror_directions[2], mirror_directions[2]) < 0
         fi = line_angle(px1[1], px1[2], px2[1], px2[2])
         if !temple_segment_intersection(temple, (px1, sqrt(squared_distance(px1, px2)), fi))

            α1 = mirror_directions[1][2] > 0 ? (acos(mirror_directions[1][1]) + fi + pi) / 2 : (-acos(mirror_directions[1][1]) + fi + pi) / 2
            α2 = mirror_directions[2][2] > 0 ? (acos(mirror_directions[2][1]) + fi) / 2 : (-acos(mirror_directions[2][1]) + fi) / 2

            (test1, point1) = adjust_mirror_pars_nopath(temple, vcat(mirrors[1:i-1], mirrors[j:end]), px1, α1)
            (test2, point2) = adjust_mirror_pars_nopath(temple, vcat(mirrors[1:i-1], mirrors[j:end]), px2, α2)
            if test1 && test2
               mirror1 = (
                  v1=point1,
                  v2=point1 + mirror_length * [cos(α1), sin(α1)],
                  α=α1,
                  e=[cos(α1), sin(α1)],
                  n=[-sin(α1), cos(α1)],
               )
               mirror2 = (
                  v1=point2,
                  v2=point2 + mirror_length * [cos(α2), sin(α2)],
                  α=α2,
                  e=[cos(α2), sin(α2)],
                  n=[-sin(α2), cos(α2)],
               )
               path_try = raytrace(temple, lamp, vcat(mirrors[1:i-1], mirror1, mirror2, mirrors[j+1:end]))
               #print_parameters(lamp, vcat(mirrors[1:i-1], mirror1, mirror2, mirrors[j+1:end]))
               #println("$i, $j")
               score_try = evaluate_path(temple, path_try, reseval)
               #print("$score_try#")

               if length(path_try.points) == p_len && score_try > score_best
                  println("		FOUND BEST: $score_try")
                  mirrors = vcat(mirrors[1:i-1], mirror1, mirror2, mirrors[j+1:end])
                  score_best = score_try
               end
               print("$step|")
               step += 1
            end
         end

         px2 += k .* mirror_directions[2]

      end
      println()
      px1 += k .* mirror_directions[1]
   end

   println("		ADJUSTMENT SUCCESSFUL")
   return mirrors
end

#MATRIX CREATORS

#lamp_pos element ... [[i,j], fi1, t1]...
#matrix element ... [v, [fi1, t1], [fi2, t2],...]
function heatgrid_matrix(temple, resgrid, fi_ε, n, lamp_n)
   println("PARAMETERS: grid resolution = $(resgrid), dα = $(fi_ε)")
   dim_M = Int(round(18 / resgrid) - 1)
   println("Excpected runs: ≈$(dim_M^2)")
   matrix = [[0.0, ([0.0, 0] for k in 1:n)...] for i in 1:dim_M, j in 1:dim_M]
   lamp_pos = [[[0, 0], 0.0, 0.0] for k in 1:lamp_n]
   max_v = 0
   s = 0
   filename = "heatgrid_" * string(resgrid) * "_" * replace(string(fi_ε), "." => ",") * "_" * string(n) * ".jld2"
   filename2 = "lamp_pos" * string(resgrid) * "_" * replace(string(fi_ε), "." => ",") * "_" * string(lamp_n) * ".jld2"

   for i in 1:dim_M
      for j in 1:dim_M
         s += 1
         println("Point ($i,$j), step $s")
         if (point_in_block(temple, (i * resgrid + 1, j * resgrid + 1)))
            continue
         end
         for α in 0:fi_ε:2*pi
            ray = (
               v=[i * resgrid + 1, j * resgrid + 1],
               e=[cos(α), sin(α)],
            )
            t = temple_ray_intersection(temple, ray)
            matrix[i, j][1] += t^2
            for k in 2:n+1
               if matrix[i, j][k][2] < t && !any(abs(angle[1] - α) < 0.08 for angle in matrix[i, j][2:k-1])
                  matrix[i, j][k] = [α, t]
               end
            end
            if i < round(dim_M / 2) && j < round(dim_M / 2)
               for l in 1:lamp_n
                  if lamp_pos[l][3] < t && !any((abs(lamp[2] - α) < 0.08) && (abs(lamp[1][1] - lamp_pos[l][1][1]) < 0.3 / resgrid || abs(lamp[1][2] - lamp_pos[l][1][2]) < 0.3 / resgrid) for lamp in lamp_pos[1:l-1])
                     lamp_pos[l] = [[i, j], α, t]
                  end
               end
            end
         end
         if matrix[i, j][1] > max_v
            max_v = matrix[i, j][1]
         end

      end
      println("ROW $i")
      @save filename matrix
      @save filename2 lamp_pos
      println("FILES SAVED FOR ROW $i")
   end
   map(x -> x[1] = x[1] / max_v, matrix)

   @save filename matrix
   @save filename2 lamp_pos
end

#lamp ... [x,y,fi,t]
function lamp_array(heatgrid_matrix, resgrid, count, dim_M=calculate_dim_M(resgrid))
   lamps = [[0.0, 0.0, 0.0, 0.0] for i in 1:count]
   half_size = ceil(Int, dim_M / 2)

   # Iterate over the first half of the rows and columns
   for i in 1:half_size
      for j in 1:half_size
         direct_upgrade = false
         point = floatpoint_gridpoint([i, j], resgrid)
         for k in 1:count
            if segment_length(point, [lamps[k][1], lamps[k][2]]) < 0.2
               if lamps[k][4] < heatgrid_matrix[i, j][2][2]
                  println("ADDED: $(heatgrid_matrix[i, j][2]) at ($i,$j) or ($(point[1]),$(point[2]))")
                  lamps[k][1] = point[1]
                  lamps[k][2] = point[2]
                  lamps[k][3] = heatgrid_matrix[i, j][2][1]
                  lamps[k][4] = heatgrid_matrix[i, j][2][2]
                  direct_upgrade = true
                  break
               else
                  direct_upgrade = true
                  break
               end
            end
         end
         if !direct_upgrade
            for k in 1:count
               if lamps[k][4] < heatgrid_matrix[i, j][2][2]
                  println("ADDED BY INSERT: $(heatgrid_matrix[i, j][2]) at ($i,$j) or ($(point[1]),$(point[2]))")
                  lamp = [point[1], point[2], heatgrid_matrix[i, j][2][1], heatgrid_matrix[i, j][2][2]]
                  insert!(lamps, k, lamp)
                  pop!(lamps)
                  break
               end
            end
         end
      end
   end

   return lamps
end

#TEST
function testFNM()
   cmc24_solution = [
      15.44 18.740000000000002 4.838999999999999;
      17.475810147392444 1.0780181506460924 0.6948146928204123;
      1.0671363111842866 6.762041629966124 4.748814692820412;
      18.78713219158912 13.345305976255643 1.6106293856408245;
      1.7565613882269626 18.999859208290427 3.865314692820412
   ]
   heatgrid_matrix = JLD2.load("heatgrid_0.06_0,005_6.jld2", "matrix")
   println(grid_indexpoint(grid_floatpoint([5.2, 7.3], 0.06), 0.06))
   println(heatgrid_matrix[grid_indexpoint(grid_floatpoint([5.2, 7.3], 0.06), 0.06)])
   println(heatgrid_matrix[70, 105])
   return
   temple = load_temple(temple_string, block_size)
   (lamp, mirrors) = load_solution_partial(cmc24_solution, mirror_length)
   find_next_mirror(heatgrid_matrix, 0.06, 20, temple, lamp, mirrors, 2, 1, 0.01)
end

function testADJMP()
   cmc24_solution = [
      15.44 18.740000000000002 4.838999999999999;
      17.475810147392444 1.0780181506460924 0.6948146928204123;
      1.0671363111842866 6.762041629966124 4.748814692820412;
      18.78713219158912 13.345305976255643 1.6106293856408245;
      1.7565613882269626 18.999859208290427 3.865314692820412
   ]
   temple = load_temple(temple_string, block_size)
   (lamp, mirrors) = load_solution_partial(cmc24_solution, mirror_length)
   path = raytrace(temple, lamp, mirrors)
   println("A:$(path.points[end-1]), B:$(path.points[end])")
   println("A:$(path.points[3]), B:$(path.points[4])")
   println(adjust_mirror_pars(temple, mirrors, path, [3.5279734608714, 7.6206597041742], 2.360907346410207))
end

function testOL2()
   cmc24_solution = [
      15.44 18.740000000000002 4.838999999999999;
      17.475810147392444 1.0780181506460924 0.6948146928204123;
      1.0671363111842866 6.762041629966124 4.748814692820412;
      18.78713219158912 13.345305976255643 1.6106293856408245;
      1.7460800934631713 18.881423874225966 3.860629385640825;
      3.77763987076741 5.879828798044856 -0.6315559215387604;
      11.554585409275912 6.689378091895963 2.847444078461238;
      10.460900076929253 13.003116753485994 2.3563706143591756;
      4.5151086859101905 11.06811738232064 -0.42529715025711257
   ]
   reseval = 0.5
   temple = load_temple(temple_string, block_size)
   (lamp, mirrors) = load_solution_partial(cmc24_solution, mirror_length)
   mirrors = last_optimization(temple, lamp, mirrors, 0.5, reseval)
   print_parameters(lamp, mirrors)
end

function tesLA()
   heatgrid_matrix = JLD2.load("heatgrid_0.06_0,005_6.jld2", "matrix")
   lamps = lamp_array(heatgrid_matrix, 0.06, 30)
   i = 1
   #println(lamps)
   for lamp in lamps
      println("lamp $i.")
      println("	(x,y) = ($(lamp[1]),$(lamp[2]))\n	fi = $(lamp[3])\n	t = $(lamp[4])")
      i += 1
   end
end

function test_ind_pars()
   cmc24_solution = [
      18.9899966959408 2.0456944803577 3.1501342250355;
      10.0120712164148 1.9689771571414 3.1415652715598;
      2.4063041114551 2.1136692357806 -0.7280893154205;
      1.0226335567006 13.5984354661066 0.388350749844;
      18.0172580316316 12.8567267177671 5.931103847365;
      17.5776224809363 5.1281172810842 0.8746528405693;
      4.3596063555755 17.2638487428852 0.1522652160651;
      8.674412896651 7.3683022170594 -0.6903934384405
   ]
   heatgrid_matrix = JLD2.load("heatgrid_0.05_0,005_30.jld2", "matrix")
   temple = load_temple(temple_string, block_size)
   lamp, mirrors = load_solution_partial(cmc24_solution, mirror_length)
   calculator = score_w
   (mirrors, _) = next_mirror(heatgrid_matrix, 0.05, 2000, temple, lamp, mirrors, 0.05, 1.5, 0.008, calculator)
   println(mirrors)
   path = raytrace(temple, lamp, mirrors)
   println(path)
   score = evaluate_path(temple, path, 0.5)
   println(score)
   print_parameters(lamp, mirrors)
end

#MAINLOOP

#best_scores, [score_try, lamp, mirrors]
function find_3(start)
   heatgrid_matrix = JLD2.load("heatgrid_0.05_0,005_30.jld2", "matrix")
   lamps = JLD2.load("lamp_pos0.05_0,005_50.jld2", "lamp_pos")
   best_scores = []
   temple = load_temple(temple_string, block_size)
   lamp_i = 1
   filename = "three_mirrors5.jld2"

   for lamp_par in lamps[start:end]
      println("SEARCHING FOR LAMP $lamp_i CONFIG")
      lamp_par = (
         v=floatpoint_gridpoint(lamp_par[1], 0.05),
         α=lamp_par[2],
         e=[cos(lamp_par[2]), sin(lamp_par[2])],
      )
      mirrors = []
      c = 0
      while c < 2

         (mirrors, _) = next_mirror(heatgrid_matrix, 0.05, 2500, temple, lamp_par, mirrors, 0.06, 1.5, 0.006, score_w)

         c += 1
         path = raytrace(temple, lamp_par, mirrors)
         println("	SCORE AFTER $(c) MIRRORS IN RES 0.5: $(evaluate_path(temple, path, 0.5))")
         print_parameters(lamp_par, mirrors)
      end

      (mirrors, _) = next_mirror(heatgrid_matrix, 0.05, 2500, temple, lamp_par, mirrors, 0.06, 1.2, 0.003, score_w)


      path = raytrace(temple, lamp_par, mirrors)
      score_try = evaluate_path(temple, path, 0.5)
      println("	SCORE AFTER $(c) MIRRORS IN RES 0.5: $(score_try)")

      push!(best_scores, [score_try, lamp_par, mirrors])
      best_scores = sort(best_scores, by=x -> -x[1])


      lamp_i += 1

      JLD2.@save filename best_scores
   end
end

function recursion(matrix, resgrid, n, temple, lamp, mirrors, resstep, reseval, resrot, calculator, depth)

   if depth != length(mirrors)
      return
   end
   if depth == 7
      path = raytrace(temple, lamp, mirrors)
      println("	SCORE AFTER $(depth) MIRRORS IN RES 0.5: $(evaluate_path(temple, path, 0.5))")
      (mirrors, _) = next_mirror(matrix, resgrid, n, temple, lamp, mirrors, resstep / 2, reseval, resrot / 2, calculator)
      path = raytrace(temple, lamp, mirrors)
      println("	SCORE AFTER 8 MIRRORS IN RES 0.5: $(evaluate_path(temple, path, 0.5))")
      print_parameters(lamp, mirrors)
      try
         mirrors = last_optimization(temple, lamp, mirrors, 0, 0.2)
      catch e
         println("Error type: ", typeof(e))  # Prints the type of the error
         return
      end
      path = raytrace(temple, lamp, mirrors)
      score_try = evaluate_path(temple, path, 0.5)
      println("	FINAL SCORE IN RES 0.5: $(evaluate_path(temple, path, 0.5))")
      print_parameters(lamp, mirrors)
      filename = string(score_try) * "-" * string(rand(1:1000)) * ".jld2"
      lamp_mirrors = vcat(lamp, mirrors)
      JLD2.@save filename lamp_mirrors
      return
   end

   path = raytrace(temple, lamp, mirrors)
   println("	SCORE AFTER $(depth) MIRRORS IN RES 0.5: $(evaluate_path(temple, path, 0.5))")
   print_parameters(lamp, mirrors)

   println("DEPTH: $depth MIRROR1")
   (mirrors1, mirrors2) = next_mirror(matrix, resgrid, n, temple, lamp, mirrors, resstep, reseval, resrot, calculator)
   path1 = raytrace(temple, lamp, mirrors1)


   recursion(matrix, resgrid, n, temple, lamp, mirrors1, resstep, reseval, resrot, calculator, depth + 1)

   if isnothing(mirrors2)
      return
   end
   for mirror in mirrors2
      if isnothing(mirror)
         return
      end
   end
   path2 = raytrace(temple, lamp, mirrors2)
   if 0.95 * evaluate_path(temple, path1, 0.5) > evaluate_path(temple, path2, 0.5)
      return
   end
   println("DEPTH: $depth MIRROR2")
   recursion(matrix, resgrid, n, temple, lamp, mirrors2, resstep, reseval, resrot, calculator, depth + 1)
end

function find_results(start)
   best_scores = JLD2.load("three_mirrors5.jld2", "best_scores")
   println(best_scores)
   heatgrid_matrix = JLD2.load("heatgrid_0.05_0,005_30.jld2", "matrix")
   take = 1



   for score in best_scores[start:end]
      println("TAKE: $take")
      lamp = score[2]
      mirrors = score[3]
      c = 3

      temple = load_temple(temple_string, block_size)
      println("ENTERING RECURSION")
      recursion(heatgrid_matrix, 0.05, 2500, temple, lamp, mirrors, 0.07, 2.0, 0.007, score_w, c)

      take += 1
   end

end

function main()
   t = @timed begin
      find_results(1)
   end

   println()
   println("Execution time: ", t.time)
   println("Memory allocated: ", t.bytes * 1e-6, " MB")
end

main()
