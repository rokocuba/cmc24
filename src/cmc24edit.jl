# Computational Modeling Challenge 2024
# CMC24 solution evaluation script
# Author: Hrvoje Abraham, hrvoje.abraham@avl.com

# using FileIO
# using ImageIO
# using Measures
# using Plots;
# gr();
# using UUIDs
# using Colors
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

"""Last 12 digits of hexadecimal format of the input integer."""
function hex12(x::Integer)
	return last(string(x, base = 16), 12)
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
					v1 = v1,  # bottom left corner
					v2 = v2,
					v3 = v3,  # up right corner
					v4 = v4, s1 = (v1, block_size, 0),
					s2 = (v2, block_size, π / 2),
					s3 = (v3, block_size, π),
					s4 = (v4, block_size, 3π / 2),
				)

				push!(temple, block)
			end
		end
	end
	# display(temple)

	println(stderr, "The temple of size $temple_shape is loaded.")

	return (
		blocks = temple,
		shape = temple_shape,
		size = block_size .* temple_shape,
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
		v = cmc24_solution[1, 1:2],
		α = α,
		e = [cos(α), sin(α)],
	)

	# preprocess the mirrors
	mirrors = []
	for m ∈ 1:8
		α = cmc24_solution[m+1, 3]

		v = cmc24_solution[m+1, 1:2]
		e = [cos(α), sin(α)]
		n = [-sin(α), cos(α)]  # normal

		mirror = (
			v1 = v,
			v2 = v + mirror_length * e,
			α = α,
			e = e,
			n = n,
		)

		push!(mirrors, mirror)
	end

	println(stderr, "The solution is loaded.")

	return (lamp, mirrors)
end

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

function line_angle(x1, y1, x2, y2)
	return atan(y2 - y1, x2 - x1)
end

function segment_length(A, B)
	(x1, y1) = A
	(x2, y2) = B
	return sqrt((x2 - x1)^2 + (y2 - y1)^2)
end

function raytrace(temple, lamp, mirrors)
	local hit_mirror

	path = (
		points = [lamp.v],
		directions = [],
	)

	ray = (
		v = lamp.v,
		e = lamp.e,
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
				v = hitting_point,
				e = ray.e - 2 * (ray.e ⋅ hit_mirror.n) * hit_mirror.n,
			)
			continue
		end

		# ray hit the temple
		break
	end

	return path
end

function point_illuminated(path, point)
	points = path.points
	return any(is_within_one_unit_of_segment(A, B, point) for (A, B) in [(points[i], points[i+1]) for i in 1:length(points)-1])
end


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
	s = 0.98
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

function squared_distance(p1, p2)
	return (p2[1] - p1[1])^2 + (p2[2] - p1[2])^2
end

# Function to check if a point is within 1 unit from a segment
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

function finalize(temple = nothing, lamp = nothing, mirrors = nothing, path = nothing)
	if temple ≠ nothing
		cmc24_plot(temple, lamp = lamp, mirrors = mirrors, path = path)
	end

	print(0)  # stdout print of a fallback result in a case of early exit

	exit()
end

function tune_α(temple, lamp, mirrors, score_best, res, fi_ε)
	improv = false
	score_temp = 0

	for (i, value) in pairs(mirrors)
		k = 0

		mirror_inc_plus = (
			v1 = value.v1,
			v2 = value.v1 + mirror_length * [cos(value.α + fi_ε), sin(value.α + fi_ε)],
			α = value.α + fi_ε,
			e = [cos(value.α + fi_ε), sin(value.α + fi_ε)],
			n = [-sin(value.α + fi_ε), cos(value.α + fi_ε)],
		)

		if !temple_segment_intersection(temple, (value.v1, mirror_length, value.α + fi_ε))
			path = raytrace(temple, lamp, vcat(mirrors[1:i-1], mirror_inc_plus, mirrors[i+1:end]))
			score_temp = evaluate(temple, path, res)
			if score_best < score_temp
				score_best = score_temp
				k = 1
				improv = true
			end

		end
		mirror_inc_minus = (
			v1 = value.v1,
			v2 = value.v1 + mirror_length * [cos(value.α - fi_ε), sin(value.α - fi_ε)],
			α = value.α - fi_ε,
			e = [cos(value.α - fi_ε), sin(value.α - fi_ε)],
			n = [-sin(value.α - fi_ε), cos(value.α - fi_ε)],
		)

		if !temple_segment_intersection(temple, (value.v1, mirror_length, value.α - fi_ε))
			path = raytrace(temple, lamp, vcat(mirrors[1:i-1], mirror_inc_minus, mirrors[i+1:end]))
			score_temp = evaluate(temple, path, res)
			if score_best < score_temp
				score_best = score_temp
				k = -1
				improv = true
			end

		end

		if k != 0

			if k == 1
				mirrors[i] = mirror_inc_plus
			else
				mirrors[i] = mirror_inc_minus
			end

		end

	end
	return (improv, mirrors, score_best)
end

function finetune1d(temple, lamp, mirrors, mirror_index, score_best, resgrid, res, fi_ε, fi, walk)
	improv = false
	mirror_best = mirrors[mirror_index]

	directions = [(0, -1), (1, 0), (0, 1), (-1, 0)]

	(improv, mirror_cur, score_cur) = finetune_α1d(temple, lamp, mirrors, mirror_index, score_best, res, fi_ε, fi)
	if improv
		score_best = score_cur
		mirror_best = mirror_cur
		println("IMPROV!")
		println(score_best)
		println(mirrors[mirror_index].α - mirror_best.α)
	end

	mirror_try = mirror_cur

	for k in 1:2:walk*2
		println("K: $(k)")

		for direction in directions[1:2]
			println("DIRECTION: $(direction)")
			for j in 1:1:k
				mirror_try = (
					v1 = mirror_try.v1 + (resgrid) .* collect(direction),
					v2 = mirror_try.v1 + resgrid .* collect(direction) + mirror_length * mirror_try.e, α = mirror_try.α,
					e = mirror_try.e,
					n = mirror_try.n,
				)
				(improv, mirror_cur, score_cur) = finetune_α1d(temple, lamp,
					vcat(mirrors[1:mirror_index-1], mirror_try, mirrors[mirror_index+1:end]),
					mirror_index, score_best, res, fi_ε, fi)
				if improv
					score_best = score_cur
					mirror_best = mirror_cur
					println("IMPROV!")
					println(score_best)
					println(mirrors[mirror_index].α - mirror_best.α)
				end
			end
		end

		for direction in directions[3:4]
			println("DIRECTION: $(direction)")
			for j in 1:1:k+1
				mirror_try = (
					v1 = mirror_try.v1 + resgrid .* collect(direction),
					v2 = mirror_try.v1 + resgrid .* collect(direction) + mirror_length * mirror_try.e, α = mirror_try.α,
					e = mirror_try.e,
					n = mirror_try.n,
				)
				(improv, mirror_cur, score_cur) = finetune_α1d(temple, lamp,
					vcat(mirrors[1:mirror_index-1], mirror_try, mirrors[mirror_index+1:end]),
					mirror_index, score_best, res, fi_ε, fi)
				if improv
					score_best = score_cur
					mirror_best = mirror_cur
					println("IMPROV!")
					println(score_best)
					println(mirrors[mirror_index].α - mirror_best.α)
				end
			end
		end

	end

	return mirror_best
end

function find_path(temple, lamp, mirrors, resgrid, res, fi_ε)
	path = raytrace(temple, lamp, mirrors)
	if length(path.points) < 2
		println("PATH IS EMPTY")
		return
	end

	path_try = (
		points = path.points,
		directions = path.directions,
	)
	score_best = 0
	point = path.points[end]
	dir_down = path.directions[end]
	dir_up = -1 .* dir_down
	α_start = dir_up[2] > 0 ? acos(dir_up[1]) : -acos(dir_up[1]) + 2 * pi

	mirror_best = nothing

	steps = segment_length(point, path.points[end-1]) / resgrid
	println("NUMBER OF STEPS: $(steps)")
	for k in 1:floor(steps)
		print("|$(k)|")
		point_dir = point .+ (((1 / 2 + k / steps) * resgrid) .* dir_up)
		for α_0 in α_start-pi:fi_ε:α_start
			mirror_exists = false
			v1 = [0, 0]
			(mirror_exists, v1) = adjust_mirror_pars(temple, point_dir, α_0)
			if !mirror_exists
				continue
			end
			mirror_try = (
				v1 = v1,
				v2 = point_dir + mirror_length * [cos(α_0), sin(α_0)],
				α = α_0,
				e = [cos(α_0), sin(α_0)],
				n = [-sin(α_0), cos(α_0)],
			)

			path_try = raytrace(temple, lamp, vcat(mirrors, mirror_try))
			score_0 = evaluate_path(temple, path_try, res)
			if score_0 > score_best
				mirror_best = mirror_try
				score_best = score_0
			end

		end

	end
	println("")
	return mirror_best
end

function finetune_lamp1d(temple, lamp, mirrors, score_best, resgrid, res, fi_ε, fi, walk)
	improv = false
	directions = [(0, -1), (1, 0), (0, 1), (-1, 0)]
	lamp_try = lamp
	lamp_best = lamp

	(improv, lamp_try, score_cur) =
		finetune_lamp_α1d(temple, lamp, mirrors, score_best, res, fi_ε, fi)
	if improv
		score_best = score_cur
		lamp_best = lamp_try
		println("IMPROV!")
		println(score_best)
	end

	for k in 1:2:walk*2
		println("K: $(k)")

		for direction in directions[1:2]

			println("DIRECTION: $(direction)")
			for j in 1:1:k
				lamp_try = (
					v = lamp_try.v .+ (resgrid) .* collect(direction),
					α = lamp_try.α,
					e = lamp_try.e,
				)

				(improv, lamp_cur, score_cur) =
					finetune_lamp_α1d(temple, lamp_try, mirrors, score_best, res, fi_ε, fi)

				if improv
					score_best = score_cur
					lamp_best = lamp_cur
					println("IMPROV!")
					println(score_best)
				end
			end
		end

		for direction in directions[3:4]
			println("DIRECTION: $(direction)")
			for j in 1:1:k+1
				lamp_try = (
					v = lamp_try.v .+ (resgrid) .* collect(direction),
					α = lamp_try.α,
					e = lamp_try.e,
				)

				(improv, lamp_cur, score_cur) =
					finetune_lamp_α1d(temple, lamp_try, mirrors, score_best, res, fi_ε, fi)

				if improv
					score_best = score_cur
					lamp_best = lamp_cur
					println("IMPROV!")
					println(score_best)
				end
			end
		end

	end

	return lamp_best
end

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
		v1 = [x1_rot, y1_rot],
		v2 = [x2_rot, y2_rot],
		α = mirror.α + α,
		e = [cos(mirror.α + α), sin(mirror.α + α)],
		n = [-sin(mirror.α + α), cos(mirror.α + α)],
	)
	return mirror_rotated
end

function translate_mirror(mirror, len)
	mirror_translated = (
		v1 = mirror.v1 + len * mirror.e,
		v2 = mirror.v2 + len * mirror.e,
		α = mirror.α,
		e = mirror.e,
		n = mirror.n,
	)
	return mirror_translated
end

function mirrors_different(mirror1, mirror2)
	mirror1_s = (mirror1.v1 .+ mirror1.v2) / 2
	mirror2_s = (mirror2.v1 .+ mirror2.v2) / 2
	return squared_distance(mirror1_s, mirror2_s) > mirror_length_half || abs(mirror1.α - mirror2.α) > 0.1
end

function rotate_lamp(lamp, α)
	lamp_rotated = (
		v = lamp.v,
		α = lamp.α + α,
		e = [cos(lamp.α + α), sin(lamp.α + α)],
	)
	return lamp_rotated
end

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

function primitive_optimization(temple, path, lamp, mirrors)
	xy_ε = 0.01
	fi_ε = 0.01
	res = 1
	score_best = evaluate(temple, path, res)
	improv = true
	while improv
		(improv, mirrors, score_best) = tune_α(temple, lamp, mirrors, score_best, res, fi_ε)
	end

	pars = [(Float64(mirror.v1[1]), Float64(mirror.v1[2]), Float64(mirror.α)) for mirror in mirrors]
	return pars
end

function load_solution_partial(cmc24_solution, mirror_length)
	cmc24_solution = float(cmc24_solution)

	# preprocess the lamp
	α = cmc24_solution[1, 3]
	lamp = (
		v = cmc24_solution[1, 1:2],
		α = α,
		e = [cos(α), sin(α)],
	)

	# preprocess the mirrors
	mirrors = []
	for m ∈ 2:size(cmc24_solution, 1)
		α = cmc24_solution[m, 3]

		v = cmc24_solution[m, 1:2]
		e = [cos(α), sin(α)]
		n = [-sin(α), cos(α)]  # normal

		mirror = (
			v1 = v,
			v2 = v + mirror_length * e, α = α,
			e = e,
			n = n,
		)

		push!(mirrors, mirror)
	end

	println(stderr, "The solution is loaded.")

	return (lamp, mirrors)
end

function assign_new_pars(lamp, pars)
	parameters = []
	lamp_vec = [Float64(lamp.v[1]), Float64(lamp.v[2]), Float64(lamp.α)]
	push!(parameters, lamp_vec)
	pars_matrix = collect([collect(p) for p in pars])
	for p in pars_matrix
		push!(parameters, p)


	end
	parameters = hcat(parameters...)'
	return parameters
end

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
]
	=#
#=
parameters = [
	5 5 0.26;
	11.5 6.5 0.9;
	11.9 16.5 0.95;
	15.2 17.6 2.45;
	13.8 12.0 0.92;
	1.6 6.2 2.53;
	2.2 14.7 0.7;
	8.5 14.2 2.325;
	8.7 3.05 2.548]=#

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

function point_value_1m(temple, res, x, y, α)

	ray = (
		v = [Float64(x), Float64(y)],
		e = [cos(α), sin(α)],
	)



	path = (
		points = [ray.v],
		directions = [ray.e],
	)

	# check where ray would hit the temple
	t = temple_ray_intersection(temple, ray)
	hitting_point = ray.v + t * ray.e
	push!(path.points, hitting_point)

	return evaluate(temple, path, res)
end

function next_mirror_optimization(temple, lamp, mirrors, score_best, depth, resstep, res, res_fine, fi_ε, fi_ε_fine, fi_fine, tolerance)
	global best_score
	mirror_len = length(mirrors)

	if mirror_len != depth
		println("WARNING: NUMBER OF MIRRORS DOES NOT CORRESPOND TO DEPTH")
	end


	if depth == 7
		print_parameters(lamp, mirrors)
		mirror1 = find_path(temple, lamp, mirrors, resstep, res, fi_ε)
		mirrors1 = optimize_last_2(temple, lamp, vcat(mirrors, mirror1), 0.01)
		path = raytrace(temple, lamp, mirrors1)
		score1 = evaluate_path(temple, path, res_fine)
		(improv1, mirrors, score1_improv) = finetune_α1d(temple, lamp, mirrors1, mirror_len + 1, score1, res_fine, fi_ε_fine, fi_fine)
		println("FINETUNE: ", improv1 ? "SUCCESSFUL" : "UNSUCCESSFUL")


		println("SCORE: $(score1_improv), DEPTH: $(depth)")
		print_parameters(lamp, mirrors)

		if score1_improv > best_score
			global n
			best_score = score1_improv
			println("\n\nTHIS IS THE BEST SCORE!\n")
			filename = "BEST_SCORE_" * string(score1_improv) * "_" * string(tolerance) * "_" * string(res_fine) * ".jld2"
			@save filename mirrors
		end

		return
	end


	println("SCORE: $(score_best), DEPTH: $(depth)")
	print_parameters(lamp, mirrors)

	if score_best > best_score
		println("\n\nTHIS SHOULD NOT HAPPEN IF NOT IN BEGINNING (BEST SCORE WIHTOUT 8 MIRRORS)\n")
	end


	depth += 1
	mirror1 = nothing
	score1 = nothing
	mirrors1 = nothing
	try
		mirror1 = find_path(temple, lamp, mirrors, resstep, res, fi_ε)
		mirrors1 = optimize_last_2(temple, lamp, vcat(mirrors, mirror1), tolerance)
		path = raytrace(temple, lamp, mirrors1)
		score1 = evaluate_path(temple, path, res_fine)
	catch e
		println("Error type: ", typeof(e))  # Prints the type of the error
		# Detailed information
		println("Stack trace: ")
		println(stacktrace(e))
		return
	end
	(improv1, mirrors1, score1) = finetune_α1d(temple, lamp, mirrors1, mirror_len + 1, score1, res_fine, fi_ε_fine, fi_fine)
	println("FINETUNE: ", improv1 ? "SUCCESSFUL" : "UNSUCCESSFUL")

	next_mirror_optimization(temple, lamp, mirrors1, score1, depth, resstep, res, res_fine, fi_ε, fi_ε_fine, fi_fine, tolerance)


end

function optimize_last_2(temple, lamp, mirrors, tolerance)
	if length(mirrors) < 2
		return mirrors
	end
	path = raytrace(temple, lamp, mirrors)
	score_best = evaluate_path(temple, path, 1)
	p_len = length(path.points)
	if p_len == 0
		println("ADJUSTMENT UNSUCCESSFUL")
		return mirrors
	end
	mirror_points = [[0.0, 0.0], [0.0, 0.0]] #p1,p2
	mirror_directions = [[0.0, 0.0], [0.0, 0.0]] #dir1,dir2,dir3

	for p in 2:p_len-1
		i = 1
		for mirror in mirrors[end-1:end]
			if is_within_epsilon_of_segment(mirror.v1, mirror.v2, path.points[p])
				if mirror_points[i] == [0.0, 0.0]
					mirror_points[i] = path.points[p]
					mirror_directions[i] = path.directions[p+i-2]
				else
					println("ADJUSTMENT UNSUCCESSFUL")
					return mirrors
				end
			end
			i += 1
		end
	end
	for mirror_point in mirror_points
		if mirror_point == [0.0, 0.0]
			println("ADJUSTMENT UNSUCCESSFUL")
			return mirrors
		end
	end
	An1 = [0.0, 0.0]
	An2 = [0.0, 0.0]
	k = 0.01
	while true
		point_dir = mirror_points[1] + k .* mirror_directions[1]
		if point_in_block(temple, point_dir)
			An1 = point_dir - tolerance .* mirror_directions[1]
			break
		end
		k += 0.01
		if k > 30
			println("ADJUSTMENT UNSUCCESSFUL")
			return mirrors
		end
	end
	k = 0.01
	while true
		point_dir = mirror_points[2] - k .* mirror_directions[2]
		if point_in_block(temple, point_dir)
			An2 = point_dir - tolerance .* mirror_directions[2]
			break
		end
		k += 0.01
		if k > 30
			println("ADJUSTMENT UNSUCCESSFUL")
			return mirrors
		end
	end
	px1 = An1
	px2 = An2
	k = 0.01
	i = 0

	while ⋅(px1 - mirror_points[1] - tolerance .* mirror_directions[1], mirror_directions[1]) > 0
		px2 = An2
		while ⋅(px2 - mirror_points[2] - tolerance .* mirror_directions[2], mirror_directions[2]) < 0
			fi = line_angle(px1[1], px1[2], px2[1], px2[2])
			if !temple_segment_intersection(temple, (px1, sqrt(squared_distance(px1, px2)), fi))
				α1 = mirror_directions[1][2] > 0 ? (acos(mirror_directions[1][1]) + fi) / 2 : (-acos(mirror_directions[1][1]) + fi) / 2
				α2 = mirror_directions[2][2] > 0 ? (acos(mirror_directions[2][1]) + fi) / 2 : (-acos(mirror_directions[2][1]) + fi) / 2
				test1 = false
				test2 = false
				point1 = nothing
				point2 = nothing
				(test1, point1) = adjust_mirror_pars(temple, px1, α1)
				(test2, point2) = adjust_mirror_pars(temple, px2, α2)
				if test1 && test2
					mirror1 = (
						v1 = point1,
						v2 = point1 + mirror_length * [cos(α1), sin(α1)],
						α = α1,
						e = [cos(α1), sin(α1)],
						n = [-sin(α1), cos(α1)],
					)
					mirror2 = (
						v1 = point2,
						v2 = point2 + mirror_length * [cos(α2), sin(α2)],
						α = α2,
						e = [cos(α2), sin(α2)],
						n = [-sin(α2), cos(α2)],
					)
					path_try = raytrace(temple, lamp, vcat(mirrors[1:end-2], mirror1, mirror2))
					score_try = evaluate_path(temple, path, 1)

					if length(path_try.points) == p_len && score_try >= score_best
						mirrors = vcat(mirrors[1:end-2], mirror1, mirror2)
						score_best = score_try
					end
				end
			end
			i += 1
			px2 += k .* mirror_directions[2]
		end
		px1 -= k .* mirror_directions[1]
		if i > 2e5
			println("ADJUSTMENT MAYBE SUCCESSFUL")
			return mirrors
		end
	end
	println("ADJUSTMENT SUCCESSFUL")
	return mirrors
end

function adjust_mirror_pars(temple, point, α)
	for k in -mirror_length+0.005:0.03:0
		if !temple_segment_intersection(temple, (point + k .* [cos(α), sin(α)], mirror_length, α))
			return (true, point + k .* [cos(α), sin(α)])
		end
	end
	return (false, [0, 0])
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

function print_parameters(lamp, mirrors)
	print(lamp.v[1], " ", lamp.v[2], " ", lamp.α)
	for mirror in mirrors
		println(";")
		print(mirror.v1[1], " ", mirror.v1[2], " ", mirror.α)
	end
	println()
end

function main()

	parameters = [5.0 5.0 0.26; 11.5 6.5 0.9; 11.9 16.5 0.95; 15.2 17.6 2.45; 13.8 12.0 0.92; 1.6 6.2 2.53; 2.2 14.7 0.7; 8.5 14.2 2.325; 8.7 3.05 2.5480000000000003]
	cmc24_solution = parameters

	res = 0.3


	# load the temple
	temple = load_temple(temple_string, block_size)

	# load the solution
	lamp, mirrors = load_solution(cmc24_solution, mirror_length)
	check_solution(temple, lamp, mirrors)

	# compute the ray path
	path = raytrace(temple, lamp, mirrors)
	score = evaluate(temple, path, res)
	println(score)
	pars = primitive_optimization(temple, path, lamp, mirrors)
	parameters = assign_new_pars(lamp, pars)
	lamp, mirrors = load_solution(parameters, mirror_length)

	check_solution(temple, lamp, mirrors)
	path = raytrace(temple, lamp, mirrors)
	println("PARAMETERS Final: ", parameters, "\n")
	# evaluate the solution
	score = evaluate(temple, path, res)
	# println(stderr, "Your CMC24 score is $(commas(score))")
	println(score)

	# create the presentation plot
	# cmc24_plot(temple, lamp=lamp, mirrors=mirrors, path=path)

end

function test()
	mirror1 = (
		v1 = [0, 0],
		v2 = [0, 0] + mirror_length * [cos(pi / 4), sin(pi / 4)],
		α = pi / 4,
		e = [cos(pi / 4), sin(pi / 4)],
		n = [-sin(pi / 4), cos(pi / 4)],
	)
	mirrors = [mirror1]
	mirror_try = mirrors[1]
	resgrid = 0.1
	directions = [(0, -1), (1, 0), (0, 1), (-1, 0)]
	println("MIRROR: $(mirror_try)")
	for k in 1:2:10

		for direction in directions[1:2]
			for j in 1:1:k
				mirror_try = (
					v1 = mirror_try.v1 + resgrid .* collect(direction),
					v2 = mirror_try.v1 + resgrid .* collect(direction) + mirror_length * mirror_try.e,
					α = mirror_try.α,
					e = mirror_try.e,
					n = mirror_try.n,
				)
				println("(j,k): $(j),$(k), DIRECTION: $(direction), MIRROR: $(mirror_try)")
			end
		end

		for direction in directions[3:4]
			for j in 1:1:k+1
				mirror_try = (
					v1 = mirror_try.v1 + resgrid .* collect(direction),
					v2 = mirror_try.v1 + resgrid .* collect(direction) + mirror_length * mirror_try.e,
					α = mirror_try.α,
					e = mirror_try.e,
					n = mirror_try.n,
				)
				println("(j,k): $(j),$(k), DIRECTION: $(direction), MIRROR: $(mirror_try)")
			end
		end

	end
end

function test2()
	temple = load_temple(temple_string, block_size)
	cmc24_solution = [
		15.45 18.91 4.84;
		17.8442373553 1.5260603781 -2.43434867321
	]
	lamp, mirrors = load_solution(cmc24_solution, mirror_length)
	println("LAMP: $(lamp)")
	println("MIRRORS: $(mirrors)")
	path = raytrace(temple, lamp, mirrors)
	score = evaluate(temple, path, 2)

end

function test3()
	c = 0
	k = 0.98
	res = 0.3
	println("0!!")
	seg = [[1.0, 3.0], [4.0, 6.0]]
	float_steps_w = 1 / res
	println("1!!")
	fi = line_angle(seg[1][1], seg[1][2], seg[2][1], seg[2][2])
	println("2!!")
	e = [cos(fi), sin(fi)]
	n = [-sin(fi), cos(fi)]
	seg[1] -= e
	seg[2] += e
	println("SEG:$(seg)")
	println("fi:$(fi)")
	gen_values = collect(res * k for k in -floor(float_steps_w):1:floor(float_steps_w))
	println("gen_values:$(gen_values)")
	starting_points = [(A, B) for (A, B) in [seg[1] + n .* value for value in gen_values]]

	for point in starting_points
		println("STARTING POINT$(c):$(point)")
		c += 1
	end
end

function test4()
	cmc24_solution = [
		5 5 0.26;
		11.5 6.5 0.9;
		11.9 16.5 0.95;
		15.2 17.6 2.45;
		13.8 12.0 0.92;
		1.6 6.2 2.53;
		2.2 14.7 0.7;
		8.5 14.2 2.325;
		8.7 3.05 2.348
	]

	temple = load_temple(temple_string, block_size)
	lamp, mirrors = load_solution_partial(cmc24_solution, mirror_length)
	path = raytrace(temple, lamp, mirrors)
	println(path)
	score = evaluate_path(temple, path, 0.1)
	println(score)
end

function heatgrid_matrix(temple, resgrid, fi_ε, n)
	println("PARAMETERS: grid resolution = $(resgrid), dα = $(fi_ε)")
	dim_M = Int(round(18 / resgrid) - 1)
	println("Excpected runs: ≈$(dim_M^2)")
	matrix = [[0.0, ([0.0, 0] for k in 1:n)...] for i in 1:dim_M, j in 1:dim_M]
	max_v = 0
	s = 0

	for i in 1:dim_M
		for j in 1:dim_M
			s += 1
			println("Point ($i,$j), step $s")
			if (point_in_block(temple, (i * resgrid + 1, j * resgrid + 1)))
				continue
			end
			for α in 0:fi_ε:2*pi
				ray = (
					v = [i * resgrid + 1, j * resgrid + 1],
					e = [cos(α), sin(α)],
				)
				t = temple_ray_intersection(temple, ray)
				matrix[i, j][1] += t^2
				for k in 2:n+1
					if matrix[i, j][k][2] < t && !any(abs(angle[1] - α) < 0.08 for angle in matrix[i, j][2:k-1])
						matrix[i, j][k] = [α, t]
					end
				end
			end
			if matrix[i, j][1] > max_v
				max_v = matrix[i, j][1]
			end

		end
	end
	map(x -> x[1] = x[1] / max_v, matrix)
	filename = "heatgrid_" * string(resgrid) * "_" * replace(string(fi_ε), "." => ",") * "_" * string(n) * ".jld2"
	@save filename matrix
end

function heatgrid_matrix_path(temple, resgrid, res, fi_ε)
	max_d = 0
	heatgrid_matrix = []
	st_coo = 1 + min(0.3, res / 3)
	p = (st_coo, st_coo)
	pass_count = 0
	println("PARAMETERS: grid resolution = $(resgrid), score_eval resolution = $(res), dα = $(fi_ε)")
	println("Excpected runs: ≈$((18/resgrid)^2)")
	for y ∈ p[2]:resgrid:(temple.size[2]-1)
		row = []
		for x ∈ p[1]:resgrid:(temple.size[1]-1)
			point = [Float64(x), Float64(y)]
			score = 0
			if point_in_block(temple, point)
				push!(point, score)
				push!(row, point)
				pass_count += 1
				print("$(pass_count)|")
				continue
			end
			for α in 0:fi_ε:2*pi
				lamp = (
					v = point,
					α = α,
					e = [cos(α), sin(α)],
				)
				path = raytrace(temple, lamp, [])
				score += evaluate_path(temple, path, res)^2
			end
			push!(point, score)
			if score > max_d
				max_d = score
			end
			push!(row, point)
			pass_count += 1
			print("$(pass_count)|")
		end
		push!(heatgrid_matrix, row)
	end

	heatgrid_matrix = hcat(heatgrid_matrix...)

	for row in eachrow(heatgrid_matrix)
		for element in row
			element[3] = element[3] / max_d
		end
	end


	@save "heatgrid_$(resgrid)_$(res)_$(fi_ε)pth" heatgrid_matrix
end

function test5()
	lamp_x = 15.44
	lamp_y = 18.740000000000002
	lamp_α = 4.838999999999999

	lamp = (
		v = [lamp_x, lamp_y],
		α = lamp_α,
		e = [cos(lamp_α), sin(lamp_α)],
	)
	temple = load_temple(temple_string, block_size)
	mirror_best1 = (
		v1 = [17.469197581800994, 1.1363885116885177],
		v2 = [17.85135612879888, 1.458808804096521],
		α = 0.7008146928204122,
		e = [0.7643170939957777, 0.6448405848160066],
		n = [-0.6448405848160066, 0.7643170939957777],
	)
	mirror_best2 = (
		v1 = [17.398355327513293, 1.7259268482117383],
		v2 = [17.77067075811628, 2.059665099143377],
		α = 0.7308146928204124,
		e = [0.7446308612059817, 0.6674765018632777],
		n = [-0.6674765018632777, 0.7446308612059817],
	)
	mirror = nothing
	score_mirrors = [[107, mirror_best1], [104, mirror_best2]]
	i = 1
	for score_mirror in score_mirrors
		improv = false
		println(score_mirror[2], "\n")
		path = raytrace(temple, lamp, [score_mirror[2]])
		score = evaluate_path(temple, path, 0.8)
		println("OLD SCORE ($(i))", score)
		(improv, mirror, score) = finetune_α1d(temple, lamp, [score_mirror[2]], 1, score_mirror[1], 0.05, 0.001, 0.004)
		if improv
			score_mirrors[i][1] = score
			score_mirrors[i][2] = mirror
		end
		println("IMPROVED: ", improv)
		println("SCORE_BEST: $(score_mirrors[i][1])")
		println("MIRROR_BEST: $(mirror)")
		i += 1
	end
	println(score_mirrors)
end

function test6()
	cmc24_solution = [
		15.44 18.740000000000002 4.838999999999999;
		17.475810147392444 1.0780181506460924 0.6948146928204123;
		1.0671363111842866 6.762041629966124 4.748814692820412;
		18.78713219158912 13.345305976255643 1.6106293856408245;
		1.7460800934631713 18.881423874225966 3.860629385640825;
		3.77763987076741 5.879828798044856 -0.6315559215387604;
		11.554585409275912 6.689378091895963 2.847444078461238
	]
	# load the solution
	temple = load_temple(temple_string, block_size)
	lamp, mirrors = load_solution_partial(cmc24_solution, mirror_length)
end

function write_mir()
	lamp_x = 15.44
	lamp_y = 18.740000000000002
	lamp_α = 4.838999999999999

	lamp = (
		v = [lamp_x, lamp_y],
		α = lamp_α,
		e = [cos(lamp_α), sin(lamp_α)],
	)

	mirrors = JLD2.load("BEST_SCORE_65770_0000110_0.05.jld2", "mirrors")
	print_parameters(lamp, mirrors)
end

function find_mirror_pars_on_path(A, B, resgrid, i, j, α)
	(x1, y1) = A
	(x2, y2) = B
	A_B_angle = line_angle(x1, y1, x2, y2)
	A_B_distance = segment_length(A, B)
	ray = (
		v = [i * resgrid + 1, j * resgrid + 1],
		e = [-cos(α), -sin(α)],
	)

	(case, t, u) = ray_segment_intersection(ray, (A, A_B_distance, A_B_angle))
end

function test_search(lamp, mirrors, depth, score)
	resstep = 0.08 #0.05
	res = 1.2 #0.8
	res_fine = 0.08 #0.05
	fi_ε = 0.01 #0.008
	fi_ε_fine = 0.002 #0.001
	fi_fine = 0.006 #0.004
	tolerance_initial = 0.3


	temple = load_temple(temple_string, block_size)
	for i in 0:0.05:1.5
		try
			next_mirror_optimization(temple, lamp, mirrors, score, depth, resstep, res, res_fine, fi_ε, fi_ε_fine, fi_fine, tolerance_initial + i)
			println(best_score)
		catch e
			println("Error type: ", typeof(e))  # Prints the type of the error
			# Detailed information
			println("Stack trace: ")
			println(stacktrace(e))
			return
		end

	end

end
n = 0
best_score = 0
t = @timed begin
	temple = load_temple(temple_string, block_size)
	heatgrid_matrix(temple, 0.06, 0.005, 6)
end

println()
println("Execution time: ", t.time)
println("Memory allocated: ", t.bytes * 1e-6, " MB")
#julia --track-allocation=user cmc24edit.jl
