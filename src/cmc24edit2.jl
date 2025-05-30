mirror_length = 0.5
function optimize_last_3(temple, lamp, mirrors)
   #if length(mirrors) > 2
   path = raytrace(temple, lamp, mirrors)
   score_best = evaluate_path(temple, path, 1)
   p_len = length(path)
   if p_len == 0
      return mirrors
   end
   mirror_points = [0, 0] #p1,p2
   mirror_directions = [0, 0] #dir1,dir2,dir3

   for p in 2:p_len-1
      i = 1
      for mirror in mirrors[end-2:end-1]
         if is_within_epsilon_of_segment(mirror.v1, mirror.v2, path.points[p])
            if mirror_points[i] == 0
               mirror_points[i] = path.points[p]
               mirror_directions[i] = path.directions[p-1]
            else
               return mirrors
            end
         end
         i += 1
      end
   end
   for mirror_point in mirror_points
      if mirror_point == 0
         return mirrors
      end
   end

   An1 = [0, 0]
   An2 = [0, 0]
   k = 0.05
   while true
      point_dir = mirror_points[1] + k .* mirror_directions[1]
      if point_in_block(temple, point_dir)
         An1 = point_dir - 0.05 .* mirror_directions[1]
         if k == 0.05
            return mirrors
         end
         break
      end
      k += k
   end

   k = 0.05
   while true
      point_dir = mirror_points[2] + k .* mirror_directions[2]
      if point_in_block(temple, point_dir)
         An2 = point_dir - 0.05 .* mirror_directions[2]
         if k == 0.05
            return mirrors
         end
         break
      end
      k += k
   end

   px1 = An1
   px2 = An2
   k = 0.05
   while ⋅(px1 - mirror_points[1], mirror_directions[1]) > 0
      while ⋅(px2 - mirror_points[2], mirror_directions[2]) > 0
         fi = line_angle(px1[1], px1[2], px2[1], px2[2])

         if !temple_segment_intersection(temple, (px1, sqrt(squared_distance(px1, px2)), fi))
            α1 = mirror_directions[1][2] > 0 ? (acos(mirror_directions[1][1]) + fi - pi) / 2 : (-acos(mirror_directions[1][1]) + pi + fi) / 2
            α2 = mirror_directions[2][2] > 0 ? (acos(mirror_directions[2][1]) + fi - pi) / 2 : (-acos(mirror_directions[2][1]) + pi + fi) / 2
            test1 = false
            test2 = false
            point1 = nothing
            point2 = nothing
            (test1, point1) = adjust_mirror_pars(px1, α1)
            (test2, point2) = adjust_mirror_pars(px2, α2)
            if test1 && test2
               mirror1 = (
                  v1=point1,
                  v2=point1 + mirror_length * [cos(α1), sin(α1)],
                  α=α1,
                  e=[cos(α1), sin(α1)],
                  n=[-sin(α1), cos(α1)]
               )
               mirror2 = (
                  v1=point2,
                  v2=point2 + mirror_length * [cos(α2), sin(α2)],
                  α=α2,
                  e=[cos(α2), sin(α2)],
                  n=[-sin(α2), cos(α2)]
               )
               path_try = raytrace(temple, lamp, vcat(mirrors[1:end-3], mirror1, mirror2, mirrors[end]))
               score_try = evaluate_path(temple, path, 1)
               if length(path_try.points) == p_len && score_try >= score_best
                  return vcat(mirrors[1:end-3], mirror1, mirror2, mirrors[end])
               end
            end
         end
         px2 -= k .* mirror_directions[2]
      end
      px1 -= k .* mirror_directions[1]
   end

   return mirrors
end

function adjust_mirror_pars(point, α)
   for k in -mirror_length+0.005:0.03:0
      if !temple_segment_intersection(temple, (point + k, mirror_length, α))
         return (true, point + k)
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

function line_angle(x1, y1, x2, y2)
   return atan(y2 - y1, x2 - x1)
end

function squared_distance(p1, p2)
   return (p2[1] - p1[1])^2 + (p2[2] - p1[2])^2
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

function ⋅(v, w)
   return v[1] * w[1] + v[2] * w[2]
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

function temple_segment_intersection(temple, segment)
   return any(segment_block_intersection(segment, block) for block ∈ temple.blocks)
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