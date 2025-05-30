using JLD2 # Za spremanje/učitavanje Julia objekata (npr. matrica, podataka)

# Funkcija: find_grid_points_of_segment
# Svrha: Pronalazi 'n' najboljih točaka na gridu duž segmenta definiranog točkama A i B.
#        "Najbolje" se vjerojatno odnosi na neki 'weighted_score' koji se računa.
# Argumenti:
#   matrix: Vjerojatno matrica koja sadrži neke predizračunate vrijednosti za svaku točku grida.
#           (npr. matrix[gx,gy][1] je neki osnovni score, matrix[gx,gy][2:end] su 'dirs' - smjerovi/kutovi?)
#   A, B: Krajnje točke segmenta (parovi koordinata).
#   resgrid: Rezolucija grida (veličina jedne ćelije grida).
#   n: Broj najboljih točaka koje treba pronaći i vratiti.
function find_grid_points_of_segment(matrix, A, B, resgrid, n)

	# dim_M: Izračunava dimenziju matrice (grida) na osnovu fiksne veličine (18) i rezolucije.
	# Pretpostavka: Hram ili područje interesa je veličine 18x18 jedinica.
	dim_M = Int(round(18 / resgrid) - 1)
	println("M_2($(dim_M ))") # Ispis dimenzije (za debug)

	# Zakomentirani kod: Inicijalizacija matrice - možda za testiranje ili ako 'matrix' nije dana.
	#matrix = [[0, ([0.0, 0] for k in 1:n)...] for i in 1:dim_M, j in 1:dim_M] #REMOVE THIS IN FINAL
	#for k in eachindex(matrix)
	#   println(k, "\n")
	#end

	# points: Lista koja će čuvati 'n' najboljih pronađenih točaka.
	# Svaki element je oblika: [weighted_score, dirs_info, [grid_x, grid_y]]
	points = [[0, [0.0, 0.0], [0, 0]] for _ in 1:n] # Inicijalizacija s nultim vrijednostima

	println(typeof(points)) # Ispis tipa (za debug)
	println(points[1][3][1]) # Ispis koordinate prve točke (za debug)

	# Određivanje početne točke na gridu blizu točke A
	(x1, y1) = A
	# pointx, pointy: koordinate na gridu najbliže A, poravnate s resgrid
	(pointx, pointy) = ((round((x1 - 1) / resgrid) * resgrid + 1), (round((y1 - 1) / resgrid) * resgrid + 1))
	println("x:$(pointx), y:$(pointy)") # Debug ispis
	# gridpointx, gridpointy: indeksi ćelije grida za pointx, pointy
	(gridpointx, gridpointy) = (Int(round((pointx - 1) / resgrid)), Int(round((pointy - 1) / resgrid)))
	# Osiguravanje da su indeksi unutar granica matrice (grida)
	gridpointx = gridpointx == 0 ? 1 : (gridpointx == dim_M + 1 ? dim_M : gridpointx)
	gridpointy = gridpointy == 0 ? 1 : (gridpointy == dim_M + 1 ? dim_M : gridpointy)

	println("i:$(gridpointx), j:$(gridpointy)") # Debug ispis indeksa


	# Određivanje krajnje točke na gridu blizu točke B (slično kao za A)
	(x2, y2) = B
	(endpointx, endpointy) = ((round((x2 - 1) / resgrid) * resgrid + 1), (round((y2 - 1) / resgrid) * resgrid + 1))
	println("x:$(endpointx), y:$(endpointy)")
	(endgridpointx, endgridpointy) = (Int(round((endpointx - 1) / resgrid)), Int(round((endpointy - 1) / resgrid)))
	endgridpointx = endgridpointx == 0 ? 1 : (endgridpointx == dim_M + 1 ? dim_M : endgridpointx)
	endgridpointy = endgridpointy == 0 ? 1 : (endgridpointy == dim_M + 1 ? dim_M : endgridpointy)

	println("i:$(endgridpointx), j:$(endgridpointy)")

	# angle_AB: Kut segmenta AB. Pretpostavka: funkcija line_angle(x1,y1,x2,y2) postoji negdje drugdje.
	angle_AB = line_angle(x1, y1, x2, y2)

	# direction: Normalizirani vektor smjera od A do B.
	direction_vector = collect(B .- A) ./ segment_length(A, B) # segment_length(A,B) također mora postojati
	println("DIRECTION_fi: ", direction_vector) # Debug ispis

	# direction_step: Određuje smjer koraka po gridu (-1 ili 1 za x i y os)
	direction_step = (x -> x < 0 ? -1 : 1).(direction_vector)
	direction_step = Int.(direction_step)
	println("DIRECTION: ", direction_step)

	#c = 0 # Brojač koraka (zakomentirano)

	# Glavna petlja: Iterira po ćelijama grida od početne do krajnje točke segmenta.
	# Koristi se neka vrsta DDA (Digital Differential Analyzer) ili sličnog algoritma za praćenje linije po gridu.
	while gridpointx != endgridpointx && gridpointy != endgridpointy
		# Trenutne koordinate centra ćelije grida
		current_cell_center_x = gridpointx * resgrid + 1
		current_cell_center_y = gridpointy * resgrid + 1
		
		# Iteriraj kroz 'dirs' (smjerove/kutove?) definirane u 'matrix' za trenutnu ćeliju grida.
		# matrix[gridpointx, gridpointy][1] je osnovni score ćelije.
		# matrix[gridpointx, gridpointy][2:end] je lista 'dirs' objekata.
		# Svaki 'dirs' je vjerojatno [kut, neki_faktor_tezine_za_taj_kut]
		for current_dirs_info in matrix[gridpointx, gridpointy][2:end]
			# Ako je kut (current_dirs_info[1]) preblizu kutu segmenta AB (unutar 0.08 radijana, tj. ~4.5 stupnjeva),
			# preskoči ovaj 'dirs'. Vjerojatno da se izbjegnu refleksije u smjeru dolaska.
			if abs(rem(current_dirs_info[1] - angle_AB, pi)) < 0.08
				continue
			end
			# weighted_score: Računa se na osnovu osnovnog scora ćelije, faktora iz 'dirs' i udaljenosti od A.
			# segment_length(A, [trenutna_tocka_centra_celije])
			weighted_score = matrix[gridpointx, gridpointy][1] * current_dirs_info[2] + segment_length(A, [current_cell_center_x, current_cell_center_y])
			
			#println("\n\n $points \n\n") # Debug ispis liste 'points'
			
			# direct_upgrade: Funkcija koja provjerava treba li ažurirati neku postojeću točku u 'points'
			# ako je nova (gridpointx, gridpointy) s ovim 'dirs' dovoljno blizu i bolja.
			# Vraća: (je_li_pronadjena_i_mozda_azurirana_slicna_tocka, treba_li_sortirati_listu_points)
			(same_point_found, req_sort) = direct_upgrade(gridpointx, gridpointy, current_dirs_info, points, n, resgrid, weighted_score)
			
			if req_sort # Ako je direct_upgrade vratio da treba sortirati
				println("UNSORTED: $points \n\n") # Debug
				points = sort(points, by = x -> -x[1]) # Sortiraj 'points' opadajuće po weighted_score
				println("SORTED: $points \n\n") # Debug
			end
			
			if same_point_found # Ako je direct_upgrade već obradio ovu točku (našla se slična i ažurirana je ili odbačena)
				continue # Preskoči dodavanje nove točke, jer je postojeća ažurirana ili je ova lošija
			end
			
			# Ako nije pronađena slična točka za ažuriranje, pokušaj umetnuti novu točku u listu 'points'
			# ako je njen weighted_score dovoljno dobar.
			for i in 1:n # Iteriraj kroz trenutno najboljih 'n' točaka
				if points[i][1] < weighted_score # Ako je nova točka bolja od i-te u listi
					insert!(points, i, [weighted_score, current_dirs_info, [gridpointx, gridpointy]]) # Umetni novu
					pop!(points) # Izbaci najgoru (sada (n+1)-tu) točku da lista ostane duljine 'n'
					#println("$gridpoints,\n $gridpointx,\n $gridpointy") # Zakomentirani debug
					#println(matrix[gridpointx, gridpointy])
					break # Prekini petlju jer je točka umetnuta
				end
			end
		end # Kraj petlje po 'dirs'
		#c += 1 # Povećaj brojač koraka (zakomentirano)

		# Određivanje sljedeće ćelije grida za provjeru.
		# Računa udaljenost od idealne linije AB do susjednih ćelija u x i y smjeru.
		# distance_from_segment(A, B, tocka_za_provjeru) - pretpostavka da ova funkcija postoji.
		d1 = distance_from_segment(A, B, [current_cell_center_x + resgrid * direction_step[1], current_cell_center_y])
		d2 = distance_from_segment(A, B, [current_cell_center_x, current_cell_center_y + resgrid * direction_step[2]])
		# Pomiče se u smjeru (x ili y) koji je bliži idealnoj liniji.
		d1 < d2 ? gridpointx += direction_step[1] : gridpointy += direction_step[2]
		#println("STEP: $c (GRIDX,GRIDY) = ($(gridpointx),$(gridpointy)) || ($(d1),$(d2))") # Debug
	end # Kraj while petlje (kretanje po segmentu)
	
	println(length(points)) # Ispis konačnog broja točaka (trebalo bi biti 'n')
	return (points) # Vraća listu 'n' najboljih točaka
end

# Funkcija: direct_upgrade
# Svrha: Provjerava postoji li u listi 'points' već točka koja je "ista" ili "vrlo blizu"
#        novoj kandidat točki (gridpointx, gridpointy) s danim 'dirs'.
#        Ako postoji i nova je bolja, ažurira postojeću. Ako postoji i nova je lošija, ne radi ništa.
# Argumenti:
#   gridpointx, gridpointy: Koordinate nove kandidat točke.
#   dirs: Informacije o smjeru/kutu za novu kandidat točku.
#   points: Trenutna lista najboljih točaka.
#   n: Duljina liste 'points'.
#   resgrid: Rezolucija grida.
#   weighted_score: Izračunati score za novu kandidat točku.
# Vraća: tuple (found_and_updated_or_discarded::Bool, needs_sort::Bool)
function direct_upgrade(gridpointx, gridpointy, dirs, points, n, resgrid, weighted_score)
	for k in 1:n # Iteriraj kroz listu 'points'
		# Uvjet sličnosti: 
		# 1. Prostorna blizina: Apsolutna razlika koordinata manja od (0.15 / resgrid).
		#    (0.15 / resgrid) je neki prag blizine u jedinicama indeksa grida.
		#    Npr. ako je resgrid 0.05, prag je 3 ćelije. Ako je resgrid 0.5, prag je 0.3 ćelije.
		#    Ovo (0.15 / resgrid) djeluje malo čudno, možda bi trebalo biti fiksni broj ćelija, npr. < 2.
		# 2. Sličnost kuta: Razlika kutova (points[k][2][1] i dirs[1]) manja od 0.08 radijana.
		if (abs(points[k][3][1] - gridpointx) < (0.15 / resgrid) || abs(points[k][3][2] - gridpointy) < (0.15 / resgrid)) && abs(rem((points[k][2][1] - dirs[1]), pi)) < 0.08
			if weighted_score > points[k][1] # Ako je nova točka bolja od postojeće slične
				points[k] = [weighted_score, dirs, [gridpointx, gridpointy]] # Ažuriraj postojeću
				return (true, true) # Pronađena, ažurirana, treba sortirati listu 'points'
			end
			return (true, false) # Pronađena, ali nova nije bolja, ne treba sortirati
		end
	end
	return (false, false) # Nije pronađena slična točka
end

# Zakomentirana funkcija: any_and_delete
# Izgleda kao pokušaj implementacije slične logike kao direct_upgrade, 
# ali s drugačijim pristupom (možda brisanje i umetanje umjesto direktnog ažuriranja).
# Sadrži dosta debug ispisa.
#=function any_and_delete(gridpointx, gridpointy, dirs, gridpoints, points, i, n, resgrid)
   if gridpointx == 6
	  println("   $gridpointx|$gridpointy|$i|$(0.15/resgrid)")
	  println("   $gridpoints\n")
	  println("   $points\n")
   end
   for k in 1:i-1
	  if (abs(gridpoints[k][1] - gridpointx) < (0.15 / resgrid) || abs(gridpoints[k][2] - gridpointy) < (0.15 / resgrid))
		 #&&
		 #abs(rem((points[k][2][1] - dirs[1]), pi)) < 0.08
		 return true  # Return true as soon as a match is found
	  end
   end
   for k in i+1:n
	  if (abs(gridpoints[k][1] - gridpointx) < (0.15 / resgrid) || abs(gridpoints[k][2] - gridpointy) < (0.15 / resgrid))
		 #&&
		 #abs(rem((points[k][2][1] - dirs[1]), pi)) < 0.08
		 for j in k:(length(points)-1)
			points[j] = points[j+1]
			gridpoints[j] = gridpoints[j+1]
		 end
		 return false
	  end

   end
   return false  # Return false if no matches were found
end=#

# Funkcija: squared_distance
# Svrha: Računa kvadrat udaljenosti između dvije točke p1 i p2.
#        Korisno za izbjegavanje sqrt operacije ako se samo uspoređuju udaljenosti.
function squared_distance(p1, p2)
	return (p2[1] - p1[1])^2 + (p2[2] - p1[2])^2
end

# Funkcija: distance_from_segment
# Svrha: Računa (najkraću) udaljenost točke P od segmenta linije AB.
#        Ova implementacija izgleda da računa udaljenost od *linije* na kojoj leži AB,
#        ili možda udaljenost od točke A ako su A i B iste. Treba provjeriti logiku.
#        Standardni algoritam za udaljenost točke od segmenta je kompleksniji i uključuje
#        provjeru je li projekcija točke P na liniju AB unutar segmenta AB.
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

	if AB_squared == 0 # A i B su ista točka
		# Udaljenost od P do A (ili B)
		# Originalno je vraćalo 0, što nije točno osim ako P==A.
		# Trebalo bi vratiti sqrt(squared_distance(A,P)) ili samo squared_distance(A,P) ako se uspoređuju kvadrati.
		# Za sada, ostavljam kako je bilo, ali ovo je potencijalni bug ako A==B.
		return 0 # Ovo je udaljenost P od segmenta A-A, što je udaljenost P od točke A.
		         # Ako je P == A, onda je 0. Inače nije. Vjerojatno bi trebalo biti sqrt(APx^2 + APy^2).
	end

	# Ostatak funkcije (nakon if AB_squared == 0) nedostaje u priloženom kodu.
	# Ovdje bi trebala ići logika za izračun udaljenosti točke od segmenta.
	# Npr. koristeći projekciju točke P na liniju AB.
	# t = ((x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1)) / AB_squared;
	# if t < 0.0: closest_point = A (x1,y1)
	# elif t > 1.0: closest_point = B (x2,y2)
	# else: closest_point = (x1 + t * ABx, y1 + t * ABy)
	# return sqrt(squared_distance(P, closest_point));
end

# Pretpostavljene funkcije koje se koriste, a nisu definirane u ovom fajlu:
# function line_angle(x1, y1, x2, y2) -> Float64
#   # Vraća kut linije od (x1,y1) do (x2,y2) u radijanima.
#   # return atan(y2-y1, x2-x1)
# end

# function segment_length(P1, P2) -> Float64
#   # Vraća duljinu segmenta između P1 i P2.
#   # return sqrt((P2[1]-P1[1])^2 + (P2[2]-P1[2])^2)
# end

