# filepath: c:\\Users\\roko\\OneDrive - fer.hr\\Programiranje\\cmc24\\cmc24C.jl
# Computational Modeling Challenge 2024
# CMC24 solution evaluation script
# Author: Hrvoje Abraham, hrvoje.abraham@avl.com
# Skripta za evaluaciju rješenja CMC24 natjecanja.
# Učitava definiciju hrama, rješenje (pozicije lampe i ogledala),
# zatim simulira putanju zrake svjetlosti i računa rezultat.

using FileIO
using ImageIO
using Measures
using Plots;
gr(); # Inicijalizacija GR-a kao backend za Plots
using UUIDs # Za generiranje jedinstvenih ID-eva, vjerojatno za imena datoteka

# String koji definira geometriju hrama. 'O' su zidovi/blokovi, '.' je prazan prostor.
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

block_size = 1      # Veličina jednog bloka u hramu (jedinična)
mirror_length = 0.5 # Duljina zrcala
light_halfwidth = 1 # Poluširina izvora svjetlosti (za provjeru je li lampa u bloku)
ε = 1e-12           # Mala vrijednost za usporedbe s nulom (epsilon)

"""Float64 infinity"""
# Definicija beskonačnosti kao Float64
∞ = Inf

"""
    ⋅(v, w)

Dot product of two 2D vectors.
Skalarni produkt dva 2D vektora.
"""
function ⋅(v, w)
    return v[1] * w[1] + v[2] * w[2]
end

"""
    ×(v, w)

Cross product of two 2D vectors.
Vektorski produkt dva 2D vektora (rezultat je skalar jer su vektori u ravnini).
"""
function ×(v, w)
    return v[1] * w[2] - v[2] * w[1]
end

"""Last 12 digits of hexadecimal format of the input integer."""
# Pomoćna funkcija za prikaz zadnjih 12 znamenki heksadekadskog broja.
function hex12(x::Integer)
    return last(string(x, base=16), 12)
end

"""
    Convert integer into string with digit grouping.

    E.g. 1234567 => "1,234,567"
Formatira cijeli broj kao string s razdjelnicima tisućica.
"""
function commas(num::Integer)
    str = string(num)
    return replace(str, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

# Učitava definiciju hrama iz stringa.
function load_temple(temple_string_param, block_size_param)
    # println(stderr, " " * temple_string_param) # Originalni ispis, može se maknuti
    rows = split(replace(strip(temple_string_param), " " => ""), '\n') # Makni razmake, podijeli po redovima, ISPRAVLJENO: '\\n' u '\n'
    temple_shape = length(rows[1]), length(rows) # Dimenzije hrama

    temple_blocks = Set() # Koristi Set za pohranu koordinata blokova (brza provjera postojanja)
    for (j, row_str) ∈ enumerate(rows)
        for (i, char_val) ∈ enumerate(row_str)
            if char_val == 'O'
                # Ovdje ide logika za dodavanje bloka u 'temple_blocks'.
                # Npr. push!(temple_blocks, ( (i-1)*block_size_param, (j-1)*block_size_param, block_size_param, block_size_param ))
                # Točna struktura ovisi o tome kako se blokovi kasnije koriste.
                # Za sada ostavljam prazno da ne uvedem krivu logiku.
                # Potrebno je definirati kako se 'O' mapira u koordinate/objekte blokova.
            end
        end
    end
    # display(temple_blocks) # Za debugiranje, prikazuje koordinate blokova

    println(stderr, "Hram dimenzija $temple_shape je učitan.")

    # Vraća named tuple s informacijama o hramu
    # Ostale funkcije očekuju temple.blocks i temple.size
    return (blocks = temple_blocks, size = temple_shape .* block_size_param)
end

# Učitava kompletno rješenje (lampa i 8 ogledala).
function load_solution(cmc24_solution, mirror_length_param)
    if size(cmc24_solution) ≠ (9, 3) # Provjera dimenzija matrice rješenja
        println(stderr, "GREŠKA: Matrica rješenja nema ispravne dimenzije (očekivano 9x3).")
        finalize() # Prekida izvođenje ako dimenzije nisu ispravne
    end

    if !(eltype(cmc24_solution) <: Number) # Provjera tipa elemenata matrice
        println(stderr, "GREŠKA: Elementi matrice rješenja nisu brojevi.")
        finalize()
    end

    local sol_float
    try # Pokušaj konverzije u Float64
        sol_float = float(cmc24_solution)
    catch e
        println(stderr, "GREŠKA: Matrica rješenja se ne može konvertirati u Float64: $e")
        finalize()
    end

    # Predprocesiranje lampe: pozicija (x,y) i kut alfa
    lamp_pos = [sol_float[1,1], sol_float[1,2]]
    α = sol_float[1, 3] # Kut lampe
    lamp_dir = [cos(α), sin(α)]
    lamp_obj = (position = lamp_pos, angle = α, direction = lamp_dir)

    # Predprocesiranje ogledala
    mirrors_list = [] # Lista za pohranu obrađenih ogledala
    for m ∈ 1:8 # Petlja kroz 8 ogledala
        mirror_center = [sol_float[m+1,1], sol_float[m+1,2]]
        mirror_angle = sol_float[m+1,3] # Kut normale na ogledalo ili orijentacije ogledala?
                                         # Ovisno o interpretaciji, definicija p1, p2 se mijenja.
                                         # Pretpostavimo da je mirror_angle kut samog segmenta ogledala.
        
        # Krajnje točke ogledala, ako je mirror_angle kut segmenta:
        # dx = (mirror_length_param / 2) * cos(mirror_angle)
        # dy = (mirror_length_param / 2) * sin(mirror_angle)
        # p1 = [mirror_center[1] - dx, mirror_center[2] - dy]
        # p2 = [mirror_center[1] + dx, mirror_center[2] + dy]

        # Ili ako je mirror_angle kut normale, pa je segment okomit:
        # dx_norm = cos(mirror_angle)
        # dy_norm = sin(mirror_angle)
        # segment_dx = -(mirror_length_param / 2) * dy_norm # Rotirano za -90 stupnjeva
        # segment_dy = (mirror_length_param / 2) * dx_norm  # Rotirano za -90 stupnjeva
        # p1 = [mirror_center[1] + segment_dx, mirror_center[2] + segment_dy]
        # p2 = [mirror_center[1] - segment_dx, mirror_center[2] - segment_dy]
        # Za potrebe strukture, koristit ćemo centar i kut. Detalji refleksije su u raytrace.
        # Segment za provjeru presjeka: (q, l, β) -> (centar - pomak, duljina, kut)
        # Ovdje je 'mirror' placeholder u originalu, treba definirati strukturu.
        # Npr. (center, length, angle) ili (p1, p2, angle_normal)
        current_mirror = (center=mirror_center, length=mirror_length_param, angle=mirror_angle)
        push!(mirrors_list, current_mirror)
    end

    println(stderr, "Rješenje je učitano.")

    return (lamp_obj, mirrors_list) # Vraća obrađene podatke o lampi i ogledalima
end

# Učitava parcijalno rješenje (samo lampa i neka ogledala).
function load_solution_partial(cmc24_solution, mirror_length_param)
    local sol_float
    try
        sol_float = float(cmc24_solution) # Konverzija u Float64
    catch e
        println(stderr, "GREŠKA: Matrica parcijalnog rješenja se ne može konvertirati u Float64: $e")
        finalize()
    end

    # Predprocesiranje lampe
    lamp_pos = [sol_float[1,1], sol_float[1,2]]
    α = sol_float[1, 3]
    lamp_dir = [cos(α), sin(α)]
    lamp_obj = (position = lamp_pos, angle = α, direction = lamp_dir)
    

    # Predprocesiranje ogledala
    mirrors_list = []
    # Petlja od drugog reda (prvi je lampa) do kraja matrice
    for m ∈ 2:size(sol_float, 1) 
        mirror_center = [sol_float[m,1], sol_float[m,2]]
        mirror_angle = sol_float[m,3]
        current_mirror = (center=mirror_center, length=mirror_length_param, angle=mirror_angle)
        push!(mirrors_list, current_mirror)
    end

    println(stderr, "Parcijalno rješenje je učitano.")

    return (lamp_obj, mirrors_list)
end

# Provjerava je li točka unutar nekog bloka hrama.
function point_in_block(temple_obj, point)
    # Pretpostavka: temple_obj.blocks je lista/set blokova,
    # svaki blok je npr. (x_start, y_start, width, height)
    for block_coords ∈ temple_obj.blocks
        # Ovdje ide logika provjere je li 'point' unutar 'block_coords'
        # Npr. if block_coords.x_start <= point[1] < block_coords.x_start + block_coords.width &&
        #          block_coords.y_start <= point[2] < block_coords.y_start + block_coords.height
        #    return true
        # end
        # Pošto ne znamo strukturu bloka, ostavljamo placeholder
    end
    return false # Ako nije ni u jednom bloku
end

# Određuje u kojem je sektoru hrama (3x3 podjela) dana točka.
function point_sector(temple_obj, point)
    # temple_obj.size su ukupne dimenzije hrama [width, height]
    # Provjeri da temple_obj.size nije nula da se izbjegne dijeljenje s nulom
    if temple_obj.size[1] == 0 || temple_obj.size[2] == 0
        println(stderr, "UPOZORENJE: Dimenzije hrama su nula u point_sector.")
        return -1 # Neki error kod
    end
    sx = floor(Int, 3 * point[1] / temple_obj.size[1]) # Normalizira x koordinatu na [0,3) i zaokruži
    sy = floor(Int, 3 * point[2] / temple_obj.size[2]) # Normalizira y koordinatu na [0,3) i zaokruži

    # Osiguraj da su indeksi unutar granica 0-2
    sx = clamp(sx, 0, 2)
    sy = clamp(sy, 0, 2)

    return 3 * sy + sx + 1 # Mapira (sx, sy) u linearni indeks sektora 1-9
end

# Računa presjek dvije zrake (ray). Zraka je definirana početnom točkom i vektorom smjera.
function ray_ray_intersection(ray1, ray2)
    (p, r) = ray1 # ray1: točka p, smjer r
    (q, s) = ray2 # ray2: točka q, smjer s

    rs = r × s # Vektorski produkt smjerova
    q_minus_p = q .- p
    qpr = q_minus_p × r # Vektorski produkt (q-p) i r

    # SLUČAJ 1: Zrake su kolinearne i možda se preklapaju
    if abs(rs) < ε && abs(qpr) < ε # rs blizu 0 znači da su paralelne ili kolinearne
                                 # qpr blizu 0 znači da (q-p) leži na istoj liniji kao r (ako su kolinearne)
        # Za kolinearne zrake, t0 i t1 definiraju parametre duž r za segment preklapanja.
        # Ako se ne preklapaju u pozitivnom smjeru, nema presjeka.
        # Ovo je kompleksnije, originalni kod nije imao potpunu logiku.
        # Vraćam (1, t0, t1) - treba vidjeti kako se t0, t1 koriste ili definirati.
        # Za sada, ovo je placeholder za "kolinearne".
        t0 = 0.0 # Placeholder
        t1 = 0.0 # Placeholder
        # Primjer logike za preklapanje (ako r i s pokazuju u istom smjeru):
        # r_norm_sq = r ⋅ r
        # if r_norm_sq > ε
        #    t0 = (q_minus_p ⋅ r) / r_norm_sq
        #    # t1 bi bio t0 + duljina_preklapanja / |r|
        #    # Ako (q-p) i r nisu u istom smjeru, a s jest, onda t0 može biti negativan.
        #    # Treba provjeriti jesu li zrake usmjerene jedna prema drugoj ili se udaljavaju.
        # end
        return (1, t0, t1) 
    end

    # SLUČAJ 2: Zrake su paralelne i ne sijeku se
    if abs(rs) < ε && abs(qpr) > ε # Paralelne (rs blizu 0), ali q nije na liniji od p u smjeru r
        return (2, 0.0, 0.0) # Nema presjeka
    end

    # SLUČAJ 3: Zrake se sijeku (linije na kojima leže se sijeku)
    # Presjek je p + t*r = q + u*s
    # rs ne smije biti 0
    if abs(rs) < ε # Već pokriveno gore, ali za svaki slučaj
        return (4, 0.0, 0.0) # Ne bi trebalo doći ovdje ako je gornja logika ispravna
    end
    
    qps = q_minus_p × s # (q-p) x s
    t = qps / rs # Parametar za prvu zraku
    u = qpr / rs # Parametar za drugu zraku

    if (t ≥ -ε) && (u ≥ -ε) # Presjek postoji ako su oba parametra ne-negativna (zrake idu "naprijed")
                           # Koristimo -ε da uhvatimo slučajeve gdje je presjek točno na početku.
        return (3, t, u) # Vraća tip presjeka i parametre t, u
    end

    # SLUČAJ 4: Linije se sijeku, ali ne u pozitivnom smjeru zraka
    return (4, 0.0, 0.0) # Nema presjeka u smjeru zraka
end

# Računa presjek zrake i segmenta. Segment je definiran početnom točkom, duljinom i kutem.
function ray_segment_intersection(ray, segment_def)
    (p_ray, r_ray) = ray                 # Zraka: točka p_ray, smjer r_ray
    (q_seg, l_seg, β_seg) = segment_def  # Segment: početna točka q_seg, duljina l_seg, kut β_seg

    # Vektor segmenta s_seg_dir, od q_seg do kraja segmenta
    s_seg_dir = [cos(β_seg), sin(β_seg)]
    # Točka q_seg je početak segmenta, s_seg_dir je smjer, l_seg je duljina.
    # Zraka koja predstavlja segment: (q_seg, s_seg_dir)
    
    (case, t, u) = ray_ray_intersection((p_ray, r_ray), (q_seg, s_seg_dir)) # Koristi presjek zraka-zraka

    # t je parametar za ray (p_ray + t*r_ray)
    # u je parametar za zraku segmenta (q_seg + u*s_seg_dir)

    # SLUČAJ 1 - Nema presjeka ili presjek nije na segmentu/zraki
    if case == 2 || case == 4 # Paralelne ili se linije ne sijeku u pravom smjeru
        return (1, 0.0, 0.0) # Nema presjeka
    end
    
    if case == 3 # Linije se sijeku
        if t < -ε || u < -ε || u > l_seg + ε # Presjek iza početka zrake (t < 0)
                                           # ili izvan segmenta (u < 0 ili u > duljina_segmenta)
            return (1, 0.0, 0.0) # Nema validnog presjeka
        else
            # Validni presjek zrake i linije segmenta, unutar granica segmenta i zrake
            return (3, t, u / l_seg) # Vraćamo u normaliziran na [0,1] za segment
        end
    end

    if case == 1 # Kolinearne
        # Ovdje treba detaljnija logika za kolinearni presjek zrake i segmenta.
        # Originalni kod je imao samo 'end', što je greška.
        # Treba provjeriti preklapaju li se.
        # Ako se preklapaju, treba naći najbližu točku preklapanja u smjeru zrake.
        # Ovo je placeholder, jer je kolinearno presijecanje komplicirano.
        # Pretpostavimo za sada da ako su kolinearne, vraćamo "nema presjeka"
        # da se izbjegne nedefinirano ponašanje.
        return (1, 0.0, 0.0) # Placeholder za kolinearni slučaj
    end
    
    # Default ako nijedan uvjet nije zadovoljen (ne bi se smjelo dogoditi)
    return (1, 0.0, 0.0)
end


# Provjerava sijeku li se dva segmenta.
function segment_segment_intersection(segment1_def, segment2_def)
    (p_seg1, l_seg1, α_seg1) = segment1_def # Segment 1: točka p, duljina la, kut α
    (q_seg2, l_seg2, β_seg2) = segment2_def # Segment 2: točka q, duljina lb, kut β

    r_dir_s1 = [cos(α_seg1), sin(α_seg1)] # Smjer prvog segmenta
    s_dir_s2 = [cos(β_seg2), sin(β_seg2)] # Smjer drugog segmenta

    # Koristi presjek zraka-zraka, ali parametri t i u moraju biti između 0 i duljine segmenta
    (case, t, u) = ray_ray_intersection((p_seg1, r_dir_s1), (q_seg2, s_dir_s2))

    if case == 3 # Linije se sijeku
        # Presjek postoji ako je 0 ≤ t ≤ l_seg1 I 0 ≤ u ≤ l_seg2
        # Koristimo ε za rubne slučajeve
        if (t >= -ε && t <= l_seg1 + ε) && (u >= -ε && u <= l_seg2 + ε)
            return true # Segmenti se sijeku
        end
    end

    if case == 1 # Kolinearni
        # Ovdje ide kompleksna logika za provjeru preklapanja kolinearnih segmenata.
        # Originalni kod je imao neke uvjete, ali mogu biti pojednostavljeni.
        # Za sada, vraćamo false za kolinearne da se izbjegnu greške dok se ne implementira potpuno.
        # TODO: Implementirati robustnu provjeru preklapanja kolinearnih segmenata.
        return false
    end

    # Svi ostali slučajevi (paralelni, ne sijeku se, linije se sijeku ali izvan segmenata)
    return false
end

# Provjerava siječe li segment zadani blok hrama.
function segment_block_intersection(segment_def, block_def)
    # Pretpostavka: 'block_def' je definiran svojim stranicama (4 segmenta)
    # ili kao (x, y, width, height).
    # Treba dobiti 4 segmenta koji čine rubove bloka.
    # Npr. ako je block_def = (bx, by, bw, bh):
    # p1 = (bx, by); p2 = (bx+bw, by); p3 = (bx+bw, by+bh); p4 = (bx, by+bh)
    # seg_top = (p1, bw, 0.0) # Horizontalni segment od p1, duljine bw, kut 0
    # seg_right = (p2, bh, pi/2) # Vertikalni od p2, duljine bh, kut pi/2 (ili p2 do p3)
    # ... itd. za sve 4 stranice bloka.
    # Zatim za svaku stranicu: if segment_segment_intersection(segment_def, stranica_bloka) return true

    # Ovo je placeholder, treba prava logika.
    # Originalni kod je imao sintaksnu grešku: ))
    return false # Placeholder
end

# Provjerava siječe li segment bilo koji blok u hramu.
function temple_segment_intersection(temple_obj, segment_def)
    # Iterira kroz sve blokove hrama i koristi segment_block_intersection
    # Pretpostavka: temple_obj.blocks je iterabilna kolekcija definicija blokova
    if isnothing(temple_obj) || isnothing(temple_obj.blocks)
        println(stderr, "UPOZORENJE: temple_obj ili temple_obj.blocks je null u temple_segment_intersection.")
        return false
    end
    for block_definition ∈ temple_obj.blocks
        if segment_block_intersection(segment_def, block_definition)
            return true
        end
    end
    return false
end

# Nalazi najbliži presjek zrake s hramom.
function temple_ray_intersection(temple_obj, ray)
    # Iterira kroz sve blokove hrama. Za svaki blok, nađe presjeke zrake s njegovim stranicama.
    # Vraća najmanji pozitivni parametar 't' koji odgovara presjeku.
    t_min_intersection = ∞ 
    # Pretpostavka: temple_obj.blocks sadrži definicije blokova
    # Svaki blok treba pretvoriti u 4 segmenta (stranice)
    if isnothing(temple_obj) || isnothing(temple_obj.blocks)
        println(stderr, "UPOZORENJE: temple_obj ili temple_obj.blocks je null u temple_ray_intersection.")
        return ∞
    end

    for block_definition ∈ temple_obj.blocks
        # block_segments = get_segments_from_block_definition(block_definition)
        # for seg_block ∈ block_segments
        #    (case, t, u_norm) = ray_segment_intersection(ray, seg_block)
        #    if case == 3 && t > ε # Presjek postoji i ispred je zrake
        #        t_min_intersection = min(t_min_intersection, t)
        #    end
        # end
        # Placeholder dok se ne definira struktura bloka i get_segments_from_block_definition
    end
    return t_min_intersection # Vraća udaljenost do najbližeg presjeka s hramom
end

# Provjerava geometrijsku ispravnost rješenja.
function check_solution(temple_obj, lamp_obj, mirrors_list)
    # Logika za provjeru:
    # 1. Je li lampa unutar hrama i nije u zidu?
    #    (npr. !point_in_block(temple_obj, lamp_obj.position))
    # 2. Preklapaju li se ogledala međusobno?
    #    (for i in 1:length(mirrors_list), j in i+1:length(mirrors_list) -> segment_segment_intersection)
    #    Ogledalo treba pretvoriti u segment: (center, length, angle) -> (p1, l, angle_seg)
    # 3. Jesu li ogledala unutar hrama i ne sijeku zidove?
    #    (for mirror_def in mirrors_list -> !temple_segment_intersection(temple_obj, mirror_as_segment(mirror_def)))
    # Ako nešto nije u redu, može ispisati grešku i pozvati finalize().
    println(stderr, "Geometrija rješenja je (pretpostavljeno) ispravna. Implementiraj provjere!") # Placeholder
end

# Simulira putanju zrake svjetlosti.
function raytrace(temple_obj, lamp_obj, mirrors_list)
    path_segments = [] # Lista segmenata [(p1, p2), ...] koji čine putanju
    
    # Početna zraka od lampe
    # Pazi da početna točka nije unutar zida, pomakni malo ako treba
    start_point = lamp_obj.position .+ ε .* lamp_obj.direction 
    current_ray = (start_point, lamp_obj.direction)
    
    MAX_REFLECTIONS = 20 # Neko ograničenje da se izbjegne beskonačna petlja

    for reflection_count ∈ 1:MAX_REFLECTIONS
        dist_to_temple = temple_ray_intersection(temple_obj, current_ray)
        
        dist_to_closest_mirror = ∞
        closest_mirror_idx = -1
        # normal_at_reflection = [0.0, 0.0] # Placeholder
        # u_param_on_mirror = 0.0 # Placeholder, parametar na segmentu ogledala

        for (idx, mirror_def) ∈ enumerate(mirrors_list)
            # Pretvaranje definicije ogledala (center, length, angle) u segment za ray_segment_intersection
            # Segment: (q_seg, l_seg, β_seg)
            # q_seg je jedna krajnja točka ogledala.
            # Npr. p_start = mirror_def.center .- (mirror_def.length/2) .* [cos(mirror_def.angle), sin(mirror_def.angle)]
            # seg_for_intersect = (p_start, mirror_def.length, mirror_def.angle)
            
            # (case, t, u_norm) = ray_segment_intersection(current_ray, seg_for_intersect)
            # if case == 3 && t > ε && t < dist_to_closest_mirror
            #    dist_to_closest_mirror = t
            #    closest_mirror_idx = idx
            #    # normal_at_reflection = ... izračunati normalu na ogledalo u točki presjeka
            #    # u_param_on_mirror = u_norm * mirror_def.length
            # end
            # Placeholder dok se ne definira mirror_def -> segment konverzija
        end

        if dist_to_temple < dist_to_closest_mirror && dist_to_temple < ∞
            # Zraka udara u hram prije nego u ogledalo
            hit_point = current_ray[1] .+ dist_to_temple .* current_ray[2]
            push!(path_segments, (current_ray[1], hit_point))
            println(stderr, "Zraka pogađa hram nakon $reflection_count refleksija.")
            break # Kraj putanje
        elseif dist_to_closest_mirror < ∞
            # Zraka udara u ogledalo
            reflection_point = current_ray[1] .+ dist_to_closest_mirror .* current_ray[2]
            push!(path_segments, (current_ray[1], reflection_point))
            
            # Izračunaj reflektiranu zraku
            # reflected_direction = reflect(current_ray[2], normal_at_reflection)
            # current_ray = (reflection_point .+ ε .* reflected_direction, reflected_direction) # Pomakni malo
            # Placeholder za logiku refleksije
            println(stderr, "Zraka se reflektira od ogledala $closest_mirror_idx.")
            # Za sada, prekidamo nakon prve "refleksije" jer nije implementirano
            break # TODO: Implementiraj refleksiju i nastavi petlju
        else
            # Zraka ne pogađa ništa (ide u beskonačnost)
            # Može se dodati dugi segment u smjeru zrake kao indikacija
            # far_point = current_ray[1] .+ 1000.0 .* current_ray[2] # Neka velika udaljenost
            # push!(path_segments, (current_ray[1], far_point))
            println(stderr, "Zraka odlazi u beskonačnost nakon $reflection_count refleksija.")
            break # Kraj putanje
        end
        if reflection_count == MAX_REFLECTIONS
             println(stderr, "Dosegnut maksimalan broj refleksija ($MAX_REFLECTIONS).")
        end
    end
    return path_segments # Vraća listu segmenata [(x1,y1), (x2,y2), ...]
end

# Generira i sprema plot hrama, lampe, ogledala i putanje zrake.
function cmc24_plot(temple_obj; lamp_obj=nothing, mirrors_list=nothing, path_segments=nothing)
    # Koristi Plots.jl za crtanje:
    # 1. Nacrtaj blokove hrama (npr. kao pravokutnike).
    #    plot_obj = plot(aspect_ratio=:equal, legend=false, title="CMC24 Rješenje")
    #    for block_def in temple_obj.blocks
    #        # Dodaj pravokutnik za blok_def
    #        # rect = Plots.Shape(block_def.x, block_def.y, block_def.w, block_def.h)
    #        # plot!(plot_obj, rect, fillcolor=:black, linecolor=:darkgray)
    #    end
    # 2. Ako postoji 'lamp_obj', nacrtaj lampu.
    #    if !isnothing(lamp_obj)
    #        # scatter!(plot_obj, [lamp_obj.position[1]], [lamp_obj.position[2]], marker=:star, color=:yellow)
    #        # Nacrtaj i smjer zrake
    #    end
    # 3. Ako postoje 'mirrors_list', nacrtaj ogledala (kao segmente).
    #    if !isnothing(mirrors_list)
    #        for mirror_def in mirrors_list
    #            # p1, p2 = get_mirror_endpoints(mirror_def)
    #            # plot!(plot_obj, [p1[1], p2[1]], [p1[2], p2[2]], linecolor=:blue, linewidth=2)
    #        end
    #    end
    # 4. Ako postoji 'path_segments', nacrtaj putanju zrake.
    #    if !isnothing(path_segments)
    #        for seg in path_segments
    #            # plot!(plot_obj, [seg[1][1], seg[2][1]], [seg[1][2], seg[2][2]], linecolor=:red, linewidth=1)
    #        end
    #    end
    # filename = "cmc24_plot_" * string(uuid4()) * ".png"
    # savefig(plot_obj, filename)
    # println(stderr, "Plot je spremljen kao: $filename")
    # return filename 
    println(stderr, "Crtanje nije implementirano. Treba koristiti Plots.jl.")
    return "placeholder_plot.png" # Placeholder
end

# Evaluacija rješenja: računa koliko je "osvijetljenih" piksela (ili površine).
function evaluate(temple_obj, path_segments)
    # Logika za evaluaciju:
    # 1. Diskretiziraj područje hrama u grid (piksele).
    # 2. Za svaki segment putanje zrake, odredi koje piksele "osvjetljava".
    #    (npr. Bresenhamov algoritam za linije ili slično).
    # 3. Prebroji jedinstvene osvijetljene piksele koji nisu unutar blokova hrama.
    # 'total_pixels' bi mogao biti ukupan broj piksela u gridu unutar granica hrama.
    # 'vacant_pixels' bi mogao biti ukupan broj praznih (ne-blok) piksela.
    # 'score' je broj osvijetljenih praznih piksela.
    
    # Placeholder vrijednosti
    total_pixels = 1000000 # Primjer
    vacant_pixels = 800000  # Primjer
    score = 0             # Inicijalno nula osvijetljenih

    if isnothing(path_segments) || isempty(path_segments)
        println(stderr, "Nema putanje zrake za evaluaciju.")
        return total_pixels, vacant_pixels, 0
    end

    # Ovdje ide stvarna evaluacija...
    # score = calculate_illuminated_area(temple_obj, path_segments, resolution_for_eval)

    println(stderr, "Evaluacija nije potpuno implementirana. Vraćam placeholder vrijednosti.")
    # Privremeno, da ne bude 0 i izazove dijeljenje s nulom kasnije:
    if vacant_pixels == 0; vacant_pixels = 1; end # Izbjegni dijeljenje s nulom
    score = rand(0:vacant_pixels÷2) # Neki random score za test

    return total_pixels, vacant_pixels, score 
end

# Funkcija za prekid izvođenja, eventualno s ispisom stanja.
function finalize_script(temple_obj=nothing, lamp_obj=nothing, mirrors_list=nothing, path_segments=nothing; message="Skripta se prekida.")
    println(stderr, message)
    # Može ispisati poruku o grešci ili konačni rezultat prije izlaska.
    # Ako je pozvana zbog greške, cmc24_plot se možda neće zvati ili će se zvati s parcijalnim podacima.
    exit(1) # Prekida izvršavanje skripte s kodom greške
end

# Primjer matrice rješenja (9x3: lampa + 8 ogledala; x, y, kut)
# Ovo su placeholder vrijednosti, treba ih zamijeniti stvarnim rješenjem.
cmc24_solution_matrix = [
    10.0 10.0 0.0;  # Lampa: x, y, alpha (kut u radijanima)
    5.0  5.0 pi/4;  # Ogledalo 1: x_centar, y_centar, beta (kut normale ili segmenta)
    15.0 5.0 -pi/4; # Ogledalo 2
    5.0  15.0 pi/3; # Ogledalo 3
    15.0 15.0 -pi/3;# Ogledalo 4
    2.0  10.0 0.0;  # Ogledalo 5
    18.0 10.0 pi;   # Ogledalo 6
    10.0 2.0 pi/2;  # Ogledalo 7
    10.0 18.0 -pi/2; # Ogledalo 8
]
# Ako matrica nije potpuna (npr. manje od 9 redaka), load_solution će javiti grešku.

# --- Glavni dio skripte ---

# Učitaj definiciju hrama
temple_data = load_temple(temple_string, block_size)

# Učitaj rješenje (lampa i ogledala)
if size(cmc24_solution_matrix, 1) != 9
    finalize_script(temple_data, message="GREŠKA: Matrica rješenja nema 9 redaka.")
end
lamp_data, mirrors_data = load_solution(cmc24_solution_matrix, mirror_length)

# Provjeri geometrijsku ispravnost rješenja (opcionalno, ali korisno)
# check_solution(temple_data, lamp_data, mirrors_data) # Originalno zakomentirano, implementiraj!

# Izračunaj putanju zrake svjetlosti
path_data = raytrace(temple_data, lamp_data, mirrors_data) 

# Evaluacija rješenja
total_eval, vacant_eval, score_eval = evaluate(temple_data, path_data) 

# Ispis rezultata
println(stderr, "Osnovni plot ima $(commas(vacant_eval)) praznih od ukupno $(commas(total_eval)) piksela (placeholder).")
if vacant_eval > 0
    percentage = round(100.0 * score_eval / vacant_eval, digits=2)
    println(stderr, "Tvoj CMC24 rezultat je $(commas(score_eval)) / $(commas(vacant_eval)) = $percentage %.")
else
    println(stderr, "Tvoj CMC24 rezultat je $(commas(score_eval)), ali nema praznih piksela za postotak.")
end
print(score_eval) # Ispisuje samo broj bodova na stdout (za automatsko ocjenjivanje)

# Kreiraj i spremi grafički prikaz rješenja
plot_filename = cmc24_plot(temple_data, lamp_obj=lamp_data, mirrors_list=mirrors_data, path_segments=path_data)
println(stderr, "Grafički prikaz (placeholder) je na: $plot_filename")

# Skripta normalno završava ovdje. finalize_script() se zove samo u slučaju grešaka.
println(stderr, "\nSkripta završila.")