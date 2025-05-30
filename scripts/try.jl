# filepath: c:\\Users\\roko\\OneDrive - fer.hr\\Programiranje\\cmc24\\try.jl
# za isprobavanje funkcije segment_length i osnovnih vektorskih operacija.

# Funkcija: segment_length
# Svrha: Računa Euklidsku udaljenost između dvije točke A i B u 2D prostoru.
# Argumenti:
#   A, B: Točke definirane kao tuple ili nizovi koordinata (npr. [x, y] ili (x, y)).
# Vraća: Duljinu segmenta (skalarnu vrijednost).
function segment_length(A, B)
   (x1, y1) = A # Destrukturiranje koordinata točke A
   (x2, y2) = B # Destrukturiranje koordinata točke B
   return sqrt((x2 - x1)^2 + (y2 - y1)^2) # Formula za udaljenost
end

# Primjer korištenja:
P = [5.0, 6.0] # Definicija točke P
M = [2.0, 2.0] # Definicija točke M

# Izračun duljine segmenta PM
d = segment_length(P, M)

# Izračun normaliziranog vektora od P prema M
# M .- P  => vektor PM (razlika koordinata)
# ./ d    => dijeljenje svakog elementa vektora PM sa skalarnom duljinom d
e = (M .- P) ./ d # e će biti jedinični vektor u smjeru od P do M

println("Točka P: ", P) # Ispis
println("Točka M: ", M)
println("Udaljenost d(P,M): ", d)
println("Normalizirani vektor od P prema M (e): ", e)