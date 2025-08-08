# CMC24 - Projekt za Natjecanje iz Računalnog Modeliranja 2024


English:
This repository contains the solution for the Computer Modeling Competition 2024 (CMC24). The project includes a simulation of a light path through a defined "temple" with the goal of illuminating as large an area as possible.

Official instructions can be found at: [https://www.fer.unizg.hr/zpm/cmc24/tehnicki_opis](https://www.fer.unizg.hr/zpm/cmc24/tehnicki_opis)

Ovaj repozitorij sadrži rješenje za Natjecanje iz Računalnog Modeliranja 2024 (CMC24). Projekt uključuje simulaciju putanje svjetlosti kroz definirani "hram" s ciljem osvjetljavanja što veće površine.

Službene upute nalaze se na: [https://www.fer.unizg.hr/zpm/cmc24/tehnicki_opis](https://www.fer.unizg.hr/zpm/cmc24/tehnicki_opis)

Solution/Rješenje:

Tree-Based Evolutionary Algorithm (Branching EA)

English:
The algorithm iteratively generates a family of mirrors by creating several new mirrors in each step. It then evaluates each mirror using a combination of a score (temple coverage) and a heuristic position metric, selecting only the two best mirrors. These two mirrors then serve as the basis for further branching—two new mirrors are generated for each of them. This process is repeated until a predefined maximum number of mirrors, N, is reached. This approach allows for a balance between exploitation (selecting the best) and exploration (generating new variants).

After this, a mirror rearrangement algorithm (a maximization algorithm) is applied to all pairs of adjacent mirrors to maximize the ray distance and the reflection angle, in a way that preserves the incoming and outgoing ray directions (relative to the pair of mirrors).

Hrvatski:
Algoritam iterativno generira familiju zrcala tako da u svakom koraku kreira više novih zrcala, zatim evaluira svako zrcalo pomoću kombinacije score-a (pokrivenost hrama) i heurističke metrike pozicije, te odabire samo dva najbolja. Ta dva zrcala zatim služe kao osnova za daljnje grananje – za svako od njih generiraju se po dva nova zrcala, što se ponavlja sve dok se ne postigne unaprijed definirani maksimalni broj zrcala N. Ovaj pristup omogućuje balans između eksploatacije (odabir najboljih) i eksploracije (generiranje novih varijanti). 

Nakon toga se algoritmom preraspoređivanja zrcala (algoritam maksimizacije) za sve parove susjednih zrcala maksimizira udaljenost zrake odnosno kut odbijanja na način da očuva ulazni i izlazni pravac zrake (u odnosu na par zrcala).

<img src="output/images/cmc24_solution.png" alt="Rješenje" width="400"/>

## Struktura Projekta

Projekt je organiziran na sljedeći način:

*   `src/`: Sadrži glavne Julia skripte (`.jl`) koje implementiraju logiku natjecanja.
    *   `cmc24.jl`: Glavna skripta.
    *   `cmc24C.jl`: Skripta za evaluaciju.
    *   Ostale `.jl` datoteke: Pomoćne skripte, editori rješenja ili specifične verzije.
*   `scripts/`: Sadrži pomoćne ili testne skripte.
    *   `points_try.jl`: Skripta za testiranje rada s točkama ili geometrijom.
    *   `try.jl`: Općenita skripta za isprobavanje manjih dijelova koda.
*   `output/`: Sadrži generirane izlazne datoteke.
    *   `images/`: Spremljene slike i vizualizacije (npr. `.png` datoteke).
    *   `heatgrids/`: Podaci vezani za "heatgrid" analize.
    *   `jld2_files/`: Spremljena rješenja i drugi podaci u Julia Data Format (`.jld2`).
*   `README.md`: Ova datoteka s opisom projekta.
*   `.gitignore`: Specificira datoteke koje Git treba ignorirati.
*   `Project.toml` i `Manifest.toml`: (Ako postoje) Definiraju ovisnosti Julia projekta.

## Pokretanje

Za pokretanje glavne simulacije ili evaluacije, vjerojatno ćete koristiti Julia interpreter.

1.  **Instalirajte Juliu:** Ako već niste, preuzmite i instalirajte Juliu s [julialang.org](https://julialang.org/).
2.  **Navigirajte do direktorija projekta**
3.  **Aktivirajte projektno okruženje (ako postoji `Project.toml`):**
    ```julia
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate() # Instalira potrebne pakete ako već nisu
    ```
4.  **Pokrenite glavnu skriptu:**
    ```bash
    julia src/cmc24.jl
    ```

## Ovisnosti

Glavne ovisnosti koje se koriste u projektu (na temelju pregledanih skripti):
*   `FileIO`
*   `ImageIO`
*   `Measures`
*   `Plots`
*   `UUIDs`
*   `JLD2`

Ako koristite `Project.toml` i `Manifest.toml`, Julia će automatski upravljati ovim ovisnostima.
