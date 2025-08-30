# Ovaj repozitorij sadrži rješenje za Natjecanje iz Računalnog Modeliranja 2024 (CMC24).
Projekt uključuje simulaciju putanje svjetlosti kroz definirani "hram" s ciljem osvjetljavanja što veće površine.

Službene upute nalaze se na: https://www.fer.unizg.hr/zpm/cmc24/tehnicki_opis

Rješenje: Evolucijski algoritam s grananjem (Branching EA)
Algoritam iterativno generira familiju zrcala. U svakom koraku kreira više novih zrcala, zatim ih evaluira pomoću kombinacije score-a (pokrivenost hrama) i heurističke metrike pozicije, te odabire samo dva najbolja. Ta dva zrcala zatim služe kao osnova za daljnje grananje – za svako od njih generiraju se po dva nova zrcala. Proces se ponavlja dok se ne postigne unaprijed definirani maksimalni broj zrcala, N, što omogućuje balans između eksploatacije i eksploracije.

Nakon toga, algoritam preraspoređivanja zrcala primjenjuje se na sve parove susjednih zrcala. Algoritam maksimizira udaljenost zrake i kut odbijanja na način da očuva ulazni i izlazni pravac zrake u odnosu na par zrcala.
