# This repository contains the solution for the Computer Modeling Competition 2024 (CMC24).
The project includes a simulation of a light path through a defined "temple" with the goal of illuminating as large an area as possible.

Official instructions can be found at: https://www.fer.unizg.hr/zpm/cmc24/tehnicki_opis

Solution: Tree-Based Evolutionary Algorithm (Branching EA)
The algorithm iteratively generates a family of mirrors. In each step, it creates several new mirrors, evaluates them based on a combined score (temple coverage) and a heuristic position metric, and selects the two best. These two mirrors then serve as the basis for further branching, with two new mirrors generated from each. This process is repeated until a predefined maximum number of mirrors, N, is reached, allowing for a balance between exploitation and exploration.

After this, a mirror rearrangement algorithm is applied to all pairs of adjacent mirrors. The algorithm maximizes the distance of the ray and the reflection angle while preserving the incoming and outgoing ray directions relative to the pair of mirrors.
