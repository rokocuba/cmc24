# CMC24: Evolutionary Algorithm for Light Path Optimization

This repository contains the solution for the 2024 Computer Modeling Competition (CMC), implemented in Julia. The project features a tree-based evolutionary algorithm designed to solve a complex light path optimization problem: illuminating the largest possible area of a simulated 2D "temple" by strategically placing a set of mirrors.

<p align="center">
  <img src="output/images/cmc24_solution.png" alt="Visualization of the final mirror placement and light path" width="500"/>
  <br>
  <em>Final optimized light path illuminating the temple.</em>
</p>

### The Challenge

The official competition task required simulating the path of a light ray through a complex 2D environment. The objective was to determine the optimal number, position, and angle of N mirrors to maximize the total illuminated area.

*Official instructions can be found here: [CMC24 Technical Description](https://www.fer.unizg.hr/zpm/cmc24/tehnicki_opis)*

### My Solution: A Tree-Based Evolutionary Algorithm

To solve this high-dimensional optimization problem, I designed and implemented a custom evolutionary algorithm from scratch. The core of the solution is a branching process that iteratively grows a "family" of mirrors.

*   **Iterative Branching:** In each step, the algorithm generates and evaluates several new mirror candidates. Evaluation is based on a combined score that measures total temple area coverage and a positional heuristic.
*   **Selective Survival:** The two best-performing candidates are selected to "survive" and serve as the basis for the next generation, each branching into two new mirrors. This process repeats until the maximum number of mirrors is reached.
*   **Local Optimization:** After the primary evolutionary phase, a mirror rearrangement algorithm is applied to adjacent pairs. This post-processing step fine-tunes the solution by maximizing ray distance and reflection angles, further improving the final score.

### Tech Stack

*   **Core Language:** Julia
*   **Key Packages:** [Plots, Measures]

### Project Context & Trade-offs

The focus of this project was on developing and testing a novel algorithmic solution under the tight deadlines of a competition. As such, the primary emphasis was on **algorithmic correctness and performance** rather than on code modularity for long-term maintenance.
