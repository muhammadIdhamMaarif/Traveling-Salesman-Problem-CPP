# Traveling Salesman Problem (TSP) Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

This repository provides a comprehensive C++ implementation of multiple algorithms to solve the **Traveling Salesman Problem (TSP)** — a classic NP-hard combinatorial optimization problem.

The repository is designed as a single C++ class `TSP` that encapsulates various algorithms, including:

- **Backtracking** (with and without fixed start city)  
- **Backtracking with mirror elimination** (reduces symmetric duplicates)  
- **Dynamic Programming (Held-Karp algorithm)**  
- **Branch and Bound**  
- **Nearest Neighbor heuristic**  
- **Greedy Edge Selection heuristic**  
- **Christofides Algorithm** (approximation algorithm with 3/2 approximation guarantee)  
- **Simulated Annealing** (metaheuristic)

Each algorithm has different time and space complexity characteristics, providing a trade-off between accuracy and performance.

---

## Repository Contents

- `tsp.hpp` — Header-only class file with all TSP algorithms implemented in a class
- `main.cpp` — Sample usage demonstrating each algorithm
- `pnr.cpp` — Paste and run code to test directly in compiler
- `LICENSE` — MIT License file  
- `README.md` — This file

---

## Algorithms Explained

### 1. Backtracking (Permutation Enumeration)

- **Methods:** `SolveWithNonFixedFirstCity()`, `SolveWithFixedFirstCity()`  
- **Idea:** Generates all permutations of city orders, computes total distance, and tracks the shortest.  
- **Complexity:**  
  - Non-fixed start: \( O(n!) \)  
  - Fixed first city: \( O((n-1)!) \) (reduces duplicates by fixing start)  
- **Reference:**  
  - Introduction to Algorithms (Cormen et al.)  
  - [TSP Wikipedia: Brute Force](https://en.wikipedia.org/wiki/Travelling_salesman_problem#Brute_force_search)

---

### 2. Mirror Elimination

- **Method:** `SolveWithEliminatingMirror()`  
- **Idea:** Avoids counting routes that are reversals of each other (which have same cost in symmetric TSP), further cutting search space roughly in half.  
- **Effect:** Reduces permutations from \( (n-1)! \) to approximately \( \frac{(n-1)!}{2} \).  
- **Reference:**  
  - Papadimitriou and Steiglitz, *Combinatorial Optimization: Algorithms and Complexity*

---

### 3. Dynamic Programming (Held-Karp)

- **Method:** `SolveWithDynamicProgramming()`  
- **Idea:** Uses bitmask DP to remember the shortest path to reach a set of visited cities ending at a particular city.  
- **Complexity:** \( O(n^2 2^n) \) time and space — still exponential but much better than brute force.  
- **Reference:**  
  - Held and Karp, "A Dynamic Programming Approach to Sequencing Problems", J. SIAM, 1962  
  - [Held-Karp Wikipedia](https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm)

---

### 4. Branch and Bound

- **Method:** `SolveWithBranchAndBound()`  
- **Idea:** Uses bounding techniques to prune large parts of the search space by estimating lower bounds on route cost.  
- **Complexity:** Varies — faster than brute force in practice but worst case still factorial.  
- **Reference:**  
  - Lawler, *Branch and Bound Methods*  
  - [TSP Wikipedia: Branch and Bound](https://en.wikipedia.org/wiki/Travelling_salesman_problem#Branch_and_bound)

---

### 5. Nearest Neighbor Heuristic

- **Method:** `SolveWithNearestNeighbor()`  
- **Idea:** Start at a city and repeatedly visit the nearest unvisited city until all visited.  
- **Complexity:** \( O(n^2) \)  
- **Approximation:** Can be arbitrarily bad in worst cases but very fast.  
- **Reference:**  
  - [Nearest Neighbor Heuristic - Wikipedia](https://en.wikipedia.org/wiki/Nearest_neighbor_algorithm)

---

### 6. Greedy Edge Selection Heuristic

- **Method:** `SolveWithGreedyEdgeSelection()`  
- **Idea:** Selects edges in order of increasing length, adding them if they don't form a cycle or cause degree > 2 until a Hamiltonian cycle forms.  
- **Complexity:** \( O(n^2 \log n) \)  
- **Reference:**  
  - [Greedy Algorithm for TSP](https://en.wikipedia.org/wiki/Travelling_salesman_problem#Greedy_algorithm)

---

### 7. Christofides Algorithm

- **Method:** `SolveWithChristofidesAlgorithm()`  
- **Idea:**  
  - Construct minimum spanning tree (MST)  
  - Find minimum-weight perfect matching on odd degree vertices  
  - Combine MST and matching to form Eulerian multigraph  
  - Shortcut Eulerian tour to Hamiltonian cycle  
- **Guarantee:** 3/2 approximation factor for metric TSP  
- **Complexity:** Polynomial  
- **Reference:**  
  - Christofides, *Worst-case analysis of a new heuristic for the travelling salesman problem*, 1976  
  - [Christofides Algorithm - Wikipedia](https://en.wikipedia.org/wiki/Christofides_algorithm)

---

### 8. Simulated Annealing (Metaheuristic)

- **Method:** `SolveWithSimulatedAnnealingAlgorithm()`  
- **Idea:** Probabilistically explores the solution space, occasionally accepting worse solutions to escape local minima, gradually reducing "temperature" to converge.  
- **Complexity:** Depends on parameters, typically much faster than exact methods for large \( n \).  
- **Reference:**  
  - Kirkpatrick, Gelatt, and Vecchi, *Optimization by Simulated Annealing*, Science, 1983  
  - [Simulated Annealing - Wikipedia](https://en.wikipedia.org/wiki/Simulated_annealing)

---

## How to Use (With Header)

1. Clone the repo:

```bash
git clone https://github.com/muhammadIdhamMaarif/Traveling-Salesman-Problem-CPP.git
cd Traveling-Salesman-Problem-CPP
```

2. Compile (requires C++11 or later):

```bash
g++ -std=c++17 main.cpp -o tsp_solver
```

3. Run:
   
```bash
./tsp_solver
```
You will see outputs of all algorithms running on the sample distance matrix.

## How to Use (Paste and Run)

1. Open pnr.cpp
2. Copy all the contents inside
3. Paste in your code editor or online compiler (Example [Programiz](https://www.programiz.com/cpp-programming/online-compiler/))


## Theoretical Notes and Proof Sketches

- Brute force correctness: Enumerates all permutations, so guaranteed to find optimal.
- Held-Karp correctness: Uses optimal substructure and overlapping subproblems, proven via DP theory.
- Branch and Bound correctness: Prunes safely by bounding costs, maintains optimality.
- Christofides approximation: Guaranteed ≤ 1.5 × OPT on metric TSP by MST and matching proofs.
- Simulated Annealing convergence: Probabilistically guaranteed to converge to global optimum as temperature → 0 and infinite iterations, practically a heuristic.

## License
This project is licensed under the MIT License — see LICENSE file for details.

## Contributions
Contributions and improvements are welcome! Feel free to open issues or submit pull requests.


