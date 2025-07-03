# Traveling Salesman Problem (TSP) Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## TL;DR
1. Copy all text in `pnr.cpp`
2. Open [Programiz](https://www.programiz.com/cpp-programming/online-compiler/)
3. Paste the code
4. Run

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

## Kompleksitas Algoritma

Berikut adalah analisis singkat kompleksitas waktu dan ruang untuk berbagai metode penyelesaian TSP yang ada di kelas ini:

| Algoritma                        | Time Complexity      | Space Complexity     | Deskripsi Singkat                                                                 |
|---------------------------------|---------------------|---------------------|-----------------------------------------------------------------------------------|
| `SolveWithNonFixedFirstCity()`   | O(n!)               | O(n)                | Mengeksplorasi semua permutasi rute, sangat mahal untuk n besar.                  |
| `SolveWithFixedFirstCity()`      | O((n-1)!)           | O(n)                | Mengurangi permutasi dengan mengunci kota awal, mempercepat sekitar faktor n.     |
| `SolveWithEliminatingMirror()`   | O((n-1)! / 2)       | O(n)                | Menghilangkan rute yang simetris terbalik, mengurangi permutasi sekitar setengahnya.|
| `SolveWithDynamicProgramming()`  | O(n² * 2^n)         | O(n * 2^n)          | Algoritma DP bitmask efisien, cocok untuk n sampai sekitar 20.                    |
| `SolveWithBranchAndBound()`      | Worst: O(n!), Avg: ? | O(n)                | Pruning mengurangi cabang pencarian, sangat tergantung pada input dan bound.      |
| `SolveWithNearestNeighbor()`     | O(n²)               | O(n)                | Heuristik cepat, hasil tidak selalu optimal.                                     |
| `SolveWithGreedyEdgeSelection()` | O(n² log n)         | O(n²)               | Membangun siklus Hamiltonian dari edge murah, menggunakan sorting dan Union-Find.|
| `SolveWithChristofidesAlgorithm()` | O(n³)            | O(n²)               | Aproksimasi dengan jaminan 1.5x solusi optimal, memakai MST dan matching.         |
| `SolveWithSimulatedAnnealingAlgorithm()` | O(iterations * n) | O(n)               | Heuristik probabilistik, iterasi tergantung parameter, cocok untuk n besar.       |

---

**Keterangan:**
- `n` = jumlah kota (numCities)
- `iterations` = parameter jumlah iterasi pada algoritma Simulated Annealing
- Pruning di Branch and Bound membuat waktu rata-rata sulit diprediksi secara pasti

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

## How to Use (Local)

1. Clone the repo:

```bash
git clone https://github.com/muhammadIdhamMaarif/Traveling-Salesman-Problem-CPP.git
cd Traveling-Salesman-Problem-CPP
```

2. Compile (requires C++11 or later):

```bash
g++ -std=c++17 pnr.cpp -o pnr
```

3. Run:
   
```bash
./pnr
```
You will see outputs of all algorithms running on the sample distance matrix.

## How to Use (Paste and Run)

1. Open `pnr.cpp`
2. Copy all the contents inside
3. Paste in your code editor or online compiler (Example [Programiz](https://www.programiz.com/cpp-programming/online-compiler/))

## How to Use (Header)

1. Copy all the contents inside `tsp.hpp`
2. Paste into your workspace
3. Use by using
```cpp
  #include "tsp.hpp"
  // your code
```


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

## References

### Core Algorithm References

#### Backtracking and Dynamic Programming
- AbdelkadirSellahi. (2024, March 2). *The-Travelling-Salesperson-Problem-TSP* [Source code]. GitHub. https://github.com/AbdelkadirSellahi/The-Travelling-Salesperson-Problem-TSP  
- ShunxiXXX. (2025, May 27). *TSP_Algorithms_Comparison: This project compares the running time of three TSP problem algorithms and visualizes them* [Source code]. GitHub. https://github.com/ShunxiXXX/TSP_Algorithms_Comparison  
- ishanjogalekar. (2020, November 25). *TSP using Dynamic Programming, Held-Karp Algorithm* [Source code]. GitHub. https://github.com/ishanjogalekar/TSP-using-Dynamic-Programming-  

#### Nearest Neighbor and Heuristics
- elifBalci. (2020, June 13). *travelling_salesman_problem_solver* [Source code]. GitHub. https://github.com/elifBalci/travelling_salesman_problem_solver  
- m3hdi-i. (2021, October 29). *tsp-with-nn: TSP with Nearest Neighbor algorithm* [Source code]. GitHub. https://github.com/m3hdi-i/tsp-with-nn  

#### Christofides Algorithm
- prakharverma. (2019, December 17). *Christofides-Algorithm* [Source code]. GitHub. https://github.com/prakharverma/Christofides-Algorithm  

#### Simulated Annealing
- tomekrzymyszkiewicz. (2021, November 20). *TSP-simulated-annealing* [Source code]. GitHub. https://github.com/tomekrzymyszkiewicz/TSP-simulated-annealing

#### Ant Colony Optimization
- nishnash54. (2019, March 31). *TSP_ACO* [Source code]. GitHub. https://github.com/nishnash54/TSP_ACO  
- EvanOman. (2015, August 13). *AntColonyOptimization-TSP* [Source code]. GitHub. https://github.com/EvanOman/AntColonyOptimization-TSP  

#### Genetic Algorithms
- GiovanniSorice. (2020, June 18). *TSPGeneticAlgorithm* [Source code]. GitHub. https://github.com/GiovanniSorice/TSPGeneticAlgorithm  
- emre-kocyigit. (2022, October 3). *genetic_algorithm_tsp* [Source code]. GitHub. https://github.com/emre-kocyigit/genetic_algorithm_tsp
  
#### Comprehensive Collections
- Piero24. (2024, February 27). *TSP_Optimization: A list of Heuristics, Metaheuristics and Matheuristic algorithms for solve the TSP* [Source code]. GitHub. https://github.com/Piero24/TSP_Optimization  
- FernandoSchett. (2024, April 4). *tsp_experiments* [Source code]. GitHub. https://github.com/FernandoSchett/tsp_experiments/  
  
### Educational Programming Resources

#### GeeksforGeeks
- GeeksforGeeks. (2024, November 26). *Travelling Salesman Problem using Dynamic Programming*. https://www.geeksforgeeks.org/dsa/travelling-salesman-problem-using-dynamic-programming/  
- GeeksforGeeks. (2023, April 30). *Traveling Salesman Problem using Branch And Bound*. https://www.geeksforgeeks.org/dsa/traveling-salesman-problem-using-branch-and-bound-2/  
- GeeksforGeeks. (2024, November 26). *Traveling Salesman Problem (TSP) Implementation*. https://www.geeksforgeeks.org/dsa/traveling-salesman-problem-tsp-implementation/  

#### W3Schools and TutorialsPoint
- W3Schools. (2025). *DSA The Traveling Salesman Problem*. https://www.w3schools.com/dsa/dsa_ref_traveling_salesman.php  
- TutorialsPoint. (n.d.). *Travelling Salesman Problem (Greedy Approach)*. https://www.tutorialspoint.com/data_structures_algorithms/travelling_salesman_problem.htm  

#### Academic and Technical Platforms
- RosettaCode. (2025, June 21). *Held–Karp algorithm*. https://rosettacode.org/wiki/Held%E2%80%93Karp_algorithm  
- Lee, S. (2025, June 11). *Dynamic Programming for TSP: A Step-by-Step Guide*. Number Analytics. https://www.numberanalytics.com/blog/dynamic-programming-for-tsp-step-by-step  

#### Python, NetworkX, and Related Tools

- NetworkX Development Team. (2025, May 29). *greedy_tsp — NetworkX 3.5 documentation*. https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.approximation.traveling_salesman.greedy_tsp.html  
- NetworkX Development Team. (2025, May 29). *simulated_annealing_tsp — NetworkX 3.5 documentation*. https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.approximation.traveling_salesman.simulated_annealing_tsp.html  
- NetworkX Development Team. (2025, May 29). *traveling_salesman_problem — NetworkX 3.5 documentation*. https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.approximation.traveling_salesman.traveling_salesman_problem.html  
- Christofides Package Maintainer. (2017, April 16). *Christofides*. PyPI. https://pypi.org/project/Christofides/  

#### MATLAB & MathWorks
- MathWorks. (2025). *Traveling Salesman Problem: Problem-Based*. https://www.mathworks.com/help/optim/ug/traveling-salesman-problem-based.html  
- MathWorks. (2025). *Traveling Salesman Problem: Solver-Based*. https://www.mathworks.com/help/optim/ug/travelling-salesman-problem.html  

#### R & CRAN
- Hahsler, M., & Hornik, K. (2025, May 27). *TSP: Infrastructure for the Traveling Salesperson Problem* [R package]. CRAN. https://cran.r-project.org/package=TSP  
- R Core Team. (2005). *Genetic Algorithm for the TSP*. gor package documentation. https://search.r-project.org/CRAN/refmans/gor/html/search_tour_genetic.html  

#### Programming Blogs and Tutorials

- Gazda, M. (n.d.). *TSP algorithms: 2-opt, 3-opt in python*. Matej Gazda Blog. http://matejgazda.com/tsp-algorithms-2-opt-3-opt-in-python/  
- The Coding Train. (2016, August 24). *Coding Challenge #35.1: Traveling Salesperson* [Video]. YouTube. https://www.youtube.com/watch?v=BAejnwN4Ccw  
- vlogize. (2025, January 27). *Solving the Traveling Salesman Problem with Dynamic Programming in C* [Video]. YouTube. https://www.youtube.com/watch?v=TWnJFgatMRU  

#### Competitive Programming Platforms

- zhongquan789. (2021, January 4). *travelling salesman problem - leetcode solutions*. GitBook. https://zhongquan789.gitbook.io/leetcode/unsensored/travelling-salesman-problem  
- KareemTahaAbdelfattah. (n.d.). *1713A - Traveling Salesman Problem.cpp* [Source code]. GitHub. https://github.com/KareemTahaAbdelfattah/Codeforces-Solutions/blob/main/1713A%20-%20Traveling%20Salesman%20Problem.cpp

### Others
- Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2009). *Introduction to algorithms* (3rd ed.). MIT Press.
- Papadimitriou, C. H., & Steiglitz, K. (1998). *Combinatorial optimization: Algorithms and complexity*. Dover Publications.
- Saller, S., Hougardy, S., & Vygen, J. (2024). Approximation algorithms for traveling salesman problems: A systematic review. *Mathematical Programming*, *204*(1), 1–89.



