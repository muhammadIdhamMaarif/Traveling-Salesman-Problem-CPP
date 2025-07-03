void RunAndDisplay(TSP tsp, const std::string& algorithmName, void (TSP::*solveFunc)()) {
    (tsp.*solveFunc)();
    std::cout << "==== " << algorithmName << " ====" << std::endl;
    std::cout << "Distance: " << tsp.GetShortestDistance() << std::endl;
    std::cout << "Route: ";
    for (int city : tsp.GetShortestRoute()) {
        std::cout << city << " ";
    }
    std::cout << "\n\n";
}

int main() {

    std::vector<std::vector<int>> distances = {
        {0, 87, 12, 56, 39, 78},
        {87, 0, 33, 64, 45, 23},
        {12, 33, 0, 97, 81, 66},
        {56, 64, 97, 0, 15, 52},
        {39, 45, 81, 15, 0, 29},
        {78, 23, 66, 52, 29, 0}
    };

    RunAndDisplay(TSP(distances), "Nearest Neighbor", &TSP::SolveWithNearestNeighbor);
    RunAndDisplay(TSP(distances), "Greedy Edge Selection", &TSP::SolveWithGreedyEdgeSelection);
    RunAndDisplay(TSP(distances), "Christofides Algorithm", &TSP::SolveWithChristofidesAlgorithm);
    RunAndDisplay(TSP(distances), "Dynamic Programming (Held-Karp)", &TSP::SolveWithDynamicProgramming);
    RunAndDisplay(TSP(distances), "Branch and Bound", &TSP::SolveWithBranchAndBound);    
    RunAndDisplay(TSP(distances), "Simulated Annealing", &TSP::SolveWithSimulatedAnnealingAlgorithm);
    RunAndDisplay(TSP(distances), "Eliminate Mirror", &TSP::SolveWithEliminatingMirror);            
    RunAndDisplay(TSP(distances), "Fixed First City", &TSP::SolveWithFixedFirstCity);    
    RunAndDisplay(TSP(distances), "Non-Fixed First City", &TSP::SolveWithNonFixedFirstCity);

    return 0;
}
