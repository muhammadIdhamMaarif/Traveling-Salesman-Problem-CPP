#ifndef TSP_HPP
#define TSP_HPP

#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <random>
#include <climits>
#include <functional>
#include <random>
#include <limits>
#include <chrono>
#include <thread>
#include <atomic>
#include <iomanip>
#include <string>

// Kelas TSP untuk menyelesaikan Traveling Salesman Problem (TSP)
class TSP {
private:
    // Matriks jarak antar kota, distanceMatrix[i][j] adalah jarak dari kota i ke kota j
    std::vector<std::vector<int>> distanceMatrix;

    // Menyimpan rute terbaik yang ditemukan selama proses pencarian solusi
    std::vector<int> bestRoute;

    // Menyimpan jarak minimum yang ditemukan selama pencarian
    int minDistance;

    // Jumlah kota yang harus dikunjungi
    int numCities;

    // Vector untuk menandai kota mana saja yang sudah dikunjungi (digunakan di beberapa algoritma)
    std::vector<bool> visited;

    // Tabel DP (Dynamic Programming)
    // dp[mask][city] menyimpan jarak minimum untuk mencapai kota 'city' dengan subset kota yang sudah dikunjungi 'mask'
    std::vector<std::vector<int>> dp;

    // Tabel untuk menyimpan parent node, berguna untuk merekonstruksi jalur rute terbaik dari hasil DP
    std::vector<std::vector<int>> parent;

    // Flag untuk menandai apakah solusi menggunakan metode DP sudah selesai dihitung
    bool solvedWithDP = false;

    // Mask bit untuk menandai semua kota sudah dikunjungi (semua bit bernilai 1 sebanyak numCities)
    int ALL_VISITED;

    // Fungsi untuk menghitung total jarak dari sebuah rute yang diberikan
    // Route di sini adalah urutan kota yang dikunjungi
    int CalculateDistance(const std::vector<int>& route) {
        int total = 0;
        // Jumlahkan jarak antar kota sepanjang rute
        for (int i = 0; i < numCities - 1; ++i) {
            total += distanceMatrix[route[i]][route[i + 1]];
        }
        // Tambahkan jarak kembali ke kota awal (rute siklik)
        total += distanceMatrix[route[numCities - 1]][route[0]];
        return total;
    }

    // Algoritma backtracking untuk mencari rute dengan jarak minimum secara brute force
    // route: vector yang menyimpan urutan kota saat ini
    // startIndex: posisi kota yang akan kita swap untuk mengeksplorasi semua kemungkinan rute
    void BacktrackAlgorithm(std::vector<int>& route, int startIndex) {
        // Jika sudah menempatkan semua kota dalam rute (base case)
        if (startIndex == numCities) {
            int currentDistance = CalculateDistance(route);
            // Jika rute ini lebih baik dari rute terbaik yang sudah ada, update bestRoute dan minDistance
            if (currentDistance < minDistance) {
                minDistance = currentDistance;
                bestRoute = route;
            }
            return;
        }

        // Loop untuk mencoba swap kota di posisi startIndex dengan kota lainnya dari startIndex ke akhir
        for (int i = startIndex; i < numCities; ++i) {
            std::swap(route[startIndex], route[i]); // tukar posisi kota
            BacktrackAlgorithm(route, startIndex + 1); // lanjut ke posisi berikutnya
            std::swap(route[startIndex], route[i]); // kembalikan ke posisi awal (backtrack)
        }
    }
    
    // Fungsi untuk mengeliminasi rute yang merupakan cerminan (mirror) dari rute lain agar tidak dihitung duplikat
    // route: vektor urutan kota, startIndex: posisi mulai swap seperti backtracking biasa
    void EliminateMirror(std::vector<int>& route, int startIndex) {
        if (startIndex == numCities) {
            // Mirror elimination condition:
            // Misal rute: 0 - 1 - 2 - 3 - 0
            // Jika elemen kedua (index 1) lebih besar dari elemen terakhir sebelum kembali ke 0,
            // berarti rute ini cerminan (mirror) dari rute lain yg sudah dihitung, jadi skip
            if (route[1] > route[numCities - 1]) return;

            int currentDistance = CalculateDistance(route);
            if (currentDistance < minDistance) {
                minDistance = currentDistance;
                bestRoute = route;
            }
            return;
        }

        // Backtracking seperti biasa untuk coba semua kemungkinan
        for (int i = startIndex; i < numCities; ++i) {
            std::swap(route[startIndex], route[i]);
            EliminateMirror(route, startIndex + 1);
            std::swap(route[startIndex], route[i]); // backtrack
        }
    }

    // Recursive helper function untuk algoritma Held-Karp (DP untuk TSP)
    // mask: bitmask untuk menandai kota yang sudah dikunjungi,
    // pos: kota saat ini
    int HeldKarpAlgorithm(int mask, int pos) {
        if (mask == ALL_VISITED) {
            // Jika semua kota sudah dikunjungi, kembalikan jarak dari posisi sekarang ke kota awal (0)
            return distanceMatrix[pos][0];
        }
        if (dp[mask][pos] != -1) {
            // Jika sudah dihitung sebelumnya (memoization), langsung return hasilnya
            return dp[mask][pos];
        }

        int ans = INT_MAX;
        int chosenCity = -1;

        // Coba semua kota yang belum dikunjungi (bit di mask = 0)
        for (int city = 0; city < numCities; city++) {
            if ((mask & (1 << city)) == 0) { // kota belum dikunjungi
                // Hitung biaya total jika pergi ke city selanjutnya
                int newCost = distanceMatrix[pos][city] + HeldKarpAlgorithm(mask | (1 << city), city);
                if (newCost < ans) {
                    ans = newCost;
                    chosenCity = city;
                }
            }
        }

        // Simpan parent untuk rekonstruksi rute nanti
        parent[mask][pos] = chosenCity;
        // Simpan hasil minimal cost ke dp dan return
        return dp[mask][pos] = ans;
    }

    // Fungsi untuk menghitung batas bawah (lower bound) dari solusi saat ini
    // Digunakan pada algoritma Branch and Bound untuk pruning
    int CalculateLowerBound(int currentCost) {
        int bound = currentCost;
        for (int city = 0; city < numCities; ++city) {
            if (!visited[city]) {
                int minEdge = INT_MAX;
                // Cari jarak minimum dari kota yang belum dikunjungi ke kota lain yang belum dikunjungi atau ke kota awal
                for (int nextCity = 0; nextCity < numCities; ++nextCity) {
                    if (city != nextCity && (!visited[nextCity] || nextCity == 0)) {
                        minEdge = std::min(minEdge, distanceMatrix[city][nextCity]);
                    }
                }
                // Tambahkan jarak minimum tersebut ke bound
                bound += (minEdge == INT_MAX) ? 0 : minEdge;
            }
        }
        return bound;
    }

    // Algoritma Branch and Bound untuk mencari solusi TSP dengan pruning menggunakan lower bound
    // currentCity: kota saat ini, count: jumlah kota yang sudah dikunjungi,
    // currentCost: total jarak yang sudah ditempuh,
    // currentRoute: rute yang sedang dibangun
    void BranchAndBound(int currentCity, int count, int currentCost, std::vector<int>& currentRoute) {
        if (count == numCities) {
            // Jika semua kota sudah dikunjungi, tambahkan jarak kembali ke kota awal
            int totalCost = currentCost + distanceMatrix[currentCity][0];
            if (totalCost < minDistance) {
                minDistance = totalCost;
                bestRoute = currentRoute;
            }
            return;
        }

        // Hitung batas bawah biaya solusi untuk rute saat ini
        int bound = CalculateLowerBound(currentCost);
        // Jika bound lebih besar atau sama dengan minDistance saat ini, cabang ini tidak akan menghasilkan solusi lebih baik, jadi prune
        if (bound >= minDistance) return;

        // Coba semua kota yang belum dikunjungi
        for (int nextCity = 0; nextCity < numCities; ++nextCity) {
            if (!visited[nextCity]) {
                visited[nextCity] = true;
                currentRoute.push_back(nextCity);

                // Rekursif lanjutkan ke kota berikutnya
                BranchAndBound(nextCity, count + 1, currentCost + distanceMatrix[currentCity][nextCity], currentRoute);

                // Backtrack: kembalikan kondisi seperti semula
                visited[nextCity] = false;
                currentRoute.pop_back();
            }
        }
    }

    // Algoritma Nearest Neighbor untuk TSP
    // Mulai dari kota startCity (default 0), selalu pilih kota terdekat yang belum dikunjungi
    void NearestNeighbor(int startCity = 0) {
        // Tandai semua kota belum dikunjungi
        std::vector<bool> visited(numCities, false);
        bestRoute.clear();  // Reset rute terbaik
        minDistance = 0;    // Reset jarak minimum

        int current = startCity;
        visited[current] = true;       // Tandai kota awal sudah dikunjungi
        bestRoute.push_back(current);  // Tambahkan ke rute

        // Loop untuk memilih kota berikutnya sebanyak numCities - 1
        for (int step = 1; step < numCities; ++step) {
            int nearest = -1;          // Menyimpan kota terdekat
            int minDist = INT_MAX;     // Jarak terkecil sementara

            // Cari kota terdekat yang belum dikunjungi
            for (int city = 0; city < numCities; ++city) {
                if (!visited[city] && distanceMatrix[current][city] < minDist) {
                    minDist = distanceMatrix[current][city];
                    nearest = city;
                }
            }

            if (nearest != -1) {
                visited[nearest] = true;         // Tandai sudah dikunjungi
                bestRoute.push_back(nearest);    // Tambahkan ke rute terbaik
                minDistance += minDist;          // Tambah jarak tempuh
                current = nearest;               // Update kota saat ini
            }
        }

        // Kembali ke kota awal untuk membentuk siklus
        minDistance += distanceMatrix[current][startCity];
        bestRoute.push_back(startCity);
    }

    // Algoritma Greedy Edge Selection untuk TSP (mirip dengan Minimum Spanning Tree tapi membentuk siklus Hamiltonian)
    // Memilih edge dengan biaya terkecil sambil menghindari siklus prematur dan menjaga derajat node <= 2
    void GreedyEdgeSelection() {
        bestRoute.clear();  // Reset rute terbaik
        minDistance = 0;    // Reset jarak minimum

        struct Edge {
            int from, to, cost;
            bool operator<(const Edge& other) const {
                return cost < other.cost;  // Sorting ascending berdasarkan cost
            }
        };

        std::vector<Edge> edges;
        // Buat list semua edge tanpa duplikat (karena graph tidak berarah)
        for (int i = 0; i < numCities; ++i) {
            for (int j = i + 1; j < numCities; ++j) {
                edges.push_back({i, j, distanceMatrix[i][j]});
            }
        }

        // Sort edge berdasarkan cost terkecil
        std::sort(edges.begin(), edges.end());

        // Struktur Union-Find untuk deteksi siklus saat menambah edge
        std::vector<int> parent(numCities);
        for (int i = 0; i < numCities; ++i) parent[i] = i;

        // Fungsi find untuk Union-Find
        std::function<int(int)> find = [&](int x) -> int {
            if (parent[x] != x)
                parent[x] = find(parent[x]);  // Path compression
            return parent[x];
        };

        // Fungsi unite untuk menggabungkan dua set
        auto unite = [&](int a, int b) {
            parent[find(a)] = find(b);
        };

        // List adjacency untuk menyimpan graph yang terbentuk
        std::vector<std::vector<int>> adj(numCities);
        int edgeCount = 0;

        // Loop untuk menambahkan edge terkecil sambil menjaga aturan:
        // - Derajat setiap node <= 2 (untuk siklus Hamiltonian)
        // - Tidak membuat siklus prematur kecuali siklus penuh
        for (const auto& edge : edges) {
            int u = edge.from;
            int v = edge.to;

            // Jika derajat node sudah 2, skip edge ini
            if (adj[u].size() >= 2 || adj[v].size() >= 2)
                continue;

            // Cek apakah menambahkan edge ini akan membuat siklus prematur
            if (find(u) == find(v)) {
                // Jika sudah hampir lengkap (edgeCount == numCities - 1), boleh tutup siklus
                if (edgeCount != numCities - 1) continue;
            }

            // Tambahkan edge ke adjacency list
            adj[u].push_back(v);
            adj[v].push_back(u);
            unite(u, v);          // Gabungkan set di Union-Find
            minDistance += edge.cost; // Tambah biaya edge ke total jarak
            edgeCount++;

            // Jika sudah cukup edge untuk membentuk siklus, selesai
            if (edgeCount == numCities) break;
        }

        // Validasi: Pastikan setiap node punya derajat tepat 2 (harus untuk siklus Hamiltonian)
        for (int i = 0; i < numCities; ++i) {
            if (adj[i].size() != 2) {
                std::cerr << "[ERROR] Invalid tour: vertex " << i << " has degree " << adj[i].size() << ".\n";
                bestRoute.clear();
                minDistance = INT_MAX;
                return;
            }
        }

        // Rekonstruksi siklus Hamiltonian dari adjacency list
        std::vector<bool> visited(numCities, false);
        std::vector<int> tour;
        int current = 0;
        tour.push_back(current);
        visited[current] = true;

        while (tour.size() < numCities) {
            int next = -1;
            // Cari tetangga yang belum dikunjungi
            for (int neighbor : adj[current]) {
                if (!visited[neighbor]) {
                    next = neighbor;
                    break;
                }
            }
            if (next == -1) {
                std::cerr << "[ERROR] Failed to complete tour. Stuck at node " << current << ".\n";
                bestRoute.clear();
                minDistance = INT_MAX;
                return;
            }
            tour.push_back(next);
            visited[next] = true;
            current = next;
        }

        // Tutup siklus dengan kembali ke kota awal
        tour.push_back(tour[0]);
        bestRoute = tour;

        // Hitung ulang jarak total dari rute akhir
        minDistance = 0;
        for (int i = 0; i < numCities; ++i) {
            minDistance += distanceMatrix[bestRoute[i]][bestRoute[i + 1]];
        }
    }

    // Fungsi Christofides untuk mencari solusi TSP dengan pendekatan approx 3/2 optimal
    // 5 langkah utama:
    // 1. Bangun MST (Minimum Spanning Tree) menggunakan Prim's Algorithm
    // 2. Cari vertex dengan derajat ganjil di MST
    // 3. Lakukan perfect matching minimum weight pada vertex ganjil
    // 4. Buat Eulerian multigraph dan cari Eulerian circuit (Hierholzer’s Algorithm)
    // 5. Buat Hamiltonian cycle dengan shortcut pada Eulerian path (menghilangkan node duplikat)
    void Christofides() {
        // Step 1: Build MST dengan Prim's algorithm
        std::vector<int> parent(numCities, -1);  // Menyimpan parent node dalam MST
        std::vector<int> key(numCities, INT_MAX); // Cost minimum untuk menghubungkan node ke MST
        std::vector<bool> inMST(numCities, false); // Penanda node sudah termasuk MST
        key[0] = 0; // Mulai dari node 0

        for (int count = 0; count < numCities - 1; ++count) {
            // Cari node u dengan key terkecil yang belum dimasukkan MST
            int u = -1;
            for (int i = 0; i < numCities; ++i) {
                if (!inMST[i] && (u == -1 || key[i] < key[u]))
                    u = i;
            }

            if (u == -1) {
                std::cerr << "[ERROR] Prim's algorithm failed to select a valid node.\n";
                return;
            }

            inMST[u] = true; // Masukkan node u ke MST

            // Update key dan parent node yang bertetangga dengan u jika jarak lebih kecil
            for (int v = 0; v < numCities; ++v) {
                if (distanceMatrix[u][v] && !inMST[v] && distanceMatrix[u][v] < key[v]) {
                    key[v] = distanceMatrix[u][v];
                    parent[v] = u;
                }
            }
        }

        // Bangun adjacency list MST dari parent array
        std::vector<std::vector<int>> mst(numCities);
        for (int i = 1; i < numCities; ++i) {
            if (parent[i] != -1) {
                mst[i].push_back(parent[i]);
                mst[parent[i]].push_back(i);
            }
        }

        // Step 2: Cari vertex dengan derajat ganjil di MST (harus genap jumlahnya)
        std::vector<int> oddVertices;
        for (int i = 0; i < numCities; ++i) {
            if (mst[i].size() % 2 != 0) {
                oddVertices.push_back(i);
            }
        }

        if (oddVertices.size() % 2 != 0) {
            std::cerr << "[ERROR] Odd number of odd-degree vertices. Cannot form perfect matching.\n";
            return;
        }

        // Step 3: Greedy Minimum Weight Perfect Matching pada vertex ganjil
        while (!oddVertices.empty()) {
            int u = oddVertices.back();
            oddVertices.pop_back();

            int minCost = INT_MAX, matchIndex = -1;
            for (int i = 0; i < oddVertices.size(); ++i) {
                int v = oddVertices[i];
                if (distanceMatrix[u][v] < minCost) {
                    minCost = distanceMatrix[u][v];
                    matchIndex = i;
                }
            }

            if (matchIndex == -1) {
                std::cerr << "[ERROR] No match found for vertex " << u << "\n";
                return;
            }

            int v = oddVertices[matchIndex];
            // Tambahkan edge matching ke MST
            mst[u].push_back(v);
            mst[v].push_back(u);
            oddVertices.erase(oddVertices.begin() + matchIndex);
        }

        // Step 4: Buat Eulerian multigraph untuk mencari Eulerian circuit
        std::vector<std::multiset<int>> graph(numCities);
        for (int u = 0; u < numCities; ++u) {
            for (int v : mst[u]) {
                graph[u].insert(v);
            }
        }

        // Validasi: semua vertex harus memiliki derajat genap (karena Eulerian graph)
        for (int i = 0; i < numCities; ++i) {
            if (graph[i].size() % 2 != 0) {
                std::cerr << "[ERROR] Vertex " << i << " has odd degree in final graph.\n";
                return;
            }
        }

        // Hierholzer's Algorithm: DFS untuk mencari Eulerian circuit
        std::vector<int> eulerPath;
        std::function<void(int)> dfs = [&](int u) {
            while (!graph[u].empty()) {
                int v = *graph[u].begin();
                graph[u].erase(graph[u].begin());
                graph[v].erase(graph[v].find(u));
                dfs(v);
            }
            eulerPath.push_back(u);
        };

        dfs(0);
        std::reverse(eulerPath.begin(), eulerPath.end());

        // Step 5: Buat Hamiltonian path dengan menghilangkan kunjungan ulang pada node yang sama
        std::unordered_set<int> visited;
        bestRoute.clear();
        minDistance = 0;

        for (int city : eulerPath) {
            if (!visited.count(city)) {
                if (!bestRoute.empty()) {
                    minDistance += distanceMatrix[bestRoute.back()][city];
                }
                visited.insert(city);
                bestRoute.push_back(city);
            }
        }

        // Kembali ke kota awal supaya siklus lengkap
        if (!bestRoute.empty()) {
            minDistance += distanceMatrix[bestRoute.back()][bestRoute[0]];
            bestRoute.push_back(bestRoute[0]);
        }
    }

    // Fungsi Simulated Annealing untuk menyelesaikan TSP secara heuristik dengan pendekatan probabilistik
    // Parameters:
    // iterations: jumlah iterasi maksimum
    // initialTemp: suhu awal (temperatur)
    // coolingRate: laju pendinginan suhu (nilai antara 0 dan 1)
    void SimulatedAnnealing(int iterations = 100000, double initialTemp = 10000.0, double coolingRate = 0.9995) {
        // Inisialisasi rute awal urut dari 0 sampai numCities-1
        std::vector<int> currentRoute(numCities);
        std::iota(currentRoute.begin(), currentRoute.end(), 0); // generate 0,1,2,...,n-1

        // Setup random generator dan distribusi
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist(0, numCities - 1); // untuk memilih indeks acak
        std::uniform_real_distribution<> prob(0.0, 1.0);        // untuk probabilitas penerimaan solusi buruk

        // Hitung jarak rute awal
        int currentDistance = CalculateDistance(currentRoute);
        bestRoute = currentRoute;
        minDistance = currentDistance;

        double temperature = initialTemp; // set suhu awal

        for (int iter = 0; iter < iterations; ++iter) {
            std::vector<int> newRoute = currentRoute;

            // Terapkan 2-opt swap: membalik urutan segmen antara dua indeks acak i dan j
            int i = dist(gen);
            int j = dist(gen);
            if (i > j) std::swap(i, j);
            std::reverse(newRoute.begin() + i, newRoute.begin() + j + 1);

            int newDistance = CalculateDistance(newRoute);
            int delta = newDistance - currentDistance; // perubahan jarak

            // Acceptance criteria:
            // Jika solusi baru lebih baik (delta < 0), terima langsung
            // Jika lebih buruk, terima dengan probabilitas exp(-delta/temperature)
            if (delta < 0 || prob(gen) < std::exp(-delta / temperature)) {
                currentRoute = newRoute;
                currentDistance = newDistance;

                // Update rute dan jarak terbaik jika solusi baru lebih baik
                if (newDistance < minDistance) {
                    bestRoute = newRoute;
                    minDistance = newDistance;
                }
            }

            // Turunkan suhu sesuai cooling rate (pendinginan)
            temperature *= coolingRate;

            // Optional early stop: jika suhu sudah sangat kecil, hentikan iterasi
            if (temperature < 1e-6) break;
        }

        // Tambahkan kota awal ke akhir rute supaya membentuk siklus lengkap
        bestRoute.push_back(bestRoute[0]);
    }

public:
    // Konstruktor kelas TSP menerima matriks jarak (distanceMatrix)
    // Inisialisasi variabel dan struktur data penting:
    // dp dan parent untuk dynamic programming (Held-Karp),
    // visited untuk penanda kota yang sudah dikunjungi,
    // ALL_VISITED bitmask untuk menandai semua kota telah dikunjungi
    TSP(const std::vector<std::vector<int>>& distances)
        : distanceMatrix(distances),
          numCities(distances.size()),
          minDistance(INT_MAX),
          bestRoute(distances.size()),
          visited(distances.size(), false)
    {
        int size = 1 << numCities;  // 2^numCities untuk bitmasking
        dp.assign(size, std::vector<int>(numCities, -1));       // Inisialisasi DP dengan -1 (belum dihitung)
        parent.assign(size, std::vector<int>(numCities, -1));   // Untuk rekonstruksi jalur
        ALL_VISITED = (1 << numCities) - 1;                     // Bitmask semua kota sudah dikunjungi (misal untuk 4 kota: 1111b = 15)
    }

    // Solusi TSP dengan backtracking tanpa mengunci kota awal
    // Mencoba semua permutasi penuh
    void SolveWithNonFixedFirstCity() {
        std::vector<int> route(numCities);
        for (int i = 0; i < numCities; ++i) {
            route[i] = i;
        }
        BacktrackAlgorithm(route, 0);
    }

    // Solusi TSP dengan backtracking, tapi mengunci kota pertama tetap di posisi 0
    // Efektif mengurangi jumlah permutasi yang dicek karena simetri rute
    void SolveWithFixedFirstCity() {
        std::vector<int> route(numCities);
        for (int i = 0; i < numCities; ++i) {
            route[i] = i;
        }
        BacktrackAlgorithm(route, 1);
    }

    // Solusi dengan backtracking dan eliminasi rute yang merupakan mirror (simetri terbalik)
    void SolveWithEliminatingMirror() {
        std::vector<int> route(numCities);
        for (int i = 0; i < numCities; ++i) {
            route[i] = i;
        }
        EliminateMirror(route, 1);
    }

    // Solusi dengan Dynamic Programming menggunakan algoritma Held-Karp
    void SolveWithDynamicProgramming() {
        HeldKarpAlgorithm(1, 0);  // Mulai dari kota 0 dengan mask sudah ter-visit
        solvedWithDP = true;
    }

    // Solusi dengan algoritma Branch and Bound untuk pruning search space
    void SolveWithBranchAndBound() {
        visited[0] = true;              // Kota awal dikunci visited
        std::vector<int> currentRoute = {0};
        BranchAndBound(0, 1, 0, currentRoute);
    }

    // Solusi dengan heuristik Nearest Neighbor
    void SolveWithNearestNeighbor() {
        NearestNeighbor();
    }

    // Solusi dengan metode Greedy Edge Selection (membangun siklus Hamiltonian dari edge murah)
    void SolveWithGreedyEdgeSelection() {
        GreedyEdgeSelection();
    }

    // Solusi menggunakan algoritma Christofides (aproksimasi dengan jaminan 3/2 dari optimal)
    void SolveWithChristofidesAlgorithm() {
        Christofides();
    }

    // Solusi dengan algoritma Simulated Annealing (heuristik berbasis probabilistik)
    void SolveWithSimulatedAnnealingAlgorithm() {
        SimulatedAnnealing();
    }

    // Mengembalikan jarak rute terpendek yang sudah ditemukan
    // Jika menggunakan DP, ambil dari tabel dp
    int GetShortestDistance() const {
        if (!solvedWithDP) return minDistance;
        else return dp[1][0];
    }

    // Mengembalikan rute terpendek yang sudah ditemukan
    // Jika menggunakan DP, rekonstruksi rute dari parent table
    std::vector<int> GetShortestRoute() const {
        if (!solvedWithDP) return bestRoute;
        else {
            std::vector<int> route;
            int mask = 1;  // Mulai dari kota 0 (visited)
            int pos = 0;   // Posisi awal kota 0

            route.push_back(pos);

            // Rekonstruksi rute dengan mengikuti parent table hingga selesai
            while (true) {
                int nextCity = parent[mask][pos];
                if (nextCity == -1) break;
                route.push_back(nextCity);
                mask |= (1 << nextCity);
                pos = nextCity;
            }

            return route;
        }
    }



    /*
        This code does not implement ILP/MIP-based TSP with MTZ constraints due to its complexity in pure C++.
        For optimal solutions using ILP, it is recommended to use solver libraries like Google OR-Tools
        that provide built-in support for TSP modeling and subtour elimination.
    */

    /*
        Note:
        Solving the Traveling Salesman Problem (TSP) using Integer Linear Programming (ILP) or Mixed-Integer Programming (MIP)
        with Miller–Tucker–Zemlin (MTZ) subtour elimination constraints is a mathematically optimal and efficient approach
        for medium-sized problems (up to ~100 cities).

        However, implementing ILP/MIP manually in pure C++ is highly complex and inefficient.
        It requires constructing a system of linear equations, managing binary integer variables,
        and integrating an optimization algorithm such as branch-and-bound or simplex,
        which is not trivial without external solver libraries.

        Therefore, this approach is not implemented directly here to maintain simplicity and code readability.

        For practical and efficient usage, it is highly recommended to use third-party libraries like **Google OR-Tools**,
        which offer a powerful modeling interface and integrate advanced MIP solvers that automatically handle
        subtour elimination, path constraints, and optimization.

        Summary:
        ILP/MIP is intentionally not implemented in native C++ here to avoid excessive complexity.
        Use OR-Tools or similar solvers if you need an efficient and production-ready ILP-based TSP solution.
    */

};

#endif // TSP_HPP