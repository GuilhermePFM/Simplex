#include "simplex/simplex.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// Simple LCG PRNG for reproducible random LPs
struct LCG {
    uint64_t state;
    explicit LCG(uint64_t seed = 42) : state(seed) {}

    // Returns a uniform double in [0, 1)
    double next()
    {
        // LCG parameters from Knuth / Numerical Recipes
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(state >> 11) / static_cast<double>(1ULL << 53);
    }
};

// Generate a random LP with m constraints and n total variables (n > m).
// Structure mirrors the Julia bench.jl random_lp: last m columns are an
// identity (slack) block, objective is only over the first (n-m) columns.
static LinearProgram random_lp(int m, int n, uint64_t seed = 42)
{
    LCG rng(seed);
    int n_orig = n - m;

    Eigen::MatrixXd A(m, n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n_orig; ++j)
            A(i, j) = rng.next();
    A.rightCols(m) = Eigen::MatrixXd::Identity(m, m);

    Eigen::VectorXd b(m);
    for (int i = 0; i < m; ++i)
        b[i] = rng.next() + 0.1;   // strictly positive

    Eigen::VectorXd c(n);
    for (int i = 0; i < n_orig; ++i) c[i] = rng.next();
    c.tail(m).setZero();

    return LinearProgram(A, b, c);
}

using Clock    = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double, std::nano>;

// Return median of sorted durations in nanoseconds
static double median_ns(std::vector<double>& v)
{
    std::sort(v.begin(), v.end());
    std::size_t n = v.size();
    return (n % 2 == 0) ? (v[n / 2 - 1] + v[n / 2]) / 2.0 : v[n / 2];
}

struct BenchCase {
    std::string    name;
    LinearProgram  lp;
    BasisState     init_basis;  // empty => two-phase solve
    bool           phase2_only;
};

static void run_bench(const BenchCase& bc, int runs = 10)
{
    std::vector<double> times;
    times.reserve(runs);

    for (int r = 0; r < runs; ++r) {
        auto t0 = Clock::now();
        if (bc.phase2_only) {
            solve_phase2(bc.lp, bc.init_basis);
        } else {
            solve(bc.lp);
        }
        auto t1 = Clock::now();
        times.push_back(Duration(t1 - t0).count());
    }

    double med = median_ns(times);
    if (med >= 1e9)
        std::printf("  %-30s  median %8.2f s\n",  bc.name.c_str(), med / 1e9);
    else if (med >= 1e6)
        std::printf("  %-30s  median %8.3f ms\n", bc.name.c_str(), med / 1e6);
    else if (med >= 1e3)
        std::printf("  %-30s  median %8.3f us\n", bc.name.c_str(), med / 1e3);
    else
        std::printf("  %-30s  median %8.1f ns\n", bc.name.c_str(), med);
}

int main()
{
    // ---- Small canonical problem ----
    Eigen::MatrixXd A_small(2, 4);
    A_small << 2, 1, 1, 0,
               1, 2, 0, 1;
    Eigen::VectorXd b_small(2); b_small << 4, 4;
    Eigen::VectorXd c_small(4); c_small << 4, 3, 0, 0;
    LinearProgram lp_small(A_small, b_small, c_small);
    BasisState    basis_small{{2, 3}, {0, 1}};

    std::vector<BenchCase> cases = {
        { "small/phase2_only", lp_small, basis_small, true  },
        { "small/two_phase",   lp_small, {},          false },
    };

    for (auto [m, n] : std::vector<std::pair<int,int>>{{10,20},{50,100},{100,200}}) {
        std::string tag = "m=" + std::to_string(m) + "_n=" + std::to_string(n);
        cases.push_back({ "random/" + tag, random_lp(m, n), {}, false });
    }

    std::cout << "=== Simplex C++ Benchmark (median over 10 runs) ===\n\n";
    for (const auto& bc : cases) run_bench(bc);
    std::cout << "\nDone.\n";
    return 0;
}
