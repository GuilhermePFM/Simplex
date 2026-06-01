#include "test_helpers.hpp"
#include "simplex/simplex.hpp"
#include <iostream>

// Julia BasisState([3,4,5],[1,2]) -> C++ basic={2,3,4}, nonbasic={0,1}
static bool test_degenerate_tie_breaking()
{
    std::cout << "  test: degenerate tie-breaking\n";
    // maximize x0+x1  s.t. x0+x1<=2, x0<=1, x1<=1
    // Standard form: add slacks x2, x3, x4
    Eigen::MatrixXd A(3, 5);
    A << 1, 1, 1, 0, 0,
         1, 0, 0, 1, 0,
         0, 1, 0, 0, 1;
    Eigen::VectorXd b(3); b << 2, 1, 1;
    Eigen::VectorXd c(5); c << 1, 1, 0, 0, 0;

    LinearProgram lp(A, b, c);
    BasisState    basis{{2, 3, 4}, {0, 1}};

    SimplexResult r = solve_phase2(lp, basis, LargestCoefficient{});

    CHECK_EQ(r.status, SolveStatus::OPTIMAL)
    CHECK_NEAR(r.z,    2.0, 1e-6)
    CHECK_NEAR(r.x[0], 1.0, 1e-6)
    CHECK_NEAR(r.x[1], 1.0, 1e-6)
    CHECK_EQ(r.iterations, 3)

    // Feasibility
    Eigen::VectorXd Ax = A * r.x;
    for (int i = 0; i < 3; ++i) CHECK_NEAR(Ax[i], b[i], 1e-6)
    for (int i = 0; i < 5; ++i) CHECK_TRUE(r.x[i] >= -1e-10)

    // Bland also finds optimal
    SimplexResult rb = solve_phase2(lp, basis, Bland{});
    CHECK_EQ(rb.status, SolveStatus::OPTIMAL)
    CHECK_NEAR(rb.z, 2.0, 1e-6)
    return true;
}

// Julia BasisState([5,6,7],[1,2,3,4]) -> C++ basic={4,5,6}, nonbasic={0,1,2,3}
static bool test_kotiah_steinberg()
{
    std::cout << "  test: Kotiah-Steinberg degenerate start\n";
    // maximize 10x0 - 57x1 - 9x2 - 24x3
    // Initial BFS: x4=x5=0, x6=1 — degenerate
    Eigen::MatrixXd A(3, 7);
    A << -0.5,  5.5,  2.5, -9, 1, 0, 0,
         -0.5,  2.5,  0.5, -1, 0, 1, 0,
          1.0,  0.0,  0.0,  0, 0, 0, 1;
    Eigen::VectorXd b(3); b << 0, 0, 1;
    Eigen::VectorXd c(7); c << 10, -57, -9, -24, 0, 0, 0;

    LinearProgram lp(A, b, c);
    BasisState    basis{{4, 5, 6}, {0, 1, 2, 3}};

    SimplexResult r = solve_phase2(lp, basis, LargestCoefficient{});

    CHECK_EQ(r.status, SolveStatus::OPTIMAL)
    CHECK_NEAR(r.z,    10.0, 1e-6)
    CHECK_NEAR(r.x[0],  1.0, 1e-6)

    Eigen::VectorXd Ax = A * r.x;
    for (int i = 0; i < 3; ++i) CHECK_NEAR(Ax[i], b[i], 1e-6)
    return true;
}

bool run_anticycling_tests()
{
    bool ok = true;
    ok = test_degenerate_tie_breaking() && ok;
    ok = test_kotiah_steinberg()        && ok;
    return ok;
}
