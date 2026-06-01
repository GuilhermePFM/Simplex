#include "test_helpers.hpp"
#include "simplex/simplex.hpp"
#include <cmath>
#include <iostream>

// Julia uses 1-based indices; C++ uses 0-based throughout.
// Julia BasisState([3,4],[1,2]) => C++ basic={2,3}, nonbasic={0,1}

static bool test_optimal_production()
{
    std::cout << "  test: optimal production problem\n";
    // maximize 4x0 + 3x1  s.t. 2x0+x1+x2=4, x0+2x1+x3=4
    Eigen::MatrixXd A(2, 4);
    A << 2, 1, 1, 0,
         1, 2, 0, 1;
    Eigen::VectorXd b(2); b << 4, 4;
    Eigen::VectorXd c(4); c << 4, 3, 0, 0;

    LinearProgram lp(A, b, c);
    BasisState    basis{{2, 3}, {0, 1}};

    SimplexResult r = solve_phase2(lp, basis);

    CHECK_EQ(r.status, SolveStatus::OPTIMAL)
    CHECK_NEAR(r.z,    28.0 / 3.0, 1e-6)
    CHECK_NEAR(r.x[0], 4.0 / 3.0, 1e-6)
    CHECK_NEAR(r.x[1], 4.0 / 3.0, 1e-6)
    CHECK_EQ(r.iterations, 3)
    return true;
}

static bool test_already_at_vertex()
{
    std::cout << "  test: already at vertex\n";
    // x0=1, x1=1 is the only feasible point; z = 3*1 + 2*1 = 5
    Eigen::MatrixXd A(2, 4);
    A << 1, 0, 1, 0,
         0, 1, 0, 1;
    Eigen::VectorXd b(2); b << 1, 1;
    Eigen::VectorXd c(4); c << 3, 2, 0, 0;

    LinearProgram lp(A, b, c);
    BasisState    basis{{0, 1}, {2, 3}};

    SimplexResult r = solve_phase2(lp, basis);

    CHECK_EQ(r.status, SolveStatus::OPTIMAL)
    CHECK_NEAR(r.z, 5.0, 1e-6)
    return true;
}

static bool test_unbounded()
{
    std::cout << "  test: unbounded\n";
    Eigen::MatrixXd A(2, 4);
    A <<  0.5, -1, 1, 0,
         -4.0,  1, 0, 1;
    Eigen::VectorXd b(2); b << 0.5, 1;
    Eigen::VectorXd c(4); c << 1, 1, 0, 0;

    LinearProgram lp(A, b, c);
    BasisState    basis{{2, 3}, {0, 1}};

    SimplexResult r = solve_phase2(lp, basis);

    CHECK_EQ(r.status, SolveStatus::UNBOUNDED)
    CHECK_TRUE(std::isinf(r.z))
    return true;
}

static bool test_bland_same_optimal()
{
    std::cout << "  test: Bland gives same optimal\n";
    Eigen::MatrixXd A(2, 4);
    A << 2, 1, 1, 0,
         1, 2, 0, 1;
    Eigen::VectorXd b(2); b << 4, 4;
    Eigen::VectorXd c(4); c << 4, 3, 0, 0;

    LinearProgram lp(A, b, c);
    BasisState    basis{{2, 3}, {0, 1}};

    SimplexResult r = solve_phase2(lp, basis, Bland{});

    CHECK_EQ(r.status, SolveStatus::OPTIMAL)
    CHECK_NEAR(r.z, 28.0 / 3.0, 1e-6)
    return true;
}

bool run_phase2_tests()
{
    bool ok = true;
    ok = test_optimal_production() && ok;
    ok = test_already_at_vertex()  && ok;
    ok = test_unbounded()          && ok;
    ok = test_bland_same_optimal() && ok;
    return ok;
}
