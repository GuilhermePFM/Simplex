#include "test_helpers.hpp"
#include "simplex/simplex.hpp"
#include <iostream>
#include <cmath>

static bool test_production_via_solve()
{
    std::cout << "  test: production problem via solve()\n";
    Eigen::MatrixXd A(2, 4);
    A << 2, 1, 1, 0,
         1, 2, 0, 1;
    Eigen::VectorXd b(2); b << 4, 4;
    Eigen::VectorXd c(4); c << 4, 3, 0, 0;

    SimplexResult r = solve(LinearProgram(A, b, c));

    CHECK_EQ(r.status, SolveStatus::OPTIMAL)
    CHECK_NEAR(r.z,    28.0 / 3.0, 1e-6)
    CHECK_NEAR(r.x[0], 4.0 / 3.0, 1e-6)
    CHECK_NEAR(r.x[1], 4.0 / 3.0, 1e-6)
    return true;
}

static bool test_negative_b_row()
{
    std::cout << "  test: problem with negative-b row\n";
    // A is 3x5, third row has b=-1 (make_b_nonnegative must flip it)
    Eigen::MatrixXd A(3, 5);
    A <<  2,  1, 1, 0, 0,
          1,  2, 0, 1, 0,
         -1, -1, 0, 0, 1;
    Eigen::VectorXd b(3); b << 4, 4, -1;
    Eigen::VectorXd c(5); c << 4, 3, 0, 0, 0;

    LinearProgram lp(A, b, c);
    SimplexResult r = solve(lp);

    CHECK_EQ(r.status, SolveStatus::OPTIMAL)
    CHECK_TRUE(r.z > 0)
    // Feasibility: Ax ≈ b
    Eigen::VectorXd Ax = A * r.x;
    for (int i = 0; i < 3; ++i) {
        CHECK_NEAR(Ax[i], b[i], 1e-6)
    }
    for (int i = 0; i < 5; ++i) {
        CHECK_TRUE(r.x[i] >= -1e-10)
    }
    return true;
}

static bool test_contradictory_infeasible()
{
    std::cout << "  test: contradictory constraints -> INFEASIBLE\n";
    // x0 + x1 = 1  AND  x0 + x1 = 0 — impossible
    Eigen::MatrixXd A(2, 2);
    A << 1, 1,
         1, 1;
    Eigen::VectorXd b(2); b << 1, 0;
    Eigen::VectorXd c(2); c << 1, 1;

    SimplexResult r = solve(LinearProgram(A, b, c));
    CHECK_EQ(r.status, SolveStatus::INFEASIBLE)
    return true;
}

static bool test_conflicting_bounds_infeasible()
{
    std::cout << "  test: conflicting bounds -> INFEASIBLE\n";
    // x0 + s0 = 1   (x0 <= 1)
    // -x0 + s1 = -2  -> after flip: x0 - s1 = 2  (x0 >= 2) — impossible
    Eigen::MatrixXd A(2, 3);
    A <<  1, 1, 0,
         -1, 0, 1;
    Eigen::VectorXd b(2); b << 1, -2;
    Eigen::VectorXd c(3); c << 1, 0, 0;

    SimplexResult r = solve(LinearProgram(A, b, c));
    CHECK_EQ(r.status, SolveStatus::INFEASIBLE)
    return true;
}

bool run_phase1_tests()
{
    bool ok = true;
    ok = test_production_via_solve()       && ok;
    ok = test_negative_b_row()             && ok;
    ok = test_contradictory_infeasible()   && ok;
    ok = test_conflicting_bounds_infeasible() && ok;
    return ok;
}
