use approx::assert_abs_diff_eq;
use nalgebra::{DMatrix, DVector};
use simplex::{solve, LargestCoefficient, LinearProgram, LogLevel, SolveStatus};

#[test]
fn feasible_production_via_solve() {
    // Same as the direct Phase 2 test, but going through solve() so Phase 1 runs
    let lp = LinearProgram::new(
        DMatrix::from_row_slice(2, 4, &[2.0, 1.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0]),
        DVector::from_vec(vec![4.0, 4.0]),
        DVector::from_vec(vec![4.0, 3.0, 0.0, 0.0]),
    );

    let r = solve(lp, &LargestCoefficient, LogLevel::Silent, None);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert_abs_diff_eq!(r.z, 28.0 / 3.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[0], 4.0 / 3.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[1], 4.0 / 3.0, epsilon = 1e-6);
}

#[test]
fn feasible_negative_b_row() {
    // Constraint row 3 has b=-1; make_b_nonnegative flips it before Phase 1
    // Original: 2x1+x2+x3=4, x1+2x2+x4=4, -x1-x2+x5=-1 (i.e., x1+x2-x5=1)
    let lp = LinearProgram::new(
        DMatrix::from_row_slice(
            3,
            5,
            &[
                2.0, 1.0, 1.0, 0.0, 0.0,
                1.0, 2.0, 0.0, 1.0, 0.0,
               -1.0,-1.0, 0.0, 0.0, 1.0,
            ],
        ),
        DVector::from_vec(vec![4.0, 4.0, -1.0]),
        DVector::from_vec(vec![4.0, 3.0, 0.0, 0.0, 0.0]),
    );

    let r = solve(lp.clone(), &LargestCoefficient, LogLevel::Silent, None);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert!(r.z > 0.0);
    // Feasibility check: Ax ≈ b
    let ax = &lp.a * &r.x;
    for i in 0..lp.b.len() {
        assert_abs_diff_eq!(ax[i], lp.b[i], epsilon = 1e-6);
    }
    assert!(r.x.iter().all(|&v| v >= -1e-10));
}

#[test]
fn infeasible_contradictory() {
    // x1 + x2 = 1  AND  x1 + x2 = 0  — impossible
    let lp = LinearProgram::new(
        DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 1.0]),
        DVector::from_vec(vec![1.0, 0.0]),
        DVector::from_vec(vec![1.0, 1.0]),
    );

    let r = solve(lp, &LargestCoefficient, LogLevel::Silent, None);

    assert_eq!(r.status, SolveStatus::Infeasible);
}

#[test]
fn infeasible_conflicting_bounds() {
    // x1 + s1 = 1  (x1 <= 1)
    // x1 - s2 = 2  (x1 >= 2, written as -x1 + s2 = -2, negated by preprocess)
    let lp = LinearProgram::new(
        DMatrix::from_row_slice(2, 3, &[1.0, 1.0, 0.0, -1.0, 0.0, 1.0]),
        DVector::from_vec(vec![1.0, -2.0]),
        DVector::from_vec(vec![1.0, 0.0, 0.0]),
    );

    let r = solve(lp, &LargestCoefficient, LogLevel::Silent, None);

    assert_eq!(r.status, SolveStatus::Infeasible);
}
