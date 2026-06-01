use approx::assert_abs_diff_eq;
use nalgebra::{DMatrix, DVector};
use simplex::{solve_phase2, Bland, BasisState, LargestCoefficient, LinearProgram, LogLevel, SolveStatus};

#[test]
fn degenerate_largest_coefficient() {
    // maximize x1 + x2
    // s.t.  x1 + x2 <= 2,  x1 <= 1,  x2 <= 1
    // Standard form: A = [[1,1,1,0,0],[1,0,0,1,0],[0,1,0,0,1]], b=[2,1,1]
    let a = DMatrix::from_row_slice(
        3,
        5,
        &[
            1.0, 1.0, 1.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 1.0,
        ],
    );
    let b = DVector::from_vec(vec![2.0, 1.0, 1.0]);
    let c = DVector::from_vec(vec![1.0, 1.0, 0.0, 0.0, 0.0]);

    let lp = LinearProgram::new(a.clone(), b.clone(), c);
    // Julia basis [3,4,5],[1,2] → 0-based: basic=[2,3,4], nonbasic=[0,1]
    let basis = BasisState::new(vec![2, 3, 4], vec![0, 1]);

    let r = solve_phase2(&lp, &basis, &LargestCoefficient, LogLevel::Silent);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert_abs_diff_eq!(r.z, 2.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[0], 1.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[1], 1.0, epsilon = 1e-6);
    assert_eq!(r.iterations, 3);

    // Feasibility certificate
    let ax = &lp.a * &r.x;
    for i in 0..lp.b.len() {
        assert_abs_diff_eq!(ax[i], lp.b[i], epsilon = 1e-6);
    }
    assert!(r.x.iter().all(|&v| v >= -1e-10));
}

#[test]
fn degenerate_bland() {
    let a = DMatrix::from_row_slice(
        3,
        5,
        &[
            1.0, 1.0, 1.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 1.0,
        ],
    );
    let b = DVector::from_vec(vec![2.0, 1.0, 1.0]);
    let c = DVector::from_vec(vec![1.0, 1.0, 0.0, 0.0, 0.0]);

    let lp = LinearProgram::new(a, b, c);
    let basis = BasisState::new(vec![2, 3, 4], vec![0, 1]);

    let r = solve_phase2(&lp, &basis, &Bland, LogLevel::Silent);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert_abs_diff_eq!(r.z, 2.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[0], 1.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[1], 1.0, epsilon = 1e-6);
}

#[test]
fn kotiah_steinberg() {
    // maximize 10x1 - 57x2 - 9x3 - 24x4
    // Degenerate start: x5=x6=0, x7=1 — optimal: x1=1, z=10
    let a = DMatrix::from_row_slice(
        3,
        7,
        &[
            -0.5,  5.5,  2.5, -9.0, 1.0, 0.0, 0.0,
            -0.5,  2.5,  0.5, -1.0, 0.0, 1.0, 0.0,
             1.0,  0.0,  0.0,  0.0, 0.0, 0.0, 1.0,
        ],
    );
    let b = DVector::from_vec(vec![0.0, 0.0, 1.0]);
    let c = DVector::from_vec(vec![10.0, -57.0, -9.0, -24.0, 0.0, 0.0, 0.0]);

    let lp = LinearProgram::new(a.clone(), b.clone(), c);
    // Julia basis [5,6,7],[1,2,3,4] → 0-based: basic=[4,5,6], nonbasic=[0,1,2,3]
    let basis = BasisState::new(vec![4, 5, 6], vec![0, 1, 2, 3]);

    let r = solve_phase2(&lp, &basis, &LargestCoefficient, LogLevel::Silent);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert_abs_diff_eq!(r.z, 10.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[0], 1.0, epsilon = 1e-6);

    let ax = &lp.a * &r.x;
    for i in 0..lp.b.len() {
        assert_abs_diff_eq!(ax[i], lp.b[i], epsilon = 1e-6);
    }
}
