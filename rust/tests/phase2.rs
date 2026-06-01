use approx::assert_abs_diff_eq;
use nalgebra::{DMatrix, DVector};
use simplex::{solve_phase2, Bland, BasisState, LargestCoefficient, LinearProgram, LogLevel, SolveStatus};

fn production_lp() -> LinearProgram {
    // maximize 4x1 + 3x2  s.t.  2x1+x2<=4, x1+2x2<=4, x1,x2>=0
    // Standard form (slacks x3, x4)
    LinearProgram::new(
        DMatrix::from_row_slice(2, 4, &[2.0, 1.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0]),
        DVector::from_vec(vec![4.0, 4.0]),
        DVector::from_vec(vec![4.0, 3.0, 0.0, 0.0]),
    )
}

#[test]
fn optimal_production_problem() {
    let lp = production_lp();
    // Julia basis [3,4],[1,2] is 1-based → 0-based: basic=[2,3], nonbasic=[0,1]
    let basis = BasisState::new(vec![2, 3], vec![0, 1]);
    let r = solve_phase2(&lp, &basis, &LargestCoefficient, LogLevel::Silent);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert_abs_diff_eq!(r.z, 28.0 / 3.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[0], 4.0 / 3.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[1], 4.0 / 3.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[2], 0.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[3], 0.0, epsilon = 1e-6);
    assert_eq!(r.iterations, 3);
}

#[test]
fn optimal_already_at_vertex() {
    // x1=1, x2=1 is the only feasible point; objective is c'x=5
    let lp = LinearProgram::new(
        DMatrix::from_row_slice(2, 4, &[1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]),
        DVector::from_vec(vec![1.0, 1.0]),
        DVector::from_vec(vec![3.0, 2.0, 0.0, 0.0]),
    );
    // Julia basis [1,2],[3,4] → 0-based: basic=[0,1], nonbasic=[2,3]
    let basis = BasisState::new(vec![0, 1], vec![2, 3]);
    let r = solve_phase2(&lp, &basis, &LargestCoefficient, LogLevel::Silent);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert_abs_diff_eq!(r.z, 5.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[0], 1.0, epsilon = 1e-6);
    assert_abs_diff_eq!(r.x[1], 1.0, epsilon = 1e-6);
}

#[test]
fn unbounded() {
    // maximize x1+x2  s.t.  0.5x1 - x2 + x3 = 0.5,  -4x1 + x2 + x4 = 1
    let lp = LinearProgram::new(
        DMatrix::from_row_slice(2, 4, &[0.5, -1.0, 1.0, 0.0, -4.0, 1.0, 0.0, 1.0]),
        DVector::from_vec(vec![0.5, 1.0]),
        DVector::from_vec(vec![1.0, 1.0, 0.0, 0.0]),
    );
    // Julia basis [3,4],[1,2] → 0-based: basic=[2,3], nonbasic=[0,1]
    let basis = BasisState::new(vec![2, 3], vec![0, 1]);
    let r = solve_phase2(&lp, &basis, &LargestCoefficient, LogLevel::Silent);

    assert_eq!(r.status, SolveStatus::Unbounded);
    assert_eq!(r.z, f64::INFINITY);
}

#[test]
fn bland_gives_same_optimal() {
    let lp = production_lp();
    let basis = BasisState::new(vec![2, 3], vec![0, 1]);
    let r = solve_phase2(&lp, &basis, &Bland, LogLevel::Silent);

    assert_eq!(r.status, SolveStatus::Optimal);
    assert_abs_diff_eq!(r.z, 28.0 / 3.0, epsilon = 1e-6);
}
