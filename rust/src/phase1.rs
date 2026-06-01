use nalgebra::{DMatrix, DVector};

use crate::logger::SimplexLogger;
use crate::pivot::{leaving_index, PivotRule};
use crate::preprocess::make_b_nonnegative;
use crate::types::{BasisState, LinearProgram, SolveStatus};

pub fn phase1(
    lp: LinearProgram,
    rule: &dyn PivotRule,
    logger: &mut SimplexLogger,
    maxit: usize,
    tol: f64,
) -> (BasisState, SolveStatus) {
    let lp = make_b_nonnegative(lp);
    let m = lp.nrows();
    let n = lp.ncols();

    // Aw = [A | I_m]
    let mut aw = DMatrix::<f64>::zeros(m, n + m);
    for i in 0..m {
        for j in 0..n {
            aw[(i, j)] = lp.a[(i, j)];
        }
        aw[(i, n + i)] = 1.0;
    }

    // cw = [0..0, -1..-1]
    let mut cw = DVector::<f64>::zeros(n + m);
    for i in n..n + m {
        cw[i] = -1.0;
    }

    let art_start = n; // 0-based: artificials are columns n..n+m
    let mut basic: Vec<usize> = (n..n + m).collect();
    let mut nonbasic: Vec<usize> = (0..n).collect();

    logger.log_phase(1);

    for it in 1..=maxit {
        let b_mat = col_submatrix(&aw, &basic);
        let n_mat = col_submatrix(&aw, &nonbasic);

        let lu = b_mat.lu();
        let xb = lu.solve(&lp.b).expect("Phase 1: B is singular");
        let binv_n = lu.solve(&n_mat).expect("Phase 1: B is singular (BinvN)");

        let cb: DVector<f64> = DVector::from_iterator(basic.len(), basic.iter().map(|&i| cw[i]));
        let cn: DVector<f64> =
            DVector::from_iterator(nonbasic.len(), nonbasic.iter().map(|&i| cw[i]));
        let cr_vec: DVector<f64> = cn - binv_n.transpose() * &cb;
        let z = cb.dot(&xb);

        logger.log_iteration(
            it,
            &BasisState::new(basic.clone(), nonbasic.clone()),
            xb.as_slice(),
            z,
        );

        let j = rule.entering_index(cr_vec.as_slice(), tol);

        if j.is_none() {
            if z < -tol {
                return (BasisState::new(vec![], vec![]), SolveStatus::Infeasible);
            }

            fix_artificials_in_basis(&mut aw, &mut basic, &mut nonbasic, art_start);
            let orig_basic: Vec<usize> = basic.iter().copied().filter(|&i| i < n).collect();
            let orig_nonbasic: Vec<usize> =
                nonbasic.iter().copied().filter(|&i| i < n).collect();
            return (BasisState::new(orig_basic, orig_nonbasic), SolveStatus::Optimal);
        }

        let j = j.unwrap();
        let d: Vec<f64> = (0..m).map(|r| binv_n[(r, j)]).collect();
        let i = leaving_index(xb.as_slice(), &d, &binv_n, tol);

        // Phase 1 auxiliary problem is always bounded (artificials provide a floor)
        let i = i.expect("Phase 1 is unexpectedly unbounded");

        std::mem::swap(&mut basic[i], &mut nonbasic[j]);
    }

    panic!("Phase 1 did not converge in {} iterations", maxit);
}

fn fix_artificials_in_basis(
    aw: &mut DMatrix<f64>,
    basic: &mut Vec<usize>,
    nonbasic: &mut Vec<usize>,
    art_start: usize,
) {
    let mut changed = true;
    while changed {
        changed = false;
        let m = basic.len();
        for pos in 0..m {
            let col = basic[pos];
            if col < art_start {
                continue;
            }

            let b_mat = col_submatrix(aw, basic);
            let lu = b_mat.lu();
            let binv_aw = lu.solve(&*aw).expect("fix_artificials: B is singular");
            let row: Vec<f64> = (0..aw.ncols()).map(|k| binv_aw[(pos, k)]).collect();

            let swap_j = nonbasic.iter().position(|&ncol| {
                ncol < art_start && row[ncol].abs() > 1e-10
            });

            if let Some(j) = swap_j {
                std::mem::swap(&mut basic[pos], &mut nonbasic[j]);
                changed = true;
                break;
            } else {
                // Zero out this row of aw — artificial is redundant
                for k in 0..aw.ncols() {
                    aw[(pos, k)] = 0.0;
                }
            }
        }
    }
}

pub fn col_submatrix(a: &DMatrix<f64>, cols: &[usize]) -> DMatrix<f64> {
    let m = a.nrows();
    let k = cols.len();
    DMatrix::from_fn(m, k, |i, j| a[(i, cols[j])])
}
