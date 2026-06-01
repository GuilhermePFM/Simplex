use nalgebra::DVector;

use crate::logger::SimplexLogger;
use crate::phase1::col_submatrix;
use crate::pivot::{leaving_index, PivotRule};
use crate::types::{BasisState, LinearProgram, SimplexResult, SolveStatus};

pub fn phase2(
    lp: &LinearProgram,
    init_basis: &BasisState,
    rule: &dyn PivotRule,
    logger: &mut SimplexLogger,
    maxit: usize,
    tol: f64,
) -> SimplexResult {
    let n = lp.ncols();

    let mut basic = init_basis.basic.clone();
    let mut nonbasic = init_basis.nonbasic.clone();

    logger.log_phase(2);

    for it in 1..=maxit {
        let b_mat = col_submatrix(&lp.a, &basic);
        let n_mat = col_submatrix(&lp.a, &nonbasic);

        let lu = b_mat.lu();
        let xb = lu.solve(&lp.b).expect("Phase 2: B is singular");
        let binv_n = lu.solve(&n_mat).expect("Phase 2: B is singular (BinvN)");

        let cb: DVector<f64> = DVector::from_iterator(basic.len(), basic.iter().map(|&i| lp.c[i]));
        let cn: DVector<f64> =
            DVector::from_iterator(nonbasic.len(), nonbasic.iter().map(|&i| lp.c[i]));
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
            let mut x = DVector::<f64>::zeros(n);
            for (k, &bi) in basic.iter().enumerate() {
                x[bi] = xb[k];
            }
            return SimplexResult { x, z, status: SolveStatus::Optimal, iterations: it };
        }

        let j = j.unwrap();
        let d: Vec<f64> = (0..basic.len()).map(|r| binv_n[(r, j)]).collect();
        let i = leaving_index(xb.as_slice(), &d, &binv_n, tol);

        if i.is_none() {
            let mut ray = DVector::<f64>::zeros(n);
            ray[nonbasic[j]] = 1.0;
            for (k, &bi) in basic.iter().enumerate() {
                ray[bi] = -d[k];
            }
            return SimplexResult { x: ray, z: f64::INFINITY, status: SolveStatus::Unbounded, iterations: it };
        }

        let i = i.unwrap();

        std::mem::swap(&mut basic[i], &mut nonbasic[j]);
    }

    panic!("Phase 2 did not converge in {} iterations — possible cycling", maxit);
}
