use crate::types::LinearProgram;

/// Return a copy of `lp` with all right-hand side values non-negative.
pub fn make_b_nonnegative(lp: LinearProgram) -> LinearProgram {
    let any_neg = lp.b.iter().any(|&v| v < 0.0);
    if !any_neg {
        return lp;
    }

    let mut a2 = lp.a.clone();
    let mut b2 = lp.b.clone();

    for i in 0..b2.len() {
        if b2[i] < 0.0 {
            b2[i] *= -1.0;
            for j in 0..a2.ncols() {
                a2[(i, j)] *= -1.0;
            }
        }
    }

    LinearProgram::new(a2, b2, lp.c)
}
