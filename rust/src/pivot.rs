use nalgebra::DMatrix;

pub trait PivotRule {
    fn entering_index(&self, cr: &[f64], tol: f64) -> Option<usize>;
}

/// Pick the non-basic variable with the largest (most positive) reduced cost.
pub struct LargestCoefficient;

/// Bland's rule: always pick the lowest-index positive reduced cost.
pub struct Bland;

impl PivotRule for LargestCoefficient {
    fn entering_index(&self, cr: &[f64], tol: f64) -> Option<usize> {
        let (j, &v) = cr.iter().enumerate().max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())?;
        if v > tol { Some(j) } else { None }
    }
}

impl PivotRule for Bland {
    fn entering_index(&self, cr: &[f64], tol: f64) -> Option<usize> {
        cr.iter().position(|&v| v > tol)
    }
}

/// Lexicographic ratio test — selects the leaving variable to prevent cycling.
///
/// Returns the row position (into `basic`) of the leaving variable, or `None`
/// if every component of `d` is non-positive (problem is unbounded).
pub fn leaving_index(xb: &[f64], d: &[f64], binv_n: &DMatrix<f64>, tol: f64) -> Option<usize> {
    let candidates: Vec<usize> = (0..d.len()).filter(|&i| d[i] > tol).collect();
    if candidates.is_empty() {
        return None;
    }

    let mut best = candidates[0];
    for &i in &candidates[1..] {
        let ratio_i = xb[i] / d[i];
        let ratio_best = xb[best] / d[best];

        if ratio_i < ratio_best - tol {
            best = i;
        } else if (ratio_i - ratio_best).abs() <= tol {
            // Lexicographic tie-break on scaled rows of BinvN
            let ncols = binv_n.ncols();
            let row_i: Vec<f64> = (0..ncols).map(|k| binv_n[(i, k)] / d[i]).collect();
            let row_best: Vec<f64> = (0..ncols).map(|k| binv_n[(best, k)] / d[best]).collect();
            if row_i < row_best {
                best = i;
            }
        }
    }

    Some(best)
}
