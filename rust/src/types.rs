use nalgebra::{DMatrix, DVector};

/// A linear program in standard equality form:
///   maximize  c'x
///   subject to Ax = b,  x >= 0
#[derive(Clone)]
pub struct LinearProgram {
    pub a: DMatrix<f64>,
    pub b: DVector<f64>,
    pub c: DVector<f64>,
}

impl LinearProgram {
    pub fn new(a: DMatrix<f64>, b: DVector<f64>, c: DVector<f64>) -> Self {
        let (m, n) = a.shape();
        assert_eq!(b.len(), m, "b has {} rows; A has {}", b.len(), m);
        assert_eq!(c.len(), n, "c has {} entries; A has {} columns", c.len(), n);
        LinearProgram { a, b, c }
    }

    pub fn nrows(&self) -> usize {
        self.a.nrows()
    }

    pub fn ncols(&self) -> usize {
        self.a.ncols()
    }
}

/// Index partition of variables into basic (in basis B) and non-basic (fixed at zero).
/// Both vectors hold 0-based column indices into the LP's A matrix.
#[derive(Clone, Debug)]
pub struct BasisState {
    pub basic: Vec<usize>,
    pub nonbasic: Vec<usize>,
}

impl BasisState {
    pub fn new(basic: Vec<usize>, nonbasic: Vec<usize>) -> Self {
        BasisState { basic, nonbasic }
    }
}

/// Termination status of a Simplex solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveStatus {
    Optimal,
    Unbounded,
    Infeasible,
}

/// Complete output of `solve` or `solve_phase2`.
pub struct SimplexResult {
    pub x: DVector<f64>,
    pub z: f64,
    pub status: SolveStatus,
    pub iterations: usize,
}
