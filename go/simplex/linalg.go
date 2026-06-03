package simplex

import (
	"gonum.org/v1/gonum/mat"
)

// matFromRowMajor wraps a row-major slice as a Gonum dense matrix.
func matFromRowMajor(A []float64, m, n int) *mat.Dense {
	M := mat.NewDense(m, n, nil)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			M.Set(i, j, A[i*n+j])
		}
	}
	return M
}

// colsFromRowMajor builds an m×len(indices) matrix from selected columns of a
// row-major m×n matrix A.
func colsFromRowMajor(A []float64, m, n int, indices []int) *mat.Dense {
	k := len(indices)
	B := mat.NewDense(m, k, nil)
	for jj, col := range indices {
		for i := 0; i < m; i++ {
			B.Set(i, jj, A[i*n+col])
		}
	}
	return B
}

// factorBasis computes the LU factorization of square basis matrix B.
func factorBasis(B *mat.Dense) mat.LU {
	var lu mat.LU
	lu.Factorize(B)
	return lu
}

// solveVec solves B·x = b using a factored LU.
func solveVec(lu mat.LU, b []float64) []float64 {
	m := len(b)
	bv := mat.NewVecDense(m, b)
	var x mat.VecDense
	if err := lu.SolveVecTo(&x, false, bv); err != nil {
		panic("singular basis in solveVec: " + err.Error())
	}
	out := make([]float64, m)
	for i := 0; i < m; i++ {
		out[i] = x.AtVec(i)
	}
	return out
}

// solveMat solves B·X = rhs using a factored LU.
func solveMat(lu mat.LU, rhs *mat.Dense) *mat.Dense {
	var out mat.Dense
	if err := lu.SolveTo(&out, false, rhs); err != nil {
		panic("singular basis in solveMat: " + err.Error())
	}
	return mat.DenseCopyOf(&out)
}

// reducedCosts returns c_N - BinvN'·cb.
func reducedCosts(cNonbasic []float64, BinvN *mat.Dense, cb []float64) []float64 {
	nb := len(cNonbasic)
	m, _ := BinvN.Dims()
	cbVec := mat.NewVecDense(m, cb)
	t := mat.Transpose{Matrix: BinvN}
	var BinvNtcb mat.VecDense
	BinvNtcb.MulVec(&t, cbVec)
	cr := make([]float64, nb)
	copy(cr, cNonbasic)
	for j := 0; j < nb; j++ {
		cr[j] -= BinvNtcb.AtVec(j)
	}
	return cr
}

// denseCol copies column j of BinvN into a new slice (for LeavingIndex).
func denseCol(BinvN *mat.Dense, j int) []float64 {
	m, _ := BinvN.Dims()
	col := make([]float64, m)
	for i := 0; i < m; i++ {
		col[i] = BinvN.At(i, j)
	}
	return col
}

// binvNToRowMajor copies BinvN into row-major layout for LeavingIndex lex tie-break.
func binvNToRowMajor(BinvN *mat.Dense) []float64 {
	m, nb := BinvN.Dims()
	out := make([]float64, m*nb)
	for i := 0; i < m; i++ {
		for j := 0; j < nb; j++ {
			out[i*nb+j] = BinvN.At(i, j)
		}
	}
	return out
}

// gather returns a new slice with elements a[indices[i]].
func gather(a []float64, indices []int) []float64 {
	out := make([]float64, len(indices))
	for i, idx := range indices {
		out[i] = a[idx]
	}
	return out
}

// dot computes the dot product of two equal-length slices.
func dot(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}
