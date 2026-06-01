package simplex

import "math"

// matGet returns element (i,j) of a row-major matrix with the given number of cols.
func matGet(A []float64, cols, i, j int) float64 {
	return A[i*cols+j]
}

// matSet sets element (i,j) of a row-major matrix with the given number of cols.
func matSet(A []float64, cols, i, j int, v float64) {
	A[i*cols+j] = v
}

// getCol returns column j of an m×cols row-major matrix as a new slice.
func getCol(A []float64, m, cols, j int) []float64 {
	col := make([]float64, m)
	for i := 0; i < m; i++ {
		col[i] = A[i*cols+j]
	}
	return col
}

// getRow returns row i of a row-major matrix as a new slice of length n.
func getRow(A []float64, cols, i, n int) []float64 {
	row := make([]float64, n)
	copy(row, A[i*cols:i*cols+n])
	return row
}

// colsOf returns a m×len(indices) row-major submatrix of selected columns.
func colsOf(A []float64, m, cols int, indices []int) []float64 {
	k := len(indices)
	out := make([]float64, m*k)
	for i := 0; i < m; i++ {
		for jj, j := range indices {
			out[i*k+jj] = A[i*cols+j]
		}
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

// luFactor computes LU factorization with partial pivoting of the m×m matrix A
// (stored row-major). Returns a copy with L and U packed in the same array, and
// a pivot index slice piv of length m such that piv[i] is the row index that
// was swapped into row i.
func luFactor(A []float64, m int) ([]float64, []int) {
	lu := make([]float64, m*m)
	copy(lu, A)
	piv := make([]int, m)
	for i := range piv {
		piv[i] = i
	}

	for k := 0; k < m; k++ {
		// find pivot row
		maxVal := math.Abs(lu[k*m+k])
		maxRow := k
		for i := k + 1; i < m; i++ {
			if v := math.Abs(lu[i*m+k]); v > maxVal {
				maxVal = v
				maxRow = i
			}
		}
		// swap rows k and maxRow
		if maxRow != k {
			piv[k], piv[maxRow] = piv[maxRow], piv[k]
			for j := 0; j < m; j++ {
				lu[k*m+j], lu[maxRow*m+j] = lu[maxRow*m+j], lu[k*m+j]
			}
		}
		pivot := lu[k*m+k]
		if math.Abs(pivot) < 1e-300 {
			continue // singular — leave as-is
		}
		for i := k + 1; i < m; i++ {
			lu[i*m+k] /= pivot
			for j := k + 1; j < m; j++ {
				lu[i*m+j] -= lu[i*m+k] * lu[k*m+j]
			}
		}
	}
	return lu, piv
}

// luSolve solves the system LU·x = b using the factored LU and pivot vector.
// Returns a new slice of length m.
func luSolve(lu []float64, piv []int, b []float64, m int) []float64 {
	x := make([]float64, m)
	// apply row permutation
	for i := 0; i < m; i++ {
		x[i] = b[piv[i]]
	}
	// forward substitution (L is unit lower triangular)
	for i := 0; i < m; i++ {
		for j := 0; j < i; j++ {
			x[i] -= lu[i*m+j] * x[j]
		}
	}
	// backward substitution (U is upper triangular)
	for i := m - 1; i >= 0; i-- {
		for j := i + 1; j < m; j++ {
			x[i] -= lu[i*m+j] * x[j]
		}
		x[i] /= lu[i*m+i]
	}
	return x
}

// luSolveMatrix solves LU·X = B where B is m×k (row-major).
// Returns the solution as an m×k row-major matrix.
func luSolveMatrix(lu []float64, piv []int, B []float64, m, k int) []float64 {
	out := make([]float64, m*k)
	colBuf := make([]float64, m)
	for j := 0; j < k; j++ {
		for i := 0; i < m; i++ {
			colBuf[i] = B[i*k+j]
		}
		sol := luSolve(lu, piv, colBuf, m)
		for i := 0; i < m; i++ {
			out[i*k+j] = sol[i]
		}
	}
	return out
}

// matVecMul computes A·x where A is m×n (row-major) and x is length n.
// Returns a new slice of length m.
func matVecMul(A []float64, m, n int, x []float64) []float64 {
	out := make([]float64, m)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			out[i] += A[i*n+j] * x[j]
		}
	}
	return out
}

// vecScale returns a new slice with each element of v multiplied by s.
func vecScale(v []float64, s float64) []float64 {
	out := make([]float64, len(v))
	for i, val := range v {
		out[i] = val * s
	}
	return out
}

// vecSub returns a - b element-wise.
func vecSub(a, b []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
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

// matTransposeVecMul computes A'·x where A is m×n (row-major) and x is length m.
// Returns a new slice of length n.
func matTransposeVecMul(A []float64, m, n int, x []float64) []float64 {
	out := make([]float64, n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			out[j] += A[i*n+j] * x[i]
		}
	}
	return out
}
