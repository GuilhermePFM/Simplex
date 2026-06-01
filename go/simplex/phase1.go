package simplex

import "math"

// Phase1 finds an initial basic feasible solution using the two-phase method.
// It introduces artificial variables, solves the auxiliary LP, then removes
// any artificial variables that remain in the basis.
//
// Returns the resulting BasisState (in original variable indices) and a status.
// Status is Optimal if a BFS exists, Infeasible otherwise.
func Phase1(lp LinearProgram, rule PivotRule, logger SimplexLogger, maxit int, tol float64) (BasisState, SolveStatus) {
	lp = MakeBNonnegative(lp)
	m, n := lp.M, lp.N

	// Build augmented matrix Aw = [A | I_m], size m × (n+m), stored row-major.
	nw := n + m
	Aw := make([]float64, m*nw)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			Aw[i*nw+j] = lp.A[i*n+j]
		}
		Aw[i*nw+n+i] = 1.0
	}

	// cw = [0...0, -1...-1], artificials have cost -1 (we maximize).
	cw := make([]float64, nw)
	for i := n; i < nw; i++ {
		cw[i] = -1.0
	}

	artStart := n // 0-based: artificials are columns n..n+m-1

	basic := make([]int, m)
	for i := range basic {
		basic[i] = n + i
	}
	nonbasic := make([]int, n)
	for i := range nonbasic {
		nonbasic[i] = i
	}

	logger.LogPhase(1)

	for it := 1; it <= maxit; it++ {
		nb := len(nonbasic)

		// B = Aw[:, basic], size m×m
		B := colsOf(Aw, m, nw, basic)
		// N = Aw[:, nonbasic], size m×nb
		N := colsOf(Aw, m, nw, nonbasic)

		lu, piv := luFactor(B, m)
		xb := luSolve(lu, piv, lp.B, m)
		BinvN := luSolveMatrix(lu, piv, N, m, nb)

		cb := gather(cw, basic)
		// cr = cw[nonbasic] - BinvN' * cb
		cr := gather(cw, nonbasic)
		BinvNtcb := matTransposeVecMul(BinvN, m, nb, cb)
		for j := 0; j < nb; j++ {
			cr[j] -= BinvNtcb[j]
		}
		z := dot(cb, xb)

		logger.LogIteration(it, BasisState{copyInts(basic), copyInts(nonbasic)}, xb, z)

		j, ok := rule.EnteringIndex(cr, tol)
		if !ok {
			// No positive reduced cost — optimal for auxiliary problem.
			if z < -tol {
				return BasisState{}, Infeasible
			}
			fixArtificialsInBasis(Aw, m, nw, basic, nonbasic, artStart)

			origBasic := filterLess(basic, n)
			origNonbasic := filterLess(nonbasic, n)
			return BasisState{origBasic, origNonbasic}, Optimal
		}

		d := getCol(BinvN, m, nb, j)
		i, ok := LeavingIndex(xb, d, BinvN, m, nb, tol)
		if !ok {
			panic("Phase 1 is unexpectedly unbounded")
		}

		basic[i], nonbasic[j] = nonbasic[j], basic[i]
	}

	panic("Phase 1 did not converge")
}

// fixArtificialsInBasis tries to pivot out any artificial variable still in the basis.
// If an artificial cannot be pivoted out (redundant row), that row of Aw is zeroed.
func fixArtificialsInBasis(Aw []float64, m, nw int, basic, nonbasic []int, artStart int) {
	changed := true
	for changed {
		changed = false
		for pos, col := range basic {
			if col < artStart {
				continue
			}
			// Compute BinvAw = B^{-1} * Aw (full m×nw matrix)
			B := colsOf(Aw, m, nw, basic)
			lu, piv := luFactor(B, m)
			BinvAw := luSolveMatrix(lu, piv, Aw, m, nw)

			// Look for a non-artificial non-basic column with non-zero entry in this row.
			swapJ := -1
			for j, ncol := range nonbasic {
				if ncol < artStart && math.Abs(BinvAw[pos*nw+ncol]) > 1e-10 {
					swapJ = j
					break
				}
			}

			if swapJ < 0 {
				// Redundant row — zero it out so it doesn't affect future factorizations.
				for j := 0; j < nw; j++ {
					Aw[pos*nw+j] = 0.0
				}
			} else {
				basic[pos], nonbasic[swapJ] = nonbasic[swapJ], basic[pos]
				changed = true
				break
			}
		}
	}
}

// copyInts returns a copy of a slice of ints.
func copyInts(s []int) []int {
	c := make([]int, len(s))
	copy(c, s)
	return c
}

// filterLess returns elements of s that are strictly less than threshold.
func filterLess(s []int, threshold int) []int {
	out := make([]int, 0, len(s))
	for _, v := range s {
		if v < threshold {
			out = append(out, v)
		}
	}
	return out
}
