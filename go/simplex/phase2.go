package simplex

import "math"

// Phase2 runs the revised simplex optimization loop starting from init_basis.
// It assumes init_basis is a valid basic feasible solution for lp.
func Phase2(lp LinearProgram, initBasis BasisState, rule PivotRule, logger SimplexLogger, maxit int, tol float64) SimplexResult {
	A, b, c := lp.A, lp.B, lp.C
	m, n := lp.M, lp.N

	basic := copyInts(initBasis.Basic)
	nonbasic := copyInts(initBasis.Nonbasic)

	logger.LogPhase(2)

	for it := 1; it <= maxit; it++ {
		nb := len(nonbasic)

		B := colsFromRowMajor(A, m, n, basic)
		N := colsFromRowMajor(A, m, n, nonbasic)

		lu := factorBasis(B)
		xb := solveVec(lu, b)
		BinvN := solveMat(lu, N)
		binvNFlat := binvNToRowMajor(BinvN)

		cr := reducedCosts(gather(c, nonbasic), BinvN, gather(c, basic))
		z := dot(gather(c, basic), xb)

		logger.LogIteration(it, BasisState{copyInts(basic), copyInts(nonbasic)}, xb, z)

		j, ok := rule.EnteringIndex(cr, tol)
		if !ok {
			// All reduced costs non-positive — optimal.
			x := make([]float64, n)
			for ii, bi := range basic {
				x[bi] = xb[ii]
			}
			return SimplexResult{X: x, Z: z, Status: Optimal, Iterations: it}
		}

		d := denseCol(BinvN, j)
		i, ok := LeavingIndex(xb, d, binvNFlat, m, nb, tol)
		if !ok {
			// Unbounded: return extreme ray.
			ray := make([]float64, n)
			ray[nonbasic[j]] = 1.0
			for ii, bi := range basic {
				ray[bi] = -d[ii]
			}
			return SimplexResult{X: ray, Z: math.Inf(1), Status: Unbounded, Iterations: it}
		}

		theta := xb[i] / d[i]

		x := make([]float64, n)
		for ii, bi := range basic {
			x[bi] = xb[ii] - theta*d[ii]
		}
		x[nonbasic[j]] = theta
		x[basic[i]] = 0.0 // numerical cleanup for the leaving variable

		basic[i], nonbasic[j] = nonbasic[j], basic[i]
		_ = x // x is recomputed at optimum; here we just track for potential ray return
	}

	panic("Phase 2 did not converge — possible cycling")
}
