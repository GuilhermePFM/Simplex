package simplex

// MakeBNonnegative returns a copy of lp with all right-hand side values non-negative.
// Rows where b[i] < 0 have their sign flipped in both A and b.
func MakeBNonnegative(lp LinearProgram) LinearProgram {
	needsFlip := false
	for _, v := range lp.B {
		if v < 0 {
			needsFlip = true
			break
		}
	}
	if !needsFlip {
		return lp
	}

	newA := make([]float64, len(lp.A))
	copy(newA, lp.A)
	newB := make([]float64, len(lp.B))
	copy(newB, lp.B)

	for i, v := range lp.B {
		if v < 0 {
			newB[i] = -v
			for j := 0; j < lp.N; j++ {
				newA[i*lp.N+j] = -lp.A[i*lp.N+j]
			}
		}
	}

	return LinearProgram{
		A: newA,
		B: newB,
		C: lp.C,
		M: lp.M,
		N: lp.N,
	}
}
