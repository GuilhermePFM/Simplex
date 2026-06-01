package simplex

import "math"

// PivotRule selects the entering variable (non-basic column with positive reduced cost).
type PivotRule interface {
	EnteringIndex(cr []float64, tol float64) (int, bool)
}

// LargestCoefficient picks the non-basic variable with the largest positive reduced cost.
type LargestCoefficient struct{}

func (LargestCoefficient) EnteringIndex(cr []float64, tol float64) (int, bool) {
	best := -1
	bestVal := tol
	for j, v := range cr {
		if v > bestVal {
			bestVal = v
			best = j
		}
	}
	if best < 0 {
		return 0, false
	}
	return best, true
}

// Bland always picks the lowest-index positive reduced cost (guarantees finite termination).
type Bland struct{}

func (Bland) EnteringIndex(cr []float64, tol float64) (int, bool) {
	for j, v := range cr {
		if v > tol {
			return j, true
		}
	}
	return 0, false
}

// LeavingIndex performs the lexicographic ratio test to select the leaving variable.
//
// Parameters:
//   - xb: values of the current basic variables (length m)
//   - d: entering column in basis coordinates, B⁻¹Aⱼ (length m)
//   - BinvN: full direction matrix B⁻¹N stored row-major (m × numNonbasic)
//   - numNonbasic: number of non-basic variables (columns of BinvN)
//
// Returns the row index (into basic) of the leaving variable, and false if
// every component of d is non-positive (problem is unbounded).
func LeavingIndex(xb, d, BinvN []float64, m, numNonbasic int, tol float64) (int, bool) {
	// collect candidates where d[i] > tol
	candidates := make([]int, 0, m)
	for i := 0; i < m; i++ {
		if d[i] > tol {
			candidates = append(candidates, i)
		}
	}
	if len(candidates) == 0 {
		return 0, false
	}

	best := candidates[0]
	for _, i := range candidates[1:] {
		ratioi := xb[i] / d[i]
		ratiob := xb[best] / d[best]

		if ratioi < ratiob-tol {
			best = i
		} else if math.Abs(ratioi-ratiob) <= tol {
			// lexicographic tie-break on scaled rows of BinvN
			for k := 0; k < numNonbasic; k++ {
				vi := BinvN[i*numNonbasic+k] / d[i]
				vb := BinvN[best*numNonbasic+k] / d[best]
				if vi < vb-1e-14 {
					best = i
					break
				} else if vi > vb+1e-14 {
					break
				}
			}
		}
	}
	return best, true
}
