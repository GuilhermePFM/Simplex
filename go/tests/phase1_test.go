package simplex_test

import (
	"math"
	"testing"

	"github.com/simplex/go/simplex"
)

func TestSolveProductionLP(t *testing.T) {
	lp := simplex.LinearProgram{
		A: []float64{
			2, 1, 1, 0,
			1, 2, 0, 1,
		},
		B: []float64{4, 4},
		C: []float64{4, 3, 0, 0},
		M: 2, N: 4,
	}

	result := simplex.Solve(lp)

	if result.Status != simplex.Optimal {
		t.Fatalf("expected Optimal, got %s", result.Status)
	}
	if !approxEq(result.Z, 28.0/3.0, tol) {
		t.Errorf("expected z≈%v, got %v", 28.0/3.0, result.Z)
	}
}

func TestSolveNegativeBRow(t *testing.T) {
	// Has one negative RHS row: -x1 - x2 + s3 = -1 (i.e., x1+x2 >= 1)
	lp := simplex.LinearProgram{
		A: []float64{
			2, 1, 1, 0, 0,
			1, 2, 0, 1, 0,
			-1, -1, 0, 0, 1,
		},
		B: []float64{4, 4, -1},
		C: []float64{4, 3, 0, 0, 0},
		M: 3, N: 5,
	}

	result := simplex.Solve(lp)

	if result.Status != simplex.Optimal {
		t.Fatalf("expected Optimal, got %s", result.Status)
	}
	// Verify feasibility: Ax ≈ b
	for i := 0; i < lp.M; i++ {
		rowSum := 0.0
		for j := 0; j < lp.N; j++ {
			rowSum += lp.A[i*lp.N+j] * result.X[j]
		}
		if math.Abs(rowSum-lp.B[i]) > tol {
			t.Errorf("row %d: Ax=%v, b=%v", i, rowSum, lp.B[i])
		}
	}
	// Non-negativity
	for j, v := range result.X {
		if v < -tol {
			t.Errorf("x[%d]=%v < 0", j, v)
		}
	}
}

func TestSolveContradictory(t *testing.T) {
	// x1 + x2 = 1 and x1 + x2 = 0: infeasible
	lp := simplex.LinearProgram{
		A: []float64{
			1, 1,
			1, 1,
		},
		B: []float64{1, 0},
		C: []float64{1, 0},
		M: 2, N: 2,
	}

	result := simplex.Solve(lp)

	if result.Status != simplex.Infeasible {
		t.Fatalf("expected Infeasible, got %s", result.Status)
	}
}

func TestSolveConflictingBounds(t *testing.T) {
	// x1 + x2 <= 1 and x1 >= 2 (via equality form):
	//   x1 + x2 + s1 = 1
	//  -x1 + s2 = -2  (i.e., x1 - s2 = 2)
	lp := simplex.LinearProgram{
		A: []float64{
			1, 1, 1, 0,
			-1, 0, 0, 1,
		},
		B: []float64{1, -2},
		C: []float64{1, 0, 0, 0},
		M: 2, N: 4,
	}

	result := simplex.Solve(lp)

	if result.Status != simplex.Infeasible {
		t.Fatalf("expected Infeasible, got %s", result.Status)
	}
}
