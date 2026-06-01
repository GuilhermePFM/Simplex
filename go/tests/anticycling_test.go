package simplex_test

import (
	"math"
	"testing"

	"github.com/simplex/go/simplex"
)

func TestAnticyclingDegenerate(t *testing.T) {
	// maximize x1 + x2
	// s.t.  x1 + x2 + s1 = 2
	//       x1      + s2 = 1
	//            x2 + s3 = 1
	// Optimal at x1=x2=1, z=2. Degenerate: all slacks zero at optimum.
	lp := simplex.LinearProgram{
		A: []float64{
			1, 1, 1, 0, 0,
			1, 0, 0, 1, 0,
			0, 1, 0, 0, 1,
		},
		B: []float64{2, 1, 1},
		C: []float64{1, 1, 0, 0, 0},
		M: 3, N: 5,
	}
	// 0-based: basic=[2,3,4], nonbasic=[0,1]
	basis := simplex.BasisState{Basic: []int{2, 3, 4}, Nonbasic: []int{0, 1}}

	result := simplex.SolvePhase2(lp, basis)

	if result.Status != simplex.Optimal {
		t.Fatalf("expected Optimal, got %s", result.Status)
	}
	if !approxEq(result.Z, 2.0, tol) {
		t.Errorf("expected z≈2.0, got %v", result.Z)
	}
	if !approxEq(result.X[0], 1.0, tol) {
		t.Errorf("expected x[0]≈1.0, got %v", result.X[0])
	}
	if !approxEq(result.X[1], 1.0, tol) {
		t.Errorf("expected x[1]≈1.0, got %v", result.X[1])
	}
	if result.Iterations != 3 {
		t.Errorf("expected 3 iterations, got %d", result.Iterations)
	}

	// Feasibility: Ax = b
	for i := 0; i < lp.M; i++ {
		rowSum := 0.0
		for j := 0; j < lp.N; j++ {
			rowSum += lp.A[i*lp.N+j] * result.X[j]
		}
		if math.Abs(rowSum-lp.B[i]) > tol {
			t.Errorf("row %d: Ax=%v, b=%v", i, rowSum, lp.B[i])
		}
	}
}

func TestAnticyclingBland(t *testing.T) {
	lp := simplex.LinearProgram{
		A: []float64{
			1, 1, 1, 0, 0,
			1, 0, 0, 1, 0,
			0, 1, 0, 0, 1,
		},
		B: []float64{2, 1, 1},
		C: []float64{1, 1, 0, 0, 0},
		M: 3, N: 5,
	}
	basis := simplex.BasisState{Basic: []int{2, 3, 4}, Nonbasic: []int{0, 1}}

	result := simplex.SolvePhase2(lp, basis, simplex.WithRule(simplex.Bland{}))

	if result.Status != simplex.Optimal {
		t.Fatalf("expected Optimal, got %s", result.Status)
	}
	if !approxEq(result.Z, 2.0, tol) {
		t.Errorf("expected z≈2.0, got %v", result.Z)
	}
	if !approxEq(result.X[0], 1.0, tol) {
		t.Errorf("expected x[0]≈1.0, got %v", result.X[0])
	}
	if !approxEq(result.X[1], 1.0, tol) {
		t.Errorf("expected x[1]≈1.0, got %v", result.X[1])
	}
}

func TestKotiahSteinberg(t *testing.T) {
	// maximize 10x1 - 57x2 - 9x3 - 24x4
	// s.t.
	//   -0.5x1 + 5.5x2 + 2.5x3 - 9x4 + x5 = 0
	//   -0.5x1 + 2.5x2 + 0.5x3 -  x4 + x6 = 0
	//    x1                              + x7 = 1
	// 0-based: basic=[4,5,6], nonbasic=[0,1,2,3]
	lp := simplex.LinearProgram{
		A: []float64{
			-0.5, 5.5, 2.5, -9, 1, 0, 0,
			-0.5, 2.5, 0.5, -1, 0, 1, 0,
			1.0, 0.0, 0.0, 0, 0, 0, 1,
		},
		B: []float64{0, 0, 1},
		C: []float64{10, -57, -9, -24, 0, 0, 0},
		M: 3, N: 7,
	}
	basis := simplex.BasisState{Basic: []int{4, 5, 6}, Nonbasic: []int{0, 1, 2, 3}}

	result := simplex.SolvePhase2(lp, basis)

	if result.Status != simplex.Optimal {
		t.Fatalf("expected Optimal, got %s", result.Status)
	}
	if !approxEq(result.Z, 10.0, tol) {
		t.Errorf("expected z≈10.0, got %v", result.Z)
	}
	if !approxEq(result.X[0], 1.0, tol) {
		t.Errorf("expected x[0]≈1.0, got %v", result.X[0])
	}

	// Feasibility: Ax = b
	for i := 0; i < lp.M; i++ {
		rowSum := 0.0
		for j := 0; j < lp.N; j++ {
			rowSum += lp.A[i*lp.N+j] * result.X[j]
		}
		if math.Abs(rowSum-lp.B[i]) > tol {
			t.Errorf("row %d: Ax=%v, b=%v", i, rowSum, lp.B[i])
		}
	}
}
