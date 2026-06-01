package benchmark_test

import (
	"testing"

	"github.com/simplex/go/simplex"
)

// randomLP generates a random LP with m constraints and n total variables (n > m).
// Appends an m×m identity block as slack columns so the LP is in standard form
// with an obvious initial BFS. Uses a simple LCG for reproducibility.
func randomLP(m, n, seed int) simplex.LinearProgram {
	state := uint64(seed)
	nextFloat := func() float64 {
		state = state*6364136223846793005 + 1442695040888963407
		return float64(state>>11) / float64(1<<53)
	}

	orig := n - m // number of original variables
	// A = [A_raw | I_m], size m×n
	A := make([]float64, m*n)
	for i := 0; i < m; i++ {
		for j := 0; j < orig; j++ {
			A[i*n+j] = nextFloat()
		}
		A[i*n+orig+i] = 1.0
	}

	b := make([]float64, m)
	for i := range b {
		b[i] = nextFloat() + 0.1
	}

	c := make([]float64, n)
	for j := 0; j < orig; j++ {
		c[j] = nextFloat()
	}

	return simplex.LinearProgram{A: A, B: b, C: c, M: m, N: n}
}

var smallLP = simplex.LinearProgram{
	A: []float64{
		2, 1, 1, 0,
		1, 2, 0, 1,
	},
	B: []float64{4, 4},
	C: []float64{4, 3, 0, 0},
	M: 2, N: 4,
}
var smallBasis = simplex.BasisState{Basic: []int{2, 3}, Nonbasic: []int{0, 1}}

func BenchmarkSmallTwoPhase(b *testing.B) {
	for i := 0; i < b.N; i++ {
		simplex.Solve(smallLP)
	}
}

func BenchmarkSmallPhase2Only(b *testing.B) {
	for i := 0; i < b.N; i++ {
		simplex.SolvePhase2(smallLP, smallBasis)
	}
}

func BenchmarkM10N20(b *testing.B) {
	lp := randomLP(10, 20, 42)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		simplex.Solve(lp)
	}
}

func BenchmarkM50N100(b *testing.B) {
	lp := randomLP(50, 100, 42)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		simplex.Solve(lp)
	}
}

func BenchmarkM100N200(b *testing.B) {
	lp := randomLP(100, 200, 42)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		simplex.Solve(lp)
	}
}
