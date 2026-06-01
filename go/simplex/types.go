package simplex

// LinearProgram represents a linear program in standard equality form:
//
//	maximize  c'x
//	subject to Ax = b,  x >= 0
type LinearProgram struct {
	A    []float64 // m×n row-major
	B    []float64 // length m
	C    []float64 // length n
	M, N int
}

// BasisState holds the index partition of variables into basic and non-basic.
// Both slices hold 0-based column indices into the LP's A matrix.
type BasisState struct {
	Basic    []int
	Nonbasic []int
}

// SimplexResult is the complete output of Solve or SolvePhase2.
type SimplexResult struct {
	X          []float64
	Z          float64
	Status     SolveStatus
	Iterations int
}

// SolveStatus is the termination status of a Simplex solve.
type SolveStatus int

const (
	Optimal    SolveStatus =  1
	Unbounded  SolveStatus = -1
	Infeasible SolveStatus = -2
)

func (s SolveStatus) String() string {
	switch s {
	case Optimal:
		return "OPTIMAL"
	case Unbounded:
		return "UNBOUNDED"
	case Infeasible:
		return "INFEASIBLE"
	default:
		return "UNKNOWN"
	}
}
