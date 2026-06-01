package simplex

import (
	"fmt"
	"io"
	"os"
)

// LogLevel controls how much the SimplexLogger outputs.
type LogLevel int

const (
	Silent LogLevel = iota
	Info
	Debug
)

// SimplexLogger writes solver progress to a stream.
type SimplexLogger struct {
	Level  LogLevel
	stream io.Writer
}

// NewSilentLogger returns a logger that produces no output.
func NewSilentLogger() SimplexLogger {
	return SimplexLogger{Level: Silent, stream: io.Discard}
}

// NewLogger returns a logger writing to stdout at the given level.
func NewLogger(level LogLevel) SimplexLogger {
	return SimplexLogger{Level: level, stream: os.Stdout}
}

// NewFileLogger returns a logger writing to the given file path at the given level.
func NewFileLogger(path string, level LogLevel) (SimplexLogger, error) {
	f, err := os.Create(path)
	if err != nil {
		return SimplexLogger{}, err
	}
	return SimplexLogger{Level: level, stream: f}, nil
}

func (l SimplexLogger) write(s string) {
	if l.stream != nil {
		fmt.Fprintln(l.stream, s)
	}
}

// LogProblem logs the LP formulation if level >= Info.
func (l SimplexLogger) LogProblem(lp LinearProgram) {
	if l.Level == Silent {
		return
	}
	l.write("=======================")
	l.write(" Revised Simplex Solver")
	l.write("=======================")
	l.write(fmt.Sprintf("A = %v", lp.A))
	l.write(fmt.Sprintf("b = %v", lp.B))
	l.write(fmt.Sprintf("c = %v", lp.C))
	l.write("")
}

// LogPhase logs a phase header if level >= Info.
func (l SimplexLogger) LogPhase(phase int) {
	if l.Level == Silent {
		return
	}
	var header string
	if phase == 1 {
		header = "Phase 1 — Finding initial BFS"
	} else {
		header = "Phase 2 — Optimizing"
	}
	l.write(header)
}

// LogIteration logs iteration details if level == Debug.
func (l SimplexLogger) LogIteration(it int, basis BasisState, xb []float64, z float64) {
	if l.Level != Debug {
		return
	}
	l.write(fmt.Sprintf("  iter %d:", it))
	l.write(fmt.Sprintf("    xb       = %v", xb))
	l.write(fmt.Sprintf("    basic    = %v", basis.Basic))
	l.write(fmt.Sprintf("    nonbasic = %v", basis.Nonbasic))
	l.write(fmt.Sprintf("    z        = %v", z))
	l.write("")
}

// LogResult logs the solve result if level >= Info.
func (l SimplexLogger) LogResult(result SimplexResult) {
	if l.Level == Silent {
		return
	}
	l.write("")
	l.write(fmt.Sprintf("Result: %s", result.Status))
	l.write(fmt.Sprintf("  x          = %v", result.X))
	l.write(fmt.Sprintf("  z          = %v", result.Z))
	l.write(fmt.Sprintf("  iterations = %d", result.Iterations))
}

// Close closes the underlying stream if it is a file.
func (l SimplexLogger) Close() {
	if f, ok := l.stream.(*os.File); ok && f != os.Stdout && f != os.Stderr {
		f.Close()
	}
}
