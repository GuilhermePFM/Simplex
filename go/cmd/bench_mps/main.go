// Usage: ./bench_mps <data.json> [N=10]
package main

import (
	"encoding/json"
	"fmt"
	"os"
	"strconv"
	"time"

	"github.com/simplex/go/simplex"
)

type lpData struct {
	Name          string      `json:"name"`
	M             int         `json:"m"`
	N             int         `json:"n"`
	A             [][]float64 `json:"A"`
	B             []float64   `json:"b"`
	C             []float64   `json:"c"`
	OptimalZ      float64     `json:"optimal_z"`
	NetlibOptimal float64     `json:"netlib_optimal"`
}

func main() {
	args := os.Args[1:]
	if len(args) < 1 {
		fmt.Fprintln(os.Stderr, "Usage: bench_mps <data.json> [N=10]")
		os.Exit(1)
	}

	dataPath := args[0]
	n := 10
	if len(args) >= 2 {
		var err error
		n, err = strconv.Atoi(args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "invalid N: %v\n", err)
			os.Exit(1)
		}
	}

	// Read and parse JSON
	raw, err := os.ReadFile(dataPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "error reading file: %v\n", err)
		os.Exit(1)
	}

	var data lpData
	if err := json.Unmarshal(raw, &data); err != nil {
		fmt.Fprintf(os.Stderr, "error parsing JSON: %v\n", err)
		os.Exit(1)
	}

	// Flatten 2-D A into row-major []float64
	aFlat := make([]float64, data.M*data.N)
	for i, row := range data.A {
		copy(aFlat[i*data.N:], row)
	}

	lp := simplex.LinearProgram{
		A: aFlat,
		B: data.B,
		C: data.C,
		M: data.M,
		N: data.N,
	}

	// Timed runs (no warmup needed for Go)
	timesMs := make([]float64, n)
	var lastResult simplex.SimplexResult

	for i := 0; i < n; i++ {
		t0 := time.Now()
		lastResult = simplex.Solve(lp)
		timesMs[i] = float64(time.Since(t0).Nanoseconds()) / 1e6
	}

	// Build JSON output manually to avoid reflection overhead and match exact format
	timesBuf := "["
	for i, ms := range timesMs {
		if i > 0 {
			timesBuf += ","
		}
		timesBuf += strconv.FormatFloat(ms, 'f', 4, 64)
	}
	timesBuf += "]"

	fmt.Printf(`{"z":%s,"iterations":%d,"times_ms":%s}`+"\n",
		strconv.FormatFloat(lastResult.Z, 'f', 4, 64),
		lastResult.Iterations,
		timesBuf,
	)
}
