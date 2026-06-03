// Two-phase Revised Simplex with hand-rolled LU factorization (no gonum).
// Usage: ./bench_mps_lu <data.json> [N=10]
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"strconv"
	"time"
)

// ---------------------------------------------------------------------------
// LU factorization with partial pivoting — flat row-major []float64
// ---------------------------------------------------------------------------

type luFact struct {
	m    []float64 // packed LU, row-major sz×sz
	perm []int
	sz   int
}

func luFactor(B []float64, sz int) luFact {
	m := make([]float64, sz*sz)
	copy(m, B)
	perm := make([]int, sz)
	for i := range perm {
		perm[i] = i
	}
	for k := 0; k < sz; k++ {
		r, maxv := k, math.Abs(m[k*sz+k])
		for i := k + 1; i < sz; i++ {
			if v := math.Abs(m[i*sz+k]); v > maxv {
				maxv = v
				r = i
			}
		}
		if r != k {
			for j := 0; j < sz; j++ {
				m[k*sz+j], m[r*sz+j] = m[r*sz+j], m[k*sz+j]
			}
			perm[k], perm[r] = perm[r], perm[k]
		}
		if math.Abs(m[k*sz+k]) < 1e-14 {
			continue
		}
		for i := k + 1; i < sz; i++ {
			m[i*sz+k] /= m[k*sz+k]
			for j := k + 1; j < sz; j++ {
				m[i*sz+j] -= m[i*sz+k] * m[k*sz+j]
			}
		}
	}
	return luFact{m, perm, sz}
}

func (lu luFact) solveVec(b []float64) []float64 {
	sz := lu.sz
	y := make([]float64, sz)
	for i := range y {
		y[i] = b[lu.perm[i]]
	}
	for i := 1; i < sz; i++ {
		for k := 0; k < i; k++ {
			y[i] -= lu.m[i*sz+k] * y[k]
		}
	}
	x := make([]float64, sz)
	copy(x, y)
	for i := sz - 1; i >= 0; i-- {
		for k := i + 1; k < sz; k++ {
			x[i] -= lu.m[i*sz+k] * x[k]
		}
		if math.Abs(lu.m[i*sz+i]) > 1e-14 {
			x[i] /= lu.m[i*sz+i]
		}
	}
	return x
}

// solveMat solves B·X = N where N is m×nb (row-major); returns m×nb row-major.
func (lu luFact) solveMat(N []float64, nb int) []float64 {
	m := lu.sz
	res := make([]float64, m*nb)
	col := make([]float64, m)
	for j := 0; j < nb; j++ {
		for i := 0; i < m; i++ {
			col[i] = N[i*nb+j]
		}
		x := lu.solveVec(col)
		for i := 0; i < m; i++ {
			res[i*nb+j] = x[i]
		}
	}
	return res
}

// ---------------------------------------------------------------------------
// Matrix helpers
// ---------------------------------------------------------------------------

// colsRM extracts columns (indices) from a row-major m×n matrix into m×k row-major.
func colsRM(A []float64, m, n int, idx []int) []float64 {
	k := len(idx)
	out := make([]float64, m*k)
	for jj, col := range idx {
		for i := 0; i < m; i++ {
			out[i*k+jj] = A[i*n+col]
		}
	}
	return out
}

func dot(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

func gather(a []float64, idx []int) []float64 {
	out := make([]float64, len(idx))
	for i, j := range idx {
		out[i] = a[j]
	}
	return out
}

// reducedCosts: cN - BinvN'·cB (BinvN is m×nb row-major)
func reducedCosts(cN, cB, BinvN []float64, m, nb int) []float64 {
	cr := make([]float64, nb)
	copy(cr, cN)
	for j := 0; j < nb; j++ {
		for i := 0; i < m; i++ {
			cr[j] -= BinvN[i*nb+j] * cB[i]
		}
	}
	return cr
}

func denseCol(BinvN []float64, m, nb, j int) []float64 {
	col := make([]float64, m)
	for i := 0; i < m; i++ {
		col[i] = BinvN[i*nb+j]
	}
	return col
}

// ---------------------------------------------------------------------------
// Pivot selection
// ---------------------------------------------------------------------------

func enteringIndex(cr []float64, tol float64) (int, bool) {
	best, bestVal := -1, tol
	for j, v := range cr {
		if v > bestVal {
			bestVal = v
			best = j
		}
	}
	return best, best >= 0
}

func leavingIndex(xb, d, BinvN []float64, m, nb int, tol float64) (int, bool) {
	var cands []int
	for i := 0; i < m; i++ {
		if d[i] > tol {
			cands = append(cands, i)
		}
	}
	if len(cands) == 0 {
		return 0, false
	}
	best := cands[0]
	for _, i := range cands[1:] {
		ri := xb[i] / d[i]
		rb := xb[best] / d[best]
		if ri < rb-tol {
			best = i
		} else if math.Abs(ri-rb) <= tol {
			for k := 0; k < nb; k++ {
				vi := BinvN[i*nb+k] / d[i]
				vb := BinvN[best*nb+k] / d[best]
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

// ---------------------------------------------------------------------------
// Phase 1
// ---------------------------------------------------------------------------

func makeBNonneg(A []float64, b []float64, m, n int) ([]float64, []float64) {
	A2 := make([]float64, len(A))
	copy(A2, A)
	b2 := make([]float64, len(b))
	copy(b2, b)
	for i := 0; i < m; i++ {
		if b2[i] < 0 {
			b2[i] = -b2[i]
			for j := 0; j < n; j++ {
				A2[i*n+j] = -A2[i*n+j]
			}
		}
	}
	return A2, b2
}

func phase1(A []float64, b []float64, m, n int, maxit int, tol float64) ([]int, []int, bool) {
	A, b = makeBNonneg(A, b, m, n)
	nw := n + m
	Aw := make([]float64, m*nw)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			Aw[i*nw+j] = A[i*n+j]
		}
		Aw[i*nw+n+i] = 1.0
	}
	cw := make([]float64, nw)
	for i := n; i < nw; i++ {
		cw[i] = -1.0
	}

	basic := make([]int, m)
	for i := range basic {
		basic[i] = n + i
	}
	nonbasic := make([]int, n)
	for i := range nonbasic {
		nonbasic[i] = i
	}

	for it := 0; it < maxit; it++ {
		nb := len(nonbasic)
		B := colsRM(Aw, m, nw, basic)
		N := colsRM(Aw, m, nw, nonbasic)
		lu := luFactor(B, m)
		xb := lu.solveVec(b)
		BinvN := lu.solveMat(N, nb)
		cb := gather(cw, basic)
		cr := reducedCosts(gather(cw, nonbasic), cb, BinvN, m, nb)
		z := dot(cb, xb)
		j, ok := enteringIndex(cr, tol)
		if !ok {
			if z < -tol {
				return nil, nil, false
			}
			fixArts(Aw, basic, nonbasic, m, nw, n)
			var bas, nonbas []int
			for _, v := range basic {
				if v < n {
					bas = append(bas, v)
				}
			}
			for _, v := range nonbasic {
				if v < n {
					nonbas = append(nonbas, v)
				}
			}
			return bas, nonbas, true
		}
		d := denseCol(BinvN, m, nb, j)
		i, ok := leavingIndex(xb, d, BinvN, m, nb, tol)
		if !ok {
			panic("Phase 1 unexpectedly unbounded")
		}
		basic[i], nonbasic[j] = nonbasic[j], basic[i]
	}
	panic("Phase 1 did not converge")
}

func fixArts(Aw []float64, basic, nonbasic []int, m, nw, artStart int) {
	for changed := true; changed; {
		changed = false
		for pos, col := range basic {
			if col < artStart {
				continue
			}
			B := colsRM(Aw, m, nw, basic)
			lu := luFactor(B, m)
			BinvAw := lu.solveMat(Aw, nw)
			swapJ := -1
			for j, nc := range nonbasic {
				if nc < artStart && math.Abs(BinvAw[pos*nw+nc]) > 1e-10 {
					swapJ = j
					break
				}
			}
			if swapJ < 0 {
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

// ---------------------------------------------------------------------------
// Phase 2
// ---------------------------------------------------------------------------

func phase2(A []float64, b, c []float64, m, n int,
	basic, nonbasic []int, maxit int, tol float64) (float64, int, bool) {
	bas := make([]int, len(basic))
	copy(bas, basic)
	nonbas := make([]int, len(nonbasic))
	copy(nonbas, nonbasic)

	for it := 1; it <= maxit; it++ {
		nb := len(nonbas)
		B := colsRM(A, m, n, bas)
		N := colsRM(A, m, n, nonbas)
		lu := luFactor(B, m)
		xb := lu.solveVec(b)
		BinvN := lu.solveMat(N, nb)
		cb := gather(c, bas)
		cr := reducedCosts(gather(c, nonbas), cb, BinvN, m, nb)
		z := dot(cb, xb)
		j, ok := enteringIndex(cr, tol)
		if !ok {
			return z, it, true
		}
		d := denseCol(BinvN, m, nb, j)
		i, ok := leavingIndex(xb, d, BinvN, m, nb, tol)
		if !ok {
			return math.Inf(1), it, false
		}
		bas[i], nonbas[j] = nonbas[j], bas[i]
	}
	panic("Phase 2 did not converge")
}

func solve(A, b, c []float64, m, n int) (float64, int, bool) {
	basic, nonbasic, feasible := phase1(A, b, m, n, 10000, 1e-8)
	if !feasible {
		return 0, 0, false
	}
	// Fill in any columns not accounted for
	inSet := make(map[int]bool)
	for _, v := range basic {
		inSet[v] = true
	}
	for _, v := range nonbasic {
		inSet[v] = true
	}
	for j := 0; j < n; j++ {
		if !inSet[j] {
			nonbasic = append(nonbasic, j)
		}
	}
	z, it, opt := phase2(A, b, c, m, n, basic, nonbasic, 10000, 1e-10)
	return z, it, opt
}

// ---------------------------------------------------------------------------
// JSON I/O
// ---------------------------------------------------------------------------

type lpData struct {
	Name string      `json:"name"`
	M    int         `json:"m"`
	N    int         `json:"n"`
	A    [][]float64 `json:"A"`
	B    []float64   `json:"b"`
	C    []float64   `json:"c"`
}

func main() {
	args := os.Args[1:]
	if len(args) < 1 {
		fmt.Fprintln(os.Stderr, "Usage: bench_mps_lu <data.json> [N=10]")
		os.Exit(1)
	}
	dataPath := args[0]
	nRuns := 10
	if len(args) >= 2 {
		if v, err := strconv.Atoi(args[1]); err == nil {
			nRuns = v
		}
	}

	raw, err := os.ReadFile(dataPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "error: %v\n", err)
		os.Exit(1)
	}
	var data lpData
	if err := json.Unmarshal(raw, &data); err != nil {
		fmt.Fprintf(os.Stderr, "error: %v\n", err)
		os.Exit(1)
	}

	m, n := data.M, data.N
	aFlat := make([]float64, m*n)
	for i, row := range data.A {
		copy(aFlat[i*n:], row)
	}
	b := data.B
	c := data.C

	timesMs := make([]float64, nRuns)
	var lastZ float64
	var lastIter int
	for i := 0; i < nRuns; i++ {
		t0 := time.Now()
		z, it, _ := solve(aFlat, b, c, m, n)
		timesMs[i] = float64(time.Since(t0).Nanoseconds()) / 1e6
		lastZ = z
		lastIter = it
	}

	timesBuf := "["
	for i, ms := range timesMs {
		if i > 0 {
			timesBuf += ","
		}
		timesBuf += strconv.FormatFloat(ms, 'f', 4, 64)
	}
	timesBuf += "]"

	fmt.Printf(`{"z":%s,"iterations":%d,"times_ms":%s}`+"\n",
		strconv.FormatFloat(lastZ, 'f', 4, 64), lastIter, timesBuf)
}
