// Two-phase Revised Simplex with hand-rolled LU factorization (no nalgebra).
// Usage: ./bench_mps_lu <data.json> [N=10]

use serde::Deserialize;
use std::fs;
use std::time::Instant;

#[derive(Deserialize)]
struct LpData {
    #[allow(dead_code)]
    name: String,
    m: usize,
    n: usize,
    #[serde(rename = "A")]
    a: Vec<Vec<f64>>,
    b: Vec<f64>,
    c: Vec<f64>,
}

// ---------------------------------------------------------------------------
// LU factorization with partial pivoting — flat row-major Vec<f64>
// ---------------------------------------------------------------------------

struct LuFact {
    m: Vec<f64>,   // packed LU sz×sz row-major
    perm: Vec<usize>,
    sz: usize,
}

fn lu_factor(b: &[f64], sz: usize) -> LuFact {
    let mut m = b.to_vec();
    let mut perm: Vec<usize> = (0..sz).collect();
    for k in 0..sz {
        let (mut r, mut maxv) = (k, m[k * sz + k].abs());
        for i in k + 1..sz {
            let v = m[i * sz + k].abs();
            if v > maxv {
                maxv = v;
                r = i;
            }
        }
        if r != k {
            for j in 0..sz {
                m.swap(k * sz + j, r * sz + j);
            }
            perm.swap(k, r);
        }
        if m[k * sz + k].abs() < 1e-14 {
            continue;
        }
        for i in k + 1..sz {
            m[i * sz + k] /= m[k * sz + k];
            for j in k + 1..sz {
                let mul = m[i * sz + k];
                m[i * sz + j] -= mul * m[k * sz + j];
            }
        }
    }
    LuFact { m, perm, sz }
}

impl LuFact {
    fn solve_vec(&self, b: &[f64]) -> Vec<f64> {
        let sz = self.sz;
        let mut y: Vec<f64> = self.perm.iter().map(|&i| b[i]).collect();
        for i in 1..sz {
            for k in 0..i {
                let v = self.m[i * sz + k] * y[k];
                y[i] -= v;
            }
        }
        let mut x = y;
        for i in (0..sz).rev() {
            for k in i + 1..sz {
                let v = self.m[i * sz + k] * x[k];
                x[i] -= v;
            }
            let d = self.m[i * sz + i];
            if d.abs() > 1e-14 {
                x[i] /= d;
            }
        }
        x
    }

    // N is m×nb row-major; returns m×nb row-major.
    fn solve_mat(&self, n_mat: &[f64], nb: usize) -> Vec<f64> {
        let m = self.sz;
        let mut res = vec![0f64; m * nb];
        let mut col = vec![0f64; m];
        for j in 0..nb {
            for i in 0..m {
                col[i] = n_mat[i * nb + j];
            }
            let x = self.solve_vec(&col);
            for i in 0..m {
                res[i * nb + j] = x[i];
            }
        }
        res
    }
}

// ---------------------------------------------------------------------------
// Matrix helpers
// ---------------------------------------------------------------------------

fn cols_rm(a: &[f64], m: usize, n: usize, idx: &[usize]) -> Vec<f64> {
    let k = idx.len();
    let mut out = vec![0f64; m * k];
    for (jj, &col) in idx.iter().enumerate() {
        for i in 0..m {
            out[i * k + jj] = a[i * n + col];
        }
    }
    out
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn gather(a: &[f64], idx: &[usize]) -> Vec<f64> {
    idx.iter().map(|&i| a[i]).collect()
}

fn reduced_costs(c_n: &[f64], c_b: &[f64], binv_n: &[f64], m: usize, nb: usize) -> Vec<f64> {
    let mut cr = c_n.to_vec();
    for j in 0..nb {
        for i in 0..m {
            cr[j] -= binv_n[i * nb + j] * c_b[i];
        }
    }
    cr
}

fn dense_col(binv_n: &[f64], m: usize, nb: usize, j: usize) -> Vec<f64> {
    (0..m).map(|i| binv_n[i * nb + j]).collect()
}

// ---------------------------------------------------------------------------
// Pivot selection
// ---------------------------------------------------------------------------

fn entering_index(cr: &[f64], tol: f64) -> Option<usize> {
    let (mut best, mut best_val) = (usize::MAX, tol);
    for (j, &v) in cr.iter().enumerate() {
        if v > best_val {
            best_val = v;
            best = j;
        }
    }
    if best == usize::MAX { None } else { Some(best) }
}

fn leaving_index(xb: &[f64], d: &[f64], binv_n: &[f64], m: usize, nb: usize, tol: f64) -> Option<usize> {
    let cands: Vec<usize> = (0..m).filter(|&i| d[i] > tol).collect();
    if cands.is_empty() {
        return None;
    }
    let mut best = cands[0];
    for &i in &cands[1..] {
        let ri = xb[i] / d[i];
        let rb = xb[best] / d[best];
        if ri < rb - tol {
            best = i;
        } else if (ri - rb).abs() <= tol {
            for k in 0..nb {
                let vi = binv_n[i * nb + k] / d[i];
                let vb = binv_n[best * nb + k] / d[best];
                if vi < vb - 1e-14 {
                    best = i;
                    break;
                } else if vi > vb + 1e-14 {
                    break;
                }
            }
        }
    }
    Some(best)
}

// ---------------------------------------------------------------------------
// Phase 1
// ---------------------------------------------------------------------------

fn make_b_nonneg(a: &mut Vec<f64>, b: &mut Vec<f64>, m: usize, n: usize) {
    for i in 0..m {
        if b[i] < 0.0 {
            b[i] = -b[i];
            for j in 0..n {
                a[i * n + j] = -a[i * n + j];
            }
        }
    }
}

fn phase1(
    a_orig: &[f64], b_orig: &[f64], m: usize, n: usize,
    maxit: usize, tol: f64,
) -> Option<(Vec<usize>, Vec<usize>)> {
    let mut a = a_orig.to_vec();
    let mut b = b_orig.to_vec();
    make_b_nonneg(&mut a, &mut b, m, n);

    let nw = n + m;
    let mut aw = vec![0f64; m * nw];
    for i in 0..m {
        for j in 0..n {
            aw[i * nw + j] = a[i * n + j];
        }
        aw[i * nw + n + i] = 1.0;
    }
    let mut cw = vec![0f64; nw];
    for i in n..nw {
        cw[i] = -1.0;
    }

    let mut basic: Vec<usize> = (n..n + m).collect();
    let mut nonbasic: Vec<usize> = (0..n).collect();

    for _ in 0..maxit {
        let nb = nonbasic.len();
        let b_mat = cols_rm(&aw, m, nw, &basic);
        let n_mat = cols_rm(&aw, m, nw, &nonbasic);
        let lu = lu_factor(&b_mat, m);
        let xb = lu.solve_vec(&b);
        let binv_n = lu.solve_mat(&n_mat, nb);
        let cb = gather(&cw, &basic);
        let cn = gather(&cw, &nonbasic);
        let cr = reduced_costs(&cn, &cb, &binv_n, m, nb);
        let z = dot(&cb, &xb);

        match entering_index(&cr, tol) {
            None => {
                if z < -tol {
                    return None;
                }
                fix_arts(&mut aw, &mut basic, &mut nonbasic, m, nw, n);
                let bas: Vec<usize> = basic.iter().copied().filter(|&v| v < n).collect();
                let nonbas: Vec<usize> = nonbasic.iter().copied().filter(|&v| v < n).collect();
                return Some((bas, nonbas));
            }
            Some(j) => {
                let d = dense_col(&binv_n, m, nb, j);
                match leaving_index(&xb, &d, &binv_n, m, nb, tol) {
                    None => panic!("Phase 1 unexpectedly unbounded"),
                    Some(i) => {
                        basic.swap(i, i); // keep borrow checker happy
                        let tmp = basic[i];
                        basic[i] = nonbasic[j];
                        nonbasic[j] = tmp;
                    }
                }
            }
        }
    }
    panic!("Phase 1 did not converge");
}

fn fix_arts(aw: &mut Vec<f64>, basic: &mut Vec<usize>, nonbasic: &mut Vec<usize>,
            m: usize, nw: usize, art_start: usize) {
    loop {
        let mut changed = false;
        for pos in 0..basic.len() {
            if basic[pos] < art_start {
                continue;
            }
            let b_mat = cols_rm(aw, m, nw, basic);
            let lu = lu_factor(&b_mat, m);
            let binv_aw = lu.solve_mat(aw, nw);
            let swap_j = nonbasic.iter().position(|&nc| {
                nc < art_start && binv_aw[pos * nw + nc].abs() > 1e-10
            });
            match swap_j {
                None => {
                    for j in 0..nw {
                        aw[pos * nw + j] = 0.0;
                    }
                }
                Some(j) => {
                    let tmp = basic[pos];
                    basic[pos] = nonbasic[j];
                    nonbasic[j] = tmp;
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// Phase 2
// ---------------------------------------------------------------------------

fn phase2(
    a: &[f64], b: &[f64], c: &[f64], m: usize, n: usize,
    basic_in: &[usize], nonbasic_in: &[usize],
    maxit: usize, tol: f64,
) -> (f64, usize, bool) {
    let mut basic = basic_in.to_vec();
    let mut nonbasic = nonbasic_in.to_vec();

    for it in 1..=maxit {
        let nb = nonbasic.len();
        let b_mat = cols_rm(a, m, n, &basic);
        let n_mat = cols_rm(a, m, n, &nonbasic);
        let lu = lu_factor(&b_mat, m);
        let xb = lu.solve_vec(b);
        let binv_n = lu.solve_mat(&n_mat, nb);
        let cb = gather(c, &basic);
        let cn = gather(c, &nonbasic);
        let cr = reduced_costs(&cn, &cb, &binv_n, m, nb);
        let z = dot(&cb, &xb);

        match entering_index(&cr, tol) {
            None => return (z, it, true),
            Some(j) => {
                let d = dense_col(&binv_n, m, nb, j);
                match leaving_index(&xb, &d, &binv_n, m, nb, tol) {
                    None => return (f64::INFINITY, it, false),
                    Some(i) => {
                        let tmp = basic[i];
                        basic[i] = nonbasic[j];
                        nonbasic[j] = tmp;
                    }
                }
            }
        }
    }
    panic!("Phase 2 did not converge");
}

fn simplex_solve(a: &[f64], b: &[f64], c: &[f64], m: usize, n: usize) -> (f64, usize) {
    match phase1(a, b, m, n, 10000, 1e-8) {
        None => (0.0, 0),
        Some((basic, mut nonbasic)) => {
            let in_set: std::collections::HashSet<usize> =
                basic.iter().chain(nonbasic.iter()).copied().collect();
            for j in 0..n {
                if !in_set.contains(&j) {
                    nonbasic.push(j);
                }
            }
            let (z, it, _) = phase2(a, b, c, m, n, &basic, &nonbasic, 10000, 1e-10);
            (z, it)
        }
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bench_mps_lu <data.json> [N]");
        std::process::exit(1);
    }
    let n_runs: usize = if args.len() >= 3 { args[2].parse().unwrap_or(10) } else { 10 };

    let raw = fs::read_to_string(&args[1]).expect("Cannot read JSON file");
    let data: LpData = serde_json::from_str(&raw).expect("Cannot parse JSON");

    let m = data.m;
    let n = data.n;
    let a_flat: Vec<f64> = (0..m).flat_map(|i| data.a[i].iter().copied()).collect();
    let b = &data.b;
    let c = &data.c;

    let mut last_z = 0.0f64;
    let mut last_it = 0usize;
    let mut times_ms: Vec<f64> = Vec::with_capacity(n_runs);

    for _ in 0..n_runs {
        let t0 = Instant::now();
        let (z, it) = simplex_solve(&a_flat, b, c, m, n);
        times_ms.push(t0.elapsed().as_secs_f64() * 1000.0);
        last_z = z;
        last_it = it;
    }

    let times_str = times_ms
        .iter()
        .map(|t| format!("{:.6}", t))
        .collect::<Vec<_>>()
        .join(", ");

    println!(
        "{{\"z\": {:.4}, \"iterations\": {}, \"times_ms\": [{}]}}",
        last_z, last_it, times_str
    );
}
