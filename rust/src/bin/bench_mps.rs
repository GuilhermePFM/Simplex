use nalgebra::{DMatrix, DVector};
use serde::Deserialize;
use simplex::{solve, LargestCoefficient, LinearProgram, LogLevel};
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bench_mps <data.json> [N]");
        std::process::exit(1);
    }
    let json_path = &args[1];
    let n_runs: usize = if args.len() >= 3 {
        args[2].parse().unwrap_or(10)
    } else {
        10
    };

    let raw = fs::read_to_string(json_path).expect("Cannot read JSON file");
    let data: LpData = serde_json::from_str(&raw).expect("Cannot parse JSON");

    let m = data.m;
    let n = data.n;

    // Build A (row-major: data.a[i][j])
    let a_flat: Vec<f64> = (0..m)
        .flat_map(|i| data.a[i].iter().copied())
        .collect();
    let a = DMatrix::from_row_slice(m, n, &a_flat);
    let b = DVector::from_vec(data.b);
    let c = DVector::from_vec(data.c);

    let rule = LargestCoefficient;
    let mut last_z = 0.0f64;
    let mut last_iterations = 0usize;
    let mut times_ms: Vec<f64> = Vec::with_capacity(n_runs);

    for _ in 0..n_runs {
        let lp = LinearProgram::new(a.clone(), b.clone(), c.clone());
        let t0 = Instant::now();
        let result = solve(lp, &rule, LogLevel::Silent, None);
        let elapsed = t0.elapsed();
        times_ms.push(elapsed.as_secs_f64() * 1000.0);
        last_z = result.z;
        last_iterations = result.iterations;
    }

    // Build JSON output manually to avoid float formatting surprises
    let times_str = times_ms
        .iter()
        .map(|t| format!("{:.6}", t))
        .collect::<Vec<_>>()
        .join(", ");

    println!(
        "{{\"z\": {:.4}, \"iterations\": {}, \"times_ms\": [{}]}}",
        last_z, last_iterations, times_str
    );
}
