use criterion::{criterion_group, criterion_main, Criterion};
use nalgebra::{DMatrix, DVector};
use simplex::{solve, solve_phase2, BasisState, LargestCoefficient, LinearProgram, LogLevel};

struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        // LCG parameters from Numerical Recipes
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn random_lp(m: usize, n: usize, seed: u64) -> LinearProgram {
    let mut rng = Lcg::new(seed);
    let orig = n - m; // number of original variables

    // A = [A_raw | I_m] where A_raw is m x orig
    let mut a_data = vec![0.0f64; m * n];
    for i in 0..m {
        for j in 0..orig {
            a_data[i * n + j] = rng.next_f64();
        }
        a_data[i * n + orig + i] = 1.0;
    }
    let a = DMatrix::from_row_slice(m, n, &a_data);

    // b: strictly positive
    let b_data: Vec<f64> = (0..m).map(|_| rng.next_f64() + 0.1).collect();
    let b = DVector::from_vec(b_data);

    // c: only original vars in objective
    let mut c_data: Vec<f64> = (0..orig).map(|_| rng.next_f64()).collect();
    c_data.extend(std::iter::repeat(0.0).take(m));
    let c = DVector::from_vec(c_data);

    LinearProgram::new(a, b, c)
}

fn small_lp() -> LinearProgram {
    LinearProgram::new(
        DMatrix::from_row_slice(2, 4, &[2.0, 1.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0]),
        DVector::from_vec(vec![4.0, 4.0]),
        DVector::from_vec(vec![4.0, 3.0, 0.0, 0.0]),
    )
}

fn bench_small_phase2(c: &mut Criterion) {
    let lp = small_lp();
    let basis = BasisState::new(vec![2, 3], vec![0, 1]);
    c.bench_function("small/phase2_only", |b| {
        b.iter(|| solve_phase2(&lp, &basis, &LargestCoefficient, LogLevel::Silent))
    });
}

fn bench_small_two_phase(c: &mut Criterion) {
    c.bench_function("small/two_phase", |b| {
        b.iter(|| {
            let lp = small_lp();
            solve(lp, &LargestCoefficient, LogLevel::Silent, None)
        })
    });
}

fn bench_random(c: &mut Criterion) {
    for &(m, n) in &[(10usize, 20usize), (50, 100), (100, 200)] {
        let lp = random_lp(m, n, 42);
        let tag = format!("random/m={}_n={}", m, n);
        c.bench_function(&tag, |b| {
            b.iter(|| {
                let lp_clone = lp.clone();
                solve(lp_clone, &LargestCoefficient, LogLevel::Silent, None)
            })
        });
    }
}

criterion_group!(benches, bench_small_phase2, bench_small_two_phase, bench_random);
criterion_main!(benches);
