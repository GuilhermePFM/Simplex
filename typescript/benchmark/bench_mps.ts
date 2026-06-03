// Usage: tsx benchmark/bench_mps.ts <data.json> [N=10]
import { readFileSync } from 'fs';
import { solve, makeLinearProgram } from '../src/index.js';

interface LPData {
  name: string;
  m: number;
  n: number;
  A: number[][];
  b: number[];
  c: number[];
  optimal_z: number;
  netlib_optimal: number;
}

const args = process.argv.slice(2);
if (args.length < 1) {
  process.stderr.write('Usage: tsx benchmark/bench_mps.ts <data.json> [N=10]\n');
  process.exit(1);
}

const dataPath = args[0];
const N = args[1] !== undefined ? parseInt(args[1], 10) : 10;

const raw: LPData = JSON.parse(readFileSync(dataPath, 'utf-8'));
const lp = makeLinearProgram(raw.A, raw.b, raw.c);

// Warmup (1 run to let V8 JIT compile)
solve(lp);

// Timed runs
const timesMs: number[] = [];
let lastResult = solve(lp); // will be overwritten

for (let i = 0; i < N; i++) {
  const t0 = performance.now();
  lastResult = solve(lp);
  timesMs.push(performance.now() - t0);
}

console.log(JSON.stringify({
  z: lastResult.z,
  iterations: lastResult.iterations,
  times_ms: timesMs,
}));
