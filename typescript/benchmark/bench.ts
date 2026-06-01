import { solve, solvePhase2 } from '../src/index.js';

/** Mulberry32 seeded PRNG — returns floats in [0, 1). */
function mulberry32(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s += 0x6d2b79f5;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t ^= t + Math.imul(t ^ (t >>> 7), 61 | t);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Generate a random LP with m constraints and n total variables (n > m).
 * Mirrors the Julia random_lp: A = [A_raw | I_m], strictly positive b, original-var objective.
 */
function randomLP(m: number, n: number, seed: number = 42) {
  const rand = mulberry32(seed);
  const nOrig = n - m;

  const A: number[][] = Array.from({ length: m }, (_, i) => {
    const row = Array.from({ length: nOrig }, () => rand());
    const slack = new Array(m).fill(0);
    slack[i] = 1;
    return [...row, ...slack];
  });
  const b: number[] = Array.from({ length: m }, () => rand() + 0.1);
  const c: number[] = [...Array.from({ length: nOrig }, () => rand()), ...new Array(m).fill(0)];

  return { A, b, c };
}

function median(arr: number[]): number {
  const sorted = arr.slice().sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0
    ? (sorted[mid - 1] + sorted[mid]) / 2
    : sorted[mid];
}

function timeit(fn: () => void, runs: number = 10): number {
  const times: number[] = [];
  for (let i = 0; i < runs; i++) {
    const t0 = performance.now();
    fn();
    times.push(performance.now() - t0);
  }
  return median(times);
}

function fmt(ms: number): string {
  if (ms < 1) return `${(ms * 1000).toFixed(1)} µs`;
  return `${ms.toFixed(3)} ms`;
}

console.log('=== Simplex TypeScript Benchmark ===\n');

// Small canonical problem
const lpSmall = {
  A: [[2,1,1,0],[1,2,0,1]],
  b: [4,4],
  c: [4,3,0,0],
};
const basisSmall = { basic: [2,3], nonbasic: [0,1] };

const t1 = timeit(() => solvePhase2(lpSmall, basisSmall));
console.log(`small/phase2_only : ${fmt(t1)}`);

const t2 = timeit(() => solve(lpSmall));
console.log(`small/two_phase   : ${fmt(t2)}`);

console.log('');

// Random LP benchmarks
for (const [m, n] of [[10, 20], [50, 100], [100, 200]] as [number, number][]) {
  const lp = randomLP(m, n);
  const tag = `m=${m}_n=${n}`;
  const t = timeit(() => solve(lp));
  console.log(`random/${tag.padEnd(12)} : ${fmt(t)}`);
}
