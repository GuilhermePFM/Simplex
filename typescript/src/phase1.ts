import type { LinearProgram, BasisState } from './types.js';
import type { SolveStatus } from './types.js';
import type { PivotRule } from './pivot.js';
import { enteringIndex, leavingIndex } from './pivot.js';
import { luFactor, luSolve, luSolveMatrix, dot } from './linalg.js';
import { getColumns } from './linalg.js';
import { makeBNonnegative } from './preprocess.js';
import type { SimplexLogger } from './logger.js';

export function phase1(
  lpIn: LinearProgram,
  rule: PivotRule,
  logger: SimplexLogger,
  maxit: number = 1000,
  tol: number = 1e-8,
): [BasisState, SolveStatus] {
  const lp = makeBNonnegative(lpIn);
  const m = lp.A.length;
  const n = lp.A[0]?.length ?? 0;

  // Augmented matrix Aw = [A | I_m], 0-based column indices
  // Artificial variables occupy columns n..n+m-1
  const artStart = n; // 0-based start of artificial columns
  const Aw: number[][] = lp.A.map((row, i) => {
    const ext = new Array(m).fill(0);
    ext[i] = 1;
    return [...row, ...ext];
  });

  // Objective for phase 1: minimize sum of artificials -> maximize cw
  // cw[0..n-1] = 0, cw[n..n+m-1] = -1
  const cw: number[] = [...new Array(n).fill(0), ...new Array(m).fill(-1)];

  // Initial basis: artificials are basic, original vars are nonbasic
  const basic: number[] = Array.from({ length: m }, (_, i) => artStart + i);
  const nonbasic: number[] = Array.from({ length: n }, (_, j) => j);

  logger.logPhase(1);

  for (let it = 1; it <= maxit; it++) {
    const B = getColumns(Aw, basic);
    const N = getColumns(Aw, nonbasic);
    const BF = luFactor(B);
    const xb = luSolve(BF, lp.b);
    const BinvN = luSolveMatrix(BF, N);

    const cb = basic.map(j => cw[j]);
    const cr = nonbasic.map((j, col) => {
      let s = cw[j];
      for (let row = 0; row < cb.length; row++) s -= BinvN[row][col] * cb[row];
      return s;
    });
    const z = dot(cb, xb);

    logger.logIteration(it, { basic: basic.slice(), nonbasic: nonbasic.slice() }, xb, z);

    const j = enteringIndex(rule, cr, tol);
    if (j === null) {
      // Phase 1 terminated; check feasibility
      if (z < -tol) return [{ basic: [], nonbasic: [] }, 'INFEASIBLE'];

      fixArtificialsInBasis(Aw, basic, nonbasic, artStart);

      const origBasic = basic.filter(i => i < artStart);
      const origNonbasic = nonbasic.filter(i => i < artStart);
      return [{ basic: origBasic, nonbasic: origNonbasic }, 'OPTIMAL'];
    }

    const d = BinvN.map(row => row[j]);
    const i = leavingIndex(xb, d, BinvN, tol);

    // Phase 1 auxiliary problem is always bounded
    if (i === null) throw new Error('Phase 1 is unexpectedly unbounded');

    const tmp = basic[i];
    basic[i] = nonbasic[j];
    nonbasic[j] = tmp;
  }

  throw new Error(`Phase 1 did not converge in ${maxit} iterations`);
}

function fixArtificialsInBasis(
  Aw: number[][],
  basic: number[],
  nonbasic: number[],
  artStart: number,
): void {
  let changed = true;
  while (changed) {
    changed = false;
    for (let pos = 0; pos < basic.length; pos++) {
      if (basic[pos] < artStart) continue;

      const B = getColumns(Aw, basic);
      const BF = luFactor(B);
      // BinvAw = B^{-1} * Aw  (solve for each column of Aw)
      const BinvAw = luSolveMatrix(BF, Aw);
      const row = BinvAw[pos];

      let swapJ: number | null = null;
      for (let j = 0; j < nonbasic.length; j++) {
        const ncol = nonbasic[j];
        if (ncol < artStart && Math.abs(row[ncol]) > 1e-10) {
          swapJ = j;
          break;
        }
      }

      if (swapJ === null) {
        // Row is linearly dependent; zero it out to remove degenerate artificial
        for (let col = 0; col < Aw[pos].length; col++) Aw[pos][col] = 0;
      } else {
        const tmp = basic[pos];
        basic[pos] = nonbasic[swapJ];
        nonbasic[swapJ] = tmp;
        changed = true;
        break;
      }
    }
  }
}
