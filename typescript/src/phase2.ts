import type { LinearProgram, BasisState, SimplexResult } from './types.js';
import type { PivotRule } from './pivot.js';
import { enteringIndex, leavingIndex } from './pivot.js';
import { luFactor, luSolve, luSolveMatrix, dot } from './linalg.js';
import { getColumns } from './linalg.js';
import type { SimplexLogger } from './logger.js';

export function phase2(
  lp: LinearProgram,
  initBasis: BasisState,
  rule: PivotRule,
  logger: SimplexLogger,
  maxit: number = 1000,
  tol: number = 1e-10,
): SimplexResult {
  const { A, b, c } = lp;
  const n = A[0]?.length ?? 0;

  const basic = initBasis.basic.slice();
  const nonbasic = initBasis.nonbasic.slice();

  logger.logPhase(2);

  for (let it = 1; it <= maxit; it++) {
    const B = getColumns(A, basic);
    const N = getColumns(A, nonbasic);
    const BF = luFactor(B);
    const xb = luSolve(BF, b);
    const BinvN = luSolveMatrix(BF, N);

    const cb = basic.map(j => c[j]);
    // cr[j] = c[nonbasic[j]] - BinvN[:,j]' * cb
    const cr = nonbasic.map((j, col) => {
      let s = c[j];
      for (let row = 0; row < cb.length; row++) s -= BinvN[row][col] * cb[row];
      return s;
    });
    const z = dot(cb, xb);

    logger.logIteration(it, { basic: basic.slice(), nonbasic: nonbasic.slice() }, xb, z);

    const j = enteringIndex(rule, cr, tol);
    if (j === null) {
      const x = new Array(n).fill(0);
      for (let i = 0; i < basic.length; i++) x[basic[i]] = xb[i];
      return { x, z, status: 'OPTIMAL', iterations: it };
    }

    const d = BinvN.map(row => row[j]);
    const i = leavingIndex(xb, d, BinvN, tol);

    if (i === null) {
      const ray = new Array(n).fill(0);
      ray[nonbasic[j]] = 1.0;
      for (let r = 0; r < basic.length; r++) ray[basic[r]] = -d[r];
      return { x: ray, z: Infinity, status: 'UNBOUNDED', iterations: it };
    }

    const theta = xb[i] / d[i];
    const x = new Array(n).fill(0);
    for (let r = 0; r < basic.length; r++) x[basic[r]] = xb[r] - theta * d[r];
    x[nonbasic[j]] = theta;
    x[basic[i]] = 0; // numerical cleanup

    // swap entering and leaving
    const tmp = basic[i];
    basic[i] = nonbasic[j];
    nonbasic[j] = tmp;
  }

  throw new Error(`Phase 2 did not converge in ${maxit} iterations — possible cycling`);
}
