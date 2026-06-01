import type { LinearProgram, BasisState, SimplexResult } from './types.js';
import { SimplexLogger } from './logger.js';
import type { LogLevel } from './logger.js';
import type { PivotRule } from './pivot.js';
import { phase1 } from './phase1.js';
import { phase2 } from './phase2.js';

export type { SolveStatus } from './types.js';
export type { LinearProgram, BasisState, SimplexResult } from './types.js';
export { makeLinearProgram } from './types.js';
export type { PivotRule } from './pivot.js';
export type { LogLevel } from './logger.js';
export { SimplexLogger } from './logger.js';

export interface SolveOptions {
  rule?: PivotRule;
  verbosity?: LogLevel;
  maxit?: number;
  tol?: number;
}

/**
 * solve — two-phase Revised Simplex solve.
 * Runs Phase 1 to find a BFS, then Phase 2 to optimize.
 */
export function solve(lp: LinearProgram, opts: SolveOptions = {}): SimplexResult {
  const rule: PivotRule = opts.rule ?? 'largest';
  const verbosity: LogLevel = opts.verbosity ?? 'SILENT';
  const logger = new SimplexLogger(verbosity);

  logger.logProblem(lp);

  const maxit = opts.maxit ?? 1000;
  const [basis, p1Status] = phase1(lp, rule, logger, maxit, opts.tol ?? 1e-8);

  if (p1Status === 'INFEASIBLE') {
    const n = lp.A[0]?.length ?? 0;
    return { x: new Array(n).fill(0), z: 0, status: 'INFEASIBLE', iterations: 0 };
  }

  const result = phase2(lp, basis, rule, logger, maxit, opts.tol ?? 1e-10);
  logger.logResult(result);
  return result;
}

/**
 * solvePhase2 — run Phase 2 only, starting from a known BFS.
 */
export function solvePhase2(
  lp: LinearProgram,
  initBasis: BasisState,
  opts: SolveOptions = {},
): SimplexResult {
  const rule: PivotRule = opts.rule ?? 'largest';
  const verbosity: LogLevel = opts.verbosity ?? 'SILENT';
  const logger = new SimplexLogger(verbosity);

  const result = phase2(lp, initBasis, rule, logger, opts.maxit ?? 1000, opts.tol ?? 1e-10);
  logger.logResult(result);
  return result;
}

