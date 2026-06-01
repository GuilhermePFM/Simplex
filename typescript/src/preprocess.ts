import type { LinearProgram } from './types.js';

/** Return a copy of lp with all right-hand side values non-negative. */
export function makeBNonnegative(lp: LinearProgram): LinearProgram {
  const mask = lp.b.map(v => v < 0);
  if (!mask.some(Boolean)) return lp;

  const A2 = lp.A.map((row, i) =>
    mask[i] ? row.map(v => -v) : row.slice()
  );
  const b2 = lp.b.map((v, i) => (mask[i] ? -v : v));
  return { A: A2, b: b2, c: lp.c };
}
