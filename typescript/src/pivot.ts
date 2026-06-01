export type PivotRule = 'largest' | 'bland';

/**
 * Return the index into cr of the entering variable, or null if optimal.
 * Indices are 0-based (index into the nonbasic array).
 */
export function enteringIndex(
  rule: PivotRule,
  cr: number[],
  tol: number = 1e-10,
): number | null {
  if (rule === 'largest') {
    let best = -1;
    let bestVal = tol;
    for (let j = 0; j < cr.length; j++) {
      if (cr[j] > bestVal) { bestVal = cr[j]; best = j; }
    }
    return best === -1 ? null : best;
  } else {
    // Bland: first positive reduced cost
    for (let j = 0; j < cr.length; j++) {
      if (cr[j] > tol) return j;
    }
    return null;
  }
}

/**
 * Lexicographic ratio test — returns the row index (into basic) of the leaving
 * variable, or null if all d[i] <= tol (unbounded direction).
 *
 * @param xb    current basic variable values (length m)
 * @param d     entering column in basis coordinates B^{-1} A_j (length m)
 * @param BinvN full B^{-1} N matrix (m x |nonbasic|), used for tie-breaking
 */
export function leavingIndex(
  xb: number[],
  d: number[],
  BinvN: number[][],
  tol: number = 1e-10,
): number | null {
  const m = d.length;
  const candidates: number[] = [];
  for (let i = 0; i < m; i++) {
    if (d[i] > tol) candidates.push(i);
  }
  if (candidates.length === 0) return null;

  let best = candidates[0];
  for (let ci = 1; ci < candidates.length; ci++) {
    const i = candidates[ci];
    const ratioI = xb[i] / d[i];
    const ratioBest = xb[best] / d[best];

    if (ratioI < ratioBest - tol) {
      best = i;
    } else if (Math.abs(ratioI - ratioBest) <= tol) {
      // Lexicographic tie-break on scaled rows of BinvN
      const k = BinvN[0]?.length ?? 0;
      for (let col = 0; col < k; col++) {
        const rowI = BinvN[i][col] / d[i];
        const rowBest = BinvN[best][col] / d[best];
        if (rowI < rowBest - tol) { best = i; break; }
        if (rowI > rowBest + tol) { break; }
      }
    }
  }

  return best;
}
