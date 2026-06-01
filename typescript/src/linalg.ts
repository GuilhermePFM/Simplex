export interface LUResult {
  L: number[][];
  U: number[][];
  piv: number[]; // row permutation applied during factorization
}

/** LU factorization with partial pivoting. Modifies a copy of A. */
export function luFactor(A: number[][]): LUResult {
  const n = A.length;
  // deep copy
  const M: number[][] = A.map(row => row.slice());
  const piv: number[] = Array.from({ length: n }, (_, i) => i);

  for (let k = 0; k < n; k++) {
    // find pivot row
    let maxVal = Math.abs(M[k][k]);
    let maxRow = k;
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(M[i][k]);
      if (v > maxVal) { maxVal = v; maxRow = i; }
    }
    // swap rows in M and piv
    if (maxRow !== k) {
      [M[k], M[maxRow]] = [M[maxRow], M[k]];
      [piv[k], piv[maxRow]] = [piv[maxRow], piv[k]];
    }
    if (Math.abs(M[k][k]) < 1e-14) continue; // singular column, skip
    for (let i = k + 1; i < n; i++) {
      M[i][k] /= M[k][k];
      for (let j = k + 1; j < n; j++) {
        M[i][j] -= M[i][k] * M[k][j];
      }
    }
  }

  // Extract L and U from packed M
  const L: number[][] = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i === j ? 1 : i > j ? M[i][j] : 0))
  );
  const U: number[][] = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i <= j ? M[i][j] : 0))
  );

  return { L, U, piv };
}

/** Solve LU * x = b where piv is the row permutation applied to A. */
export function luSolve(lu: LUResult, b: number[]): number[] {
  const { L, U, piv } = lu;
  const n = L.length;

  // Apply row permutation: pb[i] = b[piv[i]]
  const pb: number[] = piv.map(p => b[p]);

  // Forward substitution: L * y = pb
  const y: number[] = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let s = pb[i];
    for (let j = 0; j < i; j++) s -= L[i][j] * y[j];
    y[i] = s; // L[i][i] == 1
  }

  // Back substitution: U * x = y
  const x: number[] = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = y[i];
    for (let j = i + 1; j < n; j++) s -= U[i][j] * x[j];
    x[i] = Math.abs(U[i][i]) < 1e-14 ? 0 : s / U[i][i];
  }

  return x;
}

/** Solve LU * X = B for each column of B; returns matrix X (m x k). */
export function luSolveMatrix(lu: LUResult, B: number[][]): number[][] {
  const m = B.length;
  if (m === 0) return [];
  const k = B[0].length;
  // Build columns of B, solve each, then transpose back
  const result: number[][] = Array.from({ length: m }, () => new Array(k).fill(0));
  for (let col = 0; col < k; col++) {
    const bCol: number[] = B.map(row => row[col]);
    const xCol = luSolve(lu, bCol);
    for (let row = 0; row < m; row++) result[row][col] = xCol[row];
  }
  return result;
}

/** Matrix-vector multiply: A (m x n) * v (n) -> result (m). */
export function matvec(A: number[][], v: number[]): number[] {
  return A.map(row => dot(row, v));
}

/** Dot product of two vectors. */
export function dot(a: number[], b: number[]): number {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

/** Extract column col from matrix A. */
export function getColumn(A: number[][], col: number): number[] {
  return A.map(row => row[col]);
}

/** Extract submatrix: columns given by indices array. */
export function getColumns(A: number[][], indices: number[]): number[][] {
  const m = A.length;
  return A.map(row => indices.map(j => row[j]));
}
