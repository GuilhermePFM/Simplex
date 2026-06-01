export type SolveStatus = 'OPTIMAL' | 'UNBOUNDED' | 'INFEASIBLE';

export interface LinearProgram {
  A: number[][];
  b: number[];
  c: number[];
}

export interface BasisState {
  basic: number[];
  nonbasic: number[];
}

export interface SimplexResult {
  x: number[];
  z: number;
  status: SolveStatus;
  iterations: number;
}

export function makeLinearProgram(A: number[][], b: number[], c: number[]): LinearProgram {
  const m = A.length;
  const n = A[0]?.length ?? 0;
  if (b.length !== m) throw new Error(`b has ${b.length} rows; A has ${m}`);
  if (c.length !== n) throw new Error(`c has ${c.length} entries; A has ${n} columns`);
  return { A, b, c };
}
