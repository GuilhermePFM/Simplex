import { describe, it, expect } from 'vitest';
import { solvePhase2 } from '../src/index.js';

describe('Anti-cycling — degenerate tie-breaking', () => {
  const A = [
    [1,1,1,0,0],
    [1,0,0,1,0],
    [0,1,0,0,1],
  ];
  const b = [2,1,1];
  const c = [1,1,0,0,0];
  const lp = { A, b, c };
  // Julia basis [3,4,5],[1,2] -> 0-based [2,3,4],[0,1]
  const basis = { basic: [2,3,4], nonbasic: [0,1] };

  it('LargestCoefficient — terminates at tied ratio', () => {
    const r = solvePhase2(lp, basis, { rule: 'largest' });

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeCloseTo(2.0, 6);
    expect(r.x[0]).toBeCloseTo(1.0, 6);
    expect(r.x[1]).toBeCloseTo(1.0, 6);
    expect(r.iterations).toBe(3);

    // Feasibility certificate
    for (let i = 0; i < A.length; i++) {
      const ax = A[i].reduce((s, aij, j) => s + aij * r.x[j], 0);
      expect(ax).toBeCloseTo(b[i], 6);
    }
    expect(r.x.every(v => v >= -1e-10)).toBe(true);
  });

  it('Bland — same result', () => {
    const r = solvePhase2(lp, basis, { rule: 'bland' });

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeCloseTo(2.0, 6);
    expect(r.x[0]).toBeCloseTo(1.0, 6);
    expect(r.x[1]).toBeCloseTo(1.0, 6);
  });
});

describe('Anti-cycling — degenerate start (Kotiah-Steinberg)', () => {
  const A = [
    [-0.5,  5.5,  2.5, -9, 1, 0, 0],
    [-0.5,  2.5,  0.5, -1, 0, 1, 0],
    [ 1.0,  0.0,  0.0,  0, 0, 0, 1],
  ];
  const b = [0, 0, 1];
  const c = [10, -57, -9, -24, 0, 0, 0];
  const lp = { A, b, c };
  // Julia basis [5,6,7],[1,2,3,4] -> 0-based [4,5,6],[0,1,2,3]
  const basis = { basic: [4,5,6], nonbasic: [0,1,2,3] };

  it('LargestCoefficient — reaches optimal z=10, x[0]=1', () => {
    const r = solvePhase2(lp, basis, { rule: 'largest' });

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeCloseTo(10.0, 6);
    expect(r.x[0]).toBeCloseTo(1.0, 6);

    for (let i = 0; i < A.length; i++) {
      const ax = A[i].reduce((s, aij, j) => s + aij * r.x[j], 0);
      expect(ax).toBeCloseTo(b[i], 6);
    }
  });
});
