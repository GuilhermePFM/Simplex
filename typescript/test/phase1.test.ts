import { describe, it, expect } from 'vitest';
import { solve } from '../src/index.js';

describe('Phase 1 — full solve', () => {
  it('feasible: production problem with Phase 1', () => {
    const lp = {
      A: [[2,1,1,0],[1,2,0,1]],
      b: [4,4],
      c: [4,3,0,0],
    };

    const r = solve(lp);

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeCloseTo(28/3, 6);
    expect(r.x[0]).toBeCloseTo(4/3, 6);
    expect(r.x[1]).toBeCloseTo(4/3, 6);
  });

  it('feasible: problem with a negative-b row', () => {
    const lp = {
      A: [
        [ 2,  1, 1, 0, 0],
        [ 1,  2, 0, 1, 0],
        [-1, -1, 0, 0, 1],
      ],
      b: [4, 4, -1],
      c: [4, 3, 0, 0, 0],
    };

    const r = solve(lp);

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeGreaterThan(0);
    // Feasibility check: A * x ≈ b
    for (let i = 0; i < lp.A.length; i++) {
      const ax = lp.A[i].reduce((s, aij, j) => s + aij * r.x[j], 0);
      expect(ax).toBeCloseTo(lp.b[i], 5);
    }
    expect(r.x.every(v => v >= -1e-10)).toBe(true);
  });

  it('infeasible: contradictory constraints', () => {
    // x1 + x2 = 1  AND  x1 + x2 = 0 — impossible
    const lp = {
      A: [[1,1],[1,1]],
      b: [1,0],
      c: [1,1],
    };

    const r = solve(lp);

    expect(r.status).toBe('INFEASIBLE');
  });

  it('infeasible: conflicting bounds x1<=1 and x1>=2', () => {
    // x1 + s1 = 1   (x1 <= 1)
    // -x1 + s2 = -2  (x1 >= 2, negated by preprocess)
    const lp = {
      A: [[1,1,0],[-1,0,1]],
      b: [1,-2],
      c: [1,0,0],
    };

    const r = solve(lp);

    expect(r.status).toBe('INFEASIBLE');
  });
});
