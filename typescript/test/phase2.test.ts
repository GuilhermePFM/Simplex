import { describe, it, expect } from 'vitest';
import { solvePhase2 } from '../src/index.js';

// Note: all indices are 0-based (Julia uses 1-based, so Julia [3,4],[1,2] -> TS [2,3],[0,1])

describe('Phase 2 — solvePhase2', () => {
  it('optimal: production problem', () => {
    const lp = {
      A: [[2,1,1,0],[1,2,0,1]],
      b: [4,4],
      c: [4,3,0,0],
    };
    const basis = { basic: [2,3], nonbasic: [0,1] };

    const r = solvePhase2(lp, basis);

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeCloseTo(28/3, 6);
    expect(r.x[0]).toBeCloseTo(4/3, 6);
    expect(r.x[1]).toBeCloseTo(4/3, 6);
    expect(r.x[2]).toBeCloseTo(0.0, 6);
    expect(r.x[3]).toBeCloseTo(0.0, 6);
    expect(r.iterations).toBe(3);
  });

  it('optimal: already at vertex', () => {
    const lp = {
      A: [[1,0,1,0],[0,1,0,1]],
      b: [1,1],
      c: [3,2,0,0],
    };
    // Julia basis [1,2],[3,4] -> 0-based [0,1],[2,3]
    const basis = { basic: [0,1], nonbasic: [2,3] };

    const r = solvePhase2(lp, basis);

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeCloseTo(5.0, 6);
    expect(r.x[0]).toBeCloseTo(1.0, 6);
    expect(r.x[1]).toBeCloseTo(1.0, 6);
  });

  it('unbounded', () => {
    const lp = {
      A: [[0.5,-1,1,0],[-4,1,0,1]],
      b: [0.5,1],
      c: [1,1,0,0],
    };
    const basis = { basic: [2,3], nonbasic: [0,1] };

    const r = solvePhase2(lp, basis);

    expect(r.status).toBe('UNBOUNDED');
    expect(r.z).toBe(Infinity);
  });

  it('Bland rule gives same optimal', () => {
    const lp = {
      A: [[2,1,1,0],[1,2,0,1]],
      b: [4,4],
      c: [4,3,0,0],
    };
    const basis = { basic: [2,3], nonbasic: [0,1] };

    const r = solvePhase2(lp, basis, { rule: 'bland' });

    expect(r.status).toBe('OPTIMAL');
    expect(r.z).toBeCloseTo(28/3, 6);
  });
});
