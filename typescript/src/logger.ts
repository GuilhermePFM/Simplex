import type { LinearProgram, BasisState, SimplexResult } from './types.js';

export type LogLevel = 'SILENT' | 'INFO' | 'DEBUG';

export class SimplexLogger {
  readonly level: LogLevel;

  constructor(level: LogLevel = 'SILENT') {
    this.level = level;
  }

  logProblem(lp: LinearProgram): void {
    if (this.level === 'SILENT') return;
    this._write('=======================');
    this._write(' Revised Simplex Solver');
    this._write('=======================');
    this._write(`A =\n${lp.A.map(r => r.join(' ')).join('\n')}`);
    this._write(`b = ${lp.b}`);
    this._write(`c = ${lp.c}`);
    this._write('');
  }

  logPhase(phase: number): void {
    if (this.level === 'SILENT') return;
    const header = phase === 1
      ? 'Phase 1 — Finding initial BFS'
      : 'Phase 2 — Optimizing';
    this._write(header);
    this._write('─'.repeat(header.length));
  }

  logIteration(it: number, basis: BasisState, xb: number[], z: number): void {
    if (this.level !== 'DEBUG') return;
    this._write(`  iter ${it}:`);
    this._write(`    xb       = ${xb}`);
    this._write(`    basic    = ${basis.basic}`);
    this._write(`    nonbasic = ${basis.nonbasic}`);
    this._write(`    z        = ${z}`);
    this._write('');
  }

  logResult(result: SimplexResult): void {
    if (this.level === 'SILENT') return;
    this._write('');
    this._write(`Result: ${result.status}`);
    this._write(`  x          = ${result.x}`);
    this._write(`  z          = ${result.z}`);
    this._write(`  iterations = ${result.iterations}`);
  }

  private _write(s: string): void {
    console.log(s);
  }
}
