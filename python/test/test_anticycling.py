import numpy as np
import pytest
from pytest import approx

from simplex import (
    BasisState,
    Bland,
    LargestCoefficient,
    LinearProgram,
    SolveStatus,
    solve_phase2,
)


def test_degenerate_tie_breaking_largest_coefficient():
    A = np.array([[1, 1, 1, 0, 0],
                  [1, 0, 0, 1, 0],
                  [0, 1, 0, 0, 1]], dtype=float)
    b = np.array([2, 1, 1], dtype=float)
    c = np.array([1, 1, 0, 0, 0], dtype=float)

    lp    = LinearProgram(A, b, c)
    basis = BasisState(basic=[2, 3, 4], nonbasic=[0, 1])

    r = solve_phase2(lp, basis, rule=LargestCoefficient())

    assert r.status     == SolveStatus.OPTIMAL
    assert r.z          == approx(2.0, abs=1e-8)
    assert r.x[0]       == approx(1.0, abs=1e-8)
    assert r.x[1]       == approx(1.0, abs=1e-8)
    assert r.iterations == 3
    assert A @ r.x == approx(b, abs=1e-8)
    assert all(xi >= -1e-10 for xi in r.x)


def test_degenerate_tie_breaking_bland():
    A = np.array([[1, 1, 1, 0, 0],
                  [1, 0, 0, 1, 0],
                  [0, 1, 0, 0, 1]], dtype=float)
    b = np.array([2, 1, 1], dtype=float)
    c = np.array([1, 1, 0, 0, 0], dtype=float)

    lp    = LinearProgram(A, b, c)
    basis = BasisState(basic=[2, 3, 4], nonbasic=[0, 1])

    r = solve_phase2(lp, basis, rule=Bland())

    assert r.status == SolveStatus.OPTIMAL
    assert r.z      == approx(2.0, abs=1e-8)
    assert r.x[0]   == approx(1.0, abs=1e-8)
    assert r.x[1]   == approx(1.0, abs=1e-8)


def test_kotiah_steinberg():
    A = np.array([
        [-0.5,  5.5,  2.5, -9,  1,  0,  0],
        [-0.5,  2.5,  0.5, -1,  0,  1,  0],
        [ 1.0,  0.0,  0.0,  0,  0,  0,  1],
    ], dtype=float)
    b = np.array([0, 0, 1], dtype=float)
    c = np.array([10, -57, -9, -24, 0, 0, 0], dtype=float)

    lp    = LinearProgram(A, b, c)
    basis = BasisState(basic=[4, 5, 6], nonbasic=[0, 1, 2, 3])

    r = solve_phase2(lp, basis, rule=LargestCoefficient())

    assert r.status == SolveStatus.OPTIMAL
    assert r.z      == approx(10.0, abs=1e-8)
    assert r.x[0]   == approx(1.0,  abs=1e-8)
    assert A @ r.x  == approx(b,    abs=1e-8)
