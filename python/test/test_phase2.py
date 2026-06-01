import numpy as np
import pytest
from pytest import approx

from simplex import (
    BasisState,
    Bland,
    LinearProgram,
    SolveStatus,
    solve_phase2,
)


def test_optimal_production_problem():
    lp    = LinearProgram(
        np.array([[2, 1, 1, 0], [1, 2, 0, 1]], dtype=float),
        np.array([4, 4], dtype=float),
        np.array([4, 3, 0, 0], dtype=float),
    )
    basis = BasisState(basic=[2, 3], nonbasic=[0, 1])

    r = solve_phase2(lp, basis)

    assert r.status == SolveStatus.OPTIMAL
    assert r.z    == approx(28 / 3, abs=1e-8)
    assert r.x[0] == approx(4 / 3, abs=1e-8)
    assert r.x[1] == approx(4 / 3, abs=1e-8)
    assert r.x[2] == approx(0.0,   abs=1e-8)
    assert r.x[3] == approx(0.0,   abs=1e-8)
    assert r.iterations == 3


def test_optimal_already_at_vertex():
    lp    = LinearProgram(
        np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=float),
        np.array([1, 1], dtype=float),
        np.array([3, 2, 0, 0], dtype=float),
    )
    basis = BasisState(basic=[0, 1], nonbasic=[2, 3])

    r = solve_phase2(lp, basis)

    assert r.status == SolveStatus.OPTIMAL
    assert r.z    == approx(5.0,  abs=1e-8)
    assert r.x[0] == approx(1.0,  abs=1e-8)
    assert r.x[1] == approx(1.0,  abs=1e-8)


def test_unbounded():
    lp    = LinearProgram(
        np.array([[0.5, -1, 1, 0], [-4, 1, 0, 1]], dtype=float),
        np.array([0.5, 1], dtype=float),
        np.array([1, 1, 0, 0], dtype=float),
    )
    basis = BasisState(basic=[2, 3], nonbasic=[0, 1])

    r = solve_phase2(lp, basis)

    assert r.status == SolveStatus.UNBOUNDED
    assert r.z == float("inf")


def test_bland_rule_gives_same_optimal():
    lp    = LinearProgram(
        np.array([[2, 1, 1, 0], [1, 2, 0, 1]], dtype=float),
        np.array([4, 4], dtype=float),
        np.array([4, 3, 0, 0], dtype=float),
    )
    basis = BasisState(basic=[2, 3], nonbasic=[0, 1])

    r = solve_phase2(lp, basis, rule=Bland())

    assert r.status == SolveStatus.OPTIMAL
    assert r.z == approx(28 / 3, abs=1e-8)
