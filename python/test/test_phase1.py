import numpy as np
import pytest
from pytest import approx

from simplex import LinearProgram, SolveStatus, solve


def test_feasible_production_problem():
    lp = LinearProgram(
        np.array([[2, 1, 1, 0], [1, 2, 0, 1]], dtype=float),
        np.array([4, 4], dtype=float),
        np.array([4, 3, 0, 0], dtype=float),
    )

    r = solve(lp)

    assert r.status == SolveStatus.OPTIMAL
    assert r.z    == approx(28 / 3, abs=1e-8)
    assert r.x[0] == approx(4 / 3, abs=1e-8)
    assert r.x[1] == approx(4 / 3, abs=1e-8)


def test_feasible_negative_b_row():
    lp = LinearProgram(
        np.array([[2, 1, 1, 0, 0],
                  [1, 2, 0, 1, 0],
                  [-1, -1, 0, 0, 1]], dtype=float),
        np.array([4, 4, -1], dtype=float),
        np.array([4, 3, 0, 0, 0], dtype=float),
    )

    r = solve(lp)

    assert r.status == SolveStatus.OPTIMAL
    assert r.z > 0
    assert lp.A @ r.x == approx(lp.b, abs=1e-6)
    assert all(xi >= -1e-10 for xi in r.x)


def test_infeasible_contradictory_constraints():
    lp = LinearProgram(
        np.array([[1, 1], [1, 1]], dtype=float),
        np.array([1, 0], dtype=float),
        np.array([1, 1], dtype=float),
    )

    r = solve(lp)

    assert r.status == SolveStatus.INFEASIBLE


def test_infeasible_conflicting_bounds():
    lp = LinearProgram(
        np.array([[1, 1, 0], [-1, 0, 1]], dtype=float),
        np.array([1, -2], dtype=float),
        np.array([1, 0, 0], dtype=float),
    )

    r = solve(lp)

    assert r.status == SolveStatus.INFEASIBLE
