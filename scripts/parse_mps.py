"""
Parse a fixed-format MPS file into standard equality form:
    maximize  c'x
    subject to  Ax = b,  x >= 0

Netlib problems use minimization; we negate c so the solver maximises and
returns z = -netlib_optimal.

Returns a dict with keys: name, m, n, A, b, c, netlib_optimal, optimal_z.
"""

from __future__ import annotations
import re
from pathlib import Path


def parse_mps(path: str | Path) -> dict:
    path = Path(path)
    text = path.read_text()

    section = None
    obj_name = None

    row_type: dict[str, str] = {}   # row_name -> 'N'|'L'|'G'|'E'
    col_data: dict[str, dict[str, float]] = {}  # col -> {row: val}
    rhs: dict[str, float] = {}
    bounds: dict[str, tuple[float, float]] = {}  # col -> (lo, hi)

    col_order: list[str] = []
    row_order: list[str] = []
    obj_row = None

    for raw in text.splitlines():
        line = raw.rstrip()
        if not line or line.startswith("$"):
            continue

        # Section headers start in column 0 (no leading space)
        if not line[0].isspace():
            word = line.split()[0].upper()
            if word in ("NAME", "ROWS", "COLUMNS", "RHS", "BOUNDS", "RANGES", "ENDATA"):
                section = word
                continue

        tokens = line.split()

        if section == "ROWS":
            rtype, rname = tokens[0], tokens[1]
            row_type[rname] = rtype
            if rtype == "N" and obj_row is None:
                obj_row = rname
            elif rtype != "N":
                row_order.append(rname)

        elif section == "COLUMNS":
            col = tokens[0]
            if col not in col_data:
                col_data[col] = {}
                col_order.append(col)
            # Up to two (row, value) pairs per line
            for i in range(1, len(tokens) - 1, 2):
                row, val = tokens[i], float(tokens[i + 1])
                col_data[col][row] = val

        elif section == "RHS":
            # Some emps-decoded files omit the RHS vector name; detect by
            # even token count (row val [row val]) vs odd (name row val ...)
            start = 0 if len(tokens) % 2 == 0 else 1
            for i in range(start, len(tokens) - 1, 2):
                row, val = tokens[i], float(tokens[i + 1])
                rhs[row] = val

        elif section == "BOUNDS":
            btype, _, col = tokens[0], tokens[1], tokens[2]
            val = float(tokens[3]) if len(tokens) > 3 else 0.0
            lo, hi = bounds.get(col, (0.0, float("inf")))
            if btype == "UP":
                bounds[col] = (lo, val)
            elif btype == "LO":
                bounds[col] = (val, hi)
            elif btype == "FX":
                bounds[col] = (val, val)
            elif btype in ("FR", "MI"):
                bounds[col] = (float("-inf"), hi)
            elif btype == "BV":
                bounds[col] = (0.0, 1.0)

    # Original variables
    orig_cols = col_order
    n_orig = len(orig_cols)
    m_rows = len(row_order)

    # Objective coefficients for original variables (negated for maximization)
    c_orig = [-col_data[col].get(obj_row, 0.0) for col in orig_cols]

    # Build constraint matrix rows and collect slack/surplus info
    # For each structural row:
    #   L (<=): add slack  s >= 0  →  a'x + s = b
    #   G (>=): add surplus s >= 0 →  a'x - s = b
    #   E (=):  no slack
    slack_cols: list[tuple[int, float]] = []  # (row_index, sign) sign=+1 for L, -1 for G

    for i, rname in enumerate(row_order):
        rt = row_type[rname]
        if rt == "L":
            slack_cols.append((i, +1.0))
        elif rt == "G":
            slack_cols.append((i, -1.0))

    n_slack = len(slack_cols)
    n_total = n_orig + n_slack

    col_index = {col: j for j, col in enumerate(orig_cols)}

    A = [[0.0] * n_total for _ in range(m_rows)]
    b = [0.0] * m_rows

    for i, rname in enumerate(row_order):
        b[i] = rhs.get(rname, 0.0)
        for col, row_vals in col_data.items():
            if rname in row_vals:
                j = col_index[col]
                A[i][j] = row_vals[rname]

    for slack_idx, (row_i, sign) in enumerate(slack_cols):
        A[row_i][n_orig + slack_idx] = sign

    # Objective for slack variables = 0 (already)
    c = c_orig + [0.0] * n_slack

    return {
        "name": path.stem.split(".")[0],
        "m": m_rows,
        "n": n_total,
        "A": A,
        "b": b,
        "c": c,
        # netlib_optimal is the minimisation value; optimal_z is what our maximiser returns
        "netlib_optimal": None,   # filled in by preparse.py
        "optimal_z": None,
    }
