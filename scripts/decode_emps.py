"""
Decode Netlib LP data files (compressed "emps" format) to standard MPS text.

Faithful port of David M. Gay's emps.c (AT&T Bell Laboratories, 1987).
"""
from __future__ import annotations
import sys
from pathlib import Path

# 92 printable chars used by emps (no colon, no backslash)
_TRTAB = (
    '!"#$%&\'()*+,-./0123456789'   # ASCII 33-57  (skip ':' 58)
    ';<=>?@'                         # ASCII 59-64
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ'     # ASCII 65-90  (skip '\\' 92)
    '[]^_`'                          # ASCII 91, 93-96
    'abcdefghijklmnopqrstuvwxyz'     # ASCII 97-122
    '{|}~'                           # ASCII 123-126
)
assert len(_TRTAB) == 92

_INVTR = [92] * 256
for _i, _c in enumerate(_TRTAB):
    _INVTR[ord(_c)] = _i


def _tr(c: str) -> int:
    return _INVTR[ord(c)]


class _Stream:
    """Consume characters line by line, auto-advancing when a line is exhausted."""

    def __init__(self, lines: list[str]):
        self._it = iter(lines)
        self._line: str = ''
        self._pos: int = 0

    def _read_line(self) -> bool:
        """Read the next non-checksum, non-blank data line.

        Netlib emps files insert a checksum line (starts with ' ') after every
        72 data lines for integrity checking.  We skip those transparently.
        """
        while True:
            try:
                line = next(self._it).rstrip('\r\n')
            except StopIteration:
                return False
            # Skip checksum lines (start with space) and blank lines
            if line and line[0] != ' ':
                self._line = line
                self._pos = 0
                return True

    def _ensure(self):
        while self._pos >= len(self._line):
            if not self._read_line():
                raise EOFError('Premature end of emps data')

    def get(self) -> str:
        self._ensure()
        c = self._line[self._pos]
        self._pos += 1
        return c

    def peek(self) -> str:
        self._ensure()
        return self._line[self._pos]

    def consume_line(self):
        """Move position to end-of-line so next get() triggers a line read."""
        self._pos = len(self._line)

    def take_rest(self) -> str:
        """Return remaining chars on current line and mark line as consumed."""
        rest = self._line[self._pos:]
        self._pos = len(self._line)
        return rest

    def next_line(self):
        """Force read of next line (used for row names / column names)."""
        self._read_line()


def _exindx(s: _Stream) -> int:
    """Decode a super-sparse index (0-based).  0 signals a column separator."""
    k = _tr(s.get())
    if k >= 46:
        raise ValueError(f'exindx: bad char k={k}')
    if k >= 23:
        return k - 23          # single-char: values 0..22
    x = k
    while True:
        k = _tr(s.get())
        x = x * 46 + k
        if k >= 46:
            x -= 46
            break
    return x


def _exform(s: _Stream, numtab: list[str]) -> str:
    """Decode a numtab back-reference or a floating-point coefficient."""
    k = _tr(s.peek())
    if k < 46:
        # Back-reference into numtab; exindx returns 1-based index
        # (C code: ss = ss0-16; return ss + (k<<4) maps k=1..ns → ss0[0..ns-1])
        idx = _exindx(s)
        return numtab[idx - 1]

    s.get()          # consume dispatch char
    k -= 46
    negative = k >= 23
    if negative:
        k -= 23
    nelim = 11 if negative else 12

    if k >= 11:
        # Integer float: stored as "NNN."
        k -= 11
        if k >= 6:
            x = k - 6
        else:
            x = k
            while True:
                nk = _tr(s.get())
                x = x * 46 + nk
                if nk >= 46:
                    x -= 46
                    break
        mantissa = str(x) + '.'
    else:
        # General float
        ex = _tr(s.get()) - 50
        x = _tr(s.get())
        y = 0
        for _ in range(k):
            if x >= 100_000_000:
                y = x
                x = _tr(s.get())
            else:
                x = x * 92 + _tr(s.get())

        # Build decimal digit string (low-order first, then reverse)
        digits: list[str] = []
        if y:
            tmp = x
            while tmp > 1:
                digits.append(str(tmp % 10))
                tmp //= 10
            tmp = y
            while True:
                digits.append(str(tmp % 10))
                if tmp < 10:
                    break
                tmp //= 10
        elif x:
            tmp = x
            while True:
                digits.append(str(tmp % 10))
                if tmp < 10:
                    break
                tmp //= 10
        else:
            digits = ['0']
        ds = ''.join(reversed(digits))   # most-significant first
        nd = len(ds) + ex                # decimal point position from left edge

        if ex > 0:
            if nd < nelim or ex < 3:
                mantissa = ds + '0' * ex + '.'
            else:
                mantissa = _sci(ds, ex)
        elif nd >= 0:
            mantissa = (ds[:nd] or '0') + '.' + ds[nd:]
        elif ex > -nelim:
            mantissa = '.' + '0' * (-nd) + ds
        else:
            mantissa = _sci(ds, ex)

    val = ('-' if negative else '') + mantissa
    return val


def _sci(ds: str, ex: int) -> str:
    exp2 = ex + len(ds) - 1
    body = ds[0] + '.' + (ds[1:] or '0')
    return f'{body}E{exp2:+d}'


# ---------------------------------------------------------------------------
# Main decoder
# ---------------------------------------------------------------------------

def decode_emps(path: str | Path) -> str:
    raw = Path(path).read_text()
    lines = raw.split('\n')
    it = iter(lines)

    # Skip leading blank lines, find NAME
    hdr = ''
    for ln in it:
        if ln.strip():
            hdr = ln.rstrip()
            break
    if 'NAME' not in hdr:
        raise ValueError(f'Expected NAME line, got: {hdr!r}')
    out: list[str] = [hdr]

    # Statistics line 1
    s1 = next(it).split()
    nrow, ncol, _colmx, nz, _nrhs, rhsnz, _nran, ranz = [int(x) for x in s1[:8]]
    # Statistics line 2
    s2 = next(it).split()
    nbd, bdnz, ns = int(s2[0]), int(s2[1]), int(s2[2])

    strm = _Stream(list(it))

    # --- Number table (ns formatted float strings, may reference earlier entries) ---
    numtab: list[str] = []
    for _ in range(ns):
        numtab.append(_exform(strm, numtab))

    # --- ROWS section ---
    # Rows are one-per-line (type char + name)
    # Row names stored 1-based: row_names[1..nrow]
    row_names: list[str] = ['']     # index-0 unused
    out.append('ROWS')
    for _ in range(nrow):
        strm.next_line()
        ln = strm._line
        strm.consume_line()          # mark consumed so next _ensure reads next line
        rtype = ln[0]
        rname = ln[1:].strip()
        out.append(f' {rtype}  {rname}')
        row_names.append(rname)

    # --- Column name registry ---
    # Filled as we encounter column separators in COLUMNS section
    col_names: list[str] = [''] * (ncol + 2)   # 1-based
    _col_counter = [0]

    def _register_col(name: str) -> str:
        _col_counter[0] += 1
        col_names[_col_counter[0]] = name
        return name

    # --- Generic section decoder ---
    def _section(head: str, count: int, what: int) -> list[str]:
        """Decode COLUMNS (what=1), RHS (what=2), RANGES (what=3), BOUNDS (what=4)."""
        if count == 0:
            if what <= 2:
                return [head]
            return []

        sec: list[str] = [head]
        curcol = ''
        k = 0           # 0 or 1: pending pair members
        rownm = ['', '']
        rval  = ['', '']

        for _ in range(count):
            # --- Column separators (index 0) precede each entry ---
            while True:
                n = _exindx(strm)
                if n != 0:
                    break
                # Flush pending singleton
                if k:
                    sec.append(f'    {curcol:<8s}  {rownm[0]:<8s}  {rval[0]}')
                    k = 0
                # Column name = rest of current line
                raw_name = strm.take_rest().strip()[:8]
                curcol = raw_name
                if what == 1:
                    _register_col(raw_name)
                # Advance to next line of data
                strm.next_line()

            # --- Decode the nonzero ---
            if what == 4:
                # BOUNDS: n = bound type 1-6
                bound_types = ['UP', 'LO', 'FX', 'FR', 'MI', 'PL']
                if n < 1 or n > 6:
                    raise ValueError(f'Bad bound type {n}')
                bt = bound_types[n - 1]
                col_ref = _exindx(strm)
                cname = col_names[col_ref] if col_ref <= ncol else ''
                if n <= 3:
                    row_ref = _exindx(strm)
                    rname = row_names[row_ref] if row_ref < len(row_names) else ''
                    val = _exform(strm, numtab).strip()
                    sec.append(f' {bt} BND       {cname:<8s}  {rname:<8s}  {val}')
                else:
                    row_ref = _exindx(strm)
                    rname = row_names[row_ref] if row_ref < len(row_names) else ''
                    sec.append(f' {bt} BND       {cname:<8s}  {rname:<8s}')
                continue

            # COLUMNS / RHS / RANGES
            rownm[k] = row_names[n] if n < len(row_names) else f'ROW{n}'
            rval[k]  = _exform(strm, numtab).strip()

            if k == 0:
                k = 1
                continue
            # Output pair
            sec.append(
                f'    {curcol:<8s}  {rownm[0]:<8s}  {rval[0]:<15s}'
                f'{rownm[1]:<8s}  {rval[1]}'
            )
            k = 0

        # Flush trailing singleton
        if k:
            sec.append(f'    {curcol:<8s}  {rownm[0]:<8s}  {rval[0]}')
        return sec

    out.extend(_section('COLUMNS', nz, 1))
    out.extend(_section('RHS',     rhsnz, 2))
    out.extend(_section('RANGES',  ranz, 3))
    out.extend(_section('BOUNDS',  bdnz, 4))
    out.append('ENDATA')
    return '\n'.join(out) + '\n'


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.stderr.write('Usage: decode_emps.py <file.mps>\n')
        sys.exit(1)
    print(decode_emps(sys.argv[1]), end='')
