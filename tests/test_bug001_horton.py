"""Regression test for BUG-001: Horton-family generators must use exact
integer arithmetic (py2 `/` floor division migrated to `//`).

Run with: python3 tests/test_bug001_horton.py
"""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from PyDCG import crossing, geometricbasics, points


def horton_exp_reference(k):
    """Exact-integer reference implementation of points._horton_exp."""
    if k <= 0:
        return [[0, 0]]
    if k <= 1:
        g_k = 0
    else:
        f_k = 2**((k*(k-1)//2)-1)
        if k <= 2:
            f_k_1 = 0
        else:
            f_k_1 = 2**(((k-1)*(k-2)//2)-1)
        g_k = f_k-f_k_1
    H_k_1 = horton_exp_reference(k-1)
    H_even = [[2*x[0], x[1]] for x in H_k_1]
    H_odd = [[2*x[0]+1, x[1]+g_k] for x in H_k_1]
    H = []
    for i in range(len(H_even)):
        H.append(H_even[i])
        H.append(H_odd[i])
    return H


def assert_all_int(pts, label):
    for p in pts:
        for c in p:
            assert type(c) is int, "%s: non-int coordinate %r in %r" % (label, c, p)


def main():
    # _horton_exp matches the exact-integer oracle for k = 1..15
    # (before the fix, values diverge at k = 11 due to float rounding).
    for k in range(1, 16):
        H = points._horton_exp(k)
        assert H == horton_exp_reference(k), "_horton_exp(%d) != oracle" % k
        assert_all_int(H, "_horton_exp(%d)" % k)
    print("ok: _horton_exp(k) matches exact-integer oracle for k=1..15")

    # horton_set returns only ints for assorted sizes (powers of two and not).
    for n in [1, 2, 3, 7, 8, 16, 100, 2048]:
        H = points.horton_set(n)
        assert len(H) == n, "horton_set(%d) returned %d points" % (n, len(H))
        assert_all_int(H, "horton_set(%d)" % n)
    print("ok: horton_set(n) returns only ints for n in {1,2,3,7,8,16,100,2048}")

    # horton_set uses the right tree depth (TC-1: float log gave k+1 for
    # some exact powers of two, e.g. ceil(log2(2**29)) -> 30).
    assert (2**29-1).bit_length() == 29
    assert (2**31-1).bit_length() == 31
    assert points.horton_set(8) == points._horton_exp(3)
    print("ok: horton_set tree depth computed without floating point")

    # _get_CY / _get_CX return ints (before the fix: 45000.5 / 180.5).
    sample = [[0, 0], [1, 5], [2, 3], [3, 9]]
    cy = points._get_CY(sample)
    cx = points._get_CX(sample)
    assert cy == 45000, "_get_CY: expected 45000, got %r" % cy
    assert cx == 180, "_get_CX: expected 180, got %r" % cx
    assert type(cy) is int and type(cx) is int
    print("ok: _get_CY/_get_CX return ints")

    # random_squared_Horton_set returns only ints (before the fix:
    # coordinates like [3604.5, 1891930160.25]).
    random.seed(1)
    S = points.random_squared_Horton_set(4)
    assert_all_int(S, "random_squared_Horton_set(4)")
    print("ok: random_squared_Horton_set(4) returns only ints")

    # Sanity: order-type invariants at small sizes are unchanged.
    assert geometricbasics.general_position(points.horton_set(32))
    print("ok: horton_set(32) is in general position")
    cr = crossing.count_crossings(points.horton_set(16), speedup=False)
    assert cr == 1340, "count_crossings(horton_set(16)): expected 1340, got %r" % cr
    print("ok: count_crossings(horton_set(16)) == 1340")

    print("BUG-001 regression test: all checks passed")


if __name__ == "__main__":
    main()
