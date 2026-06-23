"""Rigorous Arb certificate for the R_3 trigonometric inequality.

This verifies the inequality used in the proof that dis(R_3) <= 2*pi/3:

    F_t(theta, theta') <= cos(t - 2*pi/3)

for theta, theta' in [-pi/3, pi/3] and t in [2*pi/3, pi].
The negative-t half follows by swapping theta and theta'.

The script uses python-flint (Arb balls).  On ordinary boxes it proves either
that the numerator is nonpositive, or that the squared cleared inequality is
nonnegative.  Tiny boxes around the endpoint equalities are handled by a
Taylor estimate whose constants are also bounded with Arb.
"""

from __future__ import annotations

from heapq import heappop, heappush
from itertools import product

import sympy as sp
from flint import arb, ctx


ctx.prec = 80

PI = arb.pi()
DELTA = 2 * PI / 3
ETA = 0.002
M3_BOUND = 766.0
Q2_BOUND = 5.0


def interval(lo: arb, hi: arb) -> arb:
    return lo.union(hi)


def numerator_and_cleared_gap(box: tuple[arb, arb, arb, arb, arb, arb]) -> tuple[arb, arb]:
    th = interval(box[0], box[1])
    ph = interval(box[2], box[3])
    tt = interval(box[4], box[5])

    a = arb(3) + arb(6) * (2 * th).cos()
    b = arb(3) + arb(6) * (2 * ph).cos()
    u = th - ph + tt
    c = (tt - DELTA).cos()

    numerator = a * b * u.cos() + (3 * u).cos()
    cleared_gap = c * c * (1 + a * a) * (1 + b * b) - numerator * numerator
    return numerator, cleared_gap


def widths(box: tuple[arb, arb, arb, arb, arb, arb]) -> list[float]:
    return [float(box[1] - box[0]), float(box[3] - box[2]), float(box[5] - box[4])]


def split(box: tuple[arb, arb, arb, arb, arb, arb]) -> tuple[tuple[arb, ...], tuple[arb, ...]]:
    idx = max(range(3), key=lambda i: widths(box)[i])
    data = list(box)
    mid = (data[2 * idx] + data[2 * idx + 1]) / 2
    left = data.copy()
    right = data.copy()
    left[2 * idx + 1] = mid
    right[2 * idx] = mid
    return tuple(left), tuple(right)


def is_endpoint_box(box: tuple[arb, arb, arb, arb, arb, arb]) -> bool:
    th0, th1, ph0, ph1, _t0, t1 = [float(x) for x in box]
    if t1 > float(DELTA) + ETA:
        return False
    endpoints = [-float(PI) / 3, float(PI) / 3]
    th_near = any(th0 >= e - ETA and th1 <= e + ETA for e in endpoints)
    ph_near = any(ph0 >= e - ETA and ph1 <= e + ETA for e in endpoints)
    return th_near and ph_near


def certify_ordinary_boxes() -> tuple[int, int, int]:
    root = (-PI / 3, PI / 3, -PI / 3, PI / 3, DELTA, PI)
    heap: list[tuple[float, int, tuple[arb, ...]]] = [(0.0, 0, root)]
    certified_by_numerator = 0
    certified_by_gap = 0
    endpoint_boxes = 0

    while heap:
        _priority, depth, box = heappop(heap)
        if is_endpoint_box(box):
            endpoint_boxes += 1
            continue

        numerator, gap = numerator_and_cleared_gap(box)
        if float(numerator.upper()) <= 0:
            certified_by_numerator += 1
            continue
        if float(gap.lower()) >= 0:
            certified_by_gap += 1
            continue
        if depth > 80:
            raise RuntimeError(
                "subdivision did not certify box "
                f"depth={depth}, box={[float(x) for x in box]}, "
                f"numerator={numerator}, gap={gap}"
            )

        for child in split(box):
            if is_endpoint_box(child):
                endpoint_boxes += 1
                continue
            child_numerator, child_gap = numerator_and_cleared_gap(child)
            if float(child_numerator.upper()) <= 0:
                certified_by_numerator += 1
            elif float(child_gap.lower()) >= 0:
                certified_by_gap += 1
            else:
                priority = float(child_gap.lower()) - float(child_numerator.upper())
                heappush(heap, (priority, depth + 1, child))

    return certified_by_numerator, certified_by_gap, endpoint_boxes


def endpoint_expression(sign_theta: int, sign_phi: int):
    x, y, v = sp.symbols("x y v")
    theta = sign_theta * sp.pi / 3 + (x if sign_theta < 0 else -x)
    phi = sign_phi * sp.pi / 3 + (y if sign_phi < 0 else -y)
    t = 2 * sp.pi / 3 + v
    a = 3 + 6 * sp.cos(2 * theta)
    b = 3 + 6 * sp.cos(2 * phi)
    u = theta - phi + t
    c = sp.cos(v)
    numerator = a * b * sp.cos(u) + sp.cos(3 * u)
    return c**2 * (1 + a * a) * (1 + b * b) - numerator**2


def certify_endpoint_taylor() -> float:
    x, y, v = sp.symbols("x y v")
    mods = {
        "sin": lambda z: z.sin() if hasattr(z, "sin") else arb(z).sin(),
        "cos": lambda z: z.cos() if hasattr(z, "cos") else arb(z).cos(),
    }
    box_var = arb(0).union(arb(str(ETA)))
    max_third_derivative = 0.0

    for sign_theta, sign_phi in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
        gap = endpoint_expression(sign_theta, sign_phi)
        for alpha in product(range(4), repeat=3):
            if sum(alpha) != 3:
                continue
            derivative = sp.diff(gap, x, alpha[0], y, alpha[1], v, alpha[2])
            fn = sp.lambdify((x, y, v), derivative, modules=[mods, "math"])
            val = fn(box_var, box_var, box_var)
            max_third_derivative = max(
                max_third_derivative,
                abs(float(val.lower())),
                abs(float(val.upper())),
            )

    if max_third_derivative >= M3_BOUND:
        raise RuntimeError(f"third derivative bound failed: {max_third_derivative}")

    # Taylor remainder after the quadratic term:
    # |R_3| <= M/6*(x+y+v)^3 <= (3/2)*M*ETA*(x^2+y^2+v^2).
    remainder_coefficient = 1.5 * M3_BOUND * ETA
    if remainder_coefficient >= Q2_BOUND:
        raise RuntimeError(
            f"endpoint Taylor margin failed: {remainder_coefficient} >= {Q2_BOUND}"
        )
    return max_third_derivative


def main() -> None:
    n_cert, h_cert, local = certify_ordinary_boxes()
    max_m3 = certify_endpoint_taylor()
    print("R3 inequality certified.")
    print(f"ordinary boxes certified by numerator <= 0: {n_cert}")
    print(f"ordinary boxes certified by cleared square gap >= 0: {h_cert}")
    print(f"endpoint boxes covered by Taylor estimate: {local}")
    print(f"max third derivative bound on endpoint charts: {max_m3:.6f} < {M3_BOUND}")
    print(f"Taylor remainder coefficient: {1.5 * M3_BOUND * ETA:.6f} < {Q2_BOUND}")


if __name__ == "__main__":
    main()
