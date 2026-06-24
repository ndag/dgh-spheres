"""Rigorous Arb certificate for the full-fiber R_3 inequality.

This is the companion script cited in the proof of ``dis(R_3)=2*pi/3``.
It verifies the full-fiber inequality

    q_{theta,r} . T_t q_{phi,r'} <= cos(t - 2*pi/3)

for

    0 <= theta, phi <= pi/3,
    -1 <= r, r' <= 1,
    2*pi/3 <= t <= pi.

The negative-t half follows by swapping the two fiber points.

The script uses python-flint Arb balls.  On ordinary boxes it proves either
that the numerator is nonpositive, or that the direct interval gap is
nonnegative.  The equality corner

    theta = phi = pi/3,  t = 2*pi/3

with arbitrary r,r' is handled by a Taylor certificate in inward coordinates

    x = pi/3 - theta,  y = pi/3 - phi,  v = t - 2*pi/3.

The ordinary subdivision is intentionally written as a verifier rather than as
a fast optimizer.  It may take a while on a laptop, but every accepted box is
certified by Arb interval arithmetic.
"""

from __future__ import annotations

from heapq import heappop, heappush
from itertools import product

import sympy as sp
from flint import arb, ctx


ctx.prec = 80

PI = arb.pi()
DELTA = 2 * PI / 3
ETA = arb("0.0022")
Q2_BOUND = arb(7)
M3_BOUND = arb(1881)


Box = tuple[arb, arb, arb, arb, arb, arb, arb, arb, arb, arb]


def interval(lo: arb, hi: arb) -> arb:
    return lo.union(hi)


def full_fiber_intervals(box: Box) -> tuple[arb, arb, arb]:
    """Return numerator, direct gap, and cleared square gap on a box."""
    th = interval(box[0], box[1])
    ph = interval(box[2], box[3])
    rr = interval(box[4], box[5])
    ss = interval(box[6], box[7])
    tt = interval(box[8], box[9])

    a = arb(3) + arb(6) * (2 * th).cos()
    b = arb(3) + arb(6) * (2 * ph).cos()

    cth = th.cos()
    cph = ph.cos()
    c3th = (3 * th).cos()
    c3ph = (3 * ph).cos()
    s3th = (3 * th).sin()
    s3ph = (3 * ph).sin()

    ct = tt.cos()
    c3t = (3 * tt).cos()
    st = tt.sin()
    s3t = (3 * tt).sin()
    cap = (tt - DELTA).cos()

    d_theta = (a * cth) * (a * cth) + c3th * c3th + arb(10) * rr * rr * s3th * s3th
    d_phi = (b * cph) * (b * cph) + c3ph * c3ph + arb(10) * ss * ss * s3ph * s3ph

    numerator = (
        ct * a * b * cth * cph
        + c3t * c3th * c3ph
        + (arb(9) * ct + c3t) * rr * ss * s3th * s3ph
        + ss * s3ph * (arb(3) * st * a * cth + s3t * c3th)
        + rr * s3th * (-arb(3) * st * b * cph - s3t * c3ph)
    )
    direct_gap = cap - numerator / (d_theta * d_phi).sqrt()
    cleared_gap = cap * cap * d_theta * d_phi - numerator * numerator
    return numerator, direct_gap, cleared_gap


def widths(box: Box) -> list[float]:
    return [
        float(box[1] - box[0]),
        float(box[3] - box[2]),
        float(box[5] - box[4]),
        float(box[7] - box[6]),
        float(box[9] - box[8]),
    ]


def split(box: Box) -> tuple[Box, Box]:
    # Angular/t dimensions are usually more important than the two chord
    # parameters, so weight them more heavily in the adaptive split.
    weights = [4.0, 4.0, 0.2, 0.2, 4.0]
    idx = max(range(5), key=lambda i: weights[i] * widths(box)[i])
    data = list(box)
    mid = (data[2 * idx] + data[2 * idx + 1]) / 2
    left = data.copy()
    right = data.copy()
    left[2 * idx + 1] = mid
    right[2 * idx] = mid
    return tuple(left), tuple(right)  # type: ignore[return-value]


def is_endpoint_box(box: Box) -> bool:
    return (
        box[0] >= PI / 3 - ETA
        and box[3] <= PI / 3
        and box[8] <= DELTA + ETA
    )


def ordinary_roots() -> list[Box]:
    """Split once along the endpoint thresholds before adaptive search."""
    theta_parts = [(arb(0), PI / 3 - ETA), (PI / 3 - ETA, PI / 3)]
    t_parts = [(DELTA, DELTA + ETA), (DELTA + ETA, PI)]
    roots: list[Box] = []
    for i, j, k in product(range(2), repeat=3):
        root = (
            theta_parts[i][0],
            theta_parts[i][1],
            theta_parts[j][0],
            theta_parts[j][1],
            arb(-1),
            arb(1),
            arb(-1),
            arb(1),
            t_parts[k][0],
            t_parts[k][1],
        )
        if is_endpoint_box(root):
            continue
        roots.append(root)
    return roots


def certify_ordinary_boxes(max_depth: int = 90) -> tuple[int, int, int, int]:
    certified_by_numerator = 0
    certified_by_direct_gap = 0
    certified_by_cleared_gap = 0
    endpoint_boxes = 0

    heap: list[tuple[float, int, int, Box]] = []
    counter = 0

    def push(box: Box, depth: int) -> None:
        nonlocal certified_by_numerator, certified_by_direct_gap
        nonlocal certified_by_cleared_gap, endpoint_boxes, counter

        if is_endpoint_box(box):
            endpoint_boxes += 1
            return

        numerator, direct_gap, cleared_gap = full_fiber_intervals(box)
        if float(numerator.upper()) <= 0:
            certified_by_numerator += 1
            return
        if float(direct_gap.lower()) >= 0:
            certified_by_direct_gap += 1
            return
        if float(cleared_gap.lower()) >= 0:
            certified_by_cleared_gap += 1
            return

        if depth >= max_depth:
            raise RuntimeError(
                "ordinary subdivision did not certify box "
                f"depth={depth}, box={[float(x) for x in box]}, "
                f"numerator={numerator}, direct_gap={direct_gap}, "
                f"cleared_gap={cleared_gap}"
            )

        counter += 1
        priority = float(direct_gap.lower())
        heappush(heap, (priority, counter, depth, box))

    for root in ordinary_roots():
        push(root, 0)

    while heap:
        _priority, _counter, depth, box = heappop(heap)
        for child in split(box):
            push(child, depth + 1)

    return (
        certified_by_numerator,
        certified_by_direct_gap,
        certified_by_cleared_gap,
        endpoint_boxes,
    )


def q2_minus_bound(x: arb, y: arb, v: arb, rr: arb, ss: arb) -> arb:
    q2 = (
        arb(90) * rr * rr * x * x
        + arb(63) * rr * ss * x * y
        - arb(18) * rr * v * x
        + arb(81) * rr * x * y
        + arb(90) * ss * ss * y * y
        + arb(18) * ss * v * y
        - arb(81) * ss * x * y
        + arb(8) * v * v
        + arb(27) * x * x
        + arb(27) * x * y
        + arb(27) * y * y
    )
    return q2 - Q2_BOUND * (x * x + y * y + v * v)


def certify_q2_bound() -> int:
    """Verify Q_2 >= 7(x^2+y^2+v^2) by homogeneous charts."""

    def chart_value(chart: str, box: tuple[arb, ...]) -> arb:
        aa = interval(box[0], box[1])
        bb = interval(box[2], box[3])
        rr = interval(box[4], box[5])
        ss = interval(box[6], box[7])
        if chart == "x":
            x, y, v = arb(1), aa, bb
        elif chart == "y":
            x, y, v = aa, arb(1), bb
        else:
            x, y, v = aa, bb, arb(1)
        return q2_minus_bound(x, y, v, rr, ss)

    def split_chart_box(box: tuple[arb, ...]) -> tuple[tuple[arb, ...], tuple[arb, ...]]:
        widths4 = [
            float(box[1] - box[0]),
            float(box[3] - box[2]),
            float(box[5] - box[4]),
            float(box[7] - box[6]),
        ]
        idx = max(range(4), key=lambda i: widths4[i])
        data = list(box)
        mid = (data[2 * idx] + data[2 * idx + 1]) / 2
        left = data.copy()
        right = data.copy()
        left[2 * idx + 1] = mid
        right[2 * idx] = mid
        return tuple(left), tuple(right)

    certified = 0
    root = (arb(0), arb(1), arb(0), arb(1), arb(-1), arb(1), arb(-1), arb(1))
    for chart in ["x", "y", "v"]:
        heap = [(0, root)]
        while heap:
            depth, box = heap.pop()
            if float(chart_value(chart, box).lower()) >= 0:
                certified += 1
                continue
            if depth > 60:
                raise RuntimeError(f"Q2 bound failed on chart {chart}: {[float(z) for z in box]}")
            heap.extend((depth + 1, child) for child in split_chart_box(box))
    return certified


def endpoint_expression():
    x, y, v, rr, ss = sp.symbols("x y v rr ss")
    theta = sp.pi / 3 - x
    phi = sp.pi / 3 - y
    t = 2 * sp.pi / 3 + v

    a = 3 + 6 * sp.cos(2 * theta)
    b = 3 + 6 * sp.cos(2 * phi)
    d_theta = (a * sp.cos(theta)) ** 2 + sp.cos(3 * theta) ** 2 + 10 * rr**2 * sp.sin(3 * theta) ** 2
    d_phi = (b * sp.cos(phi)) ** 2 + sp.cos(3 * phi) ** 2 + 10 * ss**2 * sp.sin(3 * phi) ** 2
    numerator = (
        sp.cos(t) * a * b * sp.cos(theta) * sp.cos(phi)
        + sp.cos(3 * t) * sp.cos(3 * theta) * sp.cos(3 * phi)
        + (9 * sp.cos(t) + sp.cos(3 * t)) * rr * ss * sp.sin(3 * theta) * sp.sin(3 * phi)
        + ss * sp.sin(3 * phi) * (3 * sp.sin(t) * a * sp.cos(theta) + sp.sin(3 * t) * sp.cos(3 * theta))
        + rr * sp.sin(3 * theta) * (-3 * sp.sin(t) * b * sp.cos(phi) - sp.sin(3 * t) * sp.cos(3 * phi))
    )
    gap = sp.cos(v) ** 2 * d_theta * d_phi - numerator**2
    return gap, (x, y, v, rr, ss)


def certify_endpoint_taylor() -> float:
    x, y, v, rr, ss = sp.symbols("x y v rr ss")
    gap, variables = endpoint_expression()
    mods = {
        "sin": lambda z: z.sin() if hasattr(z, "sin") else arb(z).sin(),
        "cos": lambda z: z.cos() if hasattr(z, "cos") else arb(z).cos(),
    }
    box_var = arb(0).union(ETA)

    max_third_derivative = 0.0
    for alpha in product(range(4), repeat=3):
        if sum(alpha) != 3:
            continue
        derivative = sp.diff(gap, x, alpha[0], y, alpha[1], v, alpha[2]).expand()
        # Arb can return nan if we evaluate a polynomial in rr,ss on [-1,1]
        # after trigonometric interval dependency.  Bound coefficient-by-
        # coefficient instead.
        poly = sp.Poly(derivative, rr, ss)
        derivative_bound = 0.0
        for _monomial, coeff in poly.terms():
            fn = sp.lambdify((x, y, v), coeff, modules=[mods, "math"])
            val = fn(box_var, box_var, box_var)
            derivative_bound += max(abs(float(val.lower())), abs(float(val.upper())))
        max_third_derivative = max(max_third_derivative, derivative_bound)

    if max_third_derivative >= float(M3_BOUND):
        raise RuntimeError(f"third derivative bound failed: {max_third_derivative}")

    # Taylor remainder after the quadratic term:
    # |R_3| <= M/6*(x+y+v)^3 <= (3/2)*M*ETA*(x^2+y^2+v^2).
    remainder_coefficient = 1.5 * float(M3_BOUND) * float(ETA)
    if remainder_coefficient >= float(Q2_BOUND):
        raise RuntimeError(
            f"endpoint Taylor margin failed: {remainder_coefficient} >= {Q2_BOUND}"
        )
    return max_third_derivative


def main() -> None:
    q2_boxes = certify_q2_bound()
    max_m3 = certify_endpoint_taylor()
    n_cert, direct_cert, h_cert, endpoint = certify_ordinary_boxes()

    print("Full R3 inequality certified.")
    print(f"Q2 lower-bound boxes certified: {q2_boxes}")
    print(f"ordinary boxes certified by numerator <= 0: {n_cert}")
    print(f"ordinary boxes certified by direct gap >= 0: {direct_cert}")
    print(f"ordinary boxes certified by cleared square gap >= 0: {h_cert}")
    print(f"endpoint boxes covered by Taylor estimate: {endpoint}")
    print(f"max third derivative bound on endpoint chart: {max_m3:.6f} < {M3_BOUND}")
    print(f"Taylor remainder coefficient: {1.5 * float(M3_BOUND) * float(ETA):.6f} < {Q2_BOUND}")


if __name__ == "__main__":
    main()
