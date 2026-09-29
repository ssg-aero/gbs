import pygbs.gbs as gbs
from math import sin, cos, pi


def distorted_grid(ni, nj):
    pts = []
    for i in range(ni):
        for j in range(nj):
            u, v = i / (ni - 1), j / (nj - 1)
            x, y = u + 0.1 * v * v, v + 0.2 * sin(pi * u)
            if 0 < i < ni - 1 and 0 < j < nj - 1:
                x += 0.3 / ni * sin(7. * i + 3. * j)
                y += 0.3 / nj * cos(5. * i - 2. * j)
            pts.append([x, y])
    return pts


def max_move(a, b):
    return max(abs(p[0] - q[0]) + abs(p[1] - q[1]) for p, q in zip(a, b))


def test_elliptic_structured_smoothing_honours_n_it_and_tol():
    ni, nj = 20, 15
    pts0 = distorted_grid(ni, nj)
    one = gbs.elliptic_structured_smoothing(distorted_grid(ni, nj), nj, 0, ni - 1, 0, nj - 1, n_it=1, tol=1e-10)
    two = gbs.elliptic_structured_smoothing(distorted_grid(ni, nj), nj, 0, ni - 1, 0, nj - 1, n_it=2, tol=1e-10)
    conv = gbs.elliptic_structured_smoothing(distorted_grid(ni, nj), nj, 0, ni - 1, 0, nj - 1, n_it=10000, tol=1e-10)
    # n_it is taken into account (it used to be ignored and tol truncated into it)
    assert max_move(one, two) > 0.
    # the fully converged grid is a fixed point of one more sweep, to ~tol
    again = gbs.elliptic_structured_smoothing(conv, nj, 0, ni - 1, 0, nj - 1, n_it=1, tol=0.)
    assert max_move(conv, again) < 1e-9
    # boundary nodes are untouched
    for i in (0, ni - 1):
        for j in range(nj):
            assert conv[j + nj * i] == pts0[j + nj * i]
