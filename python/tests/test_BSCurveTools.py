import pytest
import pygbs.gbs as gbs

def test_join():
    crv1 = gbs.BSCurve3d(
        [
            [0.,0.,0.],
            [0.,1.,0.],
            [1.,1.,0.],
            [1.,1.,1.],
            [1.,1.,2.],
            [3.,1.,1.],
            [0.,4.,1.],
        ],
        [0., 0., 0., 1, 2, 3, 4, 5., 5., 5.],
        2
    )

    crv2 = gbs.BSCurve3d(
        [
            [0.,4.,1.],
            [0.,0.,0.],
        ],
        [0., 0., 1., 1.],
        1
    )

    crv3 = gbs.join(crv1, crv2)

    assert gbs.dist(crv1(0), crv3(0)) == pytest.approx(0.)
    assert gbs.dist(crv1(5), crv3(5)) == pytest.approx(0.)
    assert gbs.dist(crv2(1), crv3(6)) == pytest.approx(0.)

def test_blossom():
    k = [0., 0., 0., 0., 0.3, 0.5, 0.5, 1.2, 2., 2., 2., 2.]
    poles = [[0., 0., 0.], [1., 2., .5], [2., -1., 1.], [3., 3., 0.], [4., 0., -1.], [5., 1., 2.], [6., -2., 0.], [7., 0., 1.]]
    crv = gbs.BSCurve3d(poles, k, 3)
    # the curve on the diagonal
    for s, u in [(3, 0.1), (4, 0.4), (6, 0.9), (7, 1.7)]:
        assert gbs.blossom(poles, k, 3, s, [u, u, u]) == pytest.approx(crv(u), abs=1e-13)
        assert gbs.blossom(crv, s, [u, u, u]) == pytest.approx(crv(u), abs=1e-13)
    # pole i at the knots after it, symmetric
    assert gbs.blossom(crv, 7, k[5:8]) == pytest.approx(poles[4], abs=1e-13)
    assert gbs.blossom(crv, 7, [0.1, 1.7, 2.6]) == pytest.approx(gbs.blossom(crv, 7, [2.6, 0.1, 1.7]), abs=1e-13)
    # a rational curve: homogeneous
    w = [1., .5, 2., 1., .7, 1.5, 1., 1.]
    rat = gbs.BSCurveRational3d([[c * wi for c in p] + [wi] for p, wi in zip(poles, w)], k, 3)
    h = gbs.blossom(rat, 6, [0.9, 0.9, 0.9])
    assert [c / h[3] for c in h[:3]] == pytest.approx(rat(0.9), abs=1e-13)
    with pytest.raises(ValueError):
        gbs.blossom(crv, 5, [1., 2., 3.])  # an empty span
