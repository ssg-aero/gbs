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