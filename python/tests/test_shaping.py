from pygbs import gbs
import pytest


def test_move_curve_point():


    crv = gbs.BSCurve2d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)

    crv.increaseDegree(4)
    crv.insertKnot(0.5,3)

    # gbs.plot_curves_2d([crv, gbs.moved_to_point(crv,[0.5,0.5],0.5)])
    assert gbs.dist(gbs.moved_to_point(crv,[0.5,0.5],0.5)(0.5),[0.5,0.5]) == pytest.approx(0.)
    gbs.move_to_point(crv,[0.5,0.5],0.5)
    assert gbs.dist(crv(0.5),[0.5,0.5]) == pytest.approx(0.)