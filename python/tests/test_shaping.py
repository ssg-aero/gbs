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


def test_move_surf_point():
    poles = [
        [0.0,0.0,0.0], [0.5,0.0,0.0],[1.0,0.0,0.0],
        [0.0,0.5,0.0], [0.5,0.5,0.0],[1.0,0.5,0.0],
        [0.0,1.0,0.0], [0.5,1.0,0.0],[1.0,1.0,0.0],
    ]
    ku = [0.,0.,0.,1.,1.,1.]
    kv = [0.,0.,0.,1.,1.,1.]

    srf = gbs.BSSurface3d(
        poles, ku, kv, 2, 2
    )

    srf1 = gbs.BSSurface3d(
        poles, ku, kv, 2, 2
    )

    gbs.move_to_point(
        srf1,[0.5,0.5,0.3], 0.5, 0.5
    )

    srf2 = gbs.moved_to_point(
        srf, 
        [0.5,0.5,0.3], 
        0.5, 
        0.5, 
        True
    )

    srf3 = gbs.moved_to_constraints_delta(
        srf, 
        [
            (0.0,0.0,[0.0,0.0,0.3],0,0),
            (0.0,0.7,[0.0,0.0,0.5],0,0)
        ],
        0,1, # clamp u max side
        0,2,
    )

    assert gbs.dist( srf1(0.5,0.5), [0.5,0.5,0.3] ) == pytest.approx(0.)

    assert gbs.dist( srf2(0.5,0.5), [0.5,0.5,0.3] ) == pytest.approx(0.)
    assert gbs.dist( srf2(0.0,0.0), srf(0.0,0.0)  ) == pytest.approx(0.)
    assert gbs.dist( srf2(0.0,1.0), srf(0.0,1.0)  ) == pytest.approx(0.)
    assert gbs.dist( srf2(1.0,1.0), srf(1.0,1.0)  ) == pytest.approx(0.)
    assert gbs.dist( srf2(1.0,0.0), srf(1.0,0.0)  ) == pytest.approx(0.)

    assert gbs.dist( srf3(0.0,0.0), [0.0,0.0,0.3] ) == pytest.approx(0.)
    assert gbs.dist( srf3(0.0,0.7), [0.0,0.7,0.5] ) == pytest.approx(0.)

    # gbs.plot_surfaces(
    #     [
    #         srf, 
    #         srf1, 
    #         srf2,
    #         srf3,
    #     ]
    # )