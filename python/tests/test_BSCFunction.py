import pygbs as gbs
import pytest
def test_interp():
    Q = [1.0,0.5,1.5]
    u = [0.0,0.5,1.0]
    p = min(len(Q)-1,5)
    f = gbs.interpolate_cn_function(
        Q = Q,
        u = u,
        p = p
    )
    for Q_,u_ in zip(Q,u):
        assert f(u_) == Q_    

    crv = f.basisCurve()
    tol = 1e-6
    u_,d_ = gbs.extrema_curve_point(crv,[0.5],tol)
    assert f(u_) == pytest.approx( 0.5 )
    assert d_ <= tol


