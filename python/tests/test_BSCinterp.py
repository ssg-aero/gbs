import unittest
import sys

# TODO: check if better using install path
sys.path.insert(1, 'build/')
import pygbs

tol = 1e-6

class TestBSCinterp(unittest.TestCase):
    def distance(self,v1,v2):
        return sum( (x-y)**2 for x , y in zip(v1,v2))
    def test_cn(self):

        pts = [
            [0.,0.,0],
            [0.,0.,1],
            [1.,0.,0.5],
            [1.,1.,1]
        ]
        constrains = []
        for p in pts:
            constrains.append([p])

        crv = pygbs.interpolate_cn_3d_d(
            constrains,
            2,
            pygbs.KnotsCalcMode.CHORD_LENGTH
        )

        self.assertLess( self.distance(crv.begin(),pts[0]),tol)
        self.assertLess( self.distance(crv.end(),pts[-1]),tol)

        for p in pts:
            result = pygbs.extrema_PC(crv,p,tol)
            self.assertLess(result.d,tol)
            self.assertLess(self.distance(crv.value(result.u),p),tol)
    