import sys
import unittest

# TODO: check if better using install path
sys.path.insert(1, 'build/python/Release')
import gbs

tol = 1e-6

class TestBSCurve(unittest.TestCase):
    def test_ctor(self):

        crv = gbs.BSCurve3d_d(
            [[0.,0.,0.],[1.,0.,0.]],
            [0.,0.,1.,1.],
            1)

        x = crv.value(0.5)

        self.assertAlmostEqual(x[0],0.5,tol)
        self.assertAlmostEqual(x[1],0.,tol)
        self.assertAlmostEqual(x[2],0.,tol)

        crv = gbs.BSCurve2d_d(
            [[0.,0.],[1.,0.]],
            [0.,0.,1.,1.],
            1)

        x = crv.value(0.5)

        self.assertAlmostEqual(x[0],0.5,tol)
        self.assertAlmostEqual(x[1],0.,tol)

        crv = gbs.BSCurve1d_d(
            [[0.],[1.]],
            [0.,0.,1.,1.],
            1)

        x = crv.value(0.5)

        self.assertAlmostEqual(x[0],0.5,tol)
