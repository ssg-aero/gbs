import sys
import unittest

# TODO: check if better using install path
sys.path.insert(1, 'build/python/Release')
import gbs

tol = 1e-6

class TestBSCinterp(unittest.TestCase):
    def test_cn(self):

        crv = gbs.interpolate_cn_3d_d(
            [
                [[0.,0.,0]],
                [[0.,0.,1]],
                [[1.,0.,0.5]],
                [[1.,1.,1]],
            ],
            2,
            gbs.KnotsCalcMode.CHORD_LENGTH
        )
    