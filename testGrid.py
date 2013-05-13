from flowTransport import FTSystem
import unittest
import numpy as np

class TestGrid(unittest.TestCase):
  """ Runs a handful of tests on a simple grid shown in the setUp method.
  """
  def setUp(self):
    xw = -75.
    xe = 50.
    ys = 0.
    yn = 100.
    nx = 3
    ny = 5
    self.S = FTSystem(xw, xe, ys, yn, nx, ny)

  def testGridSpacing(self):
    self.assertAlmostEqual(self.S._SG._dx, 125./2.)
    self.assertAlmostEqual(self.S._SG._dy, 100./4.)
    for i in range(self.S._SG._nx):
      for j in range(self.S._SG._ny):
        x = self.S._SG._X[i][j]
        y = self.S._SG._Y[i][j]
        self.assertAlmostEqual(x, self.S._SG._xw + self.S._SG._dx * i)
        self.assertAlmostEqual(y, self.S._SG._ys + self.S._SG._dy * j)

  def testHomogeneousIsotropicConductivity(self):
    self.S._SG.setConductivity(kh = 20., homogeneous = True, isotropic = True)
    for i in range(self.S._SG._nx):
      for j in range(self.S._SG._ny):
        self.assertEqual(self.S._SG._k[i][j][0], 20.)
        self.assertEqual(self.S._SG._k[i][j][1], 20.)

  def testHomogeneousThickness(self):
    self.S._SG.setThickness(H = 10., homogeneous = True)
    for i in range(self.S._SG._nx):
      for j in range(self.S._SG._ny):
        self.assertEqual(self.S._SG._H[i][j], 10.)
  
  def testHomogeneousIsotropicTransmissivity(self):
    self.assertRaises(IndexError, self.S._SG.calcTransmissivity)
    
    self.S._SG.setConductivity(kh = 20., homogeneous = True, isotropic = True)
    self.S._SG.setThickness(H = 10., homogeneous = True)
    self.S._SG.calcTransmissivity()

  def testHomogeneousIsotropicDiffusion(self):
    self.S._SG.setDiffusion(Dh = 20., homogeneous = True, isotropic = True)
    for i in range(self.S._SG._nx):
      for j in range(self.S._SG._ny):
        self.assertEqual(self.S._SG._D[i][j][0], 20.)
        self.assertEqual(self.S._SG._D[i][j][1], 20.)


if __name__ == "__main__":
  unittest.main()

