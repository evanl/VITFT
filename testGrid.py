import grid
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
    self.sg = grid.StaticGrid(xw, xe, ys, yn, nx, ny)
  def testGridSpacing(self):
    self.assertAlmostEqual(self.sg._dx, 125./2.)
    self.assertAlmostEqual(self.sg._dy, 100./4.)
    for i in range(self.sg._nx):
      for j in range(self.sg._ny):
        x = self.sg._X[i][j]
        y = self.sg._Y[i][j]
        self.assertAlmostEqual(x, self.sg._xw + self.sg._dx * i)
        self.assertAlmostEqual(y, self.sg._ys + self.sg._dy * j)

  def testHomogeneousIsotropicConductivity(self):
    self.sg.setConductivity(kh = 20., homogeneous = True, isotropic = True)
    for i in range(self.sg._nx):
      for j in range(self.sg._ny):
        self.assertEqual(self.sg._k[i][j][0], 20.)
        self.assertEqual(self.sg._k[i][j][1], 20.)

  def testHomogeneousThickness(self):
    self.sg.setThickness(H = 10., homogeneous = True)
    for i in range(self.sg._nx):
      for j in range(self.sg._ny):
        self.assertEqual(self.sg._H[i][j], 10.)
  
  def testHomogeneousIsotropicTransmissivity(self):
    self.assertRaises(IndexError, self.sg.calcTransmissivity)
    
    self.sg.setConductivity(kh = 20., homogeneous = True, isotropic = True)
    self.sg.setThickness(H = 10., homogeneous = True)
    self.sg.calcTransmissivity()





if __name__ == "__main__":
  unittest.main()

