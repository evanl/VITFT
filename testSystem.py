from flowTransport import FTSystem
import unittest
import numpy as np

class TestSystem(unittest.TestCase):
  """ Runs a handful of tests on a simple grid shown in the setUp method.
  """
  def setUp(self):
    xw = -75.
    xe = 50.
    ys = 0.
    yn = 100.
    nx = 5
    ny = 5
    self.FTS = FTSystem(xw, xe, ys, yn, nx, ny)
    self.FTS.setFlowBoundaryCondition(bctype=[0, 0, 0, 0], bcvals=[0, 0, 0, 0])
    self.FTS._SG.addSourceFlow()
    self.FTS._SG.setConductivity()
    self.FTS._SG.setThickness()
    self.FTS._SG.setStorage()
    self.FTS._SG.calcTransmissivity()
    self.FTS._SG.setPorosity()

  def testFlowHydrostatic(self):
    rhs = np.zeros(())
    self.FTS.solveFlowSystem(r =  rhs)
    for i in range(self.FTS._SG._nx):
      for j in range(self.FTS._SG._ny):
        self.assertAlmostEqual(self.FTS._FlowGrids[0]._h[i][j], 0.)
    
  def testLinearX(self):
    self.FTS.setFlowBoundaryCondition(bctype=[0, 0, 1, 1], bcvals=[10., 0., 0., 0.])
    rhs = np.zeros(())
    self.FTS.solveFlowSystem(r =  rhs)
    for i in range(self.FTS._SG._nx):
      val = 10. - 10./(self.FTS._SG._xe - self.FTS._SG._xw) * i * self.FTS._SG._dx
      for j in range(self.FTS._SG._ny):
        self.assertAlmostEqual(self.FTS._FlowGrids[0]._h[i][j], val)

  def testLinearY(self):
    self.FTS.setFlowBoundaryCondition(bctype=[1, 1, 0, 0], bcvals=[0., 0., 10., 0.])
    rhs = np.zeros(())
    self.FTS.solveFlowSystem(r =  rhs)
    for i in range(self.FTS._SG._nx):
      for j in range(self.FTS._SG._ny):
        val = 10. - 10./(self.FTS._SG._yn - self.FTS._SG._ys) * j * self.FTS._SG._dy
        self.assertAlmostEqual(self.FTS._FlowGrids[0]._h[i][j], val)

if __name__ == "__main__":
  unittest.main()

