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

  def testFlowHydrostatic(self):
    rhs = np.zeros(())
    self.FTS.solveFlowSystem(r =  rhs)
    for i in range(self.FTS._SG._nx):
      for j in range(self.FTS._SG._ny):
        self.assertAlmostEqual(self.FTS._FlowGrids[0]._h[i][j], 0.)
    


if __name__ == "__main__":
  unittest.main()

