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
    # self.FTS.plotFlow()
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

  # def testTransportDirichletZeros(self):
    # self.FTS._SG.setDiffusion(Dh = 0.1)
    # self.FTS.setTransportBoundaryCondition(bctype=[0,0,0,0], bcvals=[0.,0.,0.,0.])
    # rhs = np.zeros(())
    # self.FTS.solveFlowSystem(r =  rhs)
    # self.FTS.solveTransportSystem(self.FTS._FlowGrids[0], r = rhs)
    # for i in range(self.FTS._SG._nx):
      # for j in range(self.FTS._SG._ny):
        # self.assertAlmostEqual(0, self.FTS._FlowGrids[0]._c[i][j])

  # def testDiffusionOnly(self):
    # self.FTS._SG.setDiffusion(Dh = 0.1)
    # self.FTS.setTransportBoundaryCondition(bctype=[0 ,0 ,1,1], bcvals=[1.,0.,0.,0.])
    # self.FTS._SG.setTransportType3(kPart = 0.0, cEq = 0.0)
    # rhs = np.zeros(())
    # self.FTS.solveFlowSystem(r =  rhs)
    # self.FTS.solveTransportSystem(self.FTS._FlowGrids[0], r = rhs)
    # for i in range(self.FTS._SG._nx):
      # val = 1. - 1./(self.FTS._SG._xe - self.FTS._SG._xw) * i * self.FTS._SG._dx
      # for j in range(self.FTS._SG._ny):
        # self.assertAlmostEqual(val, self.FTS._FlowGrids[0]._c[i][j])
       
  def testAdvection(self):
    self.FTS.setFlowBoundaryCondition(bctype=[0, 0, 1, 1], bcvals=[3.75, 0., 0., 0.])
    self.FTS._SG.setDiffusion(Dh = 0.)
    self.FTS.setTransportBoundaryCondition(bctype=[0 ,0 ,1,1], bcvals=[1.,0.,0.,0.])
    self.FTS._SG.setTransportType3(kPart = 0.0, cEq = 0.0)
    rhs = np.zeros(())
    self.FTS.solveFlowSystem(r =  rhs)
    rhs = np.zeros(())
    self.FTS.solveTransportSystem(self.FTS._FlowGrids[0], r = rhs)
    print self.FTS._FlowGrids[0]._h
    print self.FTS._FlowGrids[0]._c

    # for i in range(self.FTS._SG._nx):
      # for j in range(self.FTS._SG._ny):
        # if i != self.FTS._SG._nx -1 :
          # self.assertAlmostEqual(1. , self.FTS._FlowGrids[0]._c[i][j])
        # else: 
          # self.assertAlmostEqual(0., self.FTS._FlowGrids[0]._c[i][j])


if __name__ == "__main__":
  unittest.main()

