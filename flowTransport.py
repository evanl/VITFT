from time import time, clock
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg
from grid import StaticGrid 
from grid import DynamicGrid 

class FTSystem(object):
  """ contains all of the flow and transport information. 
  Makes calls to each of the related classes such as grids and cells 
  in order to advance timesteps and plot the system
  """ 
  def __init__(self, xwest, xeast, ysouth, ynorth, nx, ny, tmax, dt):
    """ initial constructor creates the static grid of x-y locations
    """
    self._tmax = tmax
    self._dt = dt
    self._SG = StaticGrid(xwest, xeast, ysouth, ynorth, nx, ny)

  def setBoundaryCondition(self):
    return 0

  def solveFlowSystem(self):
    return 0

  def solveTransportSystem(self):
    return 0

  def flowTimeStep(self):
    return 0

  def transportTimeStep(self):
    return 0
