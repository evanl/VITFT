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
  def __init__(self, xwest, xeast, ysouth, ynorth, nx, ny, 
      tmax = 0, dt = 1):
    """ initial constructor creates the static grid of x-y locations
    """
    self._tmax = tmax
    self._dt = dt
    self._SG = StaticGrid(xwest, xeast, ysouth, ynorth, nx, ny)
    self._FlowGrids = []

  def setFlowBoundaryCondition(self, bctype = [0, 0, 0, 0], bcvals = [0, 0, 0, 0]):
    """ given in lists with the order:  [west, east, south, north]
    The integer value corresponds to the type of boundary condition:
    0: Dirichlet, 1: Neumann
    bctype lists all types. for example, bctype = [1,1,0,0] 
    Neumann at E and W; Dirichlet at N and S

    bcval is a list that describes the boundary values.
    if bctype[i] == 1: bcvals[i] = normal flux
    if bctype[i] == 0: bcvals[i] = head
    """
    self._flowbctype = bctype
    self._flowbcvals = bcvals
    return 0

  def setTransportBoundaryCondition(self, bctype=[0, 0, 0, 0], bcvals=[0, 0, 0, 0]):
    """ given in lists with the order:  [west, east, south, north]
    The integer value corresponds to the type of boundary condition:
    0: Dirichlet, 1: Neumann
    bctype lists all types. for example, bctype = [1,1,0,0] 
    Neumann at E and W; Dirichlet at N and S

    bcval is a list that describes the boundary values.
    if bctype[i] == 1: bcvals[i] = normal flux
    if bctype[i] == 0: bcvals[i] = concentration
    """
    self._transbctype = bctype
    self._transbcvals = bcvals
    return 0

  def solveFlowSystem(self, r = [] ):
    """ solves flow system for a single timestep. This system is solved implicitly
    and requires an initialization of the right hand side values from the 
    previous timestep.
    If this is an initial timestep, leave the rhs value as an empty list
    returns the solved column vector of inner head values.
    """
    print "Solving Flow Timestep"
    nx = self._SG._nx
    ny = self._SG._ny
    dx = self._SG._dx
    dy = self._SG._dy

    N = (nx - 2) * (ny - 2)
    Loff = np.zeros((N,1))
    L    = np.zeros((N,1))
    R    = np.zeros((N,1))
    Roff = np.zeros((N,1))

    C = np.zeros((N,1))

    if r.all() == False:    
      r = np.zeros((N,1))

    tinit = time()

    for i in xrange(N):
      
      # operations to determine the center node's global x and y coordinates. 
      row = ( i % (nx-2) ) + 1        
      col = i // (nx-2) + 1 
      if col == 0 :
        col = (nx-2)


      x = self._SG._xw + self._SG._dx * col
      y = self._SG._ys + self._SG._dy * row

      if r.all() != False:
        C[i] += self._SG._S[row][col] /dt

      r[i] += dx * dy * self._SG._flowSource[row][col]

      # Type 3 (leakage)
      # r[i] -= K_B * h0 * dx * dy
      # C[i] -= K_B * dx * dy

      TC = self._SG._T[row][col]
      TW = self._SG._T[row - 1][col]
      TE = self._SG._T[row + 1][col]
      TS = self._SG._T[row][col - 1]
      TN = self._SG._T[row][col + 1]

      CW =  dy /  dx * HarmAvg(TC[0], TW[0], dx, dx)
      CE =  dy /  dx * HarmAvg(TC[0], TE[0], dx, dx)
      CS =  dx /  dy * HarmAvg(TC[1], TS[1], dy, dy)
      CN =  dx /  dy * HarmAvg(TC[1], TN[1], dy, dy)

      # am i on west boundary
      if  row == 1 :
        if self._flowbctype[0] == 1:
          r[i] -= self._flowbcvals[0]
        else:
          r[i] -= CW * self._flowbcvals[0]
          C[i] -= CW
      else:
        L[i-1] = CW
        C[i] -= CW

      # east boundary
      if row == (nx-2):
        if self._flowbctype[1] == 1:
          r[i] += self._flowbcvals[1]
        else:
          r[i] -= CE * self._flowbcvals[1]
          C[i] -= CE
      else:
        R[i+1] = CE
        C[i] -= CE

      # south boundary
      if i < (nx-2):
        if self._flowbctype[2] ==  1:
          r[i] -= self._flowbcvals[2]
        else:
          r[i] -= CS * self._flowbcvals[2]
          C[i] -= CS
      else:
        Loff[i-(nx-2)] = CS
        C[i] -= CS

      # north boundary
      if i >= (N - (nx-2)):
        if self._flowbctype[3] == 1:
          r[i] += self._flowbcvals[3]
        else:
          r[i] -= CN * self._flowbcvals[3]
          C[i] -= CN
      else:
        Roff[i+(nx-2)] = CN
        C[i] -= CN

    # end of  loop
    tpop = time()
    print "matrix create time = " +str(tpop - tinit)

    # create new sparse matrix from pieces 
    data = np.hstack((Loff, L, C, R, Roff))
    data = np.transpose(data)
    diags = np.array([-(nx-2), -1, 0, 1, (nx-2)])
    Abuilt = sp.spdiags(data,diags,N,N)

    
    # print Abuilt.todense()
    # print r
    b = splinalg.spsolve(Abuilt,r)

    tmatsolve = time()
    print "linalg solve time=  " + str( tmatsolve - tpop)

    dg = DynamicGrid(self._SG._nx, self._SG._ny, self._SG._dx, self._SG._dy)
    dg.setHead(b, self._flowbctype, self._flowbcvals)
    self._FlowGrids.append(dg)
    return 0

  def solveTransportSystem(self, FG, r = []):
    """ solves flow system for a single timestep. This system is solved implicitly
    and requires an initialization of the right hand side values from the 
    previous timestep.
    If this is an initial timestep, leave the rhs value as an empty list
    returns the solved column vector of inner head values.
    requires the previous flow grid.
    """
    print "Solving Transport Timestep"
    nx = self._SG._nx
    ny = self._SG._ny
    dx = self._SG._dx
    dy = self._SG._dy
    N = (nx - 2) * (ny - 2)
    Loff = np.zeros((N,1))
    L    = np.zeros((N,1))
    R    = np.zeros((N,1))
    Roff = np.zeros((N,1))
    C = np.zeros((N,1))

    if r.all() == False:    
      r = np.zeros((N,1))

    tinit = time()

    for i in xrange(N):
      
      # specify index of row and column in grid
      row = ( i % (nx-2) ) + 1        
      col = i // (nx-2) + 1 
      if col == 0 :
        col = (nx-2)
      # determine the center node's global x and y coordinates. 
      x = self._SG._xw + self._SG._dx * col
      y = self._SG._ys + self._SG._dy * row

      if r.all() != False:
        C[i] += self._SG._S[row][col] /dt

      # Type 3 (leakage)
      r[i] -= self._SG._kPart * self._SG._cEq * dx * dy * self._SG._H[row][col]
      C[i] -= self._SG._kPart * dx * dy * self._SG._H[row][col]

      # get advective face velocities. 
      # picks concentration for upstream weighting
      #west
      vw = self._SG._k[row][col][0] / (self._SG._dx * self._SG._poro[row][col])* \
          (FG._h[row][col] - FG._h[row - 1][col])
      #east
      ve = self._SG._k[row][col][0] / (self._SG._dx * self._SG._poro[row][col])* \
          (FG._h[row + 1 ][col] - FG._h[row][col])
      #south
      vs = self._SG._k[row][col][1] / (self._SG._dy * self._SG._poro[row][col])* \
          (FG._h[row][col] - FG._h[row][col-1])
      #north
      vn = self._SG._k[row][col][1] / (self._SG._dy * self._SG._poro[row][col])* \
          (FG._h[row ][col + 1] - FG._h[row][col])

      # diffusion coefficients
      CW =  dy / dx * self._SG._D[row][col][0]
      CE =  dy / dx * self._SG._D[row][col][0]
      CS =  dx / dy * self._SG._D[row][col][1]
      CN =  dx / dy * self._SG._D[row][col][1]

      # am i on west boundary
      if  row == 1 :
        if self._transbctype[0] == 1:
          r[i] -= self._transbcvals[0]
        else:
          # adv
          r -= self._transbcvals[0] * vw
          # diff
          r[i] -= CW * self._transbcvals[0]
          C[i] -= CW
      else:
        # adv
        if vw >= 0.:
          L[i-1] -= vw
        else:
          C[i] -= vw
        # diff
        L[i-1] += CW
        C[i] -= CW

      # east boundary
      if row == (nx-2):
        if self._transbctype[1] == 1:
          r[i] += self._transbcvals[1]
        else:
          # adv
          r -= self._transbcvals[1] * ve
          # diff
          r[i] -= CE * self._transbcvals[1]
          C[i] -= CE
      else:
        if ve >= 0.:
          R[i+1] = -ve
        else:
          C[i] -= ve
        R[i+1] += CE
        C[i] -= CE

      # south boundary
      if i < (nx-2):
        if self._transbctype[2] ==  1:
          r[i] -= self._transbcvals[2]
        else:
          # adv
          r -= self._transbcvals[2] * vs
          # diff
          r[i] -= CS * self._transbcvals[2]
          C[i] -= CS
      else:
        # adv
        if vs >= 0.:
          Loff[i-(nx-2)]  = vs
        else:
          C[i] -= vs
        # diff
        Loff[i-(nx-2)] += CS
        C[i] -= CS

      # north boundary
      if i >= (N - (nx-2)):
        if self._transbctype[3] == 1:
          r[i] += self._transbcvals[3]
        else:
          # adv
          r -= self._transbcvals[3] * vw
          # diff
          r[i] -= CN * self._transbcvals[3]
          C[i] -= CN
      else:
        if vw >= 0.:
          Roff[i+(nx-2)] -= vw
        else:
          C[i] -= vw
        Roff[i+(nx-2)] += CN
        C[i] -= CN

    # end of  loop
    tpop = time()
    print "matrix create time = " +str(tpop - tinit)

    # create new sparse matrix from pieces 
    data = np.hstack((Loff, L, C, R, Roff))
    data = np.transpose(data)
    diags = np.array([-(nx-2), -1, 0, 1, (nx-2)])
    Abuilt = sp.spdiags(data,diags,N,N)

    
    # print Abuilt.todense()
    # print r 
    b = splinalg.spsolve(Abuilt,r)

    tmatsolve = time()
    print "linalg solve time=  " + str( tmatsolve - tpop)

    FG.setConcentration(b, self._transbctype, self._transbcvals)
    return 0
    return 0

  def flowTimeStep(self):
    return 0

  def transportTimeStep(self):
    return 0

  def plotFlow(self, timestep = 0):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    surf = ax.contourf(self._SG._X, self._SG_Y, self._FlowGrids[timestep]._h, 55, \
        cmap = cm.Greys, linewidth = (3,))
    CB = plt.colorbar(surf) 
    ax.set_title('Hydraulic Head Contours [m]')
    ax.set_xlabel('x-direction [m]')
    ax.set_ylabel('y-direction [m]')
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    plt.show()
    return 0
  
  def plotTransport(self, timestep = 0):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    surf = ax.contourf(self._SG._X, self._SG_Y, self._FlowGrids[timestep]._c, 55, \
        cmap = cm.Greys, linewidth = (3,))
    CB = plt.colorbar(surf) 
    ax.set_title('Concentration')
    ax.set_xlabel('x-direction [m]')
    ax.set_ylabel('y-direction [m]')
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    plt.show()
    return 0

# end FlowSystem

# convenience function for harmonic mean
def HarmAvg(T1, T2, dx1, dx2):
    return (dx1 + dx2) / ( (dx1/T1) + (dx2/T2))
