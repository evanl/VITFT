import numpy as np

class StaticGrid(object):
  def __init__(self, xwest, xeast, ysouth, ynorth, nx, ny):
    """creates grid with uniform grid spacing from xwest to xeast with nx grid cells,
    from ysouth to ynorth with ny grid cells and stores them in two numpy arrays that 
    are indexed by 
    X[i,j], Y[i,j]
    or 
    if the numpy array is not being used, they should be indexed by 
    X[i][j], Y[i][j]
    """
    # boundaries and grid
    self._xw = xwest
    self._xe = xeast
    self._ys = ysouth
    self._yn = ynorth
    self._nx = nx
    self._ny = ny
    self._dx = (xeast - xwest) /(nx - 1)
    self._dy = (ynorth - ysouth) / (ny - 1)
    X,Y = np.mgrid[xwest:xeast:nx*1j, ysouth:ynorth:ny*1j]
    self._Y = Y
    self._X = X
    # non-transient physical parameters
    self._k = []
    self._H = []
    self._T = []
    self._D = []
    self._S = []
    self._poro = []
    self._flowSource =[]
    self._flowType3 = []
    self._kPart = 0.0
    self._cEq = 0.0

  def setConductivity(self, kh = 10., homogeneous = True, isotropic = True):
    """ This function sets the hydraulic conductivity of each cell. 
    The conductivity is stored as a tuple for kx and ky so that
    k_{ij}x = k[i][j][0]
    k_{ij}y = k[i][j][1]
    """
    # clears previous hydraulic conductivity
    self._k = []
    if homogeneous == True:
      if isotropic == True:
        for i in range(self._nx):
          templist = []
          for j in range(self._ny):
            templist.append((kh, kh))
          self._k.append(templist)

    return 0

  def setThickness(self, H = 10., homogeneous = True ):
    """ This function sets the thickness of each cell. 
    the thickness is stored such that 
    H_{ij} = H[i][j]
    """
    # clears previous thickness
    self._H = []
    if homogeneous == True:
      for i in range(self._nx):
        templist = []
        for j in range(self._ny):
          templist.append(H)
        self._H.append(templist)
    return 0

  def setPorosity(self, poro = 0.30, homogeneous = True ):
    """ This function sets the porosity of each cell. 
    the porosity is stored such that 
    poro{ij} = poro[i][j]
    """
    # clears previous porosity
    self._poro = []
    if homogeneous == True:
      for i in range(self._nx):
        templist = []
        for j in range(self._ny):
          templist.append(poro)
        self._poro.append(templist)
    
  def setStorage(self, S = 0.2, homogeneous = True ):
    """ This function sets the specific storage of each cell. 
    the storage is stored such that 
    S_{ij} = S[i][j]
    """
    # clears previous Storage
    self._S = []
    if homogeneous == True:
      for i in range(self._nx):
        templist = []
        for j in range(self._ny):
          templist.append(S)
        self._S.append(templist)
    return 0

  def calcTransmissivity(self):
    """ calculates the transmissivity for the grid, only works if the conductivity
    and the transmissivity have been filled and have dimensions of nx X ny 
    """
    for i in range(self._nx):
      tempT = []
      for j in range(self._ny):
        Tx = self._k[i][j][0] * self._H[i][j]
        Ty = self._k[i][j][1] * self._H[i][j]
        tempT.append((Tx,Ty))
      self._T.append(tempT)
    return 0

  def setDiffusion(self, Dh = 0.1, homogeneous = True, isotropic = True):
    """ This function sets the Diffusion coefficient for each cell. 
    The conductivity is stored as a tuple for Dx and Dy so that
    D_{ij}x = D[i][j][0]
    D_{ij}y = D[i][j][1]
    """
    # clears previous diffusion
    self._D = []
    if homogeneous == True:
      if isotropic == True:
        for i in range(self._nx):
          templist = []
          for j in range(self._ny):
            templist.append((Dh, Dh))
          self._D.append(templist)
    return 0  

  def addSourceFlow(self, sourcefunc = 0):
    """ add a source function that is of the form:
    N = f(x,y,t)
    where N is the outflux rate, 
    x is the global x, 
    y is the global y,
    t is the time of evaluation
    
    """
    self._flowSource = []
    for i in range(self._nx):
      templist = []
      for j in range(self._ny):
        templist.append(sourcefunc)
      self._flowSource.append(templist)

    return 0

  def addType3Flow(self, type3flowfunc = 0.):
    """adds a type 3 vertical boundary condition of the form
    C,R = f(x,y,kp,h0)
    where C is the coefficient [1/T] that is in the solution matrix
    R is the right hand side term [L/T], 
    x is the global x,
    y is the global y
    kp is the conductivity of the leaky layer
    h0 is the head above or below the leaky layer
    """
    return 0

  def setTransportType3(self, kPart = 0.0, cEq = 0.0):
    self._kPart = kPart
    self._cEq = cEq
    return 0

class DynamicGrid(object):
  """ This grid must have the sae dimensions as any StaticGrid in the overall
  simulation. 

  """
  def __init__(self, nx, ny, dx, dy):
    self._nx = nx
    self._ny = ny
    self._dx = dx
    self._dy = dy
    self._h = np.zeros((ny, nx))
    self._c = np.zeros((ny, nx))

  def setHead(self, head, flowbctype, flowbcvals):
    """ This function takes in a head vector from a sparse matrix solve
    and stores it as a 2D array """
    # handles inner, sparse solved values
    for i in xrange(1,self._nx-1):
      for j in xrange(1,self._ny-1):
        self._h[i][j] = head[(i-1) + (self._nx - 2) * (j-1)]

    # handles boundary values
    
    for i in range(self._ny):
      #west
      if flowbctype[0] == 0:
        self._h[0,i] = flowbcvals[0]            
      else:
        self._h[0,i] = self._h[1,i] + 1.0/self._dx * flowbcvals[0]
      # east
      if flowbctype[1] == 0:
        self._h[-1,i] = flowbcvals[1]
      else:
        self._h[-1,i] = self._h[-2,i] - 1.0/self._dx * flowbcvals[1]

    for i in xrange(0,self._nx):
      # south
      if flowbctype[2] == 0:
        self._h[i,0] = flowbcvals[2]
      else:
        self._h[i,0] = self._h[i,1] + 1.0/self._dy * flowbcvals[2]
      # north
      if flowbctype[3] == 0:
        self._h[i,-1] = flowbcvals[3]
      else:
        self._h[i,-1] = self._h[i,-2] + 1/self._dy * flowbcvals[3]

    return 0

  def setConcentration(self, conc, transbctype, transbcvals):
    """ this function takes in a concentration vector from a sparse matrix:
    solve and stores it as a 2D array"""

    # handles inner values
    for i in xrange(1,self._nx-1):
      for j in xrange(1,self._ny-1):
        self._c[i][j] = conc[(i-1) + (self._nx - 2) * (j-1)]

    # handles boundary values
    
    for i in range(self._ny):
      #west
      if transbctype[0] == 0:
        self._c[0,i] = transbcvals[0]            
      else:
        self._c[0,i] = self._c[1,i] + 1.0/self._dx * transbcvals[0]
      # east
      if transbctype[1] == 0:
        self._c[-1,i] = transbcvals[1]
      else:
        self._c[-1,i] = self._c[-2,i] - 1.0/self._dx * transbcvals[1]

    for i in xrange(0,self._nx):
      # south
      if transbctype[2] == 0:
        self._c[i,0] = transbcvals[2]
      else:
        self._c[i,0] = self._c[i,1] + 1.0/self._dy * transbcvals[2]
      # north
      if transbctype[3] == 0:
        self._c[i,-1] = transbcvals[3]
      else:
        self._c[i,-1] = self._c[i,-2] + 1/self._dy * transbcvals[3]

    return 0
