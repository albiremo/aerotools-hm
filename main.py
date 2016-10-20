#import
#-----------------------------------------------------
import numpy             as np
from   scipy.linalg      import solve,solve_banded
import matplotlib        as mp
mp.use("Qt4Agg")
import scipy             as sp
import matplotlib.pyplot as plt
from   numpy   import pi
#-----------------------------------------------------
#-----------------------------------------------------

from geometry import *

#main
xv,yv=NACA('21015',100)
print(xv,yv)
plt.figure()
plt.plot(xv,yv,'-*')
plt.axis('equal')
plt.show()
