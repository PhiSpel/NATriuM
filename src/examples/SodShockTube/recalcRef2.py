# Libraries
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt

# Function to find the roots of!
def f(P, pL, pR, cL, cR, gg):
  a = (gg-1)*(cR/cL)*(P-1) 
  b = np.sqrt( 2*gg*(2*gg + (gg+1)*(P-1) ) )
  return P - pL/pR*( 1 - a/b )**(2.*gg/(gg-1.))

# Analtyic Sol to Sod Shock
def SodShockAnalytic(rL, uL, pL, rR, uR, pR, xs, iMid, t, gg):
  # rL, uL, pL, rR, uR, pR : Initial conditions of the Reimann problem 
  # xs: position array (e.g. xs = [0,dx,2*dx,...,(Nx-1)*dx])
  # iMid: THIS IS AN INDEX! the array index where the interface sits.
  # t: the desired solution time
  # gg: adiabatic constant 1.4=7/5 for a 3D diatomic gas
  dx = xs[1]
  Nx = len(xs)
  v_analytic = np.zeros((5,Nx),dtype='float64')

  # compute speed of sound
  cL = np.sqrt(gg*pL/rL) 
  cR = np.sqrt(gg*pR/rR)
  # compute P
  P = newton(f, 0.5, args=(pL, pR, cL, cR, gg), tol=1e-12)

  # compute region positions right to left
  # region R
  c_shock = uR + cR*np.sqrt( (gg-1+P*(gg+1)) / (2*gg) )
  iShock = iMid + int(np.floor(c_shock*t/dx))
  v_analytic[1,iShock-1:] = rR
  v_analytic[2,iShock-1:] = uR
  v_analytic[3,iShock-1:] = pR
  
  # region 2
  alpha = (gg+1)/(gg-1)
  c_contact = uL + 2*cL/(gg-1)*( 1-(P*pR/pL)**((gg-1.)/2/gg) )
  iContact = iMid + int(np.floor(c_contact*t/dx))
  v_analytic[1,iContact:iShock-1] = (1 + alpha*P)/(alpha+P)*rR
  v_analytic[2,iContact:iShock-1] = c_contact
  v_analytic[3,iContact:iShock-1] = P*pR
  
  # region 3
  r3 = rL*(P*pR/pL)**(1/gg)
  p3 = P*pR
  c_fanright = c_contact - np.sqrt(gg*p3/r3)
  iFanright = iMid + int(np.ceil(c_fanright*t/dx))
  v_analytic[1,iFanright:iContact] = r3
  v_analytic[2,iFanright:iContact] = c_contact
  v_analytic[3,iFanright:iContact] = P*pR
  
  # region 4
  c_fanleft = -cL
  iFanleft = iMid + int(np.ceil(c_fanleft*t/dx))
  u4 = 2 / (gg+1) * (cL + (xs[iFanleft:iFanright]-xs[iMid])/t )
  v_analytic[1,iFanleft:iFanright] = rL*(1 - (gg-1)/2.*u4/cL)**(2/(gg-1))
  v_analytic[2,iFanleft:iFanright] = u4
  v_analytic[3,iFanleft:iFanright] = pL*(1 - (gg-1)/2.*u4/cL)**(2*gg/(gg-1))

  # region L
  v_analytic[1,:iFanleft] = rL
  v_analytic[2,:iFanleft] = uL
  v_analytic[3,:iFanleft] = pL

  v_analytic[0] = xs
  v_analytic[4] = v_analytic[3]/v_analytic[1]  # P=rhoRT, R=1 -> T = P/rho

  return v_analytic

# Physics
gg=1.4  # gamma = C_v / C_p = 7/5 for ideal gas
# rL, uL, pL =  1.0,  0.0, 1 
# rR, uR, pR = 0.125, 0.0, .1
rL, uL, pL =  8.0,  0.0, 10
rR, uR, pR = 1, 0.0, 1

# Set Disretization
Nx = 800
X = 1.
dx = X/(Nx-1)
xs = np.linspace(0,X,Nx)
iMid = Nx//2
t = 0.2

analytic = SodShockAnalytic(rL, uL, pL, rR, uR, pR, xs, iMid, t, gg)

fig, axs = plt.subplots(1,3,figsize=(8,2), layout='constrained')
axs[0].set_title("Density")
axs[0].plot(xs,analytic[0].T)
axs[1].set_title("Velocity")
axs[1].plot(xs,analytic[1].T)
axs[1].set_yticks([0.,.2,.4,.6,.8,1.],['','','','','',''])
axs[2].set_title("Pressure")
axs[2].plot(xs,analytic[2].T)
axs[2].set_yticks([0.,.2,.4,.6,.8,1.],['','','','','',''])
# for i in range(3):
#   axs[i].set_xlim([0.,1.])
#   axs[i].set_ylim([-.05,1.05])
# plt.show()