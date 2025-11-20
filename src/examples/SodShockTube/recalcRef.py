import numpy as np
import sympy as sp
from sympy.abc import y

gamma = 1.4

def recalcRef(tmax):
  xmin=0; xmax=1; nx=1600
  PfromR = False
  R = 1

  rho1 = 8
  T1 = 1.25
  ux1 = 0
  if PfromR:
    p1 = rho1*R*T1
  else:
    p1 = (gamma-1)*rho1*T1
  a1 = np.sqrt(gamma*p1/rho1)

  rho5 = 1
  T5 = 1
  ux5 = 0
  if PfromR:
    p5 = rho5*R*T5
  else:
    p5 = (gamma-1)*rho5*T5

  # p3 = p4
  # ux3 = ux1 + 2*a1/(gamma-1)*(1-(p3/p1)**((gamma-1)/(2*gamma)))
  # ux4 = ux5 + (p4-p5)*np.sqrt(2/(p5*(gamma+1)*(p4+(gamma-1)/(gamma+1)*p5)))

  # solve for pmid by setting ux3=ux4
  # ux3ux4 = sp.Eq(ux1 + 2*a1/(gamma-1)*(1-(y/p1)**((gamma-1)/(2*gamma))),
  #              ux5 + (y-p5)*(2/(p5*(gamma+1)*(y+(gamma-1)/(gamma+1)*p5)))**0.5)
  # solveset = sp.solve(ux3ux4, y, domain=sp.Interval(rho1, rho5), simplify=False)
  # pMid = float(solveset[0])
  # pMid = 0.943160123654987
  pMid = 1.288
  p3 = pMid
  p4 = pMid
  ux3 = ux1 + 2*a1/(gamma-1)*(1-(pMid/p1)**((gamma-1)/(2*gamma)))
  # u3 = 0.613
  ux4 = ux3

  rho3 = rho1*(p3/p1)**(1/gamma)
  rho4 = rho5*(((gamma+1)*p4+(gamma-1)*p5)/((gamma-1)*p4+(gamma+1)*p5))

  a3 = np.sqrt(gamma*p3/rho3)
  a5 = np.sqrt(gamma*p5/rho5)

  x = np.linspace(xmin,xmax,nx)
  xmid = xmin+(xmax-xmin)/2

  x12 = xmid+(ux1-a1)*tmax
  x23 = xmid+(ux3-a3)*tmax
  x34 = xmid+ux3*tmax
  a5hat = a5*np.sqrt((gamma+1)/(2*gamma)*p4/p5 + (gamma-1)/(2*gamma))
  x45 = xmid + (ux5 + a5hat)*tmax
  
  x2 = np.logical_and(x>=x12, x<x23)
  x3 = np.logical_and(x>=x23, x<x34)
  x4 = np.logical_and(x>=x34, x<x45)

  chi = (x-xmid)/tmax
  uxi = 2/(gamma+1)*(a1 + chi)
  axi = a1 - (gamma-1)/2*uxi

  rho = np.zeros_like(x)
  rho[x<x12] = rho1
  rho[x2] = (rho1*(axi/a1)**(2/(gamma-1)))[x2]
  rho[x3] = rho3
  rho[x4] = rho4
  rho[x>=x45] = rho5

  ux = np.zeros_like(x)
  ux[x<x12] = ux1
  ux[x2] = (ux1*(axi/a1)**(2/(gamma-1)))[x2]
  ux[x3] = ux3
  ux[x4] = ux4
  ux[x>=x45] = ux5

  p = np.zeros_like(x)
  p[x<x12] = p1
  p[x2] = (p1*(axi/a1)**(2/(gamma-1)))[x2]
  p[x3] = p3
  p[x4] = p4
  p[x>=x45] = p5

  T = np.zeros_like(x)
  T[x<x12] = T1
  T[x2] = (T1*(axi/a1)**(2/(gamma-1)))[x2]
  T[x3] = T3
  T[x4] = T4
  T[x>=x45] = T5

  ref = np.array([x,rho,ux,p,T]).T
  return ref  # ref[:,0] = xRef, ref[:,1] = rhoRef, ref[:,2]=uxRef, ref[:,3]=pRef, ref[:,4]=Tref