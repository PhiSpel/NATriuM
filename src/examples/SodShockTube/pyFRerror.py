from glob import glob
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from getPyFRData import getPyFRData
from recalcRef2 import SodShockAnalytic
import itertools

rc('font', **{'size': 11, 'family': 'sans-serif', 'sans-serif': ['Myriad Pro', 'Arial', 'Tahoma']})
plt.rcParams['text.usetex'] = True
plt.rcParams["savefig.dpi"] = 600
plt.rcParams['markers.fillstyle'] = 'none'
plt.rcParams['figure.constrained_layout.use'] = True

sllbm_marker = 'x'
sllbm_color = 'red'
sllbm_size = 1.5*plt.rcParams['lines.markersize']
imgtype = ".png"
transparent = False

tmax = 0.15
gamma=1.4
rL, uL, pL = 8, 0, 10
rR, uR, pR = 1, 0, 1
Nx = 6400
X = 1.
dx = X/(Nx-1)
xs = np.linspace(0,X,Nx)
iMid = Nx//2
refAnalytic = SodShockAnalytic(rL, uL, pL, rR, uR, pR, xs, iMid, tmax, gamma)
# suffix = "0.15"
# suffix = "0.15-dx1.25e-3"
suffix = "0.15-dx1.25e-3-dt5e-5"
vtuFilePath = "/home/philipp/PyFR-Test-Cases/2d-viscous-shock-tube/viscous-shock-tube-"+suffix+".vtu"
refPyFR = getPyFRData(vtuFilePath)  # xi, rho, ux, p, T

imgpath = "/mnt/c/Users/phili/Desktop/sodImages/"

LList = []  # ref nu cfl L1rho L2rho L1ux L2ux L1T L2T
failedJobs = []

for dataI, dataName in zip([1,2,3,4], ["rho", "ux", "p", "T"]):
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.plot(refPyFR[:,0], refPyFR[:,dataI], label="PyFR")
  ax.plot(refAnalytic[:,0], refAnalytic[:,dataI], label="Analytic Inviscous")
  ax.set_title(f"Sod Shock Tube: {dataName} PyFRvsAnalytic")
  ax.legend()
  fig.savefig(imgpath + dataName + "_PyFRvsAnalytic" + suffix + imgtype, transparent=transparent, dpi=300)
plt.close("all")

xi = refPyFR[:,0]
rho = refPyFR[:,1]
ux = refPyFR[:,2]
T = refPyFR[:,4]
rhoRef = np.interp(xi, refAnalytic[:,0], refAnalytic[:,1])
L1rho = np.mean(np.abs(rho - rhoRef))
L2rho = np.sqrt(np.mean(np.pow(rho - rhoRef, 2)))
uxRef = np.interp(xi, refAnalytic[:,0], refAnalytic[:,2])
L1ux = np.mean(np.abs(ux - uxRef))
L2ux = np.sqrt(np.mean(np.pow(ux - uxRef, 2)))
TRef = np.interp(xi, refAnalytic[:,0], refAnalytic[:,4])
L1T = np.mean(np.abs(T - TRef))
L2T = np.sqrt(np.mean(np.pow(T - TRef, 2)))
print(f"L1rho={L1rho}, L2rho={L2rho}, L1ux={L1ux}, L2ux={L2ux}, L1T={L1T}, L2T={L2T}")
