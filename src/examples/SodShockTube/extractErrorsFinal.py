from glob import glob
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from getData import getData
from getPyFRData import getPyFRData
from NATriuM.src.examples.SodShockTube.recalcRef import SodShockAnalytic

rc('font', **{'size': 11, 'family': 'sans-serif', 'sans-serif': ['Myriad Pro', 'Arial', 'Tahoma']})
plt.rcParams['text.usetex'] = True
plt.rcParams["savefig.dpi"] = 600
plt.rcParams['markers.fillstyle'] = 'none'
plt.rcParams['figure.constrained_layout.use'] = True

sllbm_marker = 'x'
sllbm_color = 'red'
sllbm_size = 1.5*plt.rcParams['lines.markersize']
imgtype = ".pdf"
transparent = False

# === Get highest resolution data and calculate analytical inviscous solution ===
tmax = 0.15
gamma=1.4
rL, uL, pL = 8, 0, 10
rR, uR, pR = 1, 0, 1
Nx = 3200
X = 1.
dx = X/(Nx-1)
xs = np.linspace(0,X,Nx)
iMid = Nx//2
refAnalytic = SodShockAnalytic(rL, uL, pL, rR, uR, pR, xs, iMid, tmax, gamma)
suffix = "0.15-dx1.25e-3-dt5e-5"; dtPyFR=5e-5; dxPyFR=1.25e-3
vtuFilePath = "/home/philipp/PyFR-Test-Cases/2d-viscous-shock-tube/viscous-shock-tube-"+suffix+".vtu"
refPyFR = getPyFRData(vtuFilePath)  # xi, rho, ux, p, T

imgpath = "/mnt/c/Users/phili/Desktop/sodImagesFinal/"
ref = refAnalytic
if not os.path.exists(imgpath):
  os.mkdir(imgpath)
yi = 0
zi = 0

LList = []  # ref nu cfl L1rho L2rho L1ux L2ux L1T L2T

LListPath = imgpath + "LList.npy"
if os.path.exists(LListPath):
  LList = np.load(LListPath)
else:
  failedJobs = []
  for jobFolder in glob("/mnt/c/Users/phili/Desktop/sod/*nx*p4/"):
    data = getData(jobFolder)
    if type(data) == str:
      failedJobs.append(data)
      continue
    else:
      data, cs, dt, dx, jobid, cfl, p, jobName, nx, lastIteration = data

    xi = data[:,0]
    rho = data[:,1]
    ux = data[:,2]
    T = data[:,4]
    rhoRef = np.interp(xi, ref[:,0], ref[:,1])
    L1rho = np.mean(np.abs(rho - rhoRef))
    L2rho = np.sqrt(np.mean(np.pow(rho - rhoRef, 2)))
    uxRef = np.interp(xi, ref[:,0], ref[:,2])
    L1ux = np.mean(np.abs(ux - uxRef))
    L2ux = np.sqrt(np.mean(np.pow(ux - uxRef, 2)))
    TRef = np.interp(xi, ref[:,0], ref[:,4])
    L1T = np.mean(np.abs(T - TRef))
    L2T = np.sqrt(np.mean(np.pow(T - TRef, 2)))
    LList.append([nx, cfl, L1rho, L2rho, L1ux, L2ux, L1T, L2T, p, dt, dx])

  fs = "failedJobs: "
  for failedJob in failedJobs:
    fs += " " + failedJob

  LList = np.array(LList)
  np.save(LListPath, LList)

allDt = np.sort(np.unique(LList[:,9]))
allCfl = np.sort(np.unique(LList[:,1]))

data, _, _, _, _, _, _, jobName, _, lastIteration = getData("/mnt/c/Users/phili/Desktop/sod/10906738_nx400_cfl1_p4/")
fig, axs = plt.subplots(1,2,figsize=[7, 3.5])
for ax in axs:
  ax.plot(data[:,0], data[:,1], color=sllbm_color, label="SLLBM")
  ax.plot(refPyFR[:,0], refPyFR[:,1], label="4th-order Runge-Kutta")
  ax.plot(refAnalytic[:,0], refAnalytic[:,1], linestyle="--", label="Analytic Inviscous")
axs[0].add_artist(plt.Rectangle((.6,.5),.25,3.5, linestyle="--", edgecolor=".9", facecolor="none"))
axs[0].set_xlabel("$x$")
axs[0].set_ylabel(r"$\rho$")
axs[1].set_xlabel("$x$")
axs[1].set_xlim((.6,.85))
axs[1].set_ylim((.5,4))
axs[1].legend()
fig.savefig(imgpath + "rho_" + jobName + "_iT" + lastIteration + "_both" + imgtype, transparent=transparent, dpi=300)

i = 3
dataName = "rho"
dataLabel = "Density"

xOrder = np.array([1e-4, 1e-3, 1e-2, 1e-1])
# PyFR, dx = 1.25e-3, dt = 5e-5: L1rho=0.004344615895270187, L2rho=0.023828144955520027, L1ux=0.0011521760614486276, L2ux=0.01291165214905785, L1T=0.000946764660822277, L2T=0.009211594216603015
rhoL1pyFR = 0.004344615895270187


fig, ax = plt.subplots(figsize=[7, 4])
ref_m = ['o', 'v', '^', 's', 'p', 'h', 'D']
ref_c = ['royalblue', 'green', 'grey', 'black', 'cyan', 'magenta']
ref_c = ['C0', 'C1', 'C2', 'C3', "C4", "C5", "C6"]
used_dt_indices = set()
used_cfl_indices = set()
dxmin = 9e-4
LListI = LList[LList[:,10]>dxmin]
allDt = np.sort(np.unique(LListI[:,9]))
allDx = np.sort(np.unique(LListI[:,10]))
for iDt in range(len(allDt)):
  dt = allDt[iDt]
  LListI = LList[LList[:,9]==dt]
  for iCfl in range(len(allCfl)):
    cfl = allCfl[iCfl]
    LListIJ = LListI[LListI[:,1]==cfl]
    if len(LListIJ[:,0] > 0):
      for LListIJK in LListIJ[LListIJ[:,10] > dxmin]:
        ax.scatter(LListIJ[:,10], LListIJ[:,2], marker=ref_m[iCfl], color=ref_c[iDt], label=f"cfl = {cfl}")
        used_dt_indices.add(iDt)
        used_cfl_indices.add(iCfl)
line_order1 = ax.plot(xOrder, xOrder, label="Order 1", linestyle='--', color='.3')[0]
line_rk4 = ax.scatter(dxPyFR, rhoL1pyFR, color=".3", marker='D', label="4th-order Runge-Kutta,")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r"$L^1$ Norm")
ax.set_xlabel(r"$\delta x$")
ax.set_ylim((9e-4,7e-2))
ax.set_xlim((dxmin,3e-2))
ax.set_axisbelow(True)
ax.grid(which="minor", color="0.9")
ax.grid(which='major', color=".8")
legend_elements = []
# Color entries (representing dt)
for i in sorted(list(used_dt_indices)):
    legend_elements.append(
        plt.Line2D([0], [0], marker='', color=ref_c[i], label=f"$\delta t$ = {allDt[i]:.1e}", linestyle='-',
               markerfacecolor=ref_c[i], markersize=8)
    )
# Spacing
dummy_handle = plt.Line2D([0], [0], color='none', label='')
legend_elements.extend([dummy_handle])
# Marker entries (representing cfl)
for i in sorted(list(used_cfl_indices)):
    legend_elements.append(
        plt.Line2D([0], [0], marker=ref_m[i], color='k', label=f"cfl = {allCfl[i]}", linestyle='',
               markerfacecolor='k', markeredgecolor='k', markersize=8)
    )
# References
legend_elements.append(line_order1)
legend_elements.append(line_rk4)
# Spacing
dummy_handle = plt.Line2D([0], [0], color='none', label=f' $\delta t$ = {dtPyFR:.0e}')
legend_elements.extend([dummy_handle])
fig.legend(handles=legend_elements, loc='outside upper center', columnspacing=1, framealpha=.9, ncols=4,handlelength=1)
fig.savefig(imgpath + dataName + "_L1_overDxFinal" + imgtype, transparent=transparent, dpi=300)