from glob import glob
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from getData import getData
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

# === Get highest resolution data and calculate analytical inviscous solution ===
refListHighRes, _, dt, dx, _, _, _, _, nx, lastIteration = getData("/mnt/c/Users/phili/Desktop/sod/10908253_nx6400_cfl1_p4/")
tmax = float(dt)*int(lastIteration)/np.sqrt(3)  # should be 0.15
print("Getting ref from iteration", int(lastIteration), "dt=", float(dt), "tmax=", tmax)
gamma=1.4
rL, uL, pL = 8, 0, 10
rR, uR, pR = 1, 0, 1
Nx = 6400
X = 1.
dx = X/(Nx-1)
xs = np.linspace(0,X,Nx)
iMid = Nx//2
refAnalytic = SodShockAnalytic(rL, uL, pL, rR, uR, pR, xs, iMid, tmax, gamma)
refPyFR = getPyFRData()  # xi, rho, ux, p, T

for refSource in ["recalc", "PyFR"]:  # "recalc"  # "PyFR"  # "highRes"  # "txt"
  imgpath = "/mnt/c/Users/phili/Desktop/sodImages/"
  if refSource == "txt":
    imgpath += "RefOld/"
    ref = np.loadtxt("/mnt/c/Users/phili/Desktop/eval_shocktube_norm/ref14.txt")
  elif refSource == "highRes":
    imgpath += "RefHighRes/"
    ref = refListHighRes
  elif refSource == "PyFR":
    imgpath += "PyFR/"
    ref = refPyFR
  elif refSource == "recalc":
    imgpath += "Recalc/"
    ref = refAnalytic
  else:
    raise NotImplementedError("No valid reference source")
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
      
      if cfl==1:
        for dataI, dataName in zip([1,2,3,4], ["rho", "ux", "p", "T"]):
          fig, ax = plt.subplots(figsize=[7, 3.5])
          ax.scatter(data[:,0], data[:,dataI], marker=sllbm_marker, color=sllbm_color, s=sllbm_size, label="SLLBM")
          ax.plot(refPyFR[:,0], refPyFR[:,dataI], label="PyFR")
          ax.plot(refAnalytic[:,0], refAnalytic[:,dataI], label="Analytic Inviscous")
          ax.set_title(f"Sod Shock Tube: {dataName} (dx={dx}, dt={dt})")
          ax.legend()
          fig.savefig(imgpath + dataName + "_" + jobName + "_iT" + lastIteration + imgtype, transparent=transparent, dpi=300)
          # plt.show()
        plt.close("all")

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
      print("")

    fs = "failedJobs: "
    for failedJob in failedJobs:
      fs += " " + failedJob

    LList = np.array(LList)
    np.save(LListPath, LList)

  dt=5.52427e-05
  LListI = LList[(LList[:,9]-dt)<1e-10]
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("dx")
  ax.set_title(f"Density L2 Errors over dx, dt = {dt}")
  ax.scatter(LListI[:,10], LListI[:,3])
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim((8e-4, 2e-2))
  ax.set_ylim((8e-3, 1e-1))
  ax.set_aspect('equal', 'box')
  # ax.legend()
  fig.savefig(f"/mnt/c/Users/phili/Desktop/sodImages/rho_L2vs{refSource}_overDx" + imgtype, transparent=transparent, dpi=300)

  plt.close("all")