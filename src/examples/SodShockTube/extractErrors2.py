import os
from glob import glob
import vtk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import sympy as sp
from sympy.abc import y

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

imgpath = "/mnt/c/Users/phili/Desktop/sodImages/"
yi = 0
zi = 0

LList = []  # ref nu cfl L1rho L2rho L1ux L2ux L1T L2T
refOld = np.loadtxt("/mnt/c/Users/phili/Desktop/eval_shocktube_norm/ref14.txt")  # rho, ux, p, T
failedJobs = []

for jobfolder in glob("/mnt/c/Users/phili/Desktop/sod/*nx*/"):
  print("Jobfolder: ", jobfolder)
  jobName = jobfolder.split("/")[-2]
  parameters = jobName.split("_")
  logfilename = glob(jobfolder + "*.out")
  if len(logfilename) == 0:
    print("::::ERROR: No logfile found, skipping this job")
    failedJobs.append(jobName)
    continue
  logfile = open(logfilename[0], "r")
  lines = logfile.read().splitlines()
  if len(lines) < 50:
    print("::::ERROR: Logfile seems too short, skipping this job")
    failedJobs.append(jobName)
    continue
  if len([line for line in lines if ":::::NATriuM run complete." in line]) == 0:
    print("::::ERROR: Run was not completed, skipping this job")
    failedJobs.append(jobName)
    continue

  cs = float([line for line in lines if "::::Sound speed:              " in line][0].removeprefix("::::Sound speed:              "))
  dt = float([line for line in lines if "::::Actual dt:                " in line][0].removeprefix("::::Actual dt:                ").removesuffix(" s"))
  dx = float([line for line in lines if "::::dx_min:                   " in line][0].removeprefix("::::dx_min:                   "))
  jobid = int(parameters[0])
  nx = int(parameters[1].removeprefix("nx"))
  cfl = float(parameters[2].removeprefix("cfl"))
  p = float(parameters[3].removeprefix("p"))
  vtkpath = jobfolder + "output/vtk/"
  lastIteration = str(max([int(pvtuname.removesuffix(".pvtu").split("t_0.")[-1]) for pvtuname in glob(vtkpath + "*.pvtu")]))
  print("::::Last iteration: ", lastIteration)
  pvtuFilePath = vtkpath + f"t_0.{lastIteration}.pvtu"
  # xiList = np.linspace(0, 1, nx*p).tolist()
  xiList = np.linspace(0, 1, nx).tolist()

  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(pvtuFilePath)
  reader.Update()
  data = reader.GetOutput()
  if data is None:
    print(f"::::ERROR: Could not read data from {pvtuFilePath}")
    continue

  points = data.GetPoints()
  point_rho = data.GetPointData().GetArray("rho")
  point_ux = data.GetPointData().GetArray("ux")
  point_T = data.GetPointData().GetArray("T")

  rhoList = []
  uxList = []
  TList = []

  locator = vtk.vtkPointLocator()
  locator.SetDataSet(data)
  locator.BuildLocator()

  for xi in xiList:
    point_id_list = vtk.vtkIdList()
    point_id = locator.FindClosestPoint(xi, yi, zi)
    rhoList.append(point_rho.GetValue(point_id))
    uxList.append(point_ux.GetValue(point_id))
    TList.append(point_T.GetValue(point_id))
  
  # for dataList, dataName in zip([rhoList, uxList, TList], ["rho", "ux", "T"]):
  for dataList, dataName in zip([rhoList], ["rho"]):
    fig, ax = plt.subplots(figsize=[7, 3.5])
    plt.plot(refOld[:,0],refOld[:,1],'b--',label='Reference')
    ax.scatter(xiList, dataList, marker=sllbm_marker, color=sllbm_color, s=sllbm_size, label="SLLBM")
    ax.set_title(f"Sod Shock Tube: {dataName} (nx={nx}, p={p}, cfl={cfl})")
    ax.legend()
    fig.savefig(imgpath + dataName + "_" + jobName + "_iT" + lastIteration + imgtype, transparent=transparent, dpi=300)
    # plt.show()

  plt.close("all")

  rhoRef = np.interp(xiList, refOld[:,0], refOld[:,1])
  L1rho = np.mean(np.abs(np.array(rhoList) - rhoRef))
  L2rho = np.sqrt(np.mean(np.pow(np.array(rhoList) - rhoRef, 2)))
  uxRef = np.interp(xiList, refOld[:,0], refOld[:,2])
  L1ux = np.mean(np.abs(np.array(uxList) - uxRef))
  L2ux = np.sqrt(np.mean(np.pow(np.array(uxList) - uxRef, 2)))
  TRef = np.interp(xiList, refOld[:,0], refOld[:,3])
  L1T = np.mean(np.abs(np.array(TList) - TRef))
  L2T = np.sqrt(np.mean(np.pow(np.array(TList) - TRef, 2)))
  LList.append([nx, cfl, L1rho, L2rho, L1ux, L2ux, L1T, L2T, p, dt, dx])

  print("")

fs = "failedJobs: "
for failedJob in failedJobs:
  ps += " " + failedJob

LListRefLevel = np.array(LList)
LListRefLevel = LListRefLevel[LListRefLevel[:,1] == 1]  # cfl
LListRefLevel = LListRefLevel[LListRefLevel[:,8] == 4]  # p

LListCfl = np.array(LList)
LListCfl = LListCfl[LListCfl[:,8] == 4]  # p

LListCflRef3 = LListCfl[LListCfl[:,0] == 200]  # nx
LListCflRef4 = LListCfl[LListCfl[:,0] == 400]  # nx
LListCflRef5 = LListCfl[LListCfl[:,0] == 800]  # nx
LListCflRef5 = LListCfl[LListCfl[:,0] == 1600]  # nx

# for i, dataName, dataLabel in zip([[3,4],[5,6],[7,8]], ["rho", "ux", "T"], ["Density", "Velocity", "Temperature"]):
for i, dataName, dataLabel in zip([[3,4]], ["rho"], ["Density"]):
# === REF LEVEL ===
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("dx")
  ax.set_title(dataLabel + " L1 Errors over dx")
  # ax.set_xticks(LListRefLevel[:,11])
  ax.scatter(LListRefLevel[:,11], LListRefLevel[:,i[0]], label="$L^1$ Error " + dataLabel)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.legend()
  fig.savefig(imgpath + dataName + "_L1_dx" + imgtype, transparent=transparent, dpi=300)

  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("dx")
  # ax.set_xticks(LListRefLevel[:,11])
  ax.set_title(dataLabel + " L2 Errors over dx")
  ax.scatter(LListRefLevel[:,11], LListRefLevel[:,i[1]], label="$L^2$ Error " + dataLabel)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.legend()
  fig.savefig(imgpath + dataName + "_L2_dx" + imgtype, transparent=transparent, dpi=300)

# === CFL ===
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("dt")
  # ax.set_xticks(LListCflRefi[:,10])
  ax.set_title(dataLabel + f" L1 Errors over dt (refinement level all)")
  for LListCflRefi, refLevel in zip([LListCflRef3, LListCflRef4, LListCflRef5], ["3", "4", "5"]):
    ax.scatter(LListCflRefi[:,10], LListCflRefi[:,i[0]], label="$L^1$ Error " + dataLabel + f" Ref{refLevel}")
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L1_dt_refall" + imgtype, transparent=transparent, dpi=300)

  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("dt")
  # ax.set_xticks(LListCflRefi[:,10])
  ax.set_title(dataLabel + f" L2 Errors over dt (refinement level all)")
  for LListCflRefi, refLevel in zip([LListCflRef3, LListCflRef4, LListCflRef5], ["3", "4", "5"]):
    ax.scatter(LListCflRefi[:,10], LListCflRefi[:,i[0]], label="$L^2$ Error " + dataLabel + f" Ref{refLevel}")
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L2_dt_refall" + imgtype, transparent=transparent, dpi=300)

  plt.close("all")