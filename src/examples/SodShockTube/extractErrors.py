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

tmax = 0.15
gamma = 2#1.4
PfromR = False
recalcRhoRef = False
R = 1

rho1 = 8
T1 = 1.25
u1 = 0
if PfromR:
  p1 = rho1*R*T1
else:
  p1 = (gamma-1)*rho1*T1
a1 = np.sqrt(gamma*p1/rho1)

rho5 = 1
T5 = 1
u5 = 0
if PfromR:
  p5 = rho5*R*T5
else:
  p5 = (gamma-1)*rho5*T5

# p3 = p4
# u3 = u1 + 2*a1/(gamma-1)*(1-(p3/p1)**((gamma-1)/(2*gamma)))
# u4 = u5 + (p4-p5)*np.sqrt(2/(p5*(gamma+1)*(p4+(gamma-1)/(gamma+1)*p5)))

# solve for pmid by setting u3=u4
# u3u4 = sp.Eq(u1 + 2*a1/(gamma-1)*(1-(y/p1)**((gamma-1)/(2*gamma))),
#              u5 + (y-p5)*(2/(p5*(gamma+1)*(y+(gamma-1)/(gamma+1)*p5)))**0.5)
# solveset = sp.solve(u3u4, y, domain=sp.Interval(rho1, rho5), simplify=False)
# pMid = float(solveset[0])
# pMid = 0.943160123654987
pMid = 1.288
p3 = pMid
p4 = pMid
u3 = u1 + 2*a1/(gamma-1)*(1-(pMid/p1)**((gamma-1)/(2*gamma)))
# u3 = 0.613
u4 = u3

rho3 = rho1*(p3/p1)**(1/gamma)
rho4 = rho5*(((gamma+1)*p4+(gamma-1)*p5)/((gamma-1)*p4+(gamma+1)*p5))

a3 = np.sqrt(gamma*p3/rho3)
a5 = np.sqrt(gamma*p5/rho5)

LList = []  # ref nu cfl L1rho L2rho L1ux L2ux L1T L2T
refOld = np.loadtxt("/mnt/c/Users/phili/Desktop/eval_shocktube_norm/ref14.txt")  # rho, ux, p, T

for jobfolder in glob("/mnt/c/Users/phili/Desktop/sod/*/"):
  print("Jobfolder: ", jobfolder)
  jobName = jobfolder.split("/")[-2]
  parameters = jobName.split("_")
  logfilename = glob(jobfolder + "*.out")
  if len(logfilename) == 0:
    print("::::Error: No logfile found, skipping this job")
    continue
  logfile = open(logfilename[0], "r")
  lines = logfile.read().splitlines()
  if len(lines) < 50:
    print("::::Error: Logfile seems too short, skipping this job")
    continue
  tmax = float([line for line in lines if "Simulation end time will be t_max = " in line][0].removeprefix("Simulation end time will be t_max = "))
  # tmax = 0.15

  x12 = 0.5+(u1-a1)*tmax
  x23 = 0.5+(u3-a3)*tmax
  x34 = 0.5+u3*tmax
  a5hat = a5*np.sqrt((gamma+1)/(2*gamma)*p4/p5 + (gamma-1)/(2*gamma))
  x45 = 0.5 + (u5 + a5hat)*tmax
  print(f"::::Positions: {x12:.3f}, {x23:.3f}, {x34:.3f}, {x45:.3f}")

  cs = float([line for line in lines if "::::Sound speed:              " in line][0].removeprefix("::::Sound speed:              "))
  jobid = int(parameters[0])
  refLevel = int(parameters[1].removeprefix("ref"))
  viscosity = float(parameters[2].removeprefix("nu"))
  cfl = float(parameters[3].removeprefix("cfl"))
  if len(parameters) > 4:
    p = float(parameters[4].removeprefix("p"))
  else:
    p = 4
  vtkpath = jobfolder + "output/vtk/"
  lastIteration = max([int(pvtuname.removesuffix(".pvtu").split("t_0.")[-1]) for pvtuname in glob(vtkpath + "*.pvtu")])
  print("::::Last iteration: ", lastIteration)
  pvtuFilePath = vtkpath + f"t_0.{lastIteration}.pvtu"
  xiList = np.linspace(0, 1, 100*pow(2, refLevel)).tolist()

  x = np.array(xiList)
  x2 = np.logical_and(x>=x12, x<x23)
  x3 = np.logical_and(x>=x23, x<x34)
  x4 = np.logical_and(x>=x34, x<x45)

  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(pvtuFilePath)
  reader.Update()
  data = reader.GetOutput()
  if data is None:
    print(f"::::Error: Could not read data from {pvtuFilePath}")
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
    if dataName == "rho":
      ref = np.array(xiList)*0
      chi = (x-0.5)/tmax
      uxi = 2/(gamma+1)*(a1 + chi)
      axi = a1 - (gamma-1)/2*uxi
      ref[x<x12] = rho1
      ref[x2] = (rho1*(axi/a1)**(2/(gamma-1)))[x2]
      ref[x3] = rho3
      ref[x4] = rho4
      ref[x>=x45] = rho5
      ax.plot(x, ref, 'k--', label="Reference", linewidth=1)

      plt.plot(refOld[:,0],refOld[:,1],'b--',label='Reference Old')

    ax.scatter(xiList, dataList, marker=sllbm_marker, color=sllbm_color, s=sllbm_size, label="SLLBM")
    ax.set_title(f"Sod Shock Tube: {dataName} (ref={refLevel}, nu={viscosity}, cfl={cfl})")
    ax.legend()
    fig.savefig(imgpath + dataName + "_" + jobName + imgtype, transparent=transparent, dpi=300)
    # plt.show()

  plt.close("all")

  if recalcRhoRef:
    rhoRef = ref
  else:
    rhoRef = np.interp(xiList, refOld[:,0], refOld[:,1])
  L1rho = np.mean(np.abs(np.array(rhoList) - rhoRef))
  L2rho = np.sqrt(np.mean(np.pow(np.array(rhoList) - rhoRef, 2)))
  uxRef = np.interp(xiList, refOld[:,0], refOld[:,2])
  L1ux = np.mean(np.abs(np.array(uxList) - uxRef))
  L2ux = np.sqrt(np.mean(np.pow(np.array(uxList) - uxRef, 2)))
  TRef = np.interp(xiList, refOld[:,0], refOld[:,3])
  L1T = np.mean(np.abs(np.array(TList) - TRef))
  L2T = np.sqrt(np.mean(np.pow(np.array(TList) - TRef, 2)))
  LList.append([refLevel, viscosity, cfl, L1rho, L2rho, L1ux, L2ux, L1T, L2T, p])

  print("")

LListRefLevel = np.array(LList)
LListRefLevel = LListRefLevel[LListRefLevel[:,1] == 1e-3]  # viscosity
LListRefLevel = LListRefLevel[LListRefLevel[:,2] == 1]  # cfl
LListRefLevel = LListRefLevel[LListRefLevel[:,9] == 4]  # p

LListRefLevel2 = np.array(LList)
LListRefLevel2 = LListRefLevel2[LListRefLevel2[:,1] == 1e-5]  # viscosity
LListRefLevel2 = LListRefLevel2[LListRefLevel2[:,2] == 1]  # cfl
LListRefLevel2 = LListRefLevel2[LListRefLevel2[:,9] == 4]  # p

LListP = np.array(LList)
Pviscosity = 3e-4
LListP = LListP[LListP[:,0] == 0]  # refLevel
LListP = LListP[LListP[:,1] == Pviscosity]  # viscosity
LListP = LListP[LListP[:,2] == 1]  # cfl

LListCfl = np.array(LList)
Cflviscosity = 3e-4
LListCfl = LListCfl[LListCfl[:,0] == 3]  # refLevel
LListCfl = LListCfl[LListCfl[:,1] == Cflviscosity]  # viscosity
LListCfl = LListCfl[LListCfl[:,9] == 4]  # p

LListNu = np.array(LList)
LListNu = LListNu[LListNu[:,2] == 1]  # cfl
LListNu = LListNu[LListNu[:,9] == 4]  # p

# for i, dataName, dataLabel in zip([[3,4],[5,6],[7,8]], ["rho", "ux", "T"], ["Density", "Velocity", "Temperature"]):
for i, dataName, dataLabel in zip([[3,4]], ["rho"], ["Density"]):
# === REF LEVEL ===
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("Refinement Level")
  ax.set_xticks(LListRefLevel[:,0])
  ax.set_title(dataLabel + " L1 Errors over Refinement Level (nu=1e-3, cfl=1)")
  ax.scatter(LListRefLevel[:,0], LListRefLevel[:,i[0]], label="$L^1$ Error " + dataLabel)
  # for d in LListRefLevel:
  #   refLeveli = int(d[0])
  #   L1i = d[i[0]]
  #   ax.annotate(f"({refLeveli},{L1i:.2e})", (refLeveli, L1i), xytext=(refLeveli+.1, L1i*1.2), arrowprops=dict(arrowstyle="->"))
  ax.legend()
  fig.savefig(imgpath + dataName + "_L1_nu1e-3_cfl1" + imgtype, transparent=transparent, dpi=300)

  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("Refinement Level")
  ax.set_xticks(LListRefLevel[:,0])
  ax.set_title(dataLabel + " L2 Errors over Refinement Level (nu=1e-3, cfl=1)")
  ax.scatter(LListRefLevel[:,0], LListRefLevel[:,i[1]], label="$L^2$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + "_L2_nu1e-3_cfl1" + imgtype, transparent=transparent, dpi=300)

# === REF LEVEL 2 ===
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("Refinement Level")
  ax.set_xticks(LListRefLevel2[:,0])
  ax.set_title(dataLabel + " L1 Errors over Refinement Level (nu=1e-5, cfl=1)")
  ax.scatter(LListRefLevel2[:,0], LListRefLevel2[:,i[0]], label="$L^1$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + "_L1_nu1e-5_cfl1" + imgtype, transparent=transparent, dpi=300)

  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("Refinement Level")
  ax.set_xticks(LListRefLevel2[:,0])
  ax.set_title(dataLabel + " L2 Errors over Refinement Level (nu=1e-5, cfl=1)")
  ax.scatter(LListRefLevel2[:,0], LListRefLevel2[:,i[1]], label="$L^2$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + "_L2_nu1e-5_cfl1" + imgtype, transparent=transparent, dpi=300)

# === FE ORDER ===
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("FE Order")
  ax.set_xticks(LListP[:,9])
  ax.set_title(dataLabel + f"L1 Errors over FE order (cfl=1, nu={Pviscosity:.2e})")
  ax.scatter(LListP[:,9], LListP[:,i[0]], label="$L^1$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L1Perrors_nu{Pviscosity:.2e}_cfl1" + imgtype, transparent=transparent, dpi=300)

  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("FE Order")
  ax.set_xticks(LListP[:,9])
  ax.set_title(dataLabel + f"L2 Errors over FE order (cfl=1, nu={Pviscosity:.2e})")
  ax.scatter(LListP[:,9], LListP[:,i[1]], label="$L^2$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L2Perrors_nu{Pviscosity:.2e}_cfl1" + imgtype, transparent=transparent, dpi=300)

# === CFL ===
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("1/CFL")
  ax.set_xticks(1/LListCfl[:,2])
  ax.set_title(dataLabel + f"L1 Errors over 1/cfl (refLevel=3, nu={Cflviscosity:.2e}, p=4)")
  ax.scatter(1/LListCfl[:,2], LListCfl[:,i[0]], label="$L^1$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L1Cflerrors_nu{Cflviscosity:.2e}_cfl1" + imgtype, transparent=transparent, dpi=300)

  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("1/CFL")
  ax.set_xticks(1/LListCfl[:,2])
  ax.set_title(dataLabel + f"L2 Errors over 1/Cfl (refLevel=3, nu={Cflviscosity:.2e}, p=4)")
  ax.scatter(1/LListCfl[:,2], LListCfl[:,i[1]], label="$L^2$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L2Cflerrors_nu{Cflviscosity:.2e}_cfl1" + imgtype, transparent=transparent, dpi=300)

# === NU ===
  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("Nu")
  ax.set_xticks(LListNu[:,1])
  ax.set_title(dataLabel + f"L1 Errors over nu (p=4, cfl=1)")
  ax.scatter(LListNu[:,1], LListNu[:,i[0]], label="$L^1$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L1Nuerrors_cfl1" + imgtype, transparent=transparent, dpi=300)

  fig, ax = plt.subplots(figsize=[7, 3.5])
  ax.set_xlabel("Nu")
  ax.set_xticks(LListNu[:,1])
  ax.set_title(dataLabel + f"L2 Errors over nu (p=4, cfl=1)")
  ax.scatter(LListNu[:,1], LListNu[:,i[1]], label="$L^2$ Error " + dataLabel)
  ax.legend()
  fig.savefig(imgpath + dataName + f"_L2Nuerrors_cfl1" + imgtype, transparent=transparent, dpi=300)
