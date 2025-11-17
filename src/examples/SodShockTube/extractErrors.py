import os
from glob import glob
import vtk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font', **{'size': 11, 'family': 'sans-serif', 'sans-serif': ['Myriad Pro', 'Arial', 'Tahoma']})
plt.rcParams['text.usetex'] = True
plt.rcParams["savefig.dpi"] = 600
plt.rcParams['markers.fillstyle'] = 'none'
plt.rcParams['figure.constrained_layout.use'] = True

sllbm_marker = 'x'
sllbm_color = 'red'
sllbm_size = 1.5*plt.rcParams['lines.markersize']

imgpath = "/mnt/c/Users/phili/Desktop/sodImages/"
yi = 1/25/2
zi = 0

tmax = 0.15
gamma = 1.4

rho1 = 8
T1 = 1.25
u1 = 0
p1 = (gamma-1)*rho1*T1
a1 = np.sqrt(gamma*p1/rho1)

rho5 = 1
T5 = 1
u5 = 0
p5 = (gamma-1)*rho5*T5

# p3 = p4
# u3 = u1 + 2*a1/(gamma-1)*(1-(p3/p1)**((gamma-1)/(2*gamma)))
# u4 = u5 + (p4-p5)*np.sqrt(2/(p5*(gamma+1)*(p4+(gamma-1)/(gamma+1)*p5)))
# solve for pmid by setting u3=u4
# u1 + 2*a1/(gamma-1)*(1-(pMid/p1)**((gamma-1)/(2*gamma))) = u5 + (pMid-p5)*np.sqrt(2/(p5*(gamma+1)*(pMid+(gamma-1)/(gamma+1)*p5)))
p3 = 1.288
p4 = p3
u3 = 0.613
u4 = u3

rho3 = rho1*(p3/p1)**(1/gamma)
rho4 = rho5*(((gamma+1)*p4+(gamma-1)*p5)/((gamma-1)*p4+(gamma+1)*p5))

a3 = np.sqrt(gamma*p3/rho3)
a5 = np.sqrt(gamma*p5/rho5)

for jobfolder in glob("/mnt/c/Users/phili/Desktop/sod/10899319*/"):
  print("Jobfolder: ", jobfolder)
  jobName = jobfolder.split("/")[-2]
  print("::::Jobname: ", jobName)
  parameters = jobName.split("_")
  print("::::Parameters: ", parameters)
  logfile = open(jobfolder + "slurm_NATriuM_SodShock.out", "r")
  lines = logfile.read().splitlines()
  tmax = float([line for line in lines if "Simulation end time will be t_max = " in line][0].removeprefix("Simulation end time will be t_max = "))
  # tmax = 0.15

  x12 = 0.5+(u1-a1)*tmax
  x23 = 0.5+(u3-a3)*tmax
  x34 = 0.5+u3*tmax
  a5hat = a5*np.sqrt((gamma+1)/(2*gamma)*p4/p5 + (gamma-1)/(2*gamma))
  x45 = 0.5 + (u5 + a5hat)*tmax

  cs = float([line for line in lines if "::::Sound speed:              " in line][0].removeprefix("::::Sound speed:              "))
  jobid = int(parameters[0])
  refLevel = int(parameters[1].removeprefix("ref"))
  viscosity = float(parameters[2].removeprefix("nu"))
  cfl = float(parameters[3].removeprefix("cfl"))
  vtkpath = jobfolder + "output/vtk/"
  lastIteration = max([int(pvtuname.removesuffix(".pvtu").split("t_0.")[-1]) for pvtuname in glob(vtkpath + "*.pvtu")])
  print("::::Last iteration: ", lastIteration)
  pvtuFilePath = vtkpath + f"t_0.{lastIteration}.pvtu"
  xiList = np.linspace(0, 1, 100*pow(2, refLevel)).tolist()

  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(pvtuFilePath)
  reader.Update()
  data = reader.GetOutput()
  if data is None:
    print(f"::::Error: Could not read data from {pvtuFilePath}")
    exit()

  print("::::Reader set up successfully")

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

  print(f"::::Data extracted, now plotting into {imgpath}")
  
  for dataList, dataName in zip([rhoList, uxList, TList], ["rho", "ux", "T"]):
    fig, ax = plt.subplots(figsize=[7, 3.5])
    if dataName == "rho":
      x = np.array(xiList)
      ref = np.array(xiList)*0
      chi = (x-0.5)/tmax
      uxi = 2/(gamma+1)*(a1 + chi)
      axi = a1 - (gamma-1)/2*uxi
      x2 = np.logical_and(x>=x12, x<x23)
      x3 = np.logical_and(x>=x23, x<x34)
      x4 = np.logical_and(x>=x34, x<x45)
      ref[x<x12] = rho1
      ref[x2] = (rho1*(axi/a1)**(2/(gamma-1)))[x2]
      ref[x3] = rho3
      ref[x4] = rho4
      ref[x>=x45] = rho5
      ax.plot(x, ref, 'k--', label="Reference", linewidth=1)
    ax.scatter(xiList, dataList, marker=sllbm_marker, color=sllbm_color, s=sllbm_size, label="SLLBM")
    ax.set_title(f"Sod Shock Tube: {dataName} (ref={refLevel}, nu={viscosity}, cfl={cfl})")
    ax.legend()
    fig.savefig(imgpath + dataName + "_" + jobName + ".pdf", transparent=True, dpi=300)
    # plt.show()

  plt.close("all")
  print("::::Finished job ", jobid)
  print("\n")