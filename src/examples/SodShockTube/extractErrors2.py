from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import sympy as sp
from sympy.abc import y
from getData import getData

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
onlyRho = False
# refSource = "txt"
refSource = "highRes"

imgpath = "/mnt/c/Users/phili/Desktop/sodImages/"
if refSource == "txt":
  imgpath += "RefOld/"
  ref = np.loadtxt("/mnt/c/Users/phili/Desktop/eval_shocktube_norm/ref14.txt")  # rho, ux, p, T
elif refSource == "highRes":
  imgpath += "RefHighRes/"
  dataLists, cs, dt, dx, jobid, cfl, p, jobName, nx, lastIteration = getData("/mnt/c/Users/phili/Desktop/sod/10908253_nx6400_cfl1_p4/")
  ref = np.array(dataLists).T
else:
  raise NotImplementedError("No valid reference source")
yi = 0
zi = 0

LList = []  # ref nu cfl L1rho L2rho L1ux L2ux L1T L2T
failedJobs = []

for jobFolder in glob("/mnt/c/Users/phili/Desktop/sod/*nx*/"):
  data = getData(jobFolder)
  if type(data) == str:
    failedJobs.append(data)
    continue
  else:
    dataLists, cs, dt, dx, jobid, cfl, p, jobName, nx, lastIteration = data
  xiList = dataLists[0]
  rhoList = dataLists[1]
  uxList = dataLists[2]
  TList = dataLists[4]
  
  if onlyRho:
    dataIs = [1]
    dataLists = [rhoList]
    dataNames = ["rho"]
  else:
    dataIs = [1,2,4]
    dataLists = [rhoList, uxList, TList]
    dataNames = ["rho", "ux", "T"]
  for dataI, dataList, dataName in zip(dataIs, dataLists, dataNames):
    fig, ax = plt.subplots(figsize=[7, 3.5])
    plt.plot(ref[:,0],ref[:,dataI],'b--',label='Reference')
    ax.scatter(xiList, dataList, marker=sllbm_marker, color=sllbm_color, s=sllbm_size, label="SLLBM")
    ax.set_title(f"Sod Shock Tube: {dataName} (nx={nx}, p={p}, cfl={cfl})")
    ax.legend()
    fig.savefig(imgpath + dataName + "_" + jobName + "_iT" + lastIteration + imgtype, transparent=transparent, dpi=300)
    # plt.show()

  plt.close("all")

  rhoRef = np.interp(xiList, ref[:,0], ref[:,1])
  L1rho = np.mean(np.abs(np.array(rhoList) - rhoRef))
  L2rho = np.sqrt(np.mean(np.pow(np.array(rhoList) - rhoRef, 2)))
  uxRef = np.interp(xiList, ref[:,0], ref[:,2])
  L1ux = np.mean(np.abs(np.array(uxList) - uxRef))
  L2ux = np.sqrt(np.mean(np.pow(np.array(uxList) - uxRef, 2)))
  TRef = np.interp(xiList, ref[:,0], ref[:,4])
  L1T = np.mean(np.abs(np.array(TList) - TRef))
  L2T = np.sqrt(np.mean(np.pow(np.array(TList) - TRef, 2)))
  LList.append([nx, cfl, L1rho, L2rho, L1ux, L2ux, L1T, L2T, p, dt, dx])

  print("")

fs = "failedJobs: "
for failedJob in failedJobs:
  fs += " " + failedJob

LListRef = np.array(LList)
LListRef = LListRef[LListRef[:,8] == 4]  # p
LListRefCfl0125 = LListRef[LListRef[:,1] == 0.125]  # p
LListRefCfl025 = LListRef[LListRef[:,1] == 0.25]  # p
LListRefCfl05 = LListRef[LListRef[:,1] == 0.5]  # p
LListRefCfl1 = LListRef[LListRef[:,1] == 1]  # p
LListRefCfl2 = LListRef[LListRef[:,1] == 2]  # p
LListRefCfl4 = LListRef[LListRef[:,1] == 4]  # p
LListRefCfl8 = LListRef[LListRef[:,1] == 8]  # p
LListRefCfls = [LListRefCfl8,LListRefCfl4,LListRefCfl2,LListRefCfl1,LListRefCfl05,LListRefCfl025,LListRefCfl0125]
# LListRef = LListRef[LListRef[:,1] == 1]  # cfl

LListCfl = np.array(LList)
LListCfl = LListCfl[LListCfl[:,8] == 4]  # p

LListCflRef0 = LListCfl[LListCfl[:,0] == 25]  # nx
LListCflRef1 = LListCfl[LListCfl[:,0] == 50]  # nx
LListCflRef2 = LListCfl[LListCfl[:,0] == 100]  # nx
LListCflRef3 = LListCfl[LListCfl[:,0] == 200]  # nx
LListCflRef4 = LListCfl[LListCfl[:,0] == 400]  # nx
LListCflRef5 = LListCfl[LListCfl[:,0] == 800]  # nx
LListCflRef6 = LListCfl[LListCfl[:,0] == 1600]  # nx
LListCflRef7 = LListCfl[LListCfl[:,0] == 3200]  # nx
LListCflRef8 = LListCfl[LListCfl[:,0] == 6400]  # nx
LListCflRefs = [LListCflRef0, LListCflRef1, LListCflRef2, LListCflRef3, LListCflRef4, LListCflRef5, LListCflRef6, LListCflRef7, LListCflRef8]

if onlyRho:
  dataIs = [[2,3]]
  dataNames = ["rho"]
  dataLabels = ["Density"]
else:
  dataIs = [[2,3],[4,5],[6,7]]
  dataNames = ["rho", "ux", "T"]
  dataLabels = ["Density", "Velocity", "Temperature"]

for i, dataName, dataLabel in zip(dataIs, dataNames, dataLabels):
# === REF LEVEL ===
  for j, Llevel in zip([0,1], [1,2]):
    for dataLists, overLabels, k, over in zip([LListRefCfls, [LListRefCfl8], [LListRefCfl4], [LListRefCfl2], [LListRefCfl1], [LListRefCfl05], [LListRefCfl025], [LListRefCfl0125], LListCflRefs, [LListCflRef0], [LListCflRef1], [LListCflRef2], [LListCflRef3], [LListCflRef4], [LListCflRef5], [LListCflRef6], [LListCflRef7], [LListCflRef8]],
                                              [[f" (CFL {cfl})" for cfl in [8,4,2,1,0.5,0.25,0.125]], [" (CFL 8)"], [" (CFL 4)"], [" (CFL 2)"], [" (CFL 1)"], [" (CFL 0.5)"], [" (CFL 0.25)"], [" (CFL 0.125)"], [f" Ref{refLevel}" for refLevel in range(len(LListCflRefs))], [" (Ref0)"], [" (Ref1)"], [" (Ref2)"], [" (Ref3)"], [" (Ref4)"], [" (Ref5)"], [" (Ref6)"], [" (Ref7)"], [" (Ref8)"]],
                                              [10,10,10,10,10,10,10,10,9,9,9,9,9,9,9,9,9,9],
                                              ["dx","dxCfl8","dxCfl4","dxCfl2","dxCfl1","dxCfl.5","dxCfl.25","dxCfl.125","dt","dtRef0","dtRef1","dtRef2","dtRef3","dtRef4","dtRef5","dtRef6","dtRef7","dtRef8"]):
      print("Doing", overLabels)
      fig, ax = plt.subplots(figsize=[7, 3.5])
      ax.set_xlabel(f"{over}")
      ax.set_title(dataLabel + f" L{Llevel} Errors over {over}")
      for dataList, overLabel in zip(dataLists, overLabels):
        ax.scatter(dataList[:,k], dataList[:,i[j]], label=f"$L^{Llevel}$ Error " + dataLabel + overLabel)
      if ( len(dataList) > 0 ):
        ax.set_xscale('log')
        ax.set_yscale('log')
      ax.legend()
      fig.savefig(imgpath + dataName + f"_L{Llevel}_{over}" + imgtype, transparent=transparent, dpi=300)
      
      plt.close("all")
