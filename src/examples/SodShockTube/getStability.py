import matplotlib.pyplot as plt
import numpy as np
from glob import glob

stable=[]
instable=[]
allRuns=[]

for jobFolder in glob("/mnt/c/Users/phili/Desktop/sodNu/*/"):
  print("Jobfolder: ", jobFolder)
  jobName = jobFolder.split("/")[-2]
  parameters = jobName.split("_")
  jobid = int(parameters[0])
  nxIn = int(parameters[1].removeprefix("nx"))
  cfl = float(parameters[2].removeprefix("cfl"))
  p = float(parameters[3].removeprefix("p"))
  nu = float(parameters[4].removeprefix("nu"))
  logfilename = glob(jobFolder + "*.out")
  if len(logfilename) == 0:
    print("::::ERROR: No logfile found, skipping this job")
    continue
  logfile = open(logfilename[0], "r")
  lines = logfile.read().splitlines()
  if len(lines) < 50:
    print("::::ERROR: Logfile seems too short, skipping this job")
    continue

  dt = float([line for line in lines if "::::Actual dt:                " in line][0].removeprefix("::::Actual dt:                ").removesuffix(" s"))
  # dxCells = float([line for line in lines if "::::dx_min:                   " in line][0].removeprefix("::::dx_min:                   "))
  dxLBM = 1/(nxIn*8)
  
  allRuns.append([dxLBM, cfl, dt, nu])
  if len([line for line in lines if ":::::NATriuM run complete." in line]) == 0:
    print("::::ERROR: Run was not completed, probably instability")
    instable.append([dxLBM, cfl, dt, nu, nxIn])
  else:
    stable.append([dxLBM, cfl, dt, nu, nxIn])

stable = np.array(stable)
instable = np.array(instable)
allRuns = np.array(allRuns)

allDx = np.unique(allRuns[:,0])
allCfl = np.unique(allRuns[:,1])
allDt = np.unique(allRuns[:,2])

fig, ax = plt.subplots()
for cfl in allCfl:
  for dx in allDx:
    data = stable[stable[:,0]==dx]
    data = data[data[:,1]==cfl]
    if len(data) > 0:
      print(f"nuMin(dx={dx:.2e},nxIn={data[0,4]:.0f},cfl={cfl})={min(data[:,3])}")

# nuMin(dx=3.91e-05,nxIn=3200,cfl=0.5)=1e-07
# nuMin(dx=7.81e-05,nxIn=1600,cfl=0.5)=3e-07
# nuMin(dx=1.56e-04,nxIn=800,cfl=0.5)=1e-06
# nuMin(dx=3.13e-04,nxIn=400,cfl=0.5)=1e-06
# nuMin(dx=6.25e-04,nxIn=200,cfl=0.5)=3e-06
# nuMin(dx=1.25e-03,nxIn=100,cfl=0.5)=3e-05
# nuMin(dx=2.50e-03,nxIn=50,cfl=0.5)=3e-05
# nuMin(dx=5.00e-03,nxIn=25,cfl=0.5)=3e-05
# nuMin(dx=1.04e-02,nxIn=12,cfl=0.5)=3e-05
# nuMin(dx=2.08e-02,nxIn=6,cfl=0.5)=3e-05
# nuMin(dx=4.17e-02,nxIn=3,cfl=0.5)=3e-05

# nuMin(dx=3.91e-05,nxIn=3200,cfl=1.0)=3e-07
# nuMin(dx=7.81e-05,nxIn=1600,cfl=1.0)=1e-06
# nuMin(dx=1.56e-04,nxIn=800,cfl=1.0)=1e-06
# nuMin(dx=3.13e-04,nxIn=400,cfl=1.0)=3e-06
# nuMin(dx=6.25e-04,nxIn=200,cfl=1.0)=1e-05
# nuMin(dx=1.25e-03,nxIn=100,cfl=1.0)=3e-05
# nuMin(dx=2.50e-03,nxIn=50,cfl=1.0)=3e-05
# nuMin(dx=5.00e-03,nxIn=25,cfl=1.0)=3e-05
# nuMin(dx=1.04e-02,nxIn=12,cfl=1.0)=3e-05
# nuMin(dx=2.08e-02,nxIn=6,cfl=1.0)=3e-05
# nuMin(dx=4.17e-02,nxIn=3,cfl=1.0)=3e-05

# nuMin(dx=3.91e-05,nxIn=3200,cfl=2.0)=3e-06
# nuMin(dx=7.81e-05,nxIn=1600,cfl=2.0)=1e-05
# nuMin(dx=1.56e-04,nxIn=800,cfl=2.0)=1e-05
# nuMin(dx=6.25e-04,nxIn=200,cfl=2.0)=3e-05
# nuMin(dx=1.25e-03,nxIn=100,cfl=2.0)=0.0001
# nuMin(dx=2.50e-03,nxIn=50,cfl=2.0)=0.0001
# nuMin(dx=5.00e-03,nxIn=25,cfl=2.0)=0.0003
# nuMin(dx=1.04e-02,nxIn=12,cfl=2.0)=0.0003
# nuMin(dx=2.08e-02,nxIn=6,cfl=2.0)=3e-05
# nuMin(dx=4.17e-02,nxIn=3,cfl=2.0)=3e-05