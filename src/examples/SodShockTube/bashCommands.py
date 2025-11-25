import matplotlib.pyplot as plt
import numpy as np

def getCommand(nx, cfl, nu, refLevel=1, p=4):
  gPx=nx*p*pow(2, refLevel)
  gridPoints=gPx*p*pow(2, refLevel)  # 4 support points in each direction; refine once in each direction
  dx = 1/gridPoints
  dt0=5.52427e-05  # gPx=800, cfl=1, p=4
  dt=dt0*cfl*(800/gPx)
  tmax=0.259808
  iTmax=int(tmax/dt)
  iTvtk=int(iTmax/2)
  expectedTime=int(gridPoints*iTmax*1e-5/60)+5
  if nx > 500:
    expectedTime=int(gridPoints*iTmax*3e-6/60)+5
  if nx > 1000:
    expectedTime=int(gridPoints*iTmax*1e-6/60)+5
  if nx > 3000:
    expectedTime=int(gridPoints*iTmax*1e-7/60)+5
  nprocs = min(64, int(gridPoints/100))
  print(f"sbatch --time={expectedTime} --job-name=SodNx{nx}Cfl{cfl} --ntasks-per-node={nprocs} sodShockNx.sh {nx} {cfl} {iTvtk} {nu}")
  

# for nx in [3, 6, 12]:
#   for cfl in [2, 1, 0.5]:
#     for nu in [1e-2, 3e-3, 1e-1, 3e-4, 1e-4, 3e-5]:
#       getCommand(nx, cfl, nu)
# for nx in [200, 400, 800, 1600, 3200]:
#   for cfl in [2, 1, 0.5]:
#     for nu in [1e-5, 3e-6, 1e-6, 3e-7, 1e-7]:
#       getCommand(nx, cfl, nu)
for nx in [1600, 3200, 6400]:
  for cfl in [2, 1, 0.5]:
    for nu in [1e-7, 3e-8, 1e-8]:
      getCommand(nx, cfl, nu)
