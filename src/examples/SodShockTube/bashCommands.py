def getCommand(nx, cfl, refLevel=1, p=4):
  gPx=nx*p*pow(2, refLevel)
  gridPoints=gPx*p*pow(2, refLevel)  # 4 support points in each direction; refine once in each direction
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
  nprocs = min(64, int(gridPoints/100))
  print(f"sbatch --time={expectedTime} --job-name=SodNx{nx}Cfl{cfl} --ntasks-per-node={nprocs} sodShockNx.sh {nx} {cfl} {iTvtk}")


# for nx in [25, 50, 100, 200, 400, 800]:
#   for cfl in [8, 4, 2, 1, 0.5, 0.25, 0.125]:
#     getCommand(nx, cfl)
# for cfl in [8, 4, 2, 1]:
#   getCommand(1600, cfl)
# for cfl in [8, 4, 2]:
#   getCommand(3200, cfl)
# for cfl in [8]:
#   getCommand(6400, cfl)

  
for nx in [25, 50]:
  for cfl in [8, 4, 2, 1, 0.5, 0.25, 0.125]:
    getCommand(nx, cfl)
for cfl in [2, 1]:
  getCommand(3200, cfl)
for cfl in [8, 4, 2, 1]:
  getCommand(6400, cfl)