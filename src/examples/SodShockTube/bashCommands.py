# for nx in [25, 50, 100, 200, 400, 800]:
for nx in [1600, 3200, 6400]:
  # for cfl in [8, 4, 2, 1, 0.5, 0.25, 0.125]:
  for cfl in [8, 4, 2, 1]:
    gridPoints=nx*4*4*2*2  # 4 support points in each direction; refine once in each direction
    dt0=0.000220971
    dt=dt0*cfl
    tmax=0.259808
    iTmax=int(tmax/dt)
    iTvtk=int(iTmax/2)
    expectedTime=int(gridPoints*iTmax*2e-5/60)+5
    print(f"sbatch --time={expectedTime} --job-name SodNx{nx}Cfl{cfl} sodShockNx.sh {nx} {cfl} {iTvtk}")