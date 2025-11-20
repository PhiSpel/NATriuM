import vtk
import numpy as np
from glob import glob

def getData(jobFolder):
  print("Jobfolder: ", jobFolder)
  jobName = jobFolder.split("/")[-2]
  parameters = jobName.split("_")
  logfilename = glob(jobFolder + "*.out")
  if len(logfilename) == 0:
    print("::::ERROR: No logfile found, skipping this job")
    return jobName
  logfile = open(logfilename[0], "r")
  lines = logfile.read().splitlines()
  if len(lines) < 50:
    print("::::ERROR: Logfile seems too short, skipping this job")
    return jobName
  if len([line for line in lines if ":::::NATriuM run complete." in line]) == 0:
    print("::::ERROR: Run was not completed, skipping this job")
    return jobName

  cs = float([line for line in lines if "::::Sound speed:              " in line][0].removeprefix("::::Sound speed:              "))
  dt = float([line for line in lines if "::::Actual dt:                " in line][0].removeprefix("::::Actual dt:                ").removesuffix(" s"))
  dx = float([line for line in lines if "::::dx_min:                   " in line][0].removeprefix("::::dx_min:                   "))
  jobid = int(parameters[0])
  nx = int(parameters[1].removeprefix("nx"))
  cfl = float(parameters[2].removeprefix("cfl"))
  p = float(parameters[3].removeprefix("p"))
  vtkpath = jobFolder + "output/vtk/"
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
    return jobName

  # points = data.GetPoints()
  point_rho = data.GetPointData().GetArray("rho")
  point_ux = data.GetPointData().GetArray("ux")
  point_uy = data.GetPointData().GetArray("uy")
  point_T = data.GetPointData().GetArray("T")
  if point_rho is None or point_ux is None or point_uy is None or point_T is None:
    print(f"::::ERROR: Could get point data from {pvtuFilePath}")
    return jobName

  rhoList = []
  uxList = []
  uyList = []
  TList = []

  locator = vtk.vtkPointLocator()
  locator.SetDataSet(data)
  locator.BuildLocator()

  for xi in xiList:
    # point_id_list = vtk.vtkIdList()
    point_id = locator.FindClosestPoint(xi, 1/nx/2, 0)
    rhoList.append(point_rho.GetValue(point_id))
    uxList.append(point_ux.GetValue(point_id))
    uyList.append(point_uy.GetValue(point_id))
    TList.append(point_T.GetValue(point_id))

  return [xiList, rhoList, uxList, uyList, TList], cs, dt, dx, jobid, cfl, p, jobName, nx, lastIteration