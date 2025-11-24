import vtk
import numpy as np
import os

def getPyFRData(vtuFilePath):
  if not os.path.exists(vtuFilePath):
    print(f"::::ERROR: File not found at {vtuFilePath}")
    exit()
  nx = 3200
  xiList = np.linspace(0, 1, nx).tolist()

  reader = vtk.vtkXMLUnstructuredGridReader()
  reader.SetFileName(vtuFilePath)
  reader.Update()
  data = reader.GetOutput()
  if data is None:
    print(f"::::ERROR: Could not read data from {vtuFilePath}")
    exit()

  # points = data.GetPoints()
  point_rho = data.GetPointData().GetArray("Density")
  point_u = data.GetPointData().GetArray("Velocity")
  point_p = data.GetPointData().GetArray("Pressure")

  rhoList = []
  uxList = []
  pList = []
  TList = []

  locator = vtk.vtkPointLocator()
  locator.SetDataSet(data)
  locator.BuildLocator()

  for xi in xiList:
    # point_id_list = vtk.vtkIdList()
    point_id = locator.FindClosestPoint(xi, 1/nx/2, 0)
    rhoList.append(point_rho.GetValue(point_id))
    uxList.append(point_u.GetTuple(point_id)[0])
    pList.append(point_p.GetValue(point_id))

  TList = [p/rho for rho, p in zip(rhoList, pList)]

  return np.array([xiList, rhoList, uxList, pList, TList]).T