#### import the simple module from the paraview
from paraview import simple as pv
import sys
from sys import argv
import vtk.numpy_interface.dataset_adapter as dsa
import numpy as np 

solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['internalMesh']
solfoam.CellArrays = ['p']
t = np.array(solfoam.TimestepValues)
N=t.size
print(N)

pv.UpdatePipeline(time = t[1])

# if needed, evaluate the point locations used in the simulation
ofvtkdata = pv.servermanager.Fetch(solfoam)
ofdata = dsa.WrapDataObject( ofvtkdata)
ofpts = np.array(ofdata.Points.Arrays[0])
ptsmin = ofpts.min(axis=0) # minimum values of the three axes
ptsmax = ofpts.max(axis=0) # maximum values of the three axes
print(ptsmin)
print(ptsmax)

probeLoc1 = pv.ProbeLocation(Input=solfoam,ProbeType='Fixed Radius Point Source')
probeLoc1.ProbeType.Center = [0.0,0.0,ptsmin[2]]
probeLoc2 = pv.ProbeLocation(Input=solfoam,ProbeType='Fixed Radius Point Source')
probeLoc2.ProbeType.Center = [0.0,0.0,ptsmax[2]]

pd=np.zeros(N)
time=np.zeros(N)
for i in range(N):
    pv.UpdatePipeline(time=t[i], proxy=probeLoc1)
    pv.UpdatePipeline(time=t[i], proxy=probeLoc2)
    idat1 = dsa.WrapDataObject(pv.servermanager.Fetch(probeLoc1))
    idat2 = dsa.WrapDataObject(pv.servermanager.Fetch(probeLoc2))
    pd[i]=idat1.PointData['p'].item()-idat2.PointData['p'].item()
    print(t[i],pd[i])

outfile=open("pd.dat","w")
for i in range(N):
	outfile.write("%e\t%e\n"%(t[i],pd[i]))

outfile.close()
