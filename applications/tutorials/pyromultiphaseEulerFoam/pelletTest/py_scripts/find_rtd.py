# trace generated using paraview version 5.8.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 
import sys
import copy
from sys import argv
from scipy import stats


solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['internalMesh']
solfoam.CellArrays = ['U.particles']
t = np.array(solfoam.TimestepValues)
N=t.size
print(N)

nsourcepts=200
pointSource1 = pv.PointSource()
pointSource1.Center = [0.01, 0.0, 0.02]
pointSource1.NumberOfPoints = nsourcepts
pointSource1.Radius = 0.01

pt1 = pv.ParticleTracer(Input=solfoam,
    SeedSource=pointSource1)
pt1.SelectInputVectors = ['POINTS', 'U.particles']

times=np.zeros(nsourcepts)
pidold=np.arange(nsourcepts)
print(pidold)
for i in range(N):
    print(t[i])
    pv.UpdatePipeline(time=t[i], proxy=pt1)
    idat    = dsa.WrapDataObject(pv.servermanager.Fetch(pt1) )
    pidnew=idat.PointData["ParticleId"]
    diff=np.setdiff1d(pidold,pidnew)
    for ii in range(len(diff)):
        times[diff[ii]]=t[i]
    pidold=copy.deepcopy(pidnew)


allpids=np.arange(nsourcepts)
print("mean residence time:",np.mean(times))
nbins=20
(count,bedge,binnum)=stats.binned_statistic(times,allpids,statistic='count',bins=nbins,range=(0.0,22.0))
bincenters=0.5*(bedge[0:-1]+bedge[1:])
np.savetxt("rtddist.dat",np.transpose(np.vstack((bincenters,count))),delimiter="  ")
