# trace generated using paraview version 5.8.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 
import sys
from sys import argv


solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['patch/solidsinlet']
solfoam.CellArrays = ['alpha.particles', 'U.particles', 'thermo:rho.particles','alpha.gas', 'U.gas', 'thermo:rho.gas']
t = np.array(solfoam.TimestepValues)
N=t.size

withsurfnormals = pv.GenerateSurfaceNormals(Input=solfoam)
# Properties modified on calculator1
calculator1 = pv.Calculator(Input=withsurfnormals)
calculator1.AttributeType = 'Point Data'
calculator1.ResultArrayName = 'mfluxb'
calculator1.Function = '"alpha.particles"*(-"U.particles_X"*"Normals_X"-"U.particles_Y"*"Normals_Y"-"U.particles_Z"*"Normals_Z")*"thermo:rho.particles"'

Pstd=101325.0
Tstd=273.15
R=287.0
std_dens=Pstd/Tstd/R
calculator2 = pv.Calculator(Input=calculator1)
calculator2.AttributeType = 'Point Data'
calculator2.ResultArrayName = 'mfluxg'
#in slm
calculator2.Function = '"alpha.gas"*(-"U.gas_X"*"Normals_X"-"U.gas_Y"*"Normals_Y"-"U.gas_Z"*"Normals_Z")*"thermo:rho.gas"/%e*1000.0*60.0'%(std_dens)

# create a new 'Integrate Variables'
int1 = pv.IntegrateVariables(Input=calculator2)

outfile=open("mflow_bm.dat","w")
for i in range(N):
    pv.UpdatePipeline(time=t[i], proxy=int1)
    idat    = dsa.WrapDataObject(pv.servermanager.Fetch(int1) )
    area       = idat.CellData['Area'].item()
    mfluxbint   = idat.PointData['mfluxb'].item()
    alphabint   = idat.PointData['alpha.particles'].item()
    mfluxgint   = idat.PointData['mfluxg'].item()
    alphagint   = idat.PointData['alpha.gas'].item()
    print("processing time = %e\t%e\t%e\t%e\t%e\t%e" % (t[i],area,mfluxbint,alphabint/area,mfluxgint,alphagint/area))
    outfile.write("%e\t%e\t%e\t%e\t%e\t%e\n"%(t[i],area,mfluxbint,alphabint/area,mfluxgint,alphagint/area))

outfile.close()
