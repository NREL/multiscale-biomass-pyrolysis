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
solfoam.MeshRegions = ['patch/outlet']
solfoam.CellArrays = ['alpha.gas', 'U.gas' ,'SGC.gas', 'TRC.gas', 'SGH.gas', 'SGL.gas', 'TRH.gas', 
'TRL.gas', 'T.gas','thermo:rho.gas','alpha.particles', 'U.particles', 'CH.particles','thermo:rho.particles']
t = np.array(solfoam.TimestepValues)
N=t.size
print(N)
# Calculate SG species
calculator1 = pv.Calculator(Input=solfoam)
calculator1.AttributeType = 'Point Data'
calculator1.ResultArrayName = 'sgcmflux'
calculator1.Function = '"alpha.gas"*"SGC.gas"*"U.gas_Z"*"thermo:rho.gas"'

calculator2 = pv.Calculator(Input=calculator1)
calculator2.AttributeType = 'Point Data'
calculator2.ResultArrayName = 'sghmflux'
calculator2.Function = '"alpha.gas"*"SGH.gas"*"U.gas_Z"*"thermo:rho.gas"'

calculator3 = pv.Calculator(Input=calculator2)
calculator3.AttributeType = 'Point Data'
calculator3.ResultArrayName = 'sglmflux'
calculator3.Function = '"alpha.gas"*"SGL.gas"*"U.gas_Z"*"thermo:rho.gas"'

# Calculate TR species
calculator4 = pv.Calculator(Input=calculator3)
calculator4.AttributeType = 'Point Data'
calculator4.ResultArrayName = 'trcmflux'
calculator4.Function = '"alpha.gas"*"TRC.gas"*"U.gas_Z"*"thermo:rho.gas"'

calculator5 = pv.Calculator(Input=calculator4)
calculator5.AttributeType = 'Point Data'
calculator5.ResultArrayName = 'trhmflux'
calculator5.Function = '"alpha.gas"*"TRH.gas"*"U.gas_Z"*"thermo:rho.gas"'

calculator6 = pv.Calculator(Input=calculator5)
calculator6.AttributeType = 'Point Data'
calculator6.ResultArrayName = 'trlmflux'
calculator6.Function = '"alpha.gas"*"TRL.gas"*"U.gas_Z"*"thermo:rho.gas"'

calculator7 = pv.Calculator(Input=calculator6)
calculator7.AttributeType = 'Point Data'
calculator7.ResultArrayName = 'charflux'
calculator7.Function = '"alpha.particles"*"CH.particles"*"U.particles_Z"*"thermo:rho.particles"'

# create a new 'Integrate Variables'
int1 = pv.IntegrateVariables(Input=calculator7)

outfile=open("mflow_all.dat","w")
for i in range(N)[1::5]:
    pv.UpdatePipeline(time=t[i], proxy=int1)
    idat    = dsa.WrapDataObject(pv.servermanager.Fetch(int1) )
    #print('keys',idat.CellData.keys())
    #print('keys',idat.PointData.keys())
    area       = idat.CellData['Area'].item()
    sgcmfluxint   = idat.PointData['sgcmflux'].item()
    sghmfluxint   = idat.PointData['sghmflux'].item()
    sglmfluxint   = idat.PointData['sglmflux'].item()
    trcmfluxint   = idat.PointData['trcmflux'].item()
    trhmfluxint   = idat.PointData['trhmflux'].item()
    trlmfluxint   = idat.PointData['trlmflux'].item()
    alphaint   = idat.PointData['alpha.gas'].item()
    mflux_sg_tot = sgcmfluxint + sghmfluxint + sglmfluxint
    mflux_tr_tot = trcmfluxint + trhmfluxint + trlmfluxint
    chmfluxint = idat.PointData['charflux'].item()
    #print(f'SGC: {sgcmfluxint}, SGH: {sghmfluxint} SGL: {sglmfluxint} total: {mflux_sg_tot}\n')
    #print(f'TRC: {trcmfluxint}, TRH: {trhmfluxint} TRL: {trlmfluxint} total: {mflux_tr_tot}\n')
    print("processing time =",t[i],mflux_sg_tot+mflux_tr_tot+chmfluxint)
    outfile.write("%e\t%e\t%e\t%e\n"%(t[i],mflux_sg_tot,mflux_tr_tot,chmfluxint))
outfile.close()
