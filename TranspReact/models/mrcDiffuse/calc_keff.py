import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
import scipy
    
#fn=argv[1]
fn_pattern='plt*'
fn_list = sorted(glob.glob(fn_pattern), key=lambda f: int(f.split("plt")[1]))
flowdir=int(argv[1])
flowdir_char=chr(ord('x')+flowdir)
outfile='keff.csv'

with open(outfile, "w+") as f:
    f.write('pltfile,outletT,inletT,dTdX_3d,dTdx,avgJ,k_cell_avg,keff\n')
f.close()

for fn in fn_list:
    ds=yt.load(fn)
    # The gradient operator requires periodic boundaries.  This dataset has
    # open boundary conditions.
    ds.force_periodicity()
    add_concgrad = ds.add_gradient_fields(("boxlib","Temperature"))
    ncells=ds.domain_dimensions
    prob_lo=ds.domain_left_edge.d
    prob_hi=ds.domain_right_edge.d
    outlet=0.5*(prob_lo+prob_hi)
    outlet[flowdir]=0.01*prob_hi[flowdir]
    print('outlet:',outlet)
    size=prob_hi-prob_lo

    # load fields at outlet
    slc = yt.SlicePlot(ds, flowdir_char, ["solid", "Temperature"],center=[outlet[0],outlet[1],outlet[2]])
    res=np.array([ncells[(flowdir+1)%3],ncells[(flowdir+2)%3]])
    frb = slc.data_source.to_frb(((size[(flowdir+1)%3],'cm'),(size[(flowdir+2)%3],'cm')),res) 
    vfrac=np.array(frb["solid"])
    print(np.mean(vfrac))
    temp=np.array(frb["Temperature"])
    # calculate avg temperature at outlet
    T_1=np.mean(temp) # mol/m^3
    L=size[flowdir] # m
    
    # Save figures
    slc.set_log("solid",False)
    slc.set_log("Temperature",False)
    slc.save("slcdata_outlet"+"%4.4d"%(0))
    
    # Calculate temperature gradient in flow direction over center slice
    center_slice=0.5*(prob_lo+prob_hi)
    slc = yt.SlicePlot(ds, flowdir_char, ["solid", "Temperature"],center=[center_slice[0],center_slice[1],center_slice[2]])
    res=np.array([ncells[(flowdir+1)%3],ncells[(flowdir+2)%3]])
    frb = slc.data_source.to_frb(((size[(flowdir+1)%3],'cm'),(size[(flowdir+2)%3],'cm')),res) 
    vfrac=np.array(frb["solid"])
    print('center slice:',center_slice)
    dTdx=np.array(frb["Temperature_gradient_"+flowdir_char])
    area=size[(flowdir+1)%3]*size[(flowdir+2)%3] # m^2

    # calculate the thermal conductivity of each cell in slice by averaging 
    solid_dif=1.106e-5
    fluid_dif=22.39e-6
    solid_rho=3890 # kg/m^3
    fluid_rho=1.293 #kg/m^3
    # integrated over whole volume will be more accurate 3D
    k_cell=(vfrac)*solid_dif*solid_rho + (1-vfrac)*fluid_dif*fluid_rho
    alpha_cell=(vfrac)*solid_dif + (1-vfrac)*fluid_dif
    rho_eff = vfrac*solid_rho+(1-vfrac)*fluid_rho
    print('rho_eff:',np.mean(rho_eff))
    print('alpha_cell:',np.mean(alpha_cell))
    #D_avg=np.mean(D)
    print('res',res,'/n')
    print('dx',size[(flowdir+1)%3]/res[0])
    print('dy',size[(flowdir+2)%3]/res[1])
    dx=size[(flowdir+1)%3]/res[0]
    dy=size[(flowdir+2)%3]/res[0]
    # integrate over slice
    J=-np.sum(k_cell*dTdx*dx*dy)
    #J_avg=np.mean(J)
    print('J:',J) 
    print('dT/dz:',np.sum(dTdx*dx*dy))
    print('dT/dz avg:',np.sum(dTdx*dx*dy)/area)
    # Save figures
    slc.set_log("solid",False)
    slc.set_log("Temperature",False)
    slc.save("slcdata_center"+"%4.4d"%(0))



    # calculate avg concentration at inlet and dC
    inlet=0.5*(prob_lo+prob_hi)
    inlet[flowdir]=0.999*prob_hi[flowdir]
    print('inlet:',inlet)
    slc = yt.SlicePlot(ds, flowdir_char, ["solid","Temperature"],center=[inlet[0],inlet[1],inlet[2]])
    frb = slc.data_source.to_frb(((size[(flowdir+1)%3],'cm'),(size[(flowdir+2)%3],'cm')),res) 
    temp=np.array(frb["Temperature"]) #mol/m^3

    T_2=np.mean(temp)
    dT=T_1-T_2

    # save figures
    slc.set_log("Temperature",False)
    slc.set_log("solid",False)
    slc.save("slcdata_inlet"+"%4.4d"%(0))

    # calculate flux and integrate over space
    #box = ds.all_data()
    #surf = ds.surface(box, ("Temperature"),5e-27)
    #flux=surf.calculate_flux(("Temperature"), ("Temperature"),("Temperature"), ("Temperature"))
    #flux_avg = np.mean(flux)
    dTdx_3D=dT/L # mol/m^4
    k_eff=J/(dTdx_3D*area)
    k_eff_avg=np.mean(k_eff)
    print('k_eff:',k_eff,'\n')
    print('k_cell:',np.mean(k_cell),'\n')
    with open(outfile, "a+") as f:
        f.write(f'{fn},{T_1},{T_2},{dTdx_3D},{np.mean(dTdx)},{J},{np.mean(k_cell)},{k_eff}\n')
    f.close()
