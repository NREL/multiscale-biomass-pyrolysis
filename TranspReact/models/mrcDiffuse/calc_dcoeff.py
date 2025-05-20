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
outfile='dcoeff.csv'

with open(outfile, "w+") as f:
    f.write('pltfile,outletT,inletT,dTdX_3d,dTdx,avgJoutlet,avgJinlet,D,D_IP\n')
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
    outlet[flowdir]=0.001*prob_hi[flowdir]
    print('outlet:',outlet)
    size=prob_hi-prob_lo

    # load fields at outlet
    slc = yt.SlicePlot(ds, flowdir_char, ["solid", "Temperature"],center=[outlet[0],outlet[1],outlet[2]])
    res=np.array([ncells[(flowdir+1)%3],ncells[(flowdir+2)%3]])
    frb = slc.data_source.to_frb(((size[(flowdir+1)%3],'cm'),(size[(flowdir+2)%3],'cm')),res) 
    vfrac=np.array(frb["solid"])
    print(np.mean(vfrac))
    temp=np.array(frb["Temperature"])
    dTdx=np.array(frb["Temperature_gradient_"+flowdir_char])
    # conc gradient at every point on a slice
    # value of 1/((1-vfrac)/solid_dif + vfrac/fluid_dif at every point on a slice
    # solve for J

    # compare to 1/((1-vfrac)/solid_dif + vfrac/fluid_dif solved at all points on all slices with spacing dx
    solid_dif=1.106e-5
    fluid_dif=22.39e-6
    # 1 dimensional approximation on each individual cell
    # integrated over whole volume will be more accurate 3D
    k_cell=(1-vfrac)*solid_dif + vfrac*fluid_dif
    #D_avg=np.mean(D)
    J=-D*rho*dTdx

    area=size[(flowdir+1)%3]*size[(flowdir+2)%3] # m^2
    J_avg=np.mean(J)


    dX=size[flowdir] # m
    area=size[(flowdir+1)%3]*size[(flowdir+2)%3] # m^2
    print(area)
    T_1=np.mean(temp) # mol/m^3

    # Save figures
    slc.set_log("solid",False)
    slc.set_log("Temperature",False)
    slc.save("slcdata_outlet"+"%4.4d"%(0))

    # calculate avg concentration at inlet and dC
    inlet=0.5*(prob_lo+prob_hi)
    inlet[flowdir]=0.999*prob_hi[flowdir]
    print('inlet:',inlet)
    slc = yt.SlicePlot(ds, flowdir_char, ["solid","Temperature"],center=[inlet[0],inlet[1],inlet[2]])
    frb = slc.data_source.to_frb(((size[(flowdir+1)%3],'cm'),(size[(flowdir+2)%3],'cm')),res) 
    dTdx=np.array(frb["Temperature_gradient_"+flowdir_char])
    vfrac=np.array(frb["solid"])
    temp=np.array(frb["Temperature"]) #mol/m^3

    J=-D*dTdx
    J_avg_2=np.mean(J)

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
    dTdx_3D=dT/dX # mol/m^4
    flux_avg=J_avg
    D_IP=-J/dTdx_3D
    D_IP=np.mean(D_IP)
    print('D_IP:',D_IP,'\n')
    with open(outfile, "a+") as f:
        f.write(f'{fn},{T_1},{T_2},{dTdx_3D},{np.mean(dTdx)},{J_avg},{J_avg_2},{D_avg},{D_IP}\n')
    f.close()
