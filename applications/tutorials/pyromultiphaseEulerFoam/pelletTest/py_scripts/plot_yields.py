import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})

inlet_cols = ['time', 'area', 'biomassflux', 'biomassvolfrac/area', 'gasflux', 'gasflux/area']
inlet_flux = pd.read_csv('mflow_bm.dat', names=inlet_cols, sep='\t')

outlet_cols = ['time', 'syngasflux', 'tarflux', 'charflux']
outlet_flux = pd.read_csv('mflow_all.dat', names=outlet_cols, sep='\t')

tar_yield = 100*(outlet_flux['tarflux']/inlet_flux['biomassflux'].loc[0])
syngas_yield = 100*(outlet_flux['syngasflux']/inlet_flux['biomassflux'].loc[0])
char_yield = 100*(outlet_flux['charflux']/inlet_flux['biomassflux'].loc[0])

plt.figure(figsize=(7,7))
plt.plot(outlet_flux['time'], tar_yield, label='Bio-oil')
plt.plot(outlet_flux['time'], syngas_yield, label='Syngas')
plt.plot(outlet_flux['time'], char_yield, label='Char')
plt.legend(fontsize=16)
plt.xlabel('Time (s)')
plt.ylabel('Yield (%)')
plt.tight_layout()
plt.savefig('pyro_yields.png')
