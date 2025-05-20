import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
    
fn_pattern='plt*'
fn_list = sorted(glob.glob(fn_pattern), key=lambda f: int(f.split("plt")[1]))
print(fn_list)
flowdir=int(argv[1])
flowdir_char=chr(ord('x')+flowdir)
print(flowdir_char)
outfile='vol_avg_temp.csv'

with open(outfile, "w+") as f:
    f.write('pltfile,vol avg temp (K)\n')
f.close()
for i, fn in enumerate(fn_list):
    ds=yt.load(fn)
    field = ("Temperature")  # The field to average
    weight = ("cell_volume")  # The weight for the average
    ad = ds.all_data()  # This is a region describing the entire box
    # but note it doesn't read anything in yet!
    # We now use our 'quantities' call to get the average quantity
    avg_temp = ad.quantities.weighted_average_quantity(field, weight)
    avg_temp_str=str(avg_temp)
    avg_temp_num=avg_temp_str.replace("dimensionless","")
    print(avg_temp_num)
    with open(outfile, "a+") as f:
        f.write(f'{fn},{avg_temp_num}\n')
    f.close()
