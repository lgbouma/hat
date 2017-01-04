import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from astrobase import plotbase as pb

field = 'G199/'
DSP_lim = 30

head_path = '../data/DEBiL_heads/'+field[:-1]+'_DSP'+str(DSP_lim)+'.txt'
DEBiL_head = np.genfromtxt(head_path, dtype=['U15',float], names=['hatid', 'blsperiod'])

plt.ioff()
plt.close('all')
for name, period in DEBiL_head:
    LC_path = '../data/DEBiL_LCs/'+field+name+'.txt'
    LC = np.genfromtxt(LC_path, dtype=[np.float64,np.float64,np.float64], names=['time','epd000','err'])
    out_file = '../results/prelim/'+name+'.png'
    if not os.path.exists(out_file):
        pb.plot_phased_mag_series(LC['time'], LC['epd000'], period, outfile=out_file)
    else:
        continue
