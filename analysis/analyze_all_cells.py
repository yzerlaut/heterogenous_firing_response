import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator

from plotting_tools import make_3d_and_2d_figs

# modulus in SI units
from template_and_fitting import erfc_func, fitting_Vthre_then_Fout,\
    final_threshold_func, print_reduce_parameters

def loop_over_analyzed_data_for_figs(N_CELLS=30):

    FIG_LIST = []
    
    for i in range(1, N_CELLS+1):

        ##### LOADING THE DATA #####
        data = np.load('../data/cell'+str(i)+'.npz')
        
        ##### FITTING OF THE PHENOMENOLOGICAL THRESHOLD #####
        # two-steps procedure, see template_and_fitting.py
        # need SI units !!!
        P = fitting_Vthre_then_Fout(data['Fout'], 1e-3*data['muV'],\
                                    1e-3*data['sV'], data['Tv_ratio'],\
                                    data['muGn'], data['Gl'], data['Cm'],
                                    data['El'], print_things=False)
        
        ##### PLOTTING #####
        # see plotting_tools.py
        # need non SI units (electrophy units) !!!
        FIG = make_3d_and_2d_figs(P,\
                data['Fout'], data['s_Fout'], data['muV'],\
                data['sV'], data['Tv_ratio'], data['muGn'],\
                data['Gl'], data['Cm'], data['El'], 'cell'+str(i))

        FIG.savefig('../figures/cell'+str(i)+'.png', format='png')

        FIG_LIST.append(FIG)
        
    return FIG_LIST

FIG_LIST = loop_over_analyzed_data_for_figs()
