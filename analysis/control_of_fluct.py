import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('/home/yann/work/python_library/')
from my_graph import set_plot

fig, AX0 = plt.subplots(4, 3, figsize=(15,18))
AX = AX0[3,:]
plt.subplots_adjust(hspace=.4, wspace=.4)

CELL_NUMBER = 30

DEV_MUV = np.zeros(CELL_NUMBER)
DEV_SV, DEV_TV  = 0*DEV_MUV, 0*DEV_MUV
# for histograms
n_bins = 20

#### STARTING WITH THE DATA ####

muV = np.linspace(-70, -40, n_bins)
sV = np.linspace(0.5, 8, n_bins)
TvN = np.linspace(10, 110, n_bins)

ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP = [], [], []
for X, Y in zip([muV, sV, TvN],\
                [ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP]):
                for i in range(len(X)-1):
                    Y.append([])

for i in range(1, CELL_NUMBER+1):

        ##### LOADING THE DATA #####
        data = np.load('../data/cell'+str(i)+'.npz')
        for ax, X, ALL_Y, x, y in zip(AX,\
                            [muV, sV, TvN],\
                            [ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP],\
                            [data['muV'], data['sV'], 100.*data['TvN']],\
                            [data['muV_exp'], data['sV_exp'], 100.*data['TvN_exp']]):
                # ax.plot(x, y, 'o', ms=1, color='lightgray')
                for ix, xx1, xx2 in zip(range(len(X)-1), X[:-1], X[1:]):
                        ALL_Y[ix] = np.concatenate([ALL_Y[ix],y[(x>xx1) & (x<xx2) & (y<np.inf)]])
                        
        DEV_MUV[i-1] = np.sum((data['muV']-data['muV_exp'])**2)/len(data['muV'])
        DEV_SV[i-1] = np.sum((data['sV']-data['sV_exp'])**2)/len(data['sV'])
        DEV_TV[i-1] = np.sum((data['TvN']-data['TvN_exp'])**2)/len(data['TvN'])

MEAN_MUV_EXP, MEAN_SV_EXP, MEAN_TV_EXP = 0.*muV, 0.*sV, 0*TvN 
STD_MUV_EXP, STD_SV_EXP, STD_TV_EXP = 0.*muV, 0.*sV, 0*TvN 
        
for ax, X, Y, sY, ALL in zip(AX, [muV, sV, TvN],\
                [MEAN_MUV_EXP, MEAN_SV_EXP, MEAN_TV_EXP],\
                [STD_MUV_EXP, STD_SV_EXP, STD_TV_EXP],\
                [ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP]):
    for i in range(len(X)-1):
        if len(ALL[i])>1:
            Y[i] = np.mean(ALL[i])
            sY[i] = np.std(ALL[i])
                
    ax.errorbar(.5*(X[1:]+X[:-1])[Y!=0], Y[Y!=0], sY[Y!=0], color='k', lw=3, label='in vitro \n data')
    
AX[0].plot([-70, -40], [-70, -40], '--', lw=1, color='k')
AX[1].plot([0.5, 8], [0.5, 8], '--', lw=1, color='k')
AX[2].plot([10, 100], [10, 100], '--', lw=1, color='k', label='desired')

for ax, label, lim in zip(AX[:2], ['$\mu_V$', '$\sigma_V$'],\
                          [[-72,-38], [0.,9.]]):
        set_plot(ax, xlim=lim, ylim=lim,\
                 xlabel=label+' desired (mV)', ylabel=label+' observed (mV)')

                 
AX[2].legend(loc='best', frameon=False, prop={'size':'small'})

set_plot(AX[2], ylim=[0.,120.],\
     xlabel='$\\tau_V^N$ desired (%)', ylabel='$\\tau_V^N$ observed (%)')

#### Then theoretical models ####

muV = np.linspace(-70, -50, n_bins)
sV = np.linspace(0.5, 8, n_bins)
TvN = np.linspace(10, 110, n_bins)

for ii, model, color, label in zip(range(3),\
                                   ['SUBTHRE', 'LIF', 'iAdExp'],\
                                   ['b', 'k', 'r'], ['RC-circuit', 'LIF', 'iAdExp']):
                               
    AX = AX0[ii,:]
    ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP = [], [], []
    for X, Y in zip([muV, sV, TvN],\
                    [ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP]):
                    for i in range(len(X)-1):
                        Y.append([])

    ##### LOADING THE DATA #####
    data = np.load('../data/'+model+'.npz')
    for ax, X, ALL_Y, x, y in zip(AX,\
                        [muV, sV, TvN],\
                        [ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP],\
                        [data['muV'], data['sV'], 100.*data['TvN']],\
                        [data['muV_exp'], data['sV_exp'], 100.*data['TvN_exp']]):
            for ix, xx1, xx2 in zip(range(len(X)-1), X[:-1], X[1:]):
                    ALL_Y[ix] = np.concatenate([ALL_Y[ix],y[(x>xx1) & (x<xx2) & (y<np.inf)]])

    DEV_MUV[i-1] = np.sum((data['muV']-data['muV_exp'])**2)/len(data['muV'])
    DEV_SV[i-1] = np.sum((data['sV']-data['sV_exp'])**2)/len(data['sV'])
    DEV_TV[i-1] = np.sum((data['TvN']-data['TvN_exp'])**2)/len(data['TvN'])

    MEAN_MUV_EXP, MEAN_SV_EXP, MEAN_TV_EXP = 0.*muV, 0.*sV, 0*TvN 
    STD_MUV_EXP, STD_SV_EXP, STD_TV_EXP = 0.*muV, 0.*sV, 0*TvN 

    for ax, X, Y, sY, ALL in zip(AX, [muV, sV, TvN],\
                    [MEAN_MUV_EXP, MEAN_SV_EXP, MEAN_TV_EXP],\
                    [STD_MUV_EXP, STD_SV_EXP, STD_TV_EXP],\
                    [ALL_MUV_EXP, ALL_SV_EXP, ALL_TV_EXP]):
        for i in range(len(X)-1):
            if len(ALL[i])>1:
                Y[i] = np.mean(ALL[i])
                sY[i] = np.std(ALL[i])

        ax.errorbar(.5*(X[1:]+X[:-1])[Y!=0], Y[Y!=0], sY[Y!=0], color=color, lw=2, label=label)

    AX[0].plot([-70, -50], [-70, -50], '--', lw=1, color='k')
    AX[1].plot([0.5, 8], [0.5, 8], '--', lw=1, color='k')
    AX[2].plot([10, 100], [10, 100], '--', lw=1, color='k', label='desired')
    

    for ax, label, lim in zip(AX[:2], ['$\mu_V$', '$\sigma_V$'],\
                              [[-72,-48], [0.,9.]]):
            set_plot(ax, xlim=lim, ylim=lim,\
                     xlabel=label+' desired (mV)', ylabel=label+' observed (mV)')

                 
    AX[2].legend(loc='best', frameon=False, prop={'size':'small'})

    set_plot(AX[2], ylim=[0.,140.],\
         xlabel='$\\tau_V^N$ desired (%)', ylabel='$\\tau_V^N$ observed (%)')

fig.savefig('../figures/approx_allows_control_of_fluct.svg', format='svg')

plt.show()
