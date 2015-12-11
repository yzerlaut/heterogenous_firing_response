import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('/home/yann/work/python_library/')
from my_graph import set_plot
from template_and_fitting import erfc_func, fitting_Vthre_then_Fout,\
    final_threshold_func, print_reduce_parameters
from excitability_and_sensitivities import get_mean_encoding_power

fig, AX = plt.subplots(1, 3, figsize=(14,4))
plt.subplots_adjust(bottom=.25, wspace=.3, hspace=.3, right=.95)

CELL_NUMBER = 30

cells = range(1, CELL_NUMBER+1)

DEV_MUV = np.zeros(len(cells))
DEV_SV, DEV_TV  = 0*DEV_MUV, 0*DEV_MUV
E = np.zeros((len(cells), 4))


for i in range(len(cells)):

        ##### LOADING THE DATA #####
        data = np.load('../data/cell'+str(cells[i])+'.npz')

        TvN = data['Tv_exp']/data['Cm']*data['Gl']
        for ax, x, y in zip(AX,\
                            [data['muV'], data['sV'], 100.*data['Tv_ratio']],\
                            [data['muV_exp'], data['sV_exp'], 100.*TvN]):
                ax.plot(x, y, 'ko', ms=4)
                
        DEV_MUV[i] = np.sum((data['muV_exp']-data['muV']))/len(data['muV'])
        DEV_SV[i] = np.sum((data['sV_exp']-data['sV']))/len(data['sV'])
        DEV_TV[i] = np.sum(100.*(TvN-data['Tv_ratio']))/len(TvN)
        if DEV_MUV[i]>50:
                print i
        if DEV_SV[i]>1.5:
                print i
        P = fitting_Vthre_then_Fout(data['Fout'], 1e-3*data['muV'],\
                                    1e-3*data['sV'], data['Tv_ratio'],\
                                    data['muGn'], data['Gl'], data['Cm'],
                                    data['El'], print_things=False)

        E[i,:] = np.array(get_mean_encoding_power(P, data['El'], data['Gl'], data['Cm']))
                
AX[0].plot([-70, -40], [-70, -40], 'r-', lw=3)
AX[1].plot([0, 8], [0, 8], 'r-', lw=3)
AX[2].plot([10, 120], [10, 120], 'r-', lw=3)

for ax, label, lim in zip(AX[:2], ['$\mu_V$', '$\sigma_V$'],\
                          [[-75,-35], [0.5,8.5]]):
        set_plot(ax, xlim=lim, ylim=lim,\
                 xlabel=label+' desired (mV)', ylabel=label+' observed (mV)')
set_plot(AX[2], ylim=[0.,200.],\
     xlabel='$\\tau_V^N$ desired (%)', ylabel='$\\tau_V^N$ observed (%)')

     
fig2, AX = plt.subplots(4, 3, figsize=(14,12))
plt.subplots_adjust(wspace=.3, hspace=.3)
from scipy.stats.stats import pearsonr

E_LABELS = [r"$\langle V_\mathrm{thre}^\mathrm{eff} \rangle_\mathcal{D}$",\
            r"$\langle \partial \nu / \partial \mu_V \rangle_\mathcal{D}$",\
            r"$\langle \partial \nu / \partial \sigma_V \rangle_\mathcal{D}$",\
            r"$\langle \partial \nu / \partial \tau_V^{N}' \rangle_\mathcal{D}$"]
DEV_LABELS = [r"mean $\mu_V$ deviations (mV)",\
              r"mean $\sigma_V$ deviations (mV)",\
              r"mean $\tau_V^N$ deviations (%)"]

for i in range(4):
        for j, x in zip(range(3), [np.array(DEV_MUV),\
                                   np.array(DEV_SV), np.array(DEV_TV)]):
                ii = np.isfinite(x)*np.isfinite(E[:,i])
                x, y = x[ii], E[:,i][ii]
                AX[i, j].plot(x, y, 'ko')
                cc, pp = pearsonr(x, y)
                AX[i, j].annotate('c='+str(np.round(cc,1))+', p='+'%.1e' % pp,\
                         (0.2,.9), xycoords='axes fraction', fontsize=14)
                         
                if i==3 and j==0:
                        set_plot(AX[i, j], ylabel=E_LABELS[i], xlabel=DEV_LABELS[j])
                elif j==0:
                        set_plot(AX[i, j], ylabel=E_LABELS[i])
                elif i==3:
                        set_plot(AX[i, j], xlabel=DEV_LABELS[j])
                else:
                        set_plot(AX[i, j])
                        
fig.savefig('../figures/controlling_fluctuations_by_sampling.svg', format='svg')        
fig2.savefig('../figures/controlling_fluctuations_impact_on_results.svg', format='svg')

plt.show()
