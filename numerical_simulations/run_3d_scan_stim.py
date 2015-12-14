import numpy as np
from simulations import single_experiment, params_variations_calc
import models

import sys
sys.path.append('../analysis/')
from measures import measuring_subthre_dynamics

DISCRET_muV, DISCRET_sV= 4, 8
DISCRET_TvN= 4 # in those three dimensions, only the autocorrelation is varied


def make_simulation_for_model(MODEL, return_output=False,\
                              precision='low'):


    if precision is 'low':
        # discretization and seed
        SEED = np.arange(2)+1
        dt, tstop = 1e-4, 2.
    else:
        SEED = np.arange(3)+1
        dt, tstop = 1e-5, 10.

                              
    params = models.get_model_params(MODEL, {}) # paramters of the model, see models.py
        
    ### PARAMETERS OF THE EXPERIMENT
    muV_min, muV_max, sV_min1, sV_max1, sV_min2, sV_max2, Ts_ratio = params['RANGE_FOR_3D']
    muGn_min, muGn_max = 1.15, 8.

    ### cell and synaptic parameters, see models.py !!
    sim_params = {'dt':dt, 'tstop':tstop}
    t_long = np.arange(0,int(tstop/dt))*dt

    muV = np.linspace(muV_min, muV_max, DISCRET_muV, endpoint=True)

    # trying with the linear autocorrelation 
    DISCRET_muG = DISCRET_TvN
    Tv_ratio = np.linspace(1./muGn_max+Ts_ratio, 1./muGn_min+Ts_ratio, DISCRET_muG, endpoint=True)
    muGn = 1./(Tv_ratio-Ts_ratio)

    muV, sV, muGn = np.meshgrid(muV, np.zeros(DISCRET_sV), muGn)
    Tv_ratio = Ts_ratio+1./muGn

    for i in range(DISCRET_muV):
        sv1 = sV_min1+i*(sV_min2-sV_min1)/(DISCRET_muV-1)
        sv2 = sV_max1+i*(sV_max2-sV_max1)/(DISCRET_muV-1)
        for l in range(DISCRET_muG):
            sV[:,i,l] = np.linspace(sv1,sv2,DISCRET_sV,endpoint=True)

    I0, Gs, f, Q, Ts = params_variations_calc(muGn,muV,sV,Ts_ratio*np.ones(muGn.shape),params)

    Fout = np.zeros((DISCRET_sV, DISCRET_muV, DISCRET_muG, len(SEED)))

    # measuring the subthreshold fluctuations properties
    muV_exp, sV_exp, TvN_exp = 0*Fout, 0*Fout, 0*Fout

    for i_muV in range(DISCRET_muV):
        print '[[[[]]]]=====> muV : ', round(1e3*muV[0, i_muV, 0],1), 'mV'
        for i_sV in range(DISCRET_sV):
            print '[[[]]]====> sV : ', round(1e3*sV[i_sV, i_muV, 0],1), 'mV'
            for ig in range(DISCRET_muG):
                print '[]=> muGn : ', round(muGn[i_sV, i_muV, ig],1),\
                    'TvN : ', round(100*Tv_ratio[i_sV, i_muV, ig],1), '%'
                for i_s in range(len(SEED)):
                    v, spikes = single_experiment(\
                        t_long, I0[i_sV, i_muV, ig],\
                        Gs[i_sV, i_muV, ig],\
                        f[i_sV, i_muV, ig],\
                        Q[i_sV, i_muV, ig],\
                        Ts[i_sV, i_muV, ig],\
                        muV[i_sV, i_muV, ig],\
                        params, MODEL=MODEL,
                        seed=SEED[i_s]+i_muV+i_sV+ig)

                    Fout[i_sV, i_muV, ig, i_s] =\
                      len(spikes)/t_long.max()

                    muV_exp[i_sV, i_muV, ig, i_s],\
                      sV_exp[i_sV, i_muV, ig, i_s],\
                      TvN_exp[i_sV, i_muV, ig, i_s] = \
                        measuring_subthre_dynamics(v, spikes, dt,\
                                   Tm0=params['Cm']/params['Gl'])

    data_path = '../data/'+MODEL+'.npz'


    D = dict(muV=1e3*muV.flatten(), sV=1e3*sV.flatten(),\
             TvN=Ts_ratio+1./muGn.flatten(),\
             muGn=muGn.flatten(),\
             Fout=Fout.mean(axis=-1).flatten(),\
             s_Fout=Fout.std(axis=-1).flatten(),\
             muV_exp=1e3*muV_exp.mean(axis=-1).flatten(),\
             sV_exp=1e3*sV_exp.mean(axis=-1).flatten(),\
             TvN_exp=TvN_exp.mean(axis=-1).flatten(),\
             s_muV_exp=1e3*muV_exp.std(axis=-1).flatten(),\
             s_sV_exp=1e3*sV_exp.std(axis=-1).flatten(),\
             s_TvN_exp=TvN_exp.std(axis=-1).flatten(),\
             MODEL=MODEL,\
             Gl=params['Gl'], Cm=params['Cm'], El=params['El'])

    np.savez(data_path,**D)    

if __name__=='__main__':
    # for spiking properties, what model ?? see models.py
    import argparse
    parser=argparse.ArgumentParser(description=
     """ 
     Stimulate a reconstructed cell with a shotnoise and study Vm dynamics
     """
    ,formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("MODEL", help="Choose a model of NEURON")
    parser.add_argument("-t", "--WITH_TM_VARIATIONS",\
                        help="we vary tm", action="store_false")
    parser.add_argument("--precision", default='low',\
                        help="turn to 'high' for simulations as in the paper")

    args = parser.parse_args()
    
    Mlist = models.models(args.MODEL)

    if args.WITH_TM_VARIATIONS:
        def make_sim(MODEL):
            print 'by default we vary Tm !!'
            make_simulation_for_model(MODEL, precision=args.precision)
            make_simulation_for_model(MODEL+'__minus', precision=args.precision)
            make_simulation_for_model(MODEL+'__plus', precision=args.precision)
    else:
        def make_sim(MODEL):
            make_simulation_for_model(MODEL, precision=args.precision)
        
    if Mlist is None: # means it is a single model
        make_sim(args.MODEL)
    else:
        for m in Mlist:
            print m
            make_sim(m)
