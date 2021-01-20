from imports import *
# import imp
import cell_cmb
# os.remove(getattr(cell_cmb, '__cached__', 'cell_cmb.pyc'))
# reload(cell_cmb)
from cell_cmb import *

import HO02_QE as ho
# os.remove(getattr(ho, '__cached__', 'HO02_QE.pyc'))
# reload(ho)
# import HO02_QE as ho

import GMV_QE as gmv
# os.remove(getattr(gmv, '__cached__', 'TMV_QE.pyc'))
# reload(gmv)

import SQE as sqe
import OH03_flatsky as oh

import plots_compare as pltcomp

# in case need to compare HO02 estimator with GMV, use 1 otherwise 0
compare_HO02 = 0
calculate_ratio_ind_ho02_tmv = 0

# dictionary with experiment spec
HO02_ref = {"name": "HO02_ref", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
            "beam": 4., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

AdvACT = {"name": "AdvACT", "lMin": 30., "lMaxT": 6000., "lMaxP": 6000.,
         "beam": 1.4, "noise_t": 10., "noise_p": 10.*np.sqrt(2)}

SO = {"name": "SO", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1.4, "noise_t": 8., "noise_p": 8.*np.sqrt(2)}

CMBS4 = {"name": "CMBS4", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

Planck_smica = {"name": "Planck", "lMin": 100., "lMaxT": 3000., "lMaxP": 3000.,
                "beam": 5., "noise_t": 35., "noise_p": 60.}

custom = {"name": "custom", "lMin": 30., "lMaxT": 2.e3, "lMaxP": 2.e3,
         "beam": 10., "noise_t": 60., "noise_p": 100.}

time0 = time()

exp = CMBS4

cmb = Cell_cmb(exp)

est = ['TT', 'EE', 'TE', 'TB', 'EB']
# est = ['TB', 'EB']  # , 'TE', 'TB', 'EB']
# est = ['TT', 'EE', 'TE']
# est = ['TT', 'TE', 'TB', 'EB']
# est = ['TE', 'TB']


fam = "serif"
plt.rcParams["font.family"] = fam

# """
# Global minimum variance estimator (our work)
gmv_est = gmv.lensing_estimator(cmb)
print "Running GMV estimator"
gmv_est.calc_tvar()
gmv_est.interp_tvar()
# gmv_est.plot_tvar()

# """
# HO02 estimator Hu and Okamoto (2002)
ho02_est = ho.lensing_estimator(cmb)
print "Running HO02 estimator"
ho02_est.calc_var(est)
ho02_est.interp_var(est)
ho02_est.calc_cov(est)
ho02_est.interp_cov(est)
# """

if compare_HO02 == 1:
    pltcomp.plot_var_tvar(exp)
    pltcomp.plot_var_tvar_percent(exp)


if calculate_ratio_ind_ho02_tmv == 1:
    est1 = ['TT', 'EE', 'TE']
    ho02_est = ho.lensing_estimator(cmb)
    print "Running HO02 estimator for TT-EE-TE only"
    ho02_est.calc_var(est1)
    ho02_est.interp_var(est1)
    ho02_est.calc_cov(est1)
    ho02_est.interp_cov(est1)

    est2 = ['TB', 'EB']
    ho02_est = ho.lensing_estimator(cmb)
    print "Running HO02 estimator for TB-EB only"
    ho02_est.calc_var(est2)
    ho02_est.interp_var(est2)
    ho02_est.calc_cov(est2)
    ho02_est.interp_cov(est2)

    pltcomp.plotvar_tvar_ratio(exp)


# """
# Suboptimal quadratic estimator (Planck lensing 2018)
sqe_est = sqe.sqe_lensing_estimator(cmb)
print "Running SQE estimator"
sqe_est.calc_lambda()
sqe_est.interp_lambda()
sqe_est.calc_cov()
sqe_est.interp_cov()
# sqe_est.plot_cov()

# """
# OH03 estimator (Okamoto and Hu 2003)
oh_est = oh.oh03_lensing_estimator(cmb)
print "Running OH03"
oh_est.calc_var(est)
oh_est.interp_var(est)
oh_est.calc_cov(est)
oh_est.interp_cov(est)

# """

# comparison between different estimators for diff experiments
multexp = [Planck_smica, SO, CMBS4]
# pltcomp.plotvar_tvar_ratio_multexp(multexp)
# pltcomp.plotvar_var_tvar_percent_multexp(multexp)
pltcomp.plot_plancklike_var_tvar_ratio_multexp(multexp)
# pltcomp.plot_mv_multest_multexp(multexp)
# oh_est.plot_TE_comparison(multexp)
oh_est.plot_MV_comparison(multexp)

# SNT2, SNH2, cumSNT, cumSNH = tmv_est.SNR_comp(exp)

print time() - time0
