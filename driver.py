from headers import *
# import imp
import cell_cmb
# os.remove(getattr(cell_cmb, '__cached__', 'cell_cmb.pyc'))
# reload(cell_cmb)
from cell_cmb import *

import HO02_QE as ho
# os.remove(getattr(ho, '__cached__', 'HO02_QE.pyc'))
# reload(ho)
# import HO02_QE as ho

import TMV_QE as tmv
# os.remove(getattr(tmv, '__cached__', 'TMV_QE.pyc'))
# reload(tmv)
# import TMV_QE as tmv
# from HO02_QE import *
# from TMV_QE import *

# planck_like_analysis (smv): planck set C_ell(TE) = 0 in their Cfid matrix
# So if you want to compare to their results, it is
# equivalent to setting C_ell(TE) = 0 in M_1 and M_2 matrices of global MV
# but C_ell(TE) is not zero in f_TE and f_TB. This is the so called sub-optimal
# minimum variance estimator
smv = 1
compare_HO02 = 0
calculate_ratio_ind_ho02_tmv = 0

# dictionary with experiment spec
HO02_ref = {"name": "HO02_ref", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
            "beam": 4., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

AdvACT = {"name": "AdvACT", "lMin": 30., "lMaxT": 6000., "lMaxP": 6000.,
         "beam": 1.4, "noise_t": 10., "noise_p": 10.*np.sqrt(2)}

SO = {"name": "SO", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1.4, "noise_t": 5., "noise_p": 5.*np.sqrt(2)}

CMBS4 = {"name": "CMBS4", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

Planck_smica = {"name": "Planck", "lMin": 100., "lMaxT": 2000., "lMaxP": 2000.,
                "beam": 5., "noise_t": 35., "noise_p": 60.}

custom = {"name": "SO_TE_0", "lMin": 30., "lMaxT": 2.e3, "lMaxP": 2.e3,
         "beam": 1.4, "noise_t": 5., "noise_p": 5.*np.sqrt(2)}

time0 = time()

exp = SO
cmb = Cell_cmb(exp, pla)

# expind = np.array([10., 50., 500., 1500.])
# print cmb.totalTE(expind)
# cmb = Cell_cmb(beam=beam, noise_t=noise_t, noise_p=noise_p, lMin=lMin,
#                lMaxT=lMaxT, lMaxP=lMaxP)

est = ['TT', 'EE', 'TE', 'TB', 'EB']
# est = ['TB', 'EB']  # , 'TE', 'TB', 'EB']
# est = ['TT', 'EE', 'TE']
# est = ['TT', 'TE', 'TB', 'EB']
# est = ['TE', 'TB']

# cmb.plot_cell_2()

# changing the font of the text in the plots
fam = "serif"
plt.rcParams["font.family"] = fam


tmv_est = tmv.lensing_estimator(cmb)

"""
var_gmv = 'output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (exp['name'], str(exp['lMin']), str(exp['lMaxT']), str(exp['lMaxP']), str(exp['beam']), str(exp['noise_t']))

if os.path.exists(var_gmv):
    print "Global minimum variance for this configuration already calculated"
else:
    print "Global minimum variance for this configuration not yet calculated. Calculating it now."
    tmv_est.calc_tvar()
"""

tmv_est.calc_tvar()
tmv_est.interp_tvar()

if compare_HO02 == 1:
    totcov_ho02 = 'output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (exp['name'], str(exp['lMin']), str(exp['lMaxT']), str(exp['lMaxP']), str(exp['beam']), str(exp['noise_t']))
    if os.path.exists(totcov_ho02):
        # """
        ho02_est = ho.lensing_estimator(cmb)
        ho02_est.calc_var(est)
        ho02_est.interp_var(est)
        ho02_est.calc_cov(est)
        ho02_est.interp_cov(est)
        # """
        tmv_est.plot_var_tvar()
        tmv_est.plot_var_tvar_percent()
    else:
        print "HO02 files not present. Running HO02 estimator"
        ho02_est = ho.lensing_estimator(cmb)
        ho02_est.calc_var(est)
        ho02_est.interp_var(est)
        ho02_est.calc_cov(est)
        ho02_est.interp_cov(est)

        tmv_est.plot_var_tvar()
        tmv_est.plot_var_tvar_percent()

if calculate_ratio_ind_ho02_tmv == 1:

    cov1_ho02 = 'output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only.txt' % (exp['name'], str(exp['lMin']), str(exp['lMaxT']), str(exp['lMaxP']), str(exp['beam']), str(exp['noise_t']))
    cov2_ho02 = 'output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only.txt' % (exp['name'], str(exp['lMin']), str(exp['lMaxT']), str(exp['lMaxP']), str(exp['beam']), str(exp['noise_t']))

    if os.path.exists(cov1_ho02):
        print "HO02 files present for TT-EE-TE only pair"
        # """
        est = ['TT', 'EE', 'TE']
        ho02_est = ho.lensing_estimator(cmb)
        ho02_est.calc_var(est)
        ho02_est.interp_var(est)
        ho02_est.calc_cov(est)
        ho02_est.interp_cov(est)
        # """
    else:
        print "HO02 files not present for TT-EE-TE only pair. Running HO02 estimator"
        est = ['TT', 'EE', 'TE']
        ho02_est = ho.lensing_estimator(cmb)
        ho02_est.calc_var(est)
        ho02_est.interp_var(est)
        ho02_est.calc_cov(est)
        ho02_est.interp_cov(est)

    if os.path.exists(cov2_ho02):
        print "HO02 files present for TB-EB only pair"
        # """
        est = ['TB', 'EB']
        ho02_est = ho.lensing_estimator(cmb)
        ho02_est.calc_var(est)
        ho02_est.interp_var(est)
        ho02_est.calc_cov(est)
        ho02_est.interp_cov(est)
        # """
    else:
        print "HO02 files not present for TB-EB only pair. Running HO02 estimator"
        est = ['TB', 'EB']
        ho02_est = ho.lensing_estimator(cmb)
        ho02_est.calc_var(est)
        ho02_est.interp_var(est)
        ho02_est.calc_cov(est)
        ho02_est.interp_cov(est)

    tmv_est.plotvar_tvar_ratio()


tmv_est.plot_tvar()

multexp = [Planck_smica, SO, CMBS4]
# tmv_est.plotvar_tvar_ratio_multexp(multexp)
# tmv_est.plotvar_var_tvar_percent_multexp(multexp)
tmv_est.plot_plancklike_var_tvar_ratio_multexp(multexp)

"""
ho02_est = ho.lensing_estimator(cmb)
# ho02_est.calc_var(est)
ho02_est.interp_var(est)
# ho02_est.calc_cov(est)
ho02_est.interp_cov(est)
ho02_est.plot_corrcoef(est)
# """

# SNT2, SNH2, cumSNT, cumSNH = tmv_est.SNR_comp(exp)

print time() - time0
