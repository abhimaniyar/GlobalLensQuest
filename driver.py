import time
# import imp
import cell_cmb
# imp.reload(cell_cmb)
from cell_cmb import *

# imp.reload(estimator)
import HO02_QE as ho
import TMV_QE as tmv
# from HO02_QE import *
# from TMV_QE import *

compare_HO02 = 0
calculate_ratio_ind_ho02_tmv = 0

# dictionary with experiment spec
HO02_ref = {"name": "HO02_ref", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
            "beam": 4., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

AdvACT = {"name": "AdvACT", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1.4, "noise_t": 10., "noise_p": 10.*np.sqrt(2)}

SO = {"name": "SO", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1.4, "noise_t": 5., "noise_p": 5.*np.sqrt(2)}

CMBS4 = {"name": "CMBS4", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

custom = {"name": "Custom", "lMin": 30., "lMaxT": 3000., "lMaxP": 5000.,
         "beam": 1., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

time0 = time()

exp = CMBS4
cmb = Cell_cmb(exp)

# cmb = Cell_cmb(beam=beam, noise_t=noise_t, noise_p=noise_p, lMin=lMin,
#                lMaxT=lMaxT, lMaxP=lMaxP)

est = ['TT', 'EE', 'TE', 'TB', 'EB']
# est = ['TB', 'EB']  # , 'TE', 'TB', 'EB']
# est = ['TT', 'EE', 'TE']
# est = ['TT', 'TE', 'TB', 'EB']
# est = ['TE', 'TB']

# cmb.plot_cell_2(np.linspace(2, int(lMaxT), int(lMaxT)-1))
# cmb = Cell_cmb(beam=4., noise_t=1., noise_p=1.*np.sqrt(2), lMin=30.,
#                lMaxT=3.e3, lMaxP=5.e3)
# cmb = Cell_cmb(beam=7., noise_t=27., noise_p=40*np.sqrt(2.), lMin=30.,
#                lMaxT=3.e3, lMaxP=5.e3)

tmv_est = tmv.lensing_estimator(cmb)

# tmv_est.calc_tvar()

tmv_est.interp_tvar()
# tmv_est.plot_tvar()

if compare_HO02 == 1:
    totcov_ho02 = 'output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (exp['name'], str(exp['lMin']), str(exp['lMaxT']), str(exp['lMaxP']), str(exp['beam']), str(exp['noise_t']))
    if os.path.exists(totcov_ho02):
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
    else:
        print "HO02 files not present for TB-EB only pair. Running HO02 estimator"
        est = ['TB', 'EB']
        ho02_est = ho.lensing_estimator(cmb)
        ho02_est.calc_var(est)
        ho02_est.interp_var(est)
        ho02_est.calc_cov(est)
        ho02_est.interp_cov(est)

    tmv_est.plotvar_tvar_ratio()

multexp = [SO, CMBS4]
# tmv_est.plotvar_tvar_ratio_multexp(multexp)
tmv_est.plotvar_var_tvar_percent_multexp(multexp)

"""
ho02_est = ho.lensing_estimator(cmb)
# ho02_est.calc_var(est)
ho02_est.interp_var(est)
# ho02_est.calc_cov(est)
ho02_est.interp_cov(est)
ho02_est.plot_corrcoef(est)
"""

print time() - time0
