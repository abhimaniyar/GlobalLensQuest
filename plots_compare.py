from imports import *


def plot_var_tvar(exp):
    clexp = exp

    data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
    L = data[:, 0]

    data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
    L2 = data2[:, 0]

    plt.figure()
    plt.plot(L, L*(L+1)*data[:, -1]/(2*np.pi), 'r', label='TMV')

    plt.plot(L2, L2*(L2+1)*data2[:, -1]/(2*np.pi), 'b', label='HO02')
    # n1 = 1./((1./data2[:, 1]) + (1./data2[:, 2]) + (1./data2[:, 3]))
    # plt.plot(L2, L2*(L2+1)*n1/(2*np.pi), 'b', label='HO')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$L(L+1)N_L^{dd}/2\pi$', fontsize=16)
    plt.xlabel(r'$L$', fontsize=16)
    plt.legend(prop={'size':12})
    plt.tick_params(axis='both', labelsize=14)


def plot_var_tvar_percent(exp):
    clexp = exp

    data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))

    L = data[:, 0]

    data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
    L2 = data2[:, 0]

    plt.figure()
    interp_HO = interp1d(L2, data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
    interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
    L_p = np.logspace(np.log10(1.), np.log10(clexp['lMaxP']), 201, 10.)
    # L_p = np.logspace(np.log10(1.), np.log10(10000.), 51, 10.)
    # vart_var = interp_tmv(L_p)/interp_HO(L_p)  # data2[:, -1]
    vart_var = (interp_HO(L_p)-interp_tmv(L_p))*100./interp_HO(L_p)
    ind = np.argmax(vart_var)
    print("max percent improvement of %s happens at L = %s" % (vart_var[ind], L_p[ind]))
    ind2 = np.max(np.where(L_p < 100.))
    print("At L = %s percent improvement is %s" % (L_p[ind2], vart_var[ind2]))
    plt.plot(L_p, vart_var, 'b')
    plt.xscale('log')
    # plt.ylim(0.1, 1.1)
    plt.ylim(ymax=15)
    # plt.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    # plt.ylabel(r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$')
    plt.ylabel(r'$(N_{mv}^\mathrm{HO02}-N_{mv}^\mathrm{TMV}) \times 100/N_{mv}^\mathrm{HO02}$', fontsize=16)
    plt.xlabel(r'$L$', fontsize=16)
    plt.legend(prop={'size': 12})
    plt.tick_params(axis='both', labelsize=14)


def plotvar_tvar_ratio(exp):
    clexp = exp

    data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))

    L = data[:, 0]

    interp_tmv1 = interp1d(L, data[:, 1], kind='quadratic', bounds_error=False, fill_value=0.)
    interp_tmv2 = interp1d(L, data[:, 2], kind='quadratic', bounds_error=False, fill_value=0.)

    data21 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
    interp_HO_1 = interp1d(data21[:, 0], data21[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

    data22 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
    interp_HO_2 = interp1d(data22[:, 0], data22[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

    L_p = np.logspace(np.log10(1.), np.log10(clexp['lMaxP']), 201, 10.)
    plt.figure()
    vart_var1 = interp_tmv1(L_p)/interp_HO_1(L_p)
    vart_var2 = interp_tmv2(L_p)/interp_HO_2(L_p)
    plt.plot(L_p, vart_var1, 'b', label='TT-EE-TE')
    plt.plot(L_p, vart_var2, 'r', label='EB-TB')
    plt.legend()
    """
    interp2 = interp1d(L, data[:, -1], kind='linear', bounds_error=False, fill_value='extrapolate')
    vart_var2 = data2[:, -1]/interp2(L2)
    plt.plot(L2, vart_var2, 'b')
    # """
    plt.xscale('log')
    plt.ylim(0.8, 1.05)
    plt.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    plt.ylabel(r'$N_{mv}^\mathrm{TMV}/N_{mv}^\mathrm{HO02}$', fontsize=16)
    plt.xlabel(r'$L$', fontsize=16)
    plt.legend(prop={'size':12})
    plt.tick_params(axis='both', labelsize=14)
    # print vart_var2
    # print vart_var1


def plotvar_tvar_ratio_multexp(exp):
    lines = ["-", "--", "-."]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
    """
    nexp = len(exp)

    f, ax = plt.subplots(nrows=nexp, sharex=True)
    for e in range(nexp):
        clexp = exp[e]
        data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]

        interp_tmv1 = interp1d(L, data[:, 1], kind='quadratic', bounds_error=False, fill_value=0.)
        interp_tmv2 = interp1d(L, data[:, 2], kind='quadratic', bounds_error=False, fill_value=0.)
        interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data21 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO_1 = interp1d(data21[:, 0], data21[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data22 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO_2 = interp1d(data22[:, 0], data22[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO = interp1d(data2[:, 0], data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        # L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
        L_p = np.linspace(2, 5000, 4999)
        vart_var1 = interp_tmv1(L_p)/interp_HO_1(L_p)
        vart_var2 = interp_tmv2(L_p)/interp_HO_2(L_p)
        vart_var = interp_tmv(L_p)/interp_HO(L_p)

        ax[e].plot(L_p, vart_var1, 'b',
                   label='%s TT-EE-TE' % (clexp['name']))
        ax[e].plot(L_p, vart_var2, 'r',
                   label='%s TB-EB' % (clexp['name']))
        ax[e].plot(L_p, vart_var, 'k',
                   label='%s MV' % (clexp['name']))
        ax[e].set_ylim(0.8, 1.05)
        ax[e].set_ylabel(r'$N_{mv}^\mathrm{GMV}/N_{mv}^\mathrm{HO02}$',
                         fontsize=22)
        ax[e].legend(prop={'size': 17}, loc='upper left', ncol=3, frameon=False,
              labelspacing=0.2)
        ax[e].set_xscale('log')
        ax[e].set_xlim(2.0, clexp['lMaxP'])
        ax[e].set_xlabel(r'$L$', fontsize=22)
        ax[e].tick_params(axis='both', labelsize=22)
        # f.subplots_adjust(wspace=0)
        f.subplots_adjust(hspace=0)
    
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]

        interp_tmv1 = interp1d(L, data[:, 1], kind='quadratic', bounds_error=False, fill_value=0.)
        interp_tmv2 = interp1d(L, data[:, 2], kind='quadratic', bounds_error=False, fill_value=0.)
        interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data21 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO_1 = interp1d(data21[:, 0], data21[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data22 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO_2 = interp1d(data22[:, 0], data22[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO = interp1d(data2[:, 0], data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        # L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
        L_p = np.linspace(2, 5000, 4999)
        vart_var1 = interp_tmv1(L_p)/interp_HO_1(L_p)
        vart_var2 = interp_tmv2(L_p)/interp_HO_2(L_p)
        vart_var = interp_tmv(L_p)/interp_HO(L_p)

        # ax.plot(L_p, vart_var1, 'b', ls=lines[e],
        #         label='%s TT-EE-TE' % (clexp['name']))
        # ax.plot(L_p, vart_var2, 'r', ls=lines[e],
        #          label='%s TB-EB' % (clexp['name']))
        ax.plot(L_p, vart_var, 'k', ls=lines[e],
                 label='%s MV' % (clexp['name']))

    # ls = ax.get_lines()
    # leg1 = plt.legend([ls[i] for i in [0, 1]], [exp[i]['name'] for i in range(len(exp))], loc=1)
    # leg2 = plt.legend(custom_lines, ['TT-EE-TE', 'TB-EB'], loc=4)
    ax.legend(prop={'size': 17}, loc='upper left', ncol=2, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_xlim(2.0, clexp['lMaxP'])
    ax.set_ylim(0.8, 1.05)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax.set_ylabel(r'$N_{mv}^\mathrm{GMV}/N_{mv}^\mathrm{HO02}$',
                  fontsize=22)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=22)
    # ax.legend(custom_lines, ['TT-EE-TE', 'TB-EB'])
    # ax.add_artist(leg1)
    # ax.add_artist(leg2)
    # np.save('figures/multexp_vartrueind_varHO02ind.png')


def plotvar_var_tvar_percent_multexp(exp):
    lines = ["-", "--", "-."]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]
        interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L2 = data2[:, 0]
        interp_HO = interp1d(L2, data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)

        vart_var = (interp_HO(L_p)-interp_tmv(L_p))*100./interp_HO(L_p)

        ax.plot(L_p, vart_var, 'b', ls=lines[e],
                label='%s' % (clexp['name']))

    """
    plt.legend([lines[0], lines[1], lines[2]],
               [exp[0]['name'], exp[1]['name'], exp[2]['name']],
               numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               frameon=False, prop={'size': 12}, fancybox=True)
    """
    # ls = ax.get_lines()
    # leg1 = plt.legend([ls[i] for i in [0, 1]], [exp[i]['name'] for i in range(len(exp))], loc=1)
    # leg2 = plt.legend(custom_lines, ['TT-EE-TE', 'TB-EB'], loc=4)
    ax.legend(prop={'size': 17}, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_xlim(1.0, clexp['lMaxP'])
    ax.set_ylim(ymax=15)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax.set_ylabel(r'$(N_{mv}^\mathrm{HO02}-N_{mv}^\mathrm{TMV}) \times 100/N_{mv}^\mathrm{HO02}$', fontsize=22)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=22)
    # np.save('figures/multexp_percent_improvement_overHO02.png')


def plot_plancklike_var_tvar_ratio_multexp(exp):
    lines = ["-", "--", "-."]
    cl = ["b", "r", "g"]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]

        interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data3 = np.genfromtxt('output/SQE_variance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L3 = data3[:, 0]

        interp_smv = interp1d(L3, data3[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        # data2 = np.genfromtxt('output/OH03_flatsky_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO = interp1d(data2[:, 0], data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        # L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
        L_p = np.linspace(2, 5000, 3*4999)

        vart_smv = interp_tmv(L_p)/interp_smv(L_p)
        var_smv = interp_HO(L_p)/interp_smv(L_p)
        # vart_var = interp_tmv(L_p)/interp_HO(L_p)

        ax.plot(L_p, vart_smv, cl[e], ls=lines[0],
                label='%s GMV/SQE' % (clexp['name']))
        ax.plot(L_p, var_smv, cl[e], ls=lines[1],
                label='%s HO/SQE' % (clexp['name']))
        # ax.plot(L_p, vart_var, cl[e], ls=lines[2],
        #         label='%s GMV_HO' % (clexp['name']))

    # ls = ax.get_lines()
    # leg1 = plt.legend([ls[i] for i in [0, 1]], [exp[i]['name'] for i in range(len(exp))], loc=1)
    # leg2 = plt.legend(custom_lines, ['TT-EE-TE', 'TB-EB'], loc=4)
    ax.legend(prop={'size': 17}, loc='upper left', ncol=3, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_xlim(2.0, clexp['lMaxP'])
    ax.set_ylim(0.8, 1.05)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    # ax.set_ylabel(r'$N_{mv}^\mathrm{GMV, HO}/N_{mv}^\mathrm{HO, SMV}$',
    #              fontsize=22)
    ax.set_ylabel(r'$N_{mv}^\mathrm{GMV , HO}/N_{mv}^\mathrm{SQE}$',
                  fontsize=22)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=22)
    # ax.legend(custom_lines, ['TT-EE-TE', 'TB-EB'])
    # ax.add_artist(leg1)
    # ax.add_artist(leg2)
    # np.save('figures/multexp_vartrueind_varHO02ind.png')
    # """


def plot_mv_multest_multexp(exp):
    dat = np.genfromtxt("input/CAMB/Julien_lenspotentialCls.dat")

    lines = ["-", "--", "-.", ":"]
    cl = ["b", "r", "g"]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(dat[:, 0], dat[:, 5], 'k-', lw=4, label=r'signal')

    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]
        interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        interp_HO = interp1d(data2[:, 0], data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data3 = np.genfromtxt('output/SQE_variance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L3 = data3[:, 0]
        interp_smv = interp1d(L3, data3[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data4 = np.genfromtxt('output/OH03_flatsky_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L4 = data4[:, 0]
        interp_oh03 = interp1d(L4, data4[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        # L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
        L_p = np.linspace(2, 5000, 3*4999)
        Dlp = L_p*(L_p+1)/2/np.pi

        ax.plot(L_p, Dlp*interp_tmv(L_p), cl[e], ls=lines[0],
                label='%s GMV' % (clexp['name']))
        ax.plot(L_p, Dlp*interp_HO(L_p), cl[e], ls=lines[1],
                label='%s HO' % (clexp['name']))
        ax.plot(L_p, Dlp*interp_smv(L_p), cl[e], ls=lines[2],
                label='%s SMV' % (clexp['name']))
        ax.plot(L_p, Dlp*interp_oh03(L_p), cl[e], ls=lines[3],
                label='%s OH03' % (clexp['name']))

    # ls = ax.get_lines()
    # leg1 = plt.legend([ls[i] for i in [0, 1]], [exp[i]['name'] for i in range(len(exp))], loc=1)
    # leg2 = plt.legend(custom_lines, ['TT-EE-TE', 'TB-EB'], loc=4)
    ax.legend(prop={'size': 17}, loc='upper left', ncol=3, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2.0, clexp['lMaxP'])
    ax.set_ylim(ymin=1.e-9)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax.set_ylabel(r'$L(L+1)N_L^{dd}/2\pi$',
                  fontsize=22)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=22)


def plot_OH03_TE_comparison(exp):
    lines = ["-", "--", "-."]
    cl = ["b", "r", "g"]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('output/HO02_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]

        interp_HO02_TE = interp1d(L, data[:, 3], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/OH03_flatsky_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L2 = data2[:, 0]

        interp_OH03_TE = interp1d(L2, data2[:, 10], kind='quadratic', bounds_error=False, fill_value=0.)

        # L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
        L_p = np.linspace(2, 5000, 4999)

        fracdiff = (interp_OH03_TE(L_p)-interp_HO02_TE(L_p))/interp_HO02_TE(L_p)
        # print fracdiff[:10]

        ax.plot(L_p, fracdiff, cl[e], ls=lines[0],
                label='%s' % (clexp['name']))
    # ls = ax.get_lines()
    # leg1 = plt.legend([ls[i] for i in [0, 1]], [exp[i]['name'] for i in range(len(exp))], loc=1)
    # leg2 = plt.legend(custom_lines, ['TT-EE-TE', 'TB-EB'], loc=4)
    ax.legend(prop={'size': 17}, loc='upper left', ncol=1, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_xlim(2.0, clexp['lMaxP'])
    # ax.set_ylim(0.8, 1.05)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax.set_ylabel(r'$(N_{L}^\mathrm{TE, OH03} - N_{L}^\mathrm{TE, HO02})/N_{L}^\mathrm{TE, HO02}$',
                  fontsize=16)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=16)


def plot_OH03_MV_comparison(exp):
    lines = ["-", "--", "-."]
    cl = ["b", "r", "g"]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]
        interp_HO02_mv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/OH03_flatsky_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L2 = data2[:, 0]

        interp_OH03_mv = interp1d(L2, data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        # L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
        L_p = np.linspace(2, 5000, 4999)

        fracdiff = (interp_OH03_mv(L_p)-interp_HO02_mv(L_p))/interp_HO02_mv(L_p)

        ax.plot(L_p, fracdiff, cl[e], ls=lines[0],
                label='%s' % (clexp['name']))
    # ls = ax.get_lines()
    # leg1 = plt.legend([ls[i] for i in [0, 1]], [exp[i]['name'] for i in range(len(exp))], loc=1)
    # leg2 = plt.legend(custom_lines, ['TT-EE-TE', 'TB-EB'], loc=4)
    ax.hlines(0., 2.0, clexp['lMaxP'], 'k')
    ax.legend(prop={'size': 17}, loc='upper right', ncol=1, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_xlim(2.0, clexp['lMaxP'])
    # ax.set_ylim(0.8, 1.05)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax.set_ylabel(r'$(N_{L}^\mathrm{MV, OH03} - N_{L}^\mathrm{MV, HO02})/N_{L}^\mathrm{MV, HO02}$',
                  fontsize=16)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=16)


def plot_OH03_MV_comparison_per(exp):
    lines = ["-", "--", "-."]
    cl = ["b", "r", "g"]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        data200 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_200.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L = data[:, 0]
        # interp_HO02_mv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
        # interp_HO02_mv_200 = interp1d(L, data200[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/OH03_flatsky_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        data2_200 = np.genfromtxt('output/OH03_flatsky_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_200.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
        L2 = data2[:, 0]

        # interp_OH03_mv = interp1d(L2, data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
        # interp_OH03_mv_200 = interp1d(L2, data2_200[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        L_p = L
        # L_p = np.linspace(2, 5000, 4999)

        # fracdiff = (data2[:, -1]-data[:, -1])/data[:, -1]
        # fracdiff2 = (data2_200[:, -1]-data200[:, -1])/data200[:, -1]

        fracdiff = (data2[:, -1]-data2_200[:, -1])*100/data2_200[:, -1]
        fracdiff2 = (data[:, -1]-data200[:, -1])*100/data200[:, -1]

        per = (fracdiff - fracdiff2)*100/fracdiff2
        ax.plot(L_p, fracdiff, cl[e], ls=lines[0],
                label='%s' % (clexp['name']))
        ax.plot(L_p, fracdiff2, cl[e], ls=lines[1],
                label='%s' % (clexp['name']))
        # ax.plot(L_p, per, cl[e], ls=lines[0],
        #         label='%s' % (clexp['name']))
    # ls = ax.get_lines()
    # leg1 = plt.legend([ls[i] for i in [0, 1]], [exp[i]['name'] for i in range(len(exp))], loc=1)
    # leg2 = plt.legend(custom_lines, ['TT-EE-TE', 'TB-EB'], loc=4)
    ax.hlines(0., 2.0, clexp['lMaxP'], 'k')
    ax.legend(prop={'size': 17}, loc='upper right', ncol=1, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_xlim(2.0, clexp['lMaxP'])
    # ax.set_ylim(0.8, 1.05)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax.set_ylabel(r'$\% $ change between 200 and 400 steps',
                  fontsize=16)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=16)



def plot_n0_n1_multexp(exp):
    # dat = np.genfromtxt("input/CAMB/Julien_lenspotentialCls.dat")
    dat = np.genfromtxt("Julien_N1/Julien_lenspotentialCls.dat")

    lines = ["-", "--", "-.", ":"]
    cl = ["b", "r", "g"]
    # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    L = np.linspace(0, 3000, 3001)
    Dl = (L*(L+1))**2/2/np.pi

    ax.plot(dat[:, 0], dat[:, 5], 'k-', lw=4, label=r'signal')

    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('Julien_N1/N0sN1s_%s.dat' % (clexp['name']))
        # ax.plot(L, Dl*data[:, 0], cl[e], ls=lines[0],
        #         label=r'%s $N^0$ SQE' % (clexp['name']))
        ax.plot(L, Dl*data[:, 1], cl[e], ls=lines[1],
                label=r'%s $N^{(0)}_L$ GMV' % (clexp['name']))
        # ax.plot(L, Dl*data[:, 2], cl[e], ls=lines[2],
        #         label=r'%s $N^1$ SQE' % (clexp['name']))
        ax.plot(L, Dl*data[:, 3], cl[e], ls=lines[3],
                label=r'%s $N^{(1)}_L$ GMV' % (clexp['name']))

    ax.legend(prop={'size': 17}, loc='upper center', ncol=3, frameon=False,
              labelspacing=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2.0, clexp['lMaxP'])
    # ax.set_ylim(ymin=1.e-9)
    # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax.set_ylabel(r'$[L(L+1)]^2 N_L/2\pi$',
                  fontsize=22)
    ax.set_xlabel(r'$L$', fontsize=22)
    ax.tick_params(axis='both', labelsize=22)


def plot_n1_gmv_sqe_multexp(exp):
    lines = ["-", "--", "-.", ":"]
    cl = ["b", "r", "g"]

    fig = plt.figure()
    ax2 = fig.add_subplot(111)

    L = np.linspace(0, 3000, 3001)

    for e in range(len(exp)):
        clexp = exp[e]
        data = np.genfromtxt('Julien_N1/N0sN1s_%s.dat' % (clexp['name']))
        gmvn1 = data[:, 3]
        sqen1 = data[:, 2]
        ratio = sqen1/gmvn1
        # ax.plot(L, Dl*data[:, 0], cl[e], ls=lines[0],
        #         label=r'%s $N^0$ SQE' % (clexp['name']))
        ax2.plot(L, ratio, cl[e], label=r'%s' % (clexp['name']))

    ax2.legend(prop={'size': 17}, loc='upper center', ncol=3, frameon=False,
               labelspacing=0.2)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(2.0, clexp['lMaxP'])
    # ax2.set_ylim(ymin=1.e-9)
    # ax2.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    ax2.set_ylabel(r'$N^{(1) (\mathrm{SQE})}_L/N^{(1) \mathrm{GMV}}_L$', fontsize=22)
    ax2.set_xlabel(r'$L$', fontsize=22)
    ax2.tick_params(axis='both', labelsize=22)
