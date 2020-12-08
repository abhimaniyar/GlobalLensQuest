from headers import *


class Cell_cmb(object):

    def __init__(self, exp):

        self.exp = exp
        self.name = self.exp['name']
        #  beam fwhm in radians
        self.fwhm = self.exp['beam'] * (np.pi/180.)/60.
        #  detector sensitivity in muK*rad.
        self.sensitivity_t = self.exp['noise_t'] * (np.pi/180.)/60.
        self.sensitivity_p = self.exp['noise_p'] * (np.pi/180.)/60.
        # ell limits
        self.lMin = self.exp['lMin']
        self.lMaxT = self.exp['lMaxT']
        self.lMaxP = self.exp['lMaxP']

        # self.Tcmb = 2.726e6  # muK
        self.Tcmb = 1.
        #  conversion from Dl to Cl: Dl = l(l+1) Cl / 2pi
        self.dl_to_cl = lambda l: 2.*np.pi/(l*(l+1.))
        # cmb_out = 7.4311e12
        cmb_out = 1.
        # reading the power spectra
        # unlensed TT, EE, TE
        # nolens = np.loadtxt('../CAMB/qe_lenspotentialCls.dat')
        nolens = np.loadtxt('input/CAMB/qe_lenspotentialCls.dat')
        fac = self.dl_to_cl(nolens[:, 0])
        nolens[:, 1] *= fac/cmb_out
        nolens[:, 2] *= fac/cmb_out
        nolens[:, 3] *= fac/cmb_out
        nolens[:, 4] *= fac/cmb_out

        # corrul = 0.8*np.sqrt(nolens[:, 1]*nolens[:, 2])
        # interpolate
        self.unlensedTT = interp1d(nolens[:, 0], nolens[:, 1], kind='linear',
                                   bounds_error=False, fill_value=0.)
        self.unlensedEE = interp1d(nolens[:, 0], nolens[:, 2], kind='linear',
                                   bounds_error=False, fill_value=0.)
        self.unlensedBB = interp1d(nolens[:, 0], nolens[:, 3], kind='linear',
                                   bounds_error=False, fill_value=0.)
        self.unlensedTE = interp1d(nolens[:, 0], nolens[:, 4], kind='linear',
                                   bounds_error=False, fill_value=0.)
        # self.unlensedTE = interp1d(nolens[:, 0], corrul, kind='linear',
        #                            bounds_error=False, fill_value=0.)

        # lensed TT, EE, TE
        # lens = np.loadtxt('../CAMB/qe_lensedCls.dat')
        lens = np.loadtxt('input/CAMB/qe_lensedCls.dat')
        fac = self.dl_to_cl(lens[:, 0])
        lens[:, 1] *= fac/cmb_out
        lens[:, 2] *= fac/cmb_out
        lens[:, 3] *= fac/cmb_out
        lens[:, 4] *= fac/cmb_out

        # corrl = 0.8*np.sqrt(lens[:, 1]*lens[:, 2])
        # interpolate
        self.lensedTT = interp1d(lens[:, 0], lens[:, 1], kind='linear',
                                 bounds_error=False, fill_value=0.)
        self.lensedEE = interp1d(lens[:, 0], lens[:, 2], kind='linear',
                                 bounds_error=False, fill_value=0.)
        self.lensedBB = interp1d(lens[:, 0], lens[:, 3], kind='linear',
                                 bounds_error=False, fill_value=0.)
        self.lensedTE = interp1d(lens[:, 0], lens[:, 4], kind='linear',
                                 bounds_error=False, fill_value=0.)
        # self.lensedTE = interp1d(lens[:, 0], corrl, kind='linear',
        #                          bounds_error=False, fill_value=0.)

        # total lensed : lens+noise
        print 'calculating total power spectra'
        self.totalTT = lambda l: self.lensedTT(l) + self.detectorNoise(l, self.sensitivity_t) + self.artificialNoiseTT(l)
        self.totalEE = lambda l: self.lensedEE(l) + self.detectorNoise(l, self.sensitivity_p)
        self.totalBB = lambda l: self.lensedBB(l) + self.detectorNoise(l, self.sensitivity_p)
        self.totalTE = lambda l: self.lensedTE(l)

    def detectorNoise(self, l, sensitivity):
        sigma_beam = self.fwhm / np.sqrt(8.*np.log(2.))
        a = l*(l+1)*sigma_beam**2
        b = (sensitivity/self.Tcmb)**2
        return b*np.exp(a)  # *7.4311e12

    def artificialNoiseTT(self, l):
        """
        for lmaxT != lmaxP, we add artificial noise on the TT power spectra
        for ell > lmaxT. This way both TT and TE or TB and EE are calculated
        to same ell value, however, TT is dominated by this artificial noise
        and thus dies not really contrubute to the signal. This might bias
        the estimator though. So keep this in mind. For now, just for
        calculating the noise, this should be fine.
        """
        noise = np.zeros(len(l))
        if self.lMaxT != self.lMaxP:
            idx = np.where((l > self.lMaxT))[0]
            noise[idx] = 1.e6
        return noise

    def plot_cell(self):
        nolens = np.loadtxt('input/CAMB/qe_lenspotentialCls.dat')
        # np.loadtxt('../CAMB/qe_nolens_scalCls.dat')
        ell = nolens[:, 0]
        noise_t = self.detectorNoise(ell, self.sensitivity_t)
        noise_p = self.detectorNoise(ell, self.sensitivity_p)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # """
        ax.plot(ell, ell*(ell+1)*self.unlensedTT(ell)/(2*np.pi), 'k', lw=1.5,
                label=r'TT')
        ax.plot(ell, ell*(ell+1)*noise_t/(2*np.pi), 'k--', lw=1.5,
                label=r'noise TT')
        ax.plot(ell, ell*(ell+1)*self.unlensedEE(ell)/(2*np.pi), 'b', lw=1.5,
                label=r'EE')
        ax.plot(ell, ell*(ell+1)*noise_p/(2*np.pi), 'b--', lw=1.5,
                label=r'noise EE')
        ax.plot(ell, ell*(ell+1)*self.unlensedTE(ell)/(2*np.pi), 'r', lw=1.5,
                label=r'TE')
        """
        ax.plot(ell, ell*(ell+1)*self.unlensedBB(ell)/(2*np.pi), 'g', lw=1.5,
                label=r'BB')
        """
        """
        lens2 = np.loadtxt('../CAMB/manu_lenspotentialCls.dat')
        ell2 = lens2[:, 0]
        noise_t = self.detectorNoise(ell2, self.sensitivity_t)
        noise_p = self.detectorNoise(ell2, self.sensitivity_p)
        ax.plot(ell2, lens2[:, 1]/7.4311e12, 'g', label='manu TT')
        ax.plot(ell2, lens2[:, 2]/7.4311e12, 'g', label='manu EE')
        ax.plot(ell2, lens2[:, 4]/7.4311e12, 'y', label='manu TE')
        """
        ax.legend(loc=2, fontsize='8')  # , labelspacing=0.1)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        # ax.set_ylim((1.e-15, 1e-9))
        ax.set_xlabel(r'$\ell$', fontsize=16)
        ax.set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi$', fontsize=16)
        plt.show()
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(ell2, ell2*(ell2+1)*self.totalTT(ell2)/(2*np.pi), 'k', lw=1.5,
                label=r'TT')
        # ax.plot(ell2, ell2*(ell2+1)*noise_t/(2*np.pi), 'k--', lw=1.5,
        #         label=r'noise TT')
        ax.plot(ell2, ell2*(ell2+1)*self.totalEE(ell2)/(2*np.pi), 'b', lw=1.5,
                label=r'EE')
        # ax.plot(ell, ell2*(ell2+1)*noise_p/(2*np.pi), 'b--', lw=1.5,
        #         label=r'noise EE')
        ax.plot(ell2, ell2*(ell2+1)*self.totalTE(ell2)/(2*np.pi), 'r', lw=1.5,
                label=r'TE')
        """
        """
        ax.plot(ell, ell*(ell+1)*self.unlensedBB(ell)/(2*np.pi), 'g', lw=1.5,
                label=r'BB')
        """
        """
        ax.legend(loc=2, fontsize='8')  # , labelspacing=0.1)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        # ax.set_ylim((1.e-15, 1e-9))
        ax.set_xlabel(r'$\ell$', fontsize=16)
        ax.set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi$', fontsize=16)
        plt.show()
        """

    def plot_cell_2(self):
        lens = np.loadtxt('input/CAMB/qe_lensedCls.dat')
        ell = lens[:, 0]
        noise_t = self.detectorNoise(ell, self.sensitivity_t)
        noise_p = self.detectorNoise(ell, self.sensitivity_p)
        """
        fig = plt.figure(13)
        ax = fig.add_subplot(111)
        ax.plot(ell, ell*(ell+1)*self.unlensedTT(ell)/(2*np.pi), 'k', lw=1.5,
                label=r'TT')
        ax.plot(ell, ell*(ell+1)*noise_t/(2*np.pi), 'k--', lw=1.5,
                label=r'noise TT')
        ax.plot(ell, ell*(ell+1)*self.unlensedEE(ell)/(2*np.pi), 'b', lw=1.5,
                label=r'EE')
        ax.plot(ell, ell*(ell+1)*noise_p/(2*np.pi), 'b--', lw=1.5,
                label=r'noise EE')
        ax.plot(ell, ell*(ell+1)*self.unlensedTE(ell)/(2*np.pi), 'r', lw=1.5,
                label=r'TE')
        ax.legend(loc=2, fontsize='8')  # , labelspacing=0.1)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        # ax.set_ylim((1.e-15, 1e-9))
        ax.set_xlabel(r'$\ell$', fontsize=16)
        ax.set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi$', fontsize=16)
        plt.show()
        # """
        """
        ax.plot(ell, ell*(ell+1)*self.unlensedBB(ell)/(2*np.pi), 'g', lw=1.5,
                label=r'BB')
        """

        ell2 = ell
        noise_t = self.detectorNoise(ell2, self.sensitivity_t)
        noise_p = self.detectorNoise(ell2, self.sensitivity_p)
        fig = plt.figure(14)
        ax = fig.add_subplot(111)
        # """
        print self.lensedTT(np.array([100., 500., 1000., 2000., 3500., 5000.]))
        print self.totalTT(np.array([100., 500., 1000., 2000., 3500., 5000.]))
        ax.plot(ell2, ell2*(ell2+1)*self.totalTT(ell2)/(2*np.pi), 'k', lw=1.5,
                label=r'TT')
        # ax.plot(ell2, ell2*(ell2+1)*noise_t/(2*np.pi), 'k--', lw=1.5,
        #         label=r'noise TT')
        ax.plot(ell2, ell2*(ell2+1)*self.totalEE(ell2)/(2*np.pi), 'b', lw=1.5,
                label=r'EE')
        # ax.plot(ell, ell2*(ell2+1)*noise_p/(2*np.pi), 'b--', lw=1.5,
        #         label=r'noise EE')
        ax.plot(ell2, ell2*(ell2+1)*self.totalTE(ell2)/(2*np.pi), 'r', lw=1.5,
                label=r'TE')
        """
        ax.plot(ell, ell*(ell+1)*self.unlensedBB(ell)/(2*np.pi), 'g', lw=1.5,
                label=r'BB')
        """
        ax.legend(loc=2, fontsize='8')  # , labelspacing=0.1)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        # ax.set_ylim((1.e-15, 1e-9))
        ax.set_xlabel(r'$\ell$', fontsize=16)
        ax.set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi$', fontsize=16)
        plt.show()