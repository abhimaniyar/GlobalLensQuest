from imports import *


class lensing_estimator(object):

    def __init__(self, Cell_cmb):
        self.cmb = Cell_cmb
        self.name = self.cmb.name
        self.beam = self.cmb.exp['beam']
        self.noise = self.cmb.exp['noise_t']

        """
        bounds for ell integrals
        l_1 + l_2 = L
        """
        self.l1Min = self.cmb.lMin
        # max value for l1 and l2 is taken to be same
        self.l1Max = max(self.cmb.lMaxT, self.cmb.lMaxP)

        # L = l_1 + l_2. This L is for reconstructed phi field
        # self.L = np.logspace(np.log10(1.), np.log10(2.*self.l1Max+1.), 51, 10.)
        # a1 = np.logspace(np.log10(1.), np.log10(100.), 20, 10.)
        # a2 = np.logspace(np.log10(110.), np.log10(1500.), 140, 10.)
        # a3 = np.logspace(np.log10(1600.), np.log10(2*self.l1Max+1.), 51, 10.)
        # self.L = np.concatenate((a1, a2, a3))
        self.L = np.logspace(np.log10(1.), np.log10(2*self.l1Max+1.), 201, 10.)
        # self.L = np.linspace(1., 201., 1001)
        self.Nl = len(self.L)
        self.N_phi = 50  # number of steps for angular integration steps
        # reduce to 50 if you need around 0.6% max accuracy till L = 3000
        # from 200 to 400, there is just 0.03% change in the noise curves till L=3000
        self.var_out = 'output/HO02_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise), str(self.N_phi))
        self.covar_out = 'output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise), str(self.N_phi))

    """
    L = l1 + l2
    phi1 = angle betweeen vectors (L, l_1)
    phi2 = angle betweeen vectors (L, l_2)
    and phi12 = phi1 - phi2
    """

    def l2(self, L, l_1, phi1):
        """
        mod of l2 = (L-1_1) given phi1
        """
        return np.sqrt(L**2 + l_1**2 - 2*L*l_1*np.cos(phi1))

    def phi12(self, L, l_1, phi1):
        """
        phi12 = phi1 - phi2
        """
        x = L*np.cos(-phi1) - l_1
        y = L*np.sin(-phi1)
        result = -np.arctan2(y, x)
        return result

    def phi2(self, L, l_1, phi1):
        """
        phi2 = phi1 - phi12
        """
        result = phi1 - self.phi12(L, l_1, phi1)
        # result = phi1 + self.phi12(L, l_1, phi1)
        return result

    def l1max(self, XY):
        """
        max value for l1 and l2: taken to be same
        """
        if XY == 'TT':
            # return self.cmb.lMaxT
            return self.cmb.lMaxP
        elif XY == 'EE' or XY == 'BB' or XY == 'EB':
            return self.cmb.lMaxP
        elif XY == 'TE' or XY == 'TB':
            # taking the minimum of the two: this approach is suboptimal
            return max(self.cmb.lMaxT, self.cmb.lMaxP)

    def f_XY(self, L, l_1, phi1, XY):
        """
        lensing response such that
        <X_l1 Y_{L-l1}> = f_XY(l1, L-l1)*\phi_L.
        """

        l_2 = self.l2(L, l_1, phi1)
        phi12 = self.phi12(L, l_1, phi1)
        phi2 = self.phi2(L, l_1, phi1)

        Ldotl_1 = L*l_1*np.cos(phi1)
        Ldotl_2 = L*l_2*np.cos(phi2)
        # """
        if XY == 'TT':
            result = self.cmb.lensgradTT(l_1)*Ldotl_1
            result += self.cmb.lensgradTT(l_2)*Ldotl_2
            # print result
            # sys.exit()
        elif XY == 'EE':
            result = self.cmb.lensgradEE(l_1)*Ldotl_1
            result += self.cmb.lensgradEE(l_2)*Ldotl_2
            result *= np.cos(2.*phi12)
        elif XY == 'TE':
            # there is a typo in HO02!!!!!!!!!
            # instead of cos(phi12) it should be cos(2*phi12)!!!!!
            result = self.cmb.lensgradTE(l_1)*np.cos(2.*phi12)*Ldotl_1
            result += self.cmb.lensgradTE(l_2)*Ldotl_2
        elif XY == 'TB':
            result = self.cmb.lensgradTE(l_1)*np.sin(2.*phi12)*Ldotl_1
        elif XY == 'EB':
            # there is a typo in HO02!!!!!!!!!
            # instead of - it should be + between first and second term!!!!!
            result = self.cmb.lensgradEE(l_1)*Ldotl_1
            result += self.cmb.lensgradBB(l_2)*Ldotl_2
            result *= np.sin(2.*phi12)
        elif XY == 'BB':
            result = self.cmb.lensgradBB(l_1)*Ldotl_1
            result += self.cmb.lensgradBB(l_2)*Ldotl_2
            result *= np.cos(2.*phi12)
        """
        if XY == 'TT':
            result = self.cmb.lensedTT(l_1)*Ldotl_1
            result += self.cmb.lensedTT(l_2)*Ldotl_2
            # print result
            # sys.exit()
        elif XY == 'EE':
            result = self.cmb.lensedEE(l_1)*Ldotl_1
            result += self.cmb.lensedEE(l_2)*Ldotl_2
            result *= np.cos(2.*phi12)
        elif XY == 'TE':
            # there is a typo in HO02!!!!!!!!!
            # instead of cos(phi12) it should be cos(2*phi12)!!!!!
            result = self.cmb.lensedTE(l_1)*np.cos(2.*phi12)*Ldotl_1
            result += self.cmb.lensedTE(l_2)*Ldotl_2
        elif XY == 'TB':
            result = self.cmb.lensedTE(l_1)*np.sin(2.*phi12)*Ldotl_1
        elif XY == 'EB':
            # there is a typo in HO02!!!!!!!!!
            # instead of - it should be + between first and second term!!!!!
            result = self.cmb.lensedEE(l_1)*Ldotl_1
            result += self.cmb.lensedBB(l_2)*Ldotl_2
            result *= np.sin(2.*phi12)
        elif XY == 'BB':
            result = self.cmb.lensedBB(l_1)*Ldotl_1
            result += self.cmb.lensedBB(l_2)*Ldotl_2
            result *= np.cos(2.*phi12)
        # """
        # result *= 2. / L**2

        return result

    def F_XY(self, L, l_1, phi1, XY):
        """
        Weighing terms for the estimator.
        This decides the weights for a corresponding pair of multipoles for
        X and Y.
        """

        l_2 = self.l2(L, l_1, phi1)
        phi2 = self.phi2(L, l_1, phi1)

        if XY == 'TT':
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = 2.*self.cmb.totalTT(l_1)*self.cmb.totalTT(l_2)
            result = numerator/denominator
            print (result)
        elif XY == 'EE':
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = 2.*self.cmb.totalEE(l_1)*self.cmb.totalEE(l_2)
            result = numerator/denominator
        elif XY == 'TE':
            numerator = self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2)*self.f_XY(L, l_1, phi1, XY)
            numerator -= self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)*self.f_XY(L, l_2, phi2, XY)
            denominator = (self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) *
                           self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2))
            denominator -= (self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2))**2
            result = numerator/denominator
        elif XY == 'TB':
            # TB power spectrum assumed zero
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = self.cmb.totalTT(l_1)*self.cmb.totalBB(l_2)
            result = numerator/denominator
        elif XY == 'EB':
            # EB assumed zero
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2)
            result = numerator/denominator
        elif XY == 'BB':
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = 2.*self.cmb.totalBB(l_1)*self.cmb.totalBB(l_2)
            result = numerator/denominator
        return result

    def var_individual(self, L, XY):
        """
        Variance of the QE for individual XY. Same as the
        normalization of the \phi_XY estimator. Such that norm = 1/(\int F*f)
        such that \phi_XY = norm* \int F*(XY)
        """
        l1min = self.l1Min
        l1max = self.l1max(XY)
        # """
        if L > 2.*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.
        # """

        def integrand(l_1, phil):

            l_2 = self.l2(L, l_1, phil)
            """
            if l_1 < l1min or l_2 < l1min or l_1 > l1max or l_2 > l1max:
                return 0.
            """
            result = self.f_XY(L, l_1, phil, XY)*self.F_XY(L, l_1, phil, XY)

            result *= 2*l_1  # **2
            # d^2l_1 = dl_1*l_1*dphi1
            """factor of 2 above because phi integral is symmetric. Thus we've
            put instead of 0 to 2pi, 2 times 0 to pi
            Also, l_1^2 instead of l_1 if we are taking log spacing for
            l_1"""
            result /= (2.*np.pi)**2
            # """
            # idx = np.where((l_2 < l1min) | (l_2 > l1max))
            idx = np.where((l_1 < l1min) | (l_1 > l1max) | (l_2 < l1min) | (l_2 > l1max))

            result[idx] = 0.
            # """
            return result


        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        # l1 = np.logspace(np.log10(l1min), np.log10(l1max), int(l1max-l1min+1))
        phi1 = np.linspace(0., np.pi, self.N_phi)
        int_1 = np.zeros(len(phi1))
        for i in range(len(phi1)):
            intgnd = integrand(l1, phi1[i])
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')
        int_ll = integrate.simps(int_1, x=phi1, even='avg')
        result = 1./int_ll

        result *= L**2

        if not np.isfinite(result):
            result = 0.

        if result < 0.:
            print(L)

        return result

    def calc_var(self, est):
        data = np.zeros((self.Nl, 6))
        data[:, 0] = np.copy(self.L)
        pool = Pool(ncpus=4)
        n_est = len(est)
        for i_est in range(n_est):
            XY = est[i_est]
            print("Computing variance for " + XY)

            def f(l):
                return self.var_individual(l, XY)

            data[:, i_est+1] = np.array(pool.map(f, self.L))

        if est == ['TT', 'EE', 'TE']:
            np.savetxt('output/HO02_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only_%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise), str(self.N_phi)), data)
        elif est == ['TB', 'EB']:
            np.savetxt('output/HO02_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only_%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise), str(self.N_phi)), data)
        else:
            np.savetxt(self.var_out, data)

    def interp_var(self, est):

        print("Loading variances")

        self.var_d = {}
        data = np.genfromtxt(self.var_out)
        L = data[:, 0]

        n_est = len(est)
        for i_est in range(n_est):
            XY = est[i_est]
            var = data[:, i_est+1].copy()
            self.var_d[XY] = interp1d(L, var, kind='linear',
                                      bounds_error=False, fill_value=0.)

    def plot_var(self, est):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        data = np.genfromtxt("input/CAMB/Julien_lenspotentialCls.dat")
        L = data[:, 0]
        ax.plot(L, data[:, 5], 'r-', lw=1.5, label=r'signal')

        n_est = len(est)
        for i_est in range(n_est):
            XY = est[i_est]
            ax.plot(self.L, self.L*(self.L+1)*self.var_d[XY](self.L)/(2*np.pi), c=plt.cm.rainbow(i_est/6.), lw=1.5, label=XY)
            # ax.plot(self.L, self.var_d[XY](self.L), c=plt.cm.rainbow(i_est/6.), lw=1.5, label=XY)

        ax.legend(loc=2, fontsize='8')
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.set_ylabel(r'$L(L+1)C_L^{dd}/2\pi$', fontsize=16)
        ax.set_ylim((3.e-11, 0.1))
        ax.set_xlim((2., 4.e4))
        plt.show()

    def covariance(self, L, XY, AB):
        """
        Covariance of the QE for XY and AB map choice. Because the unlensed
        TB and EB power spectra are zero, they do not have any covariance
        with other map combinations.
        """
        l1min = self.l1Min
        l1max = min(self.l1max(XY), self.l1max(AB))

        if L > 2*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.

        if (XY == 'TT')*(AB == 'TB') or (XY == 'TT')*(AB == 'EB') or \
            (XY == 'EE')*(AB == 'TB') or (XY == 'EE')*(AB == 'EB') or \
        (XY == 'TE')*(AB == 'TB') or (XY == 'TE')*(AB == 'EB'):
            return 0.

        def integrand(l_1, phil):

            l_2 = self.l2(L, l_1, phil)
            phi2 = self.phi2(L, l_1, phil)

            """
            if l_1 < l1min or l_2 < l1min or l_1 > l1max or l_2 > l1max:
                return 0.
            # """

            if XY == 'TT' and AB == 'EE':
                a = self.F_XY(L, l_1, phil, AB)*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)
                a += self.F_XY(L, l_2, phi2, AB)*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

            elif XY == 'TT' and AB == 'TE':
                a = self.F_XY(L, l_1, phil, AB)*self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2)
                a += self.F_XY(L, l_2, phi2, AB)*self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2)

            elif XY == 'EE' and AB == 'TE':
                a = self.F_XY(L, l_1, phil, AB)*self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2)
                a += self.F_XY(L, l_2, phi2, AB)*self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2)

            elif XY == 'TB' and AB == 'EB':
                a = self.F_XY(L, l_1, phil, AB)*self.cmb.totalTE(l_1)*self.cmb.totalBB(l_2)

            result = a*self.F_XY(L, l_1, phil, XY)
            result *= 2*l_1  # **2
            # d^2l_1 = dl_1*l_1*dphi1
            """factor of 2 above because phi integral is symmetric. Thus we've
            put instead of 0 to 2pi, 2 times 0 to pi
            Also, l_1^2 instead of l_1 if we are taking log spacing for
            l_1"""
            result /= (2.*np.pi)**2

            # idx = np.where((l_2 < l1min) | (l_2 > l1max))
            idx = np.where((l_1 < l1min) | (l_1 > l1max) | (l_2 < l1min) | (l_2 > l1max))

            result[idx] = 0.
            return result

        # """        
        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        phi1 = np.linspace(0., np.pi, self.N_phi)
        int_1 = np.zeros(len(phi1))
        for i in range(len(phi1)):
            intgnd = integrand(l1, phi1[i])
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')
        int_l1 = integrate.simps(int_1, x=phi1, even='avg')
        result = int_l1
        # """

        result *= self.var_d[XY](L)*self.var_d[AB](L)
        result *= 1./L**2

        if not np.isfinite(result):
            result = 0.
        return result

    def calc_cov(self, est):
        data = np.zeros((self.Nl, 12))
        data[:, 0] = np.copy(self.L)
        pool = Pool(ncpus=4)

        cov_XY_AB = {}

        n_est = len(est)
        counter = 1
        for i_est in range(n_est):
            XY = est[i_est]
            cov_XY_AB[XY+XY] = self.var_d[XY](self.L)
            for i2 in range(i_est+1, n_est):
                AB = est[i2]
                print("Computing covariance for " + XY + "-" + AB)

                def f(l):
                    return self.covariance(l, XY, AB)

                cov_XY_AB[XY+AB] = np.array(pool.map(f, self.L))
                data[:, counter] = cov_XY_AB[XY+AB]
                counter += 1

        # min variance estimator noise
        n_mv = np.zeros(self.Nl)
        for el in range(self.Nl):
            covmat = np.zeros((n_est, n_est))
            for i_est in range(n_est):
                XY = est[i_est]
                covmat[i_est, i_est] = self.var_d[XY](self.L[el])
                for i2 in range(i_est+1, n_est):
                    AB = est[i2]
                    covmat[i_est, i2] = covmat[i2, i_est] = cov_XY_AB[XY+AB][el]
            # invert the matrix
            try:
                invcov = np.linalg.inv(covmat)
                n_mv[el] = 1./np.sum(invcov)
                # np.savetxt('covmat.txt', covmat)
            except:
                print("exception while inverting the covariance matrix at L = %s !" % str(el))
                pass

        data[:, -1] = n_mv

        if est == ['TT', 'EE', 'TE']:
            np.savetxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only_%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise), str(self.N_phi)), data)
        elif est == ['TB', 'EB']:
            np.savetxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only_%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise), str(self.N_phi)), data)
        else:
            np.savetxt(self.covar_out, data)

    def interp_cov(self, est):
        print("Interpolating covariances")

        self.cov_d = {}
        data = np.genfromtxt(self.covar_out)
        L = data[:, 0]

        n_est = len(est)
        counter = 1
        for i_est in range(n_est):
            XY = est[i_est]
            for i2 in range(i_est+1, n_est):
                AB = est[i2]
                norm = data[:, counter].copy()
                self.cov_d[XY+AB] = interp1d(L, norm, kind='linear', bounds_error=False, fill_value=0.)
                counter += 1

        nmv = data[:, -1]
        self.var_d['mv'] = interp1d(L, nmv, kind='linear', bounds_error=False, fill_value=0.)

    def plot_cov(self, est):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        data2 = np.genfromtxt("input/CAMB/Julien_lenspotentialCls.dat")
        L = data2[:, 0]
        ax.plot(L, data2[:, 5], 'r-', lw=1.5, label=r'signal')

        n_est = len(est)
        for i_est in range(n_est):
            XY = est[i_est]
            ax.plot(self.L, self.L*(self.L+1)*self.var_d[XY](self.L)/(2*np.pi), c=plt.cm.rainbow(i_est/6.), lw=1.5, label=XY)

        ax.plot(self.L, self.L*(self.L+1)*self.var_d['mv'](self.L)/(2*np.pi), 'k', lw=1.5, label='min var')

        ax.legend(loc=2, fontsize='12')  # , labelspacing=0.1)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.set_ylabel(r'$L(L+1)C_L^{dd}/2\pi$', fontsize=16)
        ax.set_ylim((3.e-11, 0.1))
        ax.set_xlim((2., 4.e4))
        ax.tick_params(axis='both', labelsize=12)
        plt.show()

    def plot_cov_pair(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.L, self.L*(self.L+1)*np.abs(self.cov_d['TTEE'](self.L))/(2*np.pi), label='TT-EE')
        ax.plot(self.L, self.L*(self.L+1)*np.abs(self.cov_d['TTTE'](self.L))/(2*np.pi), label='TT-TE')
        ax.plot(self.L, self.L*(self.L+1)*np.abs(self.cov_d['EETE'](self.L))/(2*np.pi), label='EE-TE')
        ax.plot(self.L, self.L*(self.L+1)*np.abs(self.cov_d['TBEB'](self.L))/(2*np.pi), label='TB-EB')
        ax.legend(loc=2, fontsize='12')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.set_ylabel(r'$L(L+1)|{N_L}^\mathrm{cov}|/2\pi$', fontsize=16)
        # ax.set_ylim((1.e-15))
        ax.set_xlim((2., 4.e4))
        ax.tick_params(axis='both', labelsize=12)
        plt.show()

    def corr_coef(self, est):
        self.pear_coeff = {}
        n_est = len(est)
        for i_est in range(n_est):
            XY = est[i_est]
            for i2 in range(i_est+1, n_est):
                AB = est[i2]
                if self.cov_d[XY+AB] == 0:
                    break
                else:
                    # np.seterr(divide='raise', invalid='raise')
                    # print XY, AB, self.N_d[XY](self.L), self.N_d[AB](self.L)
                    numerator = self.cov_d[XY+AB](self.L)
                    denominator = np.sqrt(self.var_d[XY](self.L)*self.var_d[AB](self.L))
                    # result[:, counter] = numerator/denominator
                    self.pear_coeff[XY+AB] = numerator/denominator
                    # print counter, XY+AB
        return self.pear_coeff

    def plot_corrcoef(self, est):
        pear = self.corr_coef(est)
        pairs = ['TT-EE', 'TT-TE', 'EE-TE', 'TB-EB']
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.L, np.abs(pear['TTEE']), label='TT-EE')
        ax.plot(self.L, np.abs(pear['TTTE']), label='TT-TE')
        ax.plot(self.L, np.abs(pear['EETE']), label='EE-TE')
        ax.plot(self.L, np.abs(pear['TBEB']), label='TB-EB')
        # for i in range(4):
        #     ax.plot(self.L, pear[:, i], label=pairs[i])
        ax.legend(prop={'size': 14}, frameon=False)
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.set_ylabel(r'$|{N_{XY}}(L)| / \sqrt{{N_{XX}(L)} \times {N_{YY}(L)}}$', fontsize=16)
        ax.set_ylim(ymin=-0.001)
        ax.set_xlim((2., self.l1Max))
        ax.set_xscale('log')
        ax.tick_params(axis='both', labelsize=14)
        plt.show()


if __name__ == '__main__':
    import time
    # import imp
    # import cell_cmb
    # imp.reload(cell_cmb)
    from cell_cmb import *

    time0 = time()
    
    """
    # S4
    beam = 1. # 7.
    noise_t = 1.  # 27.
    noise_p = 1.*np.sqrt(2)  # 40.*np.sqrt(2)
    lMinT = 30.
    lMaxT = 3.5e3
    lMaxP = 3.5e3
    # """

    HO02_ref = {"name": "HO02_ref", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
            "beam": 4., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

    AdvACT = {"name": "AdvACT", "lMin": 30., "lMaxT": 6000., "lMaxP": 6000.,
              "beam": 1.4, "noise_t": 10., "noise_p": 10.*np.sqrt(2)}

    SO = {"name": "SO", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
          "beam": 1.4, "noise_t": 8., "noise_p": 8.*np.sqrt(2)}

    CMBS4 = {"name": "CMBS4", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
             "beam": 1., "noise_t": 1., "noise_p": 1.*np.sqrt(2)}

    Planck_smica = {"name": "Planck", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
                    "beam": 5., "noise_t": 35., "noise_p": 60.}

    custom = {"name": "Ideal", "lMin": 30., "lMaxT": 3.e3, "lMaxP": 3.e3,
              "beam": 0., "noise_t": 0., "noise_p": 0.}

    time0 = time()

    exp = Planck_smica
    cmb = Cell_cmb(exp)

    l_est = lensing_estimator(cmb)
    # l_est.cmb.plot_cell()
    # ell = np.linspace(2, int(lMaxP), int(lMaxP)-1)
    l_est.cmb.plot_cell_2()
    # est = ['TT']
    est = ['TT', 'EE', 'TE', 'TB', 'EB']
    # est = ['TB', 'EB']
    # est = ['TT', 'EE', 'TE']
    l_est.calc_var(est)
    # l_est.interp_var(est)
    # l_est.plot_var(est)
    # l_est.comp_manu()
    # l_est.calc_cov(est)

    # l_est.interp_cov(est)
    # l_est.plot_cov(est)
    # l_est.plot_cov_pair()

    # l_est.plot_corrcoef(est)
    # l_est.comp_manu()
    # """
    print(time()-time0)
