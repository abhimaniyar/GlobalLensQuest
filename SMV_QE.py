# SQE_MV

from headers import *


class sqe_lensing_estimator(object):

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
        self.var_out = 'output/SQE_variance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise))
        self.lambda_out = 'output/SQE_lambda_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise))

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
        """
        if XY == 'TT':
            result = self.cmb.unlensedTT(l_1)*Ldotl_1
            result += self.cmb.unlensedTT(l_2)*Ldotl_2
            # print result
            # sys.exit()
        elif XY == 'EE':
            result = self.cmb.unlensedEE(l_1)*Ldotl_1
            result += self.cmb.unlensedEE(l_2)*Ldotl_2
            result *= np.cos(2.*phi12)
        elif XY == 'TE':
            # there is a typo in HO02!!!!!!!!!
            # instead of cos(phi12) it should be cos(2*phi12)!!!!!
            result = self.cmb.unlensedTE(l_1)*np.cos(2.*phi12)*Ldotl_1
            result += self.cmb.unlensedTE(l_2)*Ldotl_2
        elif XY == 'TB':
            result = self.cmb.unlensedTE(l_1)*np.sin(2.*phi12)*Ldotl_1
        elif XY == 'EB':
            # there is a typo in HO02!!!!!!!!!
            # instead of - it should be + between first and second term!!!!!
            result = self.cmb.unlensedEE(l_1)*Ldotl_1
            result += self.cmb.unlensedBB(l_2)*Ldotl_2
            result *= np.sin(2.*phi12)
        elif XY == 'BB':
            result = self.cmb.unlensedBB(l_1)*Ldotl_1
            result += self.cmb.unlensedBB(l_2)*Ldotl_2
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

    def lambda_ind(self, L, l_1, phi1, XY):
        l_2 = self.l2(L, l_1, phi1)
        if XY == 'TT':
            numerator = self.f_XY(L, l_1, phi1, XY)**2
            denominator = 2.*self.cmb.totalTT(l_1)*self.cmb.totalTT(l_2)
            result = numerator/denominator
        elif XY == 'EE':
            numerator = self.f_XY(L, l_1, phi1, XY)**2
            denominator = 2.*self.cmb.totalEE(l_1)*self.cmb.totalEE(l_2)
            result = numerator/denominator
        elif XY == 'TE':
            numerator = self.f_XY(L, l_1, phi1, XY)**2
            denominator = self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2)
            result = numerator/denominator
        elif XY == 'TB':
            # TB power spectrum assumed zero
            numerator = self.f_XY(L, l_1, phi1, XY)**2
            denominator = self.cmb.totalTT(l_1)*self.cmb.totalBB(l_2)
            result = numerator/denominator
        elif XY == 'EB':
            # EB assumed zero
            numerator = self.f_XY(L, l_1, phi1, XY)**2
            denominator = self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2)
            result = numerator/denominator
        elif XY == 'BB':
            numerator = self.f_XY(L, l_1, phi1, XY)**2
            denominator = 2.*self.cmb.totalBB(l_1)*self.cmb.totalBB(l_2)
            result = numerator/denominator

        return result

    def lambda_sqe(self, L):
        l1min = self.l1Min
        def integrand(l_1, phi1):
            l_2 = self.l2(L, l_1, phi1)
            est = ['TT', 'EE', 'TE', 'TB', 'EB']
            n_est = len(est)
            res = np.zeros((len(l_1), n_est))
            for i_est in range(n_est):
                XY = est[i_est]
                el1max = self.l1max(XY)
                res[:, i_est] = self.lambda_ind(L, l_1, phi1, XY)

            res1 = np.sum(res, 1)
            res1 *= 2*l_1  # **2
            # d^2l_1 = dl_1*l_1*dphi1
            """factor of 2 above because phi integral is symmetric. Thus we've
            put instead of 0 to 2pi, 2 times 0 to pi
            """
            res1 /= (2.*np.pi)**2
            # """
            # idx = np.where((l_2 < l1min) | (l_2 > l1max))
            idx = np.where((l_1 < l1min) | (l_1 > el1max) | (l_2 < l1min) | (l_2 > el1max))

            res1[idx] = 0.
            # """
            return res1

        l1max = self.l1Max
        # """
        if L > 2.*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.
        # """
        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        phil = np.linspace(0., np.pi, 30)
        int_1 = np.zeros(len(phil))
        for i in range(len(phil)):
            intgnd = integrand(l1, phil[i])
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')
        result = integrate.simps(int_1, x=phil, even='avg')

        if not np.isfinite(result):
            result = 0.

        if result < 0.:
            print L

        return result

    def calc_lambda(self):

        data = np.zeros((self.Nl, 2))
        data[:, 0] = np.copy(self.L)

        pool = Pool(ncpus=4)

        print "Computing Lambda(L)"

        def f(l):
            return self.lambda_sqe(l)

        intarray = np.array(pool.map(f, self.L))

        result = 1./intarray
        result *= self.L**2
        data[:, -1] = result
        np.savetxt(self.lambda_out, data)

    def interp_lambda(self):
        print "Interpolating Lambda(L)"
        data = np.genfromtxt(self.lambda_out)
        L = data[:, 0]
        self.lambda_L = interp1d(L, data[:, -1], kind='linear',
                                 bounds_error=False, fill_value=0.)

    def F_XY(self, L, l_1, phi1, XY):
        """
        Weighing terms for the estimator.
        This decides the weights for a corresponding pair of multipoles for
        X and Y.
        """

        l_2 = self.l2(L, l_1, phi1)

        if XY == 'TT':
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = 2.*self.cmb.totalTT(l_1)*self.cmb.totalTT(l_2)
            result = numerator/denominator
        elif XY == 'EE':
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = 2.*self.cmb.totalEE(l_1)*self.cmb.totalEE(l_2)
            result = numerator/denominator
        elif XY == 'TE':
            numerator = self.f_XY(L, l_1, phi1, XY)
            denominator = self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2)
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

        lambda_l = self.lambda_L(L)
        result *= lambda_l
        return result

    def cov_ind(self, L, l_1, phi1, XY, AB):
        """
        Covariance of the QE for XY and AB map choice. Because the unlensed
        TB and EB power spectra are zero, they do not have any covariance
        with other map combinations.
        """
        # print "here"
        l_2 = self.l2(L, l_1, phi1)
        # print "here"
        phi2 = self.phi2(L, l_1, phi1)

        # print "here"

        if (XY == 'TT')*(AB == 'TB') or (XY == 'TT')*(AB == 'EB') or \
            (XY == 'EE')*(AB == 'TB') or (XY == 'EE')*(AB == 'EB') or \
           (XY == 'TE')*(AB == 'TB') or (XY == 'TE')*(AB == 'EB'):
            return 0.

        if XY == 'TT' and AB == 'TT':
            a = 2*self.F_XY(L, l_1, phi1, XY)*self.cmb.totalTT(l_1)*self.cmb.totalTT(l_2)

        if XY == 'TT' and AB == 'EE':
            a = self.F_XY(L, l_1, phi1, AB)*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)
            a += self.F_XY(L, l_2, phi2, AB)*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        elif XY == 'TT' and AB == 'TE':
            a = self.F_XY(L, l_1, phi1, AB)*self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2)
            a += self.F_XY(L, l_2, phi2, AB)*self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2)

        elif XY == 'EE' and AB == 'EE':
            a = 2*self.F_XY(L, l_1, phi1, XY)*self.cmb.totalEE(l_1)*self.cmb.totalEE(l_2)

        elif XY == 'EE' and AB == 'TE':
            a = self.F_XY(L, l_1, phi1, AB)*self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2)
            a += self.F_XY(L, l_2, phi2, AB)*self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2)

        elif XY == 'TE' and AB == 'TE':
            a = self.F_XY(L, l_1, phi1, XY)*self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2)
            a += self.F_XY(L, l_2, phi2, XY)*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        elif XY == 'TB' and AB == 'TB':
            a = self.F_XY(L, l_1, phi1, XY)*self.cmb.totalTT(l_1)*self.cmb.totalBB(l_2)

        elif XY == 'TB' and AB == 'EB':
            a = self.F_XY(L, l_1, phi1, AB)*self.cmb.totalTE(l_1)*self.cmb.totalBB(l_2)

        elif XY == 'EB' and AB == 'EB':
            a = self.F_XY(L, l_1, phi1, XY)*self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2)

        result = a*self.F_XY(L, l_1, phi1, XY)

        return result

    def cov_sqe(self, L):
        l1min = self.l1Min
        l1max = self.l1Max

        if L > 2*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.

        def integrand(l_1, phi1):
            l_2 = self.l2(L, l_1, phi1)
            # print "here"
            est = ['TT', 'EE', 'TE', 'TB', 'EB']
            n_est = len(est)
            res = np.zeros((len(l_1), 16))
            counter = 0
            for i_est in range(n_est):
                XY = est[i_est]
                for i2 in range(i_est, n_est):
                    AB = est[i2]
                    # print "here"
                    res[:, counter] = self.cov_ind(L, l_1, phi1, XY, AB)
                    counter += 1
            # print counter
            res1 = np.sum(res, 1)
            res1 *= 2*l_1  # **2
            # print "here"
            res1 /= (2.*np.pi)**2
            idx = np.where((l_1 < l1min) | (l_1 > l1max) | (l_2 < l1min) | (l_2 > l1max))
            res1[idx] = 0.
            return res1

        # """        
        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        phil = np.linspace(0., np.pi, 30)
        int_1 = np.zeros(len(phil))
        for i in range(len(phil)):
            intgnd = integrand(l1, phil[i])
            # print "here"
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')
        int_l1 = integrate.simps(int_1, x=phil, even='avg')
        result = int_l1
        # """
        result *= 1./L**2

        if not np.isfinite(result):
            result = 0.
        return result

    def calc_cov(self):
        data = np.zeros((self.Nl, 2))
        data[:, 0] = np.copy(self.L)
        pool = Pool(ncpus=4)

        print "Computing covariance"

        def f(l):
            return self.cov_sqe(l)

        data[:, -1] = np.array(pool.map(f, self.L))
        np.savetxt(self.var_out, data)

    def interp_cov(self):
        print "Interpolating covariances"

        data = np.genfromtxt(self.var_out)
        L = data[:, 0]
        self.cov_d = interp1d(L, data[:, -1], kind='linear', bounds_error=False, fill_value=0.)

    def plot_cov(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        data2 = np.genfromtxt("input/CAMB/qe_lenspotentialCls.dat")
        L = data2[:, 0]
        ax.plot(L, data2[:, 5], 'r-', lw=1.5, label=r'signal')

        ax.plot(self.L, self.L*(self.L+1)*self.cov_d(self.L)/(2*np.pi), 'k', lw=1.5, label='SQE min var')

        ax.legend(loc=2, fontsize='12')  # , labelspacing=0.1)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.set_ylabel(r'$L(L+1)C_L^{dd}/2\pi$', fontsize=16)
        ax.set_ylim((4.e-9, 3.e-7))
        ax.set_xlim((2., 3.e3))
        ax.tick_params(axis='both', labelsize=12)
        plt.show()


if __name__ == '__main__':
    import time
    import imp
    import cell_cmb
    imp.reload(cell_cmb)
    from cell_cmb import *

    time0 = time()
    pla = 1
    SO = {"name": "SO", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1.4, "noise_t": 5., "noise_p": 5.*np.sqrt(2)}
    exp = SO
    cmb = Cell_cmb(exp, pla)
    smv_est = sqe_lensing_estimator(cmb)
    smv_est.calc_lambda()
    smv_est.interp_lambda()
    smv_est.calc_cov()
    smv_est.interp_cov()
    smv_est.plot_cov()

    print time() - time0
