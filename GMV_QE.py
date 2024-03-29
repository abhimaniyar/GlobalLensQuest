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
        # a1 = np.logspace(np.log10(1.), np.log10(100.), 20, 10.)
        # a2 = np.logspace(np.log10(110.), np.log10(1500.), 140, 10.)
        # a3 = np.logspace(np.log10(1600.), np.log10(2*self.l1Max+1.), 51, 10.)
        # self.L = np.concatenate((a1, a2, a3))
        self.L = np.logspace(np.log10(1.), np.log10(2*self.l1Max+1.), 201, 10.)
        # self.L = np.logspace(np.log10(1.), np.log10(2*self.l1Max+1.), 51, 10.)
        self.Nl = len(self.L)
        self.N_phi = 50  # number of steps for angular integration steps
        # reduce to 50 if you need around 0.6% max accuracy till L = 3000
        # from 200 to 400, there is just 0.03% change in the noise curves till L=3000
        self.var_out = 'output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise), str(self.N_phi))


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
        # result = np.arctan2(y, x)
        result = -np.arctan2(y, x)
        # - sign because we want phi1 - phi2.
        return result

    def phi2(self, L, l_1, phi1):
        """
        phi2 = phi1 - phi12
        """
        result = phi1 - self.phi12(L, l_1, phi1)
        # result = self.phi12(L, l_1, phi1) + phi1
        return result

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
        return result

    # """
    def M_1(self, L, l_1, phi1):
        l_2 = self.l2(L, l_1, phi1)
        # m1 = np.zeros((4, 4))
        m1 = np.zeros((len(l_1), 4, 4))
        m1[:, 0, 0] = 2.*self.cmb.totalTT(l_1)*self.cmb.totalTT(l_2)

        m1[:, 1, 1] = 2.*self.cmb.totalEE(l_1)*self.cmb.totalEE(l_2)

        m1[:, 2, 2] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) +
                           self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2)) + \
                           self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        m1[:, 3, 3] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) +
                           self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2)) - \
                           self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        m1[:, 0, 1] = m1[:, 1, 0] = 2.*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        # ###############################

        m1[:, 0, 2] = m1[:, 2, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2) +
                                     self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2))

        m1[:, 0, 3] = m1[:, 3, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2) -
                                     self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2))

        m1[:, 1, 2] = m1[:, 2, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2) +
                                     self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2))

        m1[:, 1, 3] = m1[:, 3, 1] = -(self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2) -
                                     self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2))

        m1[:, 2, 3] = m1[:, 3, 2] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) -
                                         self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2))

        return m1  # [np.ix_(np.arange(len(l_1)), [2, 3], [2, 3])]

    def f_1(self, L, l_1, phi1):

        l_2 = self.l2(L, l_1, phi1)
        phi2 = self.phi2(L, l_1, phi1)

        f_TE_sym = (self.f_XY(L, l_1, phi1, 'TE')+self.f_XY(L, l_2, phi2, 'TE'))/2.
        f_TE_asym = (self.f_XY(L, l_1, phi1, 'TE')-self.f_XY(L, l_2, phi2, 'TE'))/2.

        est1 = ['TT', 'EE', 'TE', 'TE']
        n1 = len(est1)
        # f1 = np.zeros((n1, 1))
        f1 = np.zeros((len(l_1), n1))

        for i in range(2):
            # f1[i, 0] = self.f_XY(L, l_1, phi1, est1[i])
            f1[:, i] = self.f_XY(L, l_1, phi1, est1[i])

        # f1[2, 0] = f_TE_sym
        # f1[3, 0] = f_TE_asym
        f1[:, 2] = f_TE_sym
        f1[:, 3] = f_TE_asym
        return f1  # [np.ix_(np.arange(len(l_1)), [2, 3])]

    def M1_inv(self, L, l_1, phi1):
        # nl1 = len(l_1[:, 0])
        l_2 = self.l2(L, l_1, phi1)
        nl2 = len(l_2)
        inv_m1 = np.zeros((nl2, 4, 4))
        det = self.cmb.totalTT(l_1)*self.cmb.totalEE(l_1)-self.cmb.totalTE(l_1)**2
        det *= self.cmb.totalTT(l_2)*self.cmb.totalEE(l_2)-self.cmb.totalTE(l_2)**2
        # determinant = 1./det
        inv_m1[:, 0, 0] = 0.5*self.cmb.totalEE(l_1)*self.cmb.totalEE(l_2)

        inv_m1[:, 1, 1] = 0.5*self.cmb.totalTT(l_1)*self.cmb.totalTT(l_2)

        inv_m1[:, 2, 2] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) +
                            self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2)) + \
                            self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        inv_m1[:, 3, 3] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) +
                            self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2)) - \
                            self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        inv_m1[:, 0, 1] = inv_m1[:, 1, 0] = 0.5*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        inv_m1[:, 0, 2] = inv_m1[:, 2, 0] = -0.5*(self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2) +
                                            self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2))

        inv_m1[:, 0, 3] = inv_m1[:, 3, 0] = 0.5*(self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2) -
                                           self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2))

        inv_m1[:, 1, 2] = inv_m1[:, 2, 1] = -0.5*(self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2) +
                                            self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2))

        inv_m1[:, 1, 3] = inv_m1[:, 3, 1] = -0.5*(self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2) -
                                           self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2))

        inv_m1[:, 2, 3] = inv_m1[:, 3, 2] = 0.5*(self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2) -
                                           self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2))
        return inv_m1/det[:, None, None]

    def F1prime(self, L, l_1, phi1):

        """
        F1 = A1(L)*M1^{-1}*f1
        F1prime = M1^{-1}*f1
        """
        f_1 = self.f_1(L, l_1, phi1)
        """
        M_1 = self.M_1(L, l_1, phi1)
        M1invf1 = np.linalg.solve(M_1, f_1)
        """
        M1_inv = self.M1_inv(L, l_1, phi1)
        M1invf1 = np.einsum('ijk, ij -> ik', M1_inv, f_1)
        # """
        # M1_inv = np.linalg.inv(M_1)
        # M1_inv = pinvh(M_1)
        # M1_inv = inv(M_1)

        # F1 = np.matmul(M1_inv, f_1)
        # M1invf1 = np.matmul(M1_inv, f_1)
        return M1invf1

    def A_1(self, L):

        # print "calculating A_1"
        l1min = self.l1Min
        l1max = max(self.cmb.lMaxT, self.cmb.lMaxP)
        # """
        if L > 2.*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.
        # """
        def integrand(l_1, phil):

            l_2 = self.l2(L, l_1, phil)

            M1invf1 = self.F1prime(L, l_1, phil)
            f_1 = self.f_1(L, l_1, phil)
            Fdotf = np.sum(M1invf1*f_1, -1)
            result = Fdotf
            result *= 2*l_1  # **2
            """factor of 2 above because phi integral is symmetric. Thus we've
            put instead of 0 to 2pi, 2 times 0 to pi
            Also, l_1^2 instead of l_1 if we are taking log spacing for
            l_1"""
            result /= (2.*np.pi)**2
            idx = np.where((l_1 < l1min) | (l_1 > l1max) | (l_2 < l1min) | (l_2 > l1max))[0]
            # idx = np.where((l_2 < l1min) | (l_2 > l1max))
            result[idx] = 0.
            return result

        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        # l1 = np.logspace(np.log10(l1min), np.log10(l1max), int(l1max-l1min+1))
        phi1 = np.linspace(0., np.pi, self.N_phi)
        int_1 = np.zeros(len(phi1))
        for i in range(len(phi1)):
            intgnd = integrand(l1, phi1[i])
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')

        int_ll = integrate.simps(int_1, x=phi1, even='avg')
        # int_ll = np.trapz(int_1, x=phi1)

        result = 1./int_ll
        result *= L**2  # factor of L**2 here means we are basically
        # calculating the reconstruction noise for d field instead of the
        # phi field.

        if not np.isfinite(result):
            result = 0.

        if result < 0.:
            print(L)
        return result

    def M_2(self, L, l_1, phi1):
        m2 = np.zeros((len(l_1), 2, 2))
        l_2 = self.l2(L, l_1, phi1)
        m2[:, 0, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalBB(l_2))

        m2[:, 1, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2))

        m2[:, 0, 1] = m2[:, 1, 0] = (self.cmb.totalTE(l_1)*self.cmb.totalBB(l_2))
        return m2

    def f_2(self, L, l_1, phi1):

        est2 = ['TB', 'EB']
        n2 = len(est2)
        # f2 = np.zeros((n2, 1))
        f2 = np.zeros((len(l_1), n2))

        for i in range(n2):
            f2[:, i] = self.f_XY(L, l_1, phi1, est2[i])
            # f2[i, 0] = self.f_XY(L, l_1, phi1, est2[i])

        return f2

    def M2_inv(self, L, l_1, phi1):
        l_2 = self.l2(L, l_1, phi1)
        nl2 = len(l_2)
        inv_m2 = np.zeros((nl2, 2, 2))
        det = self.cmb.totalTT(l_1)*self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2)**2
        det -= self.cmb.totalTE(l_1)**2*self.cmb.totalBB(l_2)**2

        inv_m2[:, 0, 0] = self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2)
        inv_m2[:, 1, 1] = self.cmb.totalTT(l_1)*self.cmb.totalBB(l_2)
        inv_m2[:, 0, 1] = inv_m2[:, 1, 0] = -self.cmb.totalTE(l_1)*self.cmb.totalBB(l_2)
        return inv_m2/det[:, None, None]

    def F2prime(self, L, l_1, phi1):
        """
        F2 = A2(L)*M2^{-1}*f2
        F2prime = M2^{-1}*f2
        """

        f_2 = self.f_2(L, l_1, phi1)
        """
        M_2 = self.M_2(L, l_1, phi1)
        M2invf2 = np.linalg.solve(M_2, f_2)
        """
        M2_inv = self.M2_inv(L, l_1, phi1)
        M2invf2 = np.einsum('ijk, ij -> ik', M2_inv, f_2)
        # """
        # M2_inv = np.linalg.inv(M_2)
        # M2_inv = pinvh(M_2)
        # M2_inv = inv(M_2)
        # det = M_2[0, 0]*M_2[1, 1] - M_2[0, 1]**2
        # M2_inv = np.array([[M_2[1, 1], -M_2[0, 1]],[-M_2[0, 1], M_2[0, 0]]])/det

        # M2invf2 = np.matmul(M2_inv, f_2)
        return M2invf2

    def A_2(self, L):

        # print "calculating A_2"
        l1min = self.l1Min
        l1max = max(self.cmb.lMaxT, self.cmb.lMaxP)
        # """
        if L > 2.*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.
        # """

        def integrand(l_1, phil):
            l_2 = self.l2(L, l_1, phil)
            # """

            F2p = self.F2prime(L, l_1, phil)
            f_2 = self.f_2(L, l_1, phil)
            Fdotf = np.sum(F2p*f_2, -1)
            result = Fdotf
            result *= 2*l_1  # **2
            """factor of 2 above because phi integral is symmetric. Thus we've
            put instead of 0 to 2pi, 2 times 0 to pi
            Also, l_1^2 instead of l_1 because if are taking log spacing for
            l_1"""
            result /= (2.*np.pi)**2
            # result *= 2.
            # idx = np.where((l_2 < l1min) | (l_2 > l1max))
            """
            # check integration bounds
            # """
            idx = np.where((l_1 < l1min) | (l_1 > l1max) | (l_2 < l1min) | (l_2 > l1max))
            # just as a precaution
            result[idx] = 0.
            return result

        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        # l1 = np.logspace(np.log10(l1min), np.log10(l1max), int(l1max-l1min+1))
        phi1 = np.linspace(0., np.pi, self.N_phi)
        int_1 = np.zeros(len(phi1))
        for i in range(len(phi1)):
            intgnd = integrand(l1, phi1[i])
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')
        int_ll = integrate.simps(int_1, x=phi1, even='avg')
        # int_ll = np.trapz(int_1, x=phi1)

        result = 1./int_ll
        result *= L**2  # factor of L**2 here means we are basically
        # calculating the reconstruction noise for d field instead of the
        # phi field.

        if not np.isfinite(result):
            result = 0.

        if result < 0.:
            print(L)
        return result

    def var_d(self, var_d1, var_d2):
        inv_vard1 = 1./var_d1
        inv_vard2 = 1./var_d2
        vard = 1./(inv_vard1+inv_vard2)
        return vard

    def calc_tvar(self):
        data = np.zeros((self.Nl, 4))
        data[:, 0] = np.copy(self.L)
        pool = Pool(ncpus=4)

        def f1(l):
            return self.A_1(l)

        def f2(l):
            return self.A_2(l)

        """
        for i in range(len(self.L)):
            # print self.L[i]
            data[i, 1] = f1(self.L[i])
        """
        print("Computing variance for d1")
        data[:, 1] = np.array(pool.map(f1, self.L))

        print("Computing variance for d2")
        data[:, 2] = np.array(pool.map(f2, self.L))

        print("Computing variance for d")
        data[:, 3] = self.var_d(data[:, 1], data[:, 2])

        # data[:, 3] = self.var_d(2.*data[:, 1], 2.*data[:, 1])

        np.savetxt(self.var_out, data)

    def interp_tvar(self):
        print("Interpolating variances")

        self.N_d = {}
        data = np.genfromtxt(self.var_out)

        L = data[:, 0]

        norm1 = data[:, 1].copy()
        self.N_d['d1'] = interp1d(L, norm1, kind='linear', bounds_error=False, fill_value=0.)

        norm2 = data[:, 2].copy()
        self.N_d['d2'] = interp1d(L, norm2, kind='linear', bounds_error=False, fill_value=0.)

        norm = data[:, 3].copy()
        self.N_d['d'] = interp1d(L, norm, kind='quadratic', bounds_error=False, fill_value=0.)

    def plot_tvar(self):
        data = np.genfromtxt("input/CAMB/Julien_lenspotentialCls.dat")
        L = data[:, 0]  # data[:, 5] = l(l+1)C^dd_ell/2 pi = (l(l+1))**2*C^phiphi_ell/2 pi
        clphiphi = data[:, 5]*(2.*np.pi)/(L*(L+1))**2

        # data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise)))
        # L2 = data2[:, 0]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(L, data[:, 5], 'k-', lw=2.5, label=r'signal')
        # ax.plot(L, (L*(L+1))**2*clphiphi/2./np.pi, 'k-', lw=2.5, label=r'signal')

        est = ['d1', 'd2', 'd']
        lbl = ['GMV TT-EE-TE', 'GMV TB-EB', 'GMV combined']

        nest = len(est)
        for iEst in range(nest):
            XY = est[iEst]
            ax.plot(self.L, (self.L*(self.L+1))*self.N_d[XY](self.L)/(2*np.pi),
                    lw=1.5, label=lbl[iEst])
            # ax.plot(self.L, self.N_d[XY](self.L)/(2*np.pi),
            #     lw=1.5, label=lbl[iEst])

        # ax.plot(L2, (L2*(L2+1))*data2[:, -1]/(2*np.pi), 'k--', lw=1.5,
        #         label='HO02 MV')
        # ax.plot(L2, data2[:, -1]/(2*np.pi), 'k--', lw=1.5,
        #       label='HO02 MV')

        ax.legend(prop={'size': 17}, loc='upper left', ncol=2, frameon=False,
                  labelspacing=0.2)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        ax.set_xlabel(r'$L$', fontsize=22)
        ax.set_ylabel(r'$L(L+1)C_L^{dd}/2\pi$', fontsize=22)
        ax.set_ylim(4.e-9, 3.e-7)
        ax.set_xlim(2., self.cmb.lMaxT)
        ax.tick_params(axis='both', labelsize=22)
        plt.show()

    def SNR_comp(self, exp):
        data1 = np.genfromtxt("input/CAMB/Julien_lenspotentialCls.dat")
        L1 = data1[:, 0]
        cldd = data1[:, 5]*2*np.pi/(L1*(L1+1))
        clphiphi = 2*np.pi*cldd/(L1*(L1+1))**2

        clddint = interp1d(L1, cldd, kind='quadratic', bounds_error=False, fill_value=0.)
        # clppint = interp1d(L1, clphiphi, kind='quadratic', bounds_error=False, fill_value=0.)

        data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (exp['name'], str(exp['lMin']), str(exp['lMaxT']), str(exp['lMaxP']), str(exp['beam']), str(exp['noise_t'])))
        L = data[:, 0]
        interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (exp['name'], str(exp['lMin']), str(exp['lMaxT']), str(exp['lMaxP']), str(exp['beam']), str(exp['noise_t'])))
        L2 = data2[:, 0]
        interp_HO = interp1d(L2, data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        L_p = np.logspace(np.log10(2.), np.log10(6000.), 201, 10.)

        SNT2 = (clddint(L_p)/(clddint(L_p)+interp_tmv(L_p)))**2
        # SNT2 = (clppint(L_p)/(clppint(L_p)+interp_tmv(L_p)))**2
        SNT2 *= (2*L_p+1)
        cumSNT2 = np.cumsum(SNT2)
        cumSNT = np.sqrt(cumSNT2)

        SNH2 = (clddint(L_p)/(clddint(L_p)+interp_HO(L_p)))**2
        # SNH2 = (clppint(L_p)/(clppint(L_p)+interp_HO(L_p)))**2
        SNH2 *= (2*L_p+1)
        cumSNH2 = np.cumsum(SNH2)
        cumSNH = np.sqrt(cumSNH2)

        SNT = np.sqrt(SNT2)
        SNH = np.sqrt(SNH2)
        # ratio = SNT/SNH
        # cum_ratio = cumSNT/cumSNH

        # percent = (SNT - SNH)*100/SNH
        cum_percent = (cumSNT - cumSNH)*100/cumSNH

        plt.figure()
        plt.loglog(L_p, cumSNH, 'r')
        plt.loglog(L_p, cumSNT, 'b')
        return SNT2, SNH2, cumSNT, cumSNH


if __name__ == '__main__':
    import time
    import imp
    import cell_cmb
    imp.reload(cell_cmb)
    from cell_cmb import *

    time0 = time()

    SO = {"name": "SO", "lMin": 30., "lMaxT": 3000., "lMaxP": 3000.,
         "beam": 1.4, "noise_t": 5., "noise_p": 5.*np.sqrt(2)}

    exp = SO
    cmb = Cell_cmb(exp)

    l_est = lensing_estimator(cmb)
    l_est.cmb.plot_cell()
    # """

    print(time()-time0)

    # """
    # l_est.calc_tvar()
    l_est.interp_tvar()
    l_est.plot_tvar()
    plt.figure()
    data = np.genfromtxt('true_variance_lmin%s_lmaxT%s_lmaxP%s_fin_ownintegrator_simps_simps.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    # data = np.genfromtxt('true_variance_lmin%s_lmaxT%s_lmaxP%s_fin_ownintegrator_trapz_trapz.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    L = data[:, 0]
    plt.plot(L, L*(L+1)*data[:, -1]/(2*np.pi), 'r', label='TMV')
    data2 = np.genfromtxt('covariance_minvar_lmin%s_lmaxT%s_lmaxP%s_own.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    # data2 = np.genfromtxt('variance_ind_lmin%s_lmaxT%s_lmaxP%s.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    L2 = data2[:, 0]
    plt.plot(L2, L2*(L2+1)*data2[:, -1]/(2*np.pi), 'b', label='HO')
    # n1 = 1./((1./data2[:, 1]) + (1./data2[:, 2]) + (1./data2[:, 3]))
    # plt.plot(L2, L2*(L2+1)*n1/(2*np.pi), 'b', label='HO')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$L(L+1)N_L^{dd}/2\pi$', fontsize=16)
    plt.xlabel(r'$L$', fontsize=16)
    plt.legend(prop={'size':12})
    plt.tick_params(axis='both', labelsize=14)
    # """
    # """
    plt.figure()
    interp_HO = interp1d(L2, data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
    interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
    L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
    # L_p = np.logspace(np.log10(1.), np.log10(10000.), 51, 10.)
    # vart_var = interp_tmv(L_p)/interp_HO(L_p)  # data2[:, -1]
    vart_var = (interp_HO(L_p)-interp_tmv(L_p))*100./interp_HO(L_p)  # data2[:, -1]
    plt.plot(L_p, vart_var, 'b')  # , label=r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$')
    plt.xscale('log')
    # plt.ylim(0.1, 1.1)
    plt.ylim(ymax=15)
    # plt.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    # plt.ylabel(r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$')
    plt.ylabel(r'$(N_{mv}^\mathrm{HO02}-N_{mv}^\mathrm{true}) \times 100/N_{mv}^\mathrm{HO02}$', fontsize=16)
    plt.xlabel(r'$L$', fontsize=16)
    plt.legend(prop={'size':12})
    plt.tick_params(axis='both', labelsize=14)

    interp_tmv1 = interp1d(L, data[:, 1], kind='linear', bounds_error=False, fill_value=0.)
    interp_tmv2 = interp1d(L, data[:, 2], kind='linear', bounds_error=False, fill_value=0.)
    data21 = np.genfromtxt('covariance_minvar_lmin%s_lmaxT%s_lmaxP%s_TT_EE_TE_only_own.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    interp_HO_1 = interp1d(data21[:, 0], data21[:, -1], kind='linear', bounds_error=False, fill_value=0.)
    # data21 = np.genfromtxt('variance_ind_lmin%s_lmaxT%s_lmaxP%s_TT_EE_TE_only.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    # interp_HO_1 = interp1d(L2, data21[:, 3], kind='linear', bounds_error=False, fill_value=0.)

    data22 = np.genfromtxt('covariance_minvar_lmin%s_lmaxT%s_lmaxP%s_TB_EB_only_own.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    interp_HO_2 = interp1d(data22[:, 0], data22[:, -1], kind='linear', bounds_error=False, fill_value=0.)

    plt.figure()
    vart_var1 = interp_tmv1(L_p)/interp_HO_1(L_p)  # data2[:, -1]
    vart_var2 = interp_tmv2(L_p)/interp_HO_2(L_p)  # data2[:, -1]
    plt.plot(L_p, vart_var1, 'b', label='TT-EE-TE')  # , label=r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$')
    plt.plot(L_p, vart_var2, 'r', label='EB-TB')  # , label=r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$')
    plt.legend()
    plt.xscale('log')
    plt.ylim(0.8, 1.05)
    plt.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    plt.ylabel(r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$', fontsize=16)
    plt.xlabel(r'$L$', fontsize=16)
    plt.legend(prop={'size':12})
    plt.tick_params(axis='both', labelsize=14)
    # """
    """
    interp2 = interp1d(L, data[:, -1], kind='linear', bounds_error=False, fill_value='extrapolate')
    vart_var2 = data2[:, -1]/interp2(L2)
    plt.plot(L2, vart_var2, 'b')
    # """

    """
    plt.figure()
    data1 = np.genfromtxt('true_variance_lmin%s_lmaxT%s_lmaxP%s_newbasis.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    ratio = data[:, -1]/data1[:, -1]
    plt.semilogx(L, ratio)
    # """
    # print time()-time0
