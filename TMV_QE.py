from headers import *


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
        # self.L = np.linspace(1., 201., 1001)
        self.Nl = len(self.L)
        self.var_out = 'output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise))

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

    """
    def l1max(self, XY):
        #
        # max value for l1 and l2: taken to be same
        #
        if XY == 'TT':
            return self.cmb.lMaxT
        elif XY == 'EE' or XY == 'BB' or XY == 'EB':
            return self.cmb.lMaxP
        elif XY == 'TE' or XY == 'TB':
            # taking the minimum of the two: this approach is suboptimal
            return min(self.cmb.lMaxT, self.cmb.lMaxP)
    """

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
        # print l_2, phi2, np.cos(phi2)
        # sys.exit()
        if XY == 'TT':
            result = self.cmb.unlensedTT(l_1)*Ldotl_1
            result += self.cmb.unlensedTT(l_2)*Ldotl_2
            # print Ldotl_1, Ldotl_2, self.cmb.unlensedTT(l_1), self.cmb.unlensedTT(l_2)
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
            result = self.cmb.unlensedEE(l_1)*Ldotl_1
            result -= self.cmb.unlensedBB(l_2)*Ldotl_2
            result *= np.sin(2.*phi12)
        elif XY == 'BB':
            result = self.cmb.unlensedBB(l_1)*Ldotl_1
            result += self.cmb.unlensedBB(l_2)*Ldotl_2
            result *= np.cos(2.*phi12)
        return result

    # """
    def M_1(self, L, l_1, phi1):
        l_2 = self.l2(L, l_1, phi1)
        # m1 = np.zeros((4, 4))
        m1 = np.zeros((len(l_1), 4, 4))
        # print len(l_1), np.shape(m1)
        """
        m1[0, 0] = 2.*self.cmb.totalTT(l_1)*self.cmb.totalTT(l_2)

        m1[1, 1] = 2.*self.cmb.totalEE(l_1)*self.cmb.totalEE(l_2)

        m1[2, 2] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) +
                        self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2)) + \
                        self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        m1[3, 3] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) +
                        self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2)) - \
                        self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        m1[0, 1] = m1[1, 0] = 2.*self.cmb.totalTE(l_1)*self.cmb.totalTE(l_2)

        # ###############################
        m1[0, 2] = m1[2, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2) +
                               self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2))

        m1[0, 3] = m1[3, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2) -
                               self.cmb.totalTE(l_1)*self.cmb.totalTT(l_2))

        m1[1, 2] = m1[2, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2) +
                               self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2))

        m1[1, 3] = m1[3, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2) -
                               self.cmb.totalTE(l_1)*self.cmb.totalEE(l_2))
        m1[2, 3] = m1[3, 2] = 0.5*(self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2) -
                                   self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2))
        """
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
        # """

        # ##############################

        # m1[0, 2] = m1[2, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2))

        # m1[0, 3] = m1[3, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalTE(l_2))

        # m1[1, 2] = m1[2, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2))

        # m1[1, 3] = m1[3, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalTE(l_2))
        # ###############################

        # a = np.array([1., 1., 1., 1.])
        # d = np.diag(a)
        # m1 *= d
        # print m1
        # print np.shape(m1)
        # inv(l_1)
        # print l_1, l_2
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
        # print np.shape(f1), f1
        return f1  # [np.ix_(np.arange(len(l_1)), [2, 3])]

    def F1prime(self, L, l_1, phi1):

        """
        F1 = A1(L)*M1^{-1}*f1
        F1prime = M1^{-1}*f1
        """
        f_1 = self.f_1(L, l_1, phi1)
        # print np.shape(f_1)
        # print f_1[:, 0]
        M_1 = self.M_1(L, l_1, phi1)
        # print np.shape(M_1)
        # M1_inv = np.linalg.inv(M_1)
        # M1_inv = pinvh(M_1)
        # M1_inv = inv(M_1)

        # F1 = np.matmul(M1_inv, f_1)
        # M1invf1 = np.matmul(M1_inv, f_1)
        M1invf1 = np.linalg.solve(M_1, f_1)
        # print M1invf1
        # print np.allclose(np.matmul(M_1[0], M1invf1[0]), f_1)
        # inv(l_1)
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
            
            # """
            """
            if l_1 < l1min or l_2 < l1min or l_1 > l1max or l_2 > l1max:
                return 0.
            # """

            M1invf1 = self.F1prime(L, l_1, phil)
            # print ('here 1')
            f_1 = self.f_1(L, l_1, phil)
            # print ('here 2')
            Fdotf = np.sum(M1invf1*f_1, -1)
            # Fdotf = np.diag(np.matmul(M1invf1, f_1.transpose()))
            # print np.shape(M1invf1), np.shape(f_1), np.shape(M1invf1*f_1), np.shape(Fdotf)
            # Fdotf = np.sum(F1p*f_1)
            # print ('here 3')
            result = Fdotf
            # print ('here 4')
            result *= 2*l_1  # **2
            # print ('here 5')
            """factor of 2 above because phi integral is symmetric. Thus we've
            put instead of 0 to 2pi, 2 times 0 to pi
            Also, l_1^2 instead of l_1 because we are taking log spacing for
            l_1"""
            # print "this again"
            result /= (2.*np.pi)**2
            # result *= 2.
            idx = np.where((l_1 < l1min) | (l_1 > l1max) | (l_2 < l1min) | (l_2 > l1max))[0]
            # print np.shape(idx)
            # idx = np.where((l_2 < l1min) | (l_2 > l1max))
            # print idx
            # print np.shape(result)
            result[idx] = 0.
            # print ('here 6')
            return result

        def phi_integral(ll):
            res = integrate.quad(integrand, 0., 2*np.pi, args=(ll))[0]
            return res

        def ll_integral(phil):
            res = integrate.quad(integrand, 0., 3*l1max, args=(phil))[0]
            return res

        # int_ll = integrate.quad(phi_integral, 0., lmax)[0]
        # """
        # ll = 100.
        # int_ll = phi_integral(ll)
        # """
        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        # l1 = np.logspace(np.log10(l1min), np.log10(l1max), int(l1max-l1min+1))
        phi1 = np.linspace(0., np.pi, 30)
        int_1 = np.zeros(len(phi1))
        for i in range(len(phi1)):
            intgnd = integrand(l1, phi1[i])
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')
            # int_1[i] = np.trapz(intgnd, x=l1)
            # print i
            # int_1[i] = ll_integral(phi1[i])
        int_ll = integrate.simps(int_1, x=phi1, even='avg')
        # int_ll = np.trapz(int_1, x=phi1)

        result = 1./int_ll
        result *= L**2
        # print ('here 7')
        if not np.isfinite(result):
            result = 0.

        if result < 0.:
            print L
        return result

    """
    def var_d1(self, L):

        # print "var d1 integral"
        l1min = self.l1Min
        l1max = max(self.cmb.lMaxT, self.cmb.lMaxP)

        if L > 2.*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.

        def integrand(x):
            phi1 = x[1]
            l_1 = np.exp(x[0])
            # check integration bounds
            l_2 = self.l2(L, l_1, phi1)

            if l_1 < l1min or l_2 < l1min or l_1 > l1max or l_2 > l1max:
                return 0.

            # print np.shape(l_1), np.shape(l_2)# , len(l_1)
            F_1 = self.F_1(L, l_1, phi1)
            # print F_1
            # sys.exit()
            M_1 = self.M_1(L, l_1, phi1)
            # print M_1
            # sys.exit()
            # print np.shape(F_1), np.shape(M_1)

            matmult = np.matmul(M_1, F_1)
            # print np.shape(M_1), np.shape(F_1), np.shape(matmult)
            result = np.matmul(F_1.transpose(), matmult)[0][0]
            # print np.shape(np.matmul(F_1.transpose(), matmult))
            # print len(l_1)
            result *= 2*l_1**2

            # print np.shape(result[0]), np.shape(l_1)
            # factor of 2 above because phi integral is symmetric. Thus we've
            # put instead of 0 to 2pi, 2 times 0 to pi
            # Also, l_1^2 instead of l_1 because we are taking log spacing for
            # l_1

            result /= (2.*np.pi)**2
            result *= 2.

            # factor of 2 here because we have var(d1) = 2*F_1*M_1*F_1.
            # That is how we fix the coefficients for M_1

            return result

        intg = vegas.Integrator([[np.log(l1min), np.log(l1max)], [0., np.pi]])
        intg(integrand, nitn=8, neval=1000)
        result = intg(integrand, nitn=1, neval=5000)
        result = result.mean
        # print type(result)
        A_1 = self.A_1(L)
        # print "type is %s" % type(A_1)
        result *= A_1**2/L**2  # /4
        # print "this again"

        if not np.isfinite(result):
            result = 0.

        if result < 0:
            print L
            # print len(L)
        return result
    """

    def M_2(self, L, l_1, phi1):
        """
        m2 = np.zeros((2, 2))

        l_2 = self.l2(L, l_1, phi1)

        m2[0, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalBB(l_2))

        m2[1, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2))

        m2[0, 1] = m2[1, 0] = (self.cmb.totalTE(l_1)*self.cmb.totalBB(l_2))
        """
        m2 = np.zeros((len(l_1), 2, 2))
        l_2 = self.l2(L, l_1, phi1)
        m2[:, 0, 0] = (self.cmb.totalTT(l_1)*self.cmb.totalBB(l_2))

        m2[:, 1, 1] = (self.cmb.totalEE(l_1)*self.cmb.totalBB(l_2))

        m2[:, 0, 1] = m2[:, 1, 0] = (self.cmb.totalTE(l_1)*self.cmb.totalBB(l_2))
        # """
        # a = np.array([1., 1.])
        # d = np.diag(a)
        # m2 *= d
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

    def F2prime(self, L, l_1, phi1):
        """
        F2 = A2(L)*M2^{-1}*f2
        F2prime = M2^{-1}*f2
        """

        f_2 = self.f_2(L, l_1, phi1)
        M_2 = self.M_2(L, l_1, phi1)
        # print (np.shape(f_2), np.shape(M_2))
        # M2_inv = np.linalg.inv(M_2)
        # M2_inv = pinvh(M_2)
        # M2_inv = inv(M_2)
        # det = M_2[0, 0]*M_2[1, 1] - M_2[0, 1]**2
        # M2_inv = np.array([[M_2[1, 1], -M_2[0, 1]],[-M_2[0, 1], M_2[0, 0]]])/det

        # M2invf2 = np.matmul(M2_inv, f_2)
        # print np.shape(F2)
        M2invf2 = np.linalg.solve(M_2, f_2)
        # print np.shape(M2invf2)
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
            """
            # check integration bounds
            # """
            l_2 = self.l2(L, l_1, phil)
            # """
            """
            if l_1 < l1min or l_2 < l1min or l_1 > l1max or l_2 > l1max:
                return 0.
            # """

            F2p = self.F2prime(L, l_1, phil)
            f_2 = self.f_2(L, l_1, phil)
            # print np.shape(F2p*f_2)
            Fdotf = np.sum(F2p*f_2, -1)
            result = Fdotf
            result *= 2*l_1  # **2
            """factor of 2 above because phi integral is symmetric. Thus we've
            put instead of 0 to 2pi, 2 times 0 to pi
            Also, l_1^2 instead of l_1 because we are taking log spacing for
            l_1"""
            result /= (2.*np.pi)**2
            # result *= 2.
            # print (np.shape(l_1), np.shape(l_2))
            # idx = np.where((l_2 < l1min) | (l_2 > l1max))
            idx = np.where((l_1 < l1min) | (l_1 > l1max) | (l_2 < l1min) | (l_2 > l1max))
            # print idx
            # print np.shape(result)
            result[idx] = 0.
            return result

        l1 = np.linspace(l1min, l1max, int(l1max-l1min+1))
        # l1 = np.logspace(np.log10(l1min), np.log10(l1max), int(l1max-l1min+1))
        phi1 = np.linspace(0., np.pi, 30)
        int_1 = np.zeros(len(phi1))
        for i in range(len(phi1)):
            intgnd = integrand(l1, phi1[i])
            int_1[i] = integrate.simps(intgnd, x=l1, even='avg')
            # int_1[i] = np.trapz(intgnd, x=l1)
            # print i
            # int_1[i] = ll_integral(phi1[i])
        int_ll = integrate.simps(int_1, x=phi1, even='avg')
        # int_ll = np.trapz(int_1, x=phi1)

        result = 1./int_ll
        result *= L**2
        # print ('here 7')
        if not np.isfinite(result):
            result = 0.

        if result < 0.:
            print L
        return result

    """
    def var_d2(self, L):

        # print "var d2 integral"
        l1min = self.l1Min
        l1max = max(self.cmb.lMaxT, self.cmb.lMaxP)

        if L > 2.*l1max:  # L = l1 + l2 thus max L = 2*l1
            return 0.

        def integrand(x):
            phi1 = x[1]
            l_1 = np.exp(x[0])
            # check integration bounds
            l_2 = self.l2(L, l_1, phi1)

            if l_1 < l1min or l_2 < l1min or l_1 > l1max or l_2 > l1max:
                # print "this again"
                return 0.

            F_2 = self.F_2(L, l_1, phi1)
            M_2 = self.M_2(L, l_1, phi1)

            matmul = np.matmul(M_2, F_2)
            result = np.matmul(F_2.transpose(), matmul)[0][0]
            result *= 2*l_1**2

            # factor of 2 above because phi integral is symmetric. Thus we've
            # put instead of 0 to 2pi, 2 times 0 to pi
            # Also, l_1^2 instead of l_1 because we are taking log spacing for
            # l_1 

            result /= (2.*np.pi)**2
            result *= 2.

            # factor of 2 here because we have var(d2) = 2*F_2*M_2*F_2.
            # That is how we fix the coefficients for M_1

            return result

        A_2 = self.A_2(L)
        integrator = vegas.Integrator([[np.log(l1min), np.log(l1max)], [0., np.pi]])
        integrator(integrand, nitn=8, neval=1000)
        result = integrator(integrand, nitn=1, neval=5000)
        result = result.mean
        result *= A_2**2/L**2  # /4

        if not np.isfinite(result):
            # print "this again"
            result = 0.

        return result
    """

    def var_d(self, var_d1, var_d2):
        inv_vard1 = 1./var_d1
        inv_vard2 = 1./var_d2
        vard = 1./(inv_vard1+inv_vard2)
        # vard = 1./(inv_vard2)
        return vard

    def calc_tvar(self):
        data = np.zeros((self.Nl, 4))
        data[:, 0] = np.copy(self.L)
        pool = Pool(ncpus=4)
        # pool = Pool(4)

        def f1(l):
            # print XY
            return self.A_1(l)

        def f2(l):
            # print XY
            return self.A_2(l)

        """
        for i in range(len(self.L)):
            # print self.L[i]
            data[i, 1] = f1(self.L[i])
        """
        print "Computing variance for d1"
        data[:, 1] = np.array(pool.map(f1, self.L))
        # """
        print "Computing variance for d2"
        data[:, 2] = np.array(pool.map(f2, self.L))

        print "Computing variance for d"
        data[:, 3] = self.var_d(data[:, 1], data[:, 2])
        # data[:, 3] = self.var_d(data[:, 1], data[:, 1])
        # data[:, 3] = self.var_d(2.*data[:, 1], 2.*data[:, 1])

        np.savetxt(self.var_out, data)
        # np.savetxt('true_variance_lmin%s_lmaxT%s_lmaxP%s_fin_ownintegrator_trapz_trapz.txt' % (str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP)), data)

    def interp_tvar(self):
        print "Interpolating variances"

        self.N_d = {}
        data = np.genfromtxt(self.var_out)
        # data = np.genfromtxt('true_variance_lmin%s_lmaxT%s_lmaxP%s_fin_ownintegrator_trapz_trapz.txt' % (str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP)))
        L = data[:, 0]

        norm1 = data[:, 1].copy()
        self.N_d['d1'] = interp1d(L, norm1, kind='linear', bounds_error=False, fill_value=0.)

        norm2 = data[:, 2].copy()
        self.N_d['d2'] = interp1d(L, norm2, kind='linear', bounds_error=False, fill_value=0.)

        norm = data[:, 3].copy()
        self.N_d['d'] = interp1d(L, norm, kind='linear', bounds_error=False, fill_value=0.)

    def plot_tvar(self):
        data = np.genfromtxt("../CAMB/qe_lens_lenspotentialCls.dat")
        L = data[:, 0]

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise)))
        L2 = data2[:, 0]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(L, data[:, 5], 'r-', lw=1.5, label=r'signal')

        est = ['d1', 'd2', 'd']
        lbl = ['TMV TT-EE-TE', 'TMV TB-EB', 'TMV combined']

        nest = len(est)
        for iEst in range(nest):
            XY = est[iEst]
            ax.plot(self.L, self.L*(self.L+1)*self.N_d[XY](self.L)/(2*np.pi),
                    c=plt.cm.rainbow(iEst/6.), lw=1.5, label=lbl[iEst])

        ax.plot(L2, L2*(L2+1)*data2[:, -1]/(2*np.pi), 'k--', label='HO02 MV')

        ax.legend(prop={'size': 14}, loc='upper right', frameon=False)
        ax.set_xscale('log')
        ax.set_yscale('log', nonposy='mask')
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.set_ylabel(r'$L(L+1)C_L^{dd}/2\pi$', fontsize=16)
        ax.set_ylim(1.e-9, 5.e-7)
        ax.set_xlim(2., self.cmb.lMaxT)
        ax.tick_params(axis='both', labelsize=14)
        plt.show()

    def plot_var_tvar(self):
        data = np.genfromtxt(self.var_out)
        L = data[:, 0]

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise)))
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

    def plot_var_tvar_percent(self):
        data = np.genfromtxt(self.var_out)
        L = data[:, 0]

        data2 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise)))
        L2 = data2[:, 0]

        plt.figure()
        interp_HO = interp1d(L2, data2[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
        interp_tmv = interp1d(L, data[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
        L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
        # L_p = np.logspace(np.log10(1.), np.log10(10000.), 51, 10.)
        # vart_var = interp_tmv(L_p)/interp_HO(L_p)  # data2[:, -1]
        vart_var = (interp_HO(L_p)-interp_tmv(L_p))*100./interp_HO(L_p)
        plt.plot(L_p, vart_var, 'b')
        plt.xscale('log')
        # plt.ylim(0.1, 1.1)
        plt.ylim(ymax=15)
        # plt.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
        # plt.ylabel(r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$')
        plt.ylabel(r'$(N_{mv}^\mathrm{HO02}-N_{mv}^\mathrm{TMV}) \times 100/N_{mv}^\mathrm{HO02}$', fontsize=16)
        plt.xlabel(r'$L$', fontsize=16)
        plt.legend(prop={'size':12})
        plt.tick_params(axis='both', labelsize=14)

    def plotvar_tvar_ratio(self):
        data = np.genfromtxt(self.var_out)
        L = data[:, 0]

        interp_tmv1 = interp1d(L, data[:, 1], kind='quadratic', bounds_error=False, fill_value=0.)
        interp_tmv2 = interp1d(L, data[:, 2], kind='quadratic', bounds_error=False, fill_value=0.)

        data21 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise)))
        interp_HO_1 = interp1d(data21[:, 0], data21[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)
    
        data22 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only.txt' % (self.name, str(self.cmb.lMin), str(self.cmb.lMaxT), str(self.cmb.lMaxP), str(self.beam), str(self.noise)))
        interp_HO_2 = interp1d(data22[:, 0], data22[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

        L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
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

    def plotvar_tvar_ratio_multexp(self, exp):
        lines = ["-", "--", "-."]
        # custom_lines = [Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for e in range(len(exp)):
            clexp = exp[e]
            data = np.genfromtxt('output/True_variance_individual_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
            L = data[:, 0]

            interp_tmv1 = interp1d(L, data[:, 1], kind='quadratic', bounds_error=False, fill_value=0.)
            interp_tmv2 = interp1d(L, data[:, 2], kind='quadratic', bounds_error=False, fill_value=0.)

            data21 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TT_EE_TE_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
            interp_HO_1 = interp1d(data21[:, 0], data21[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

            data22 = np.genfromtxt('output/HO02_covariance_%s_lmin%s_lmaxT%s_lmaxP%s_beam%s_noise%s_TB_EB_only.txt' % (clexp['name'], str(clexp['lMin']), str(clexp['lMaxT']), str(clexp['lMaxP']), str(clexp['beam']), str(clexp['noise_t'])))
            interp_HO_2 = interp1d(data22[:, 0], data22[:, -1], kind='quadratic', bounds_error=False, fill_value=0.)

            L_p = np.logspace(np.log10(1.), np.log10(5000.), 201, 10.)
            vart_var1 = interp_tmv1(L_p)/interp_HO_1(L_p)
            vart_var2 = interp_tmv2(L_p)/interp_HO_2(L_p)

            ax.plot(L_p, vart_var1, 'b', ls=lines[e],
                    label='%s TT-EE-TE' % (clexp['name']))
            ax.plot(L_p, vart_var2, 'r', ls=lines[e],
                     label='%s TB-EB' % (clexp['name']))
                
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
        ax.legend(prop={'size': 14}, frameon=False)
        ax.set_xscale('log')
        ax.set_xlim(1.0, clexp['lMaxP'])
        ax.set_ylim(0.8, 1.05)
        # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
        ax.set_ylabel(r'$N_{mv}^\mathrm{TMV}/N_{mv}^\mathrm{HO02}$',
                      fontsize=16)
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.tick_params(axis='both', labelsize=14)
        # ax.legend(custom_lines, ['TT-EE-TE', 'TB-EB'])
        # ax.add_artist(leg1)
        # ax.add_artist(leg2)
        # np.save('figures/multexp_vartrueind_varHO02ind.png')

    def plotvar_var_tvar_percent_multexp(self, exp):
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
        ax.legend(prop={'size': 14}, frameon=False)
        ax.set_xscale('log')
        ax.set_xlim(1.0, clexp['lMaxP'])
        ax.set_ylim(ymax=15)
        # ax.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
        ax.set_ylabel(r'$(N_{mv}^\mathrm{HO02}-N_{mv}^\mathrm{TMV}) \times 100/N_{mv}^\mathrm{HO02}$', fontsize=16)
        ax.set_xlabel(r'$L$', fontsize=16)
        ax.tick_params(axis='both', labelsize=14)
        # np.save('figures/multexp_percent_improvement_overHO02.png')


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
    # l_est.cmb.plot_cell()
    # """
    """
    LL = np.logspace(np.log10(1.), np.log10(4000+1.), 3, 10.)
    A1 = np.zeros(len(LL))
    for i in range(len(LL)):
        A1[i] = l_est.A_1(LL[i])
    """
    l_est.calc_tvar()

    print time()-time0

    # """
    l_est.interp_tvar()
    l_est.plot_tvar()
    # """
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
    """
    interp2 = interp1d(L, data[:, -1], kind='linear', bounds_error=False, fill_value='extrapolate')
    vart_var2 = data2[:, -1]/interp2(L2)
    plt.plot(L2, vart_var2, 'b')
    # """
    plt.xscale('log')
    plt.ylim(0.8, 1.05)
    plt.hlines(y=1, xmin=min(L_p), xmax=max(L_p))  # , color='k--')
    plt.ylabel(r'$N_{mv}^\mathrm{true}/N_{mv}^\mathrm{HO02}$', fontsize=16)
    plt.xlabel(r'$L$', fontsize=16)
    plt.legend(prop={'size':12})
    plt.tick_params(axis='both', labelsize=14)

    """
    plt.figure()
    data1 = np.genfromtxt('true_variance_lmin%s_lmaxT%s_lmaxP%s_newbasis.txt' % (str(cmb.lMin), str(cmb.lMaxT), str(cmb.lMaxP)))
    ratio = data[:, -1]/data1[:, -1]
    plt.semilogx(L, ratio)
    # """
    # print time()-time0
