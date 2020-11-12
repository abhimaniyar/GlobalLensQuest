# evaluating the inverse of M1 and M2 mattrices at several points to speed
# up the process


import numpy as np


class sample_matrix_inversions(object):

    def __init__(self, Cell_cmb):
        self.cmb = Cell_cmb

    def M1_inv(self, l_1, l_2):
        # nl1 = len(l_1[:, 0])
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

        inv_m1[:, 2, 3] = inv_m1[:, 3, 1] = 0.5*(self.cmb.totalEE(l_1)*self.cmb.totalTT(l_2) -
                                           self.cmb.totalTT(l_1)*self.cmb.totalEE(l_2))
        return inv_m1/det[:, None, None]

    def eval_inv(self, lmin, lmax):
        # lens = np.loadtxt('../CAMB/manu_lensedCls.dat')
        # lmin = 30.
        # lmax = 1e4
        l_1 = np.linspace(lmin, lmax, int(lmax-lmin+1))
        l_2 = np.linspace(lmin, lmax, int(lmax-lmin+1))
        # m1inv = np.zeros((len(l_2), 4, 4))
        M1inv = np.zeros((len(l_1), len(l_2), 4, 4))
        for i in range(len(l_1)):
            m1inv = self.M1_inv(l_1[i], l_2)
            M1inv[i, :] = m1inv
        np.save('m1_inverse_precalculated.npy', M1inv)
        np.savez_compressed('m1_inverse_precalculated.npz', M1inv)


if __name__ == '__main__':
    import time
    from cell_cmb import *
    beam = 1.  # 7.
    noise_t = 1.  # 27.
    noise_p = 1.*np.sqrt(2)  # 40.*np.sqrt(2)
    lMinT = 30.
    lMaxT = 10.e3
    lMaxP = 10.e3

    time0 = time()
    cmb = Cell_cmb(beam=beam, noise_t=noise_t, noise_p=noise_p, lMin=lMinT,
                   lMaxT=lMaxT, lMaxP=lMaxP)
    # cmb = Cell_cmb(beam=7., noise_t=27., noise_p=40*np.sqrt(2.), lMin=30.,
    #                lMaxT=3.e3, lMaxP=5.e3)
    l_est = sample_matrix_inversions(cmb)

    l_est.eval_inv(lMinT, lMaxT)
    print time()-time0
