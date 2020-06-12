import numpy as np

class prometheusChemistry():

    def __init__(self):

        self.mm = 4;
        self.kk = 6;
        self.ii = 2;
        self.one_atm   = 1.01325e5;
        self.one_third = 1.0 / 3.0;
        self.gas_constant = 8314.4621;
        self.big_number   = 1.0e300;

        self.wts = np.array( [ 2.805376e+01, 3.199880e+01, 4.400980e+01, 2.801040e+01, 1.801528e+01, 2.801348e+01 ] )

        return

    def get_mixture_specific_heat_mass(self,T,Y):

        nz,ny,nx = T.shape
        cp = np.zeros( [ nz,ny,nx ] )

        cp0_R = self.get_species_specific_heats_R( T )
        cp += Y[0] * cp0_R[0] / self.wts[0]
        cp += Y[1] * cp0_R[1] / self.wts[1]
        cp += Y[2] * cp0_R[2] / self.wts[2]
        cp += Y[3] * cp0_R[3] / self.wts[3]
        cp += Y[4] * cp0_R[4] / self.wts[4]
        cp += Y[5] * cp0_R[5] / self.wts[5]
        cp *= self.gas_constant

        return cp

    def get_species_specific_heats_R(self,T):

        tt0 = T
        tt1 = T * tt0
        tt2 = T * tt1
        tt3 = T * tt2
        tt4 = np.power( T, -1.0 )
        tt5 = tt4 * tt4

        nz,ny,nx = T.shape
        cp0_R = np.zeros( [ self.kk,nz,ny,nx ] )

        cp_high  = 2.036111e+00 + 1.464542e-02 * tt0 - 6.710779e-06 * tt1 + 1.472229e-09 * tt2 - 1.257061e-13 * tt3
        cp_low   = 3.959201e+00 - 7.570522e-03 * tt0 + 5.709903e-05 * tt1 - 6.915888e-08 * tt2 + 2.698844e-11 * tt3
        cp0_R[0] = np.where( tt0 < 1.000000e+03, cp_low, cp_high )

        cp_high  = 3.282538e+00 + 1.483088e-03 * tt0 - 7.579667e-07 * tt1 + 2.094706e-10 * tt2 - 2.167178e-14 * tt3
        cp_low   = 3.782456e+00 - 2.996734e-03 * tt0 + 9.847302e-06 * tt1 - 9.681295e-09 * tt2 + 3.243728e-12 * tt3
        cp0_R[1] = np.where( tt0 < 1.000000e+03, cp_low, cp_high )

        cp_high  = 3.857460e+00 + 4.414370e-03 * tt0 - 2.214814e-06 * tt1 + 5.234902e-10 * tt2 - 4.720842e-14 * tt3
        cp_low   = 2.356774e+00 + 8.984597e-03 * tt0 - 7.123563e-06 * tt1 + 2.459190e-09 * tt2 - 1.436995e-13 * tt3
        cp0_R[2] = np.where( tt0 < 1.000000e+03, cp_low, cp_high )

        cp_high  = 2.715186e+00 + 2.062527e-03 * tt0 - 9.988258e-07 * tt1 + 2.300530e-10 * tt2 - 2.036477e-14 * tt3
        cp_low   = 3.579533e+00 - 6.103537e-04 * tt0 + 1.016814e-06 * tt1 + 9.070059e-10 * tt2 - 9.044245e-13 * tt3
        cp0_R[3] = np.where( tt0 < 1.000000e+03, cp_low, cp_high )

        cp_high  = 3.033992e+00 + 2.176918e-03 * tt0 - 1.640725e-07 * tt1 - 9.704199e-11 * tt2 + 1.682010e-14 * tt3
        cp_low   = 4.198641e+00 - 2.036434e-03 * tt0 + 6.520402e-06 * tt1 - 5.487971e-09 * tt2 + 1.771978e-12 * tt3
        cp0_R[4] = np.where( tt0 < 1.000000e+03, cp_low, cp_high )

        cp_high  = 2.926640e+00 + 1.487977e-03 * tt0 - 5.684760e-07 * tt1 + 1.009704e-10 * tt2 - 6.753351e-15 * tt3
        cp_low   = 3.298677e+00 + 1.408240e-03 * tt0 - 3.963222e-06 * tt1 + 5.641515e-09 * tt2 - 2.444854e-12 * tt3
        cp0_R[5] = np.where( tt0 < 1.000000e+03, cp_low, cp_high )

        return cp0_R

    def get_species_enthalpies_RT(self,T):

        tt0 = T
        tt1 = T * tt0
        tt2 = T * tt1
        tt3 = T * tt2
        tt4 = np.power( T, -1.0 )
        tt5 = tt4 * tt4
        tt6 = np.log(tt0) * tt4

        nz,ny,nx = T.shape
        h0_RT = np.zeros( [ self.kk,nz,ny,nx ] )

        h_high  = 2.036111e+00 + 1.464542e-02 * 0.50 * tt0 - 6.710779e-06 * self.one_third * tt1 + 1.472229e-09 * 0.25 * tt2 - 1.257061e-13 * 0.20 * tt3 + 4.939886e+03 * tt4
        h_low   = 3.959201e+00 - 7.570522e-03 * 0.50 * tt0 + 5.709903e-05 * self.one_third * tt1 - 6.915888e-08 * 0.25 * tt2 + 2.698844e-11 * 0.20 * tt3 + 5.089776e+03 * tt4
        h_RT[0] = np.where( tt0 < 1.000000e+03, h_low, h_high )

        h_high  = 3.282538e+00 + 1.483088e-03 * 0.50 * tt0 - 7.579667e-07 * self.one_third * tt1 + 2.094706e-10 * 0.25 * tt2 - 2.167178e-14 * 0.20 * tt3 - 1.088458e+03 * tt4
        h_low   = 3.782456e+00 - 2.996734e-03 * 0.50 * tt0 + 9.847302e-06 * self.one_third * tt1 - 9.681295e-09 * 0.25 * tt2 + 3.243728e-12 * 0.20 * tt3 - 1.063944e+03 * tt4
        h_RT[1] = np.where( tt0 < 1.000000e+03, h_low, h_high )

        h_high  = 3.857460e+00 + 4.414370e-03 * 0.50 * tt0 - 2.214814e-06 * self.one_third * tt1 + 5.234902e-10 * 0.25 * tt2 - 4.720842e-14 * 0.20 * tt3 - 4.875917e+04 * tt4
        h_low   = 2.356774e+00 + 8.984597e-03 * 0.50 * tt0 - 7.123563e-06 * self.one_third * tt1 + 2.459190e-09 * 0.25 * tt2 - 1.436995e-13 * 0.20 * tt3 - 4.837197e+04 * tt4
        h_RT[2] = np.where( tt0 < 1.000000e+03, h_low, h_high )

        h_high  = 2.715186e+00 + 2.062527e-03 * 0.50 * tt0 - 9.988258e-07 * self.one_third * tt1 + 2.300530e-10 * 0.25 * tt2 - 2.036477e-14 * 0.20 * tt3 - 1.415187e+04 * tt4
        h_low   = 3.579533e+00 - 6.103537e-04 * 0.50 * tt0 + 1.016814e-06 * self.one_third * tt1 + 9.070059e-10 * 0.25 * tt2 - 9.044245e-13 * 0.20 * tt3 - 1.434409e+04 * tt4
        h_RT[3] = np.where( tt0 < 1.000000e+03, h_low, h_high )

        h_high  = 3.033992e+00 + 2.176918e-03 * 0.50 * tt0 - 1.640725e-07 * self.one_third * tt1 - 9.704199e-11 * 0.25 * tt2 + 1.682010e-14 * 0.20 * tt3 - 3.000430e+04 * tt4
        h_low   = 4.198641e+00 - 2.036434e-03 * 0.50 * tt0 + 6.520402e-06 * self.one_third * tt1 - 5.487971e-09 * 0.25 * tt2 + 1.771978e-12 * 0.20 * tt3 - 3.029373e+04 * tt4
        h_RT[4] = np.where( tt0 < 1.000000e+03, h_low, h_high )

        h_high  = 2.926640e+00 + 1.487977e-03 * 0.50 * tt0 - 5.684760e-07 * self.one_third * tt1 + 1.009704e-10 * 0.25 * tt2 - 6.753351e-15 * 0.20 * tt3 - 9.227977e+02 * tt4
        h_low   = 3.298677e+00 + 1.408240e-03 * 0.50 * tt0 - 3.963222e-06 * self.one_third * tt1 + 5.641515e-09 * 0.25 * tt2 - 2.444854e-12 * 0.20 * tt3 - 1.020900e+03 * tt4
        h_RT[5] = np.where( tt0 < 1.000000e+03, h_low, h_high )

        return h0_RT

    def get_species_entropies_R(self,T):

        tt0 = T
        tt1 = T * tt0
        tt2 = T * tt1
        tt3 = T * tt2
        tt4 = np.power( T, -1.0 )
        tt5 = tt4 * tt4
        tt6 = np.log(tt0)

        nz,ny,nx = T.shape
        s0_R = np.zeros( [ self.kk,nz,ny,nx ] )

        s_high  = 2.036111e+00 * tt6 + 1.464542e-02 * tt0 - 6.710779e-06 * 0.50 * tt1 + 1.472229e-09 * self.one_third * tt2 - 1.257061e-13 * 0.25 * tt3 + 1.030537e+01
        s_low   = 3.959201e+00 * tt6 - 7.570522e-03 * tt0 + 5.709903e-05 * 0.50 * tt1 - 6.915888e-08 * self.one_third * tt2 + 2.698844e-11 * 0.25 * tt3 + 4.097331e+00
        s0_R[0] = np.where( tt0 < 1.000000e+03, s_low, s_high )

        s_high  = 3.282538e+00 * tt6 + 1.483088e-03 * tt0 - 7.579667e-07 * 0.50 * tt1 + 2.094706e-10 * self.one_third * tt2 - 2.167178e-14 * 0.25 * tt3 + 5.453231e+00
        s_low   = 3.782456e+00 * tt6 - 2.996734e-03 * tt0 + 9.847302e-06 * 0.50 * tt1 - 9.681295e-09 * self.one_third * tt2 + 3.243728e-12 * 0.25 * tt3 + 3.657676e+00
        s0_R[1] = np.where( tt0 < 1.000000e+03, s_low, s_high )

        s_high  = 3.857460e+00 * tt6 + 4.414370e-03 * tt0 - 2.214814e-06 * 0.50 * tt1 + 5.234902e-10 * self.one_third * tt2 - 4.720842e-14 * 0.25 * tt3 + 2.271638e+00
        s_low   = 2.356774e+00 * tt6 + 8.984597e-03 * tt0 - 7.123563e-06 * 0.50 * tt1 + 2.459190e-09 * self.one_third * tt2 - 1.436995e-13 * 0.25 * tt3 + 9.901052e+00
        s0_R[2] = np.where( tt0 < 1.000000e+03, s_low, s_high )

        s_high  = 2.715186e+00 * tt6 + 2.062527e-03 * tt0 - 9.988258e-07 * 0.50 * tt1 + 2.300530e-10 * self.one_third * tt2 - 2.036477e-14 * 0.25 * tt3 + 7.818688e+00
        s_low   = 3.579533e+00 * tt6 - 6.103537e-04 * tt0 + 1.016814e-06 * 0.50 * tt1 + 9.070059e-10 * self.one_third * tt2 - 9.044245e-13 * 0.25 * tt3 + 3.508409e+00
        s0_R[3] = np.where( tt0 < 1.000000e+03, s_low, s_high )

        s_high  = 3.033992e+00 * tt6 + 2.176918e-03 * tt0 - 1.640725e-07 * 0.50 * tt1 - 9.704199e-11 * self.one_third * tt2 + 1.682010e-14 * 0.25 * tt3 + 4.966770e+00
        s_low   = 4.198641e+00 * tt6 - 2.036434e-03 * tt0 + 6.520402e-06 * 0.50 * tt1 - 5.487971e-09 * self.one_third * tt2 + 1.771978e-12 * 0.25 * tt3 - 8.490322e-01
        s0_R[4] = np.where( tt0 < 1.000000e+03, s_low, s_high )

        s_high  = 2.926640e+00 * tt6 + 1.487977e-03 * tt0 - 5.684760e-07 * 0.50 * tt1 + 1.009704e-10 * self.one_third * tt2 - 6.753351e-15 * 0.25 * tt3 + 5.980528e+00
        s_low   = 3.298677e+00 * tt6 + 1.408240e-03 * tt0 - 3.963222e-06 * 0.50 * tt1 + 5.641515e-09 * self.one_third * tt2 - 2.444854e-12 * 0.25 * tt3 + 3.950372e+00
        s0_R[5] = np.where( tt0 < 1.000000e+03, s_low, s_high )

        return s0_R

    def get_species_gibbs_RT(self,T):

        h0_RT = self.get_enthalpies_RT( T )
        s0_R  = self.get_entropies_R( T )
        g0_RT = h0_RT - s0_R

        return g0_RT

    def get_equilibrium_constants(self,T):

        RT = self.gas_constant * T
        C0 = self.one_atm * np.power( RT, -1.0 )

        nz,ny,nx = T.shape
        keq = np.zeros( [ self.ii,nz,ny,nx ] )

        g0_RT = self.get_gibbs_RT( T )

        keq[1] =  C0 + ( g0_RT[2] ) - ( g0_RT[3] + g0_RT[1] )

        return keq

