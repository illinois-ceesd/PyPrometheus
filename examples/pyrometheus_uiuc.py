import numpy as np

class Thermochemistry:
    model_name    = 'uiuc.cti'
    num_elements  = 4
    num_species   = 7
    num_reactions = 3
    num_falloff   = 0

    one_atm = 101325.0
    gas_constant = 8314.46261815324
    big_number = 1.0e300

    species_names = ['C2H4', 'O2', 'CO2', 'CO', 'H2O', 'H2', 'N2']
    species_indices = {'C2H4': 0, 'O2': 1, 'CO2': 2, 'CO': 3, 'H2O': 4, 'H2': 5, 'N2': 6}

    wts = np.array([28.054, 31.998, 44.009, 28.009999999999998, 18.015, 2.016, 28.014])
    iwts = 1/wts

    def species_name(self, species_index):
        return self.species_name[species_index]

    def species_index(self, species_name):
        return self.species_indices[species_name]

    def get_specific_gas_constant(self, Y):
        return self.gas_constant * np.dot( self.iwts, Y )

    def get_density(self, p, T, Y):
        mmw = self.get_mix_molecular_weight( Y )
        RT  = self.gas_constant * T
        return p * mmw / RT

    def get_pressure(self, rho, T, Y):
        mmw = self.get_mix_molecular_weight( Y )
        RT  = self.gas_constant * T
        return rho * RT / mmw

    def get_mix_molecular_weight(self, Y):
        return 1/np.dot( self.iwts, Y )

    def get_concentrations(self, rho, Y):
        return self.iwts * rho * Y

    def get_mixture_specific_heat_cp_mass(self, T, Y):
        return self.gas_constant * np.sum(
            self.get_species_specific_heats_R(T)* Y * self.iwts)

    def get_mixture_specific_heat_cv_mass(self, T, Y):
        cp0_R = self.get_species_specific_heats_R( T ) - 1.0
        return self.gas_constant * np.sum(Y * cp0_R * self.iwts)

    def get_mixture_enthalpy_mass(self, T, Y):
        h0_RT = self.get_species_enthalpies_RT( T )
        return self.gas_constant * T * np.sum(Y * h0_RT * self.iwts)

    def get_mixture_internal_energy_mass(self, T, Y):
        e0_RT = self.get_species_enthalpies_RT( T ) - 1.0
        return self.gas_constant * T * np.sum(Y * e0_RT * self.iwts)

    def get_species_specific_heats_R(self, T):
        return np.array([
            np.where(T > 1000.0, 2.03611116 + 0.0146454151*T + -6.71077915e-06*T**2 + 1.47222923e-09*T**3 + -1.25706061e-13*T**4, 3.95920148 + -0.00757052247*T + 5.70990292e-05*T**2 + -6.91588753e-08*T**3 + 2.69884373e-11*T**4),
            np.where(T > 1000.0, 3.28253784 + 0.00148308754*T + -7.57966669e-07*T**2 + 2.09470555e-10*T**3 + -2.16717794e-14*T**4, 3.78245636 + -0.00299673416*T + 9.84730201e-06*T**2 + -9.68129509e-09*T**3 + 3.24372837e-12*T**4),
            np.where(T > 1000.0, 3.85746029 + 0.00441437026*T + -2.21481404e-06*T**2 + 5.23490188e-10*T**3 + -4.72084164e-14*T**4, 2.35677352 + 0.00898459677*T + -7.12356269e-06*T**2 + 2.45919022e-09*T**3 + -1.43699548e-13*T**4),
            np.where(T > 1000.0, 2.71518561 + 0.00206252743*T + -9.98825771e-07*T**2 + 2.30053008e-10*T**3 + -2.03647716e-14*T**4, 3.57953347 + -0.00061035368*T + 1.01681433e-06*T**2 + 9.07005884e-10*T**3 + -9.04424499e-13*T**4),
            np.where(T > 1000.0, 3.03399249 + 0.00217691804*T + -1.64072518e-07*T**2 + -9.7041987e-11*T**3 + 1.68200992e-14*T**4, 4.19864056 + -0.0020364341*T + 6.52040211e-06*T**2 + -5.48797062e-09*T**3 + 1.77197817e-12*T**4),
            np.where(T > 1000.0, 3.3372792 + -4.94024731e-05*T + 4.99456778e-07*T**2 + -1.79566394e-10*T**3 + 2.00255376e-14*T**4, 2.34433112 + 0.00798052075*T + -1.9478151e-05*T**2 + 2.01572094e-08*T**3 + -7.37611761e-12*T**4),
            np.where(T > 1000.0, 2.92664 + 0.0014879768*T + -5.68476e-07*T**2 + 1.0097038e-10*T**3 + -6.753351e-15*T**4, 3.298677 + 0.0014082404*T + -3.963222e-06*T**2 + 5.641515e-09*T**3 + -2.444854e-12*T**4),
            ])

    def get_species_enthalpies_RT(self, T):
        return np.array([
            np.where(T > 1000.0, 2.03611116 + 0.00732270755*T + -2.2369263833333335e-06*T**2 + 3.680573075e-10*T**3 + -2.51412122e-14*T**4 + 4939.88614 / T, 3.95920148 + -0.003785261235*T + 1.9033009733333333e-05*T**2 + -1.7289718825e-08*T**3 + 5.3976874600000004e-12*T**4 + 5089.77593 / T),
            np.where(T > 1000.0, 3.28253784 + 0.00074154377*T + -2.526555563333333e-07*T**2 + 5.236763875e-11*T**3 + -4.33435588e-15*T**4 + -1088.45772 / T, 3.78245636 + -0.00149836708*T + 3.282434003333333e-06*T**2 + -2.4203237725e-09*T**3 + 6.48745674e-13*T**4 + -1063.94356 / T),
            np.where(T > 1000.0, 3.85746029 + 0.00220718513*T + -7.382713466666667e-07*T**2 + 1.30872547e-10*T**3 + -9.44168328e-15*T**4 + -48759.166 / T, 2.35677352 + 0.004492298385*T + -2.3745208966666665e-06*T**2 + 6.14797555e-10*T**3 + -2.8739909599999997e-14*T**4 + -48371.9697 / T),
            np.where(T > 1000.0, 2.71518561 + 0.001031263715*T + -3.329419236666667e-07*T**2 + 5.7513252e-11*T**3 + -4.07295432e-15*T**4 + -14151.8724 / T, 3.57953347 + -0.00030517684*T + 3.3893811e-07*T**2 + 2.26751471e-10*T**3 + -1.808848998e-13*T**4 + -14344.086 / T),
            np.where(T > 1000.0, 3.03399249 + 0.00108845902*T + -5.469083933333333e-08*T**2 + -2.426049675e-11*T**3 + 3.36401984e-15*T**4 + -30004.2971 / T, 4.19864056 + -0.00101821705*T + 2.17346737e-06*T**2 + -1.371992655e-09*T**3 + 3.54395634e-13*T**4 + -30293.7267 / T),
            np.where(T > 1000.0, 3.3372792 + -2.470123655e-05*T + 1.6648559266666665e-07*T**2 + -4.48915985e-11*T**3 + 4.00510752e-15*T**4 + -950.158922 / T, 2.34433112 + 0.003990260375*T + -6.4927169999999995e-06*T**2 + 5.03930235e-09*T**3 + -1.4752235220000002e-12*T**4 + -917.935173 / T),
            np.where(T > 1000.0, 2.92664 + 0.0007439884*T + -1.8949200000000001e-07*T**2 + 2.5242595e-11*T**3 + -1.3506701999999999e-15*T**4 + -922.7977 / T, 3.298677 + 0.0007041202*T + -1.3210739999999999e-06*T**2 + 1.41037875e-09*T**3 + -4.889707999999999e-13*T**4 + -1020.8999 / T),
            ])

    def get_species_entropies_R(self, T):
        return np.array([
                np.where(T > 1000.0, 2.03611116*np.log(T) + 0.0146454151*T + -3.355389575e-06*T**2 + 4.907430766666667e-10*T**3 + -3.142651525e-14*T**4 + 10.3053693, 3.95920148*np.log(T) + -0.00757052247*T + 2.85495146e-05*T**2 + -2.3052958433333332e-08*T**3 + 6.747109325e-12*T**4 + 4.09733096),
                np.where(T > 1000.0, 3.28253784*np.log(T) + 0.00148308754*T + -3.789833345e-07*T**2 + 6.982351833333333e-11*T**3 + -5.41794485e-15*T**4 + 5.45323129, 3.78245636*np.log(T) + -0.00299673416*T + 4.923651005e-06*T**2 + -3.2270983633333334e-09*T**3 + 8.109320925e-13*T**4 + 3.65767573),
                np.where(T > 1000.0, 3.85746029*np.log(T) + 0.00441437026*T + -1.10740702e-06*T**2 + 1.7449672933333335e-10*T**3 + -1.18021041e-14*T**4 + 2.27163806, 2.35677352*np.log(T) + 0.00898459677*T + -3.561781345e-06*T**2 + 8.197300733333333e-10*T**3 + -3.5924887e-14*T**4 + 9.90105222),
                np.where(T > 1000.0, 2.71518561*np.log(T) + 0.00206252743*T + -4.994128855e-07*T**2 + 7.6684336e-11*T**3 + -5.0911929e-15*T**4 + 7.81868772, 3.57953347*np.log(T) + -0.00061035368*T + 5.08407165e-07*T**2 + 3.023352946666667e-10*T**3 + -2.2610612475e-13*T**4 + 3.50840928),
                np.where(T > 1000.0, 3.03399249*np.log(T) + 0.00217691804*T + -8.2036259e-08*T**2 + -3.2347329e-11*T**3 + 4.2050248e-15*T**4 + 4.9667701, 4.19864056*np.log(T) + -0.0020364341*T + 3.260201055e-06*T**2 + -1.82932354e-09*T**3 + 4.429945425e-13*T**4 + -0.849032208),
                np.where(T > 1000.0, 3.3372792*np.log(T) + -4.94024731e-05*T + 2.49728389e-07*T**2 + -5.985546466666667e-11*T**3 + 5.0063844e-15*T**4 + -3.20502331, 2.34433112*np.log(T) + 0.00798052075*T + -9.7390755e-06*T**2 + 6.7190698e-09*T**3 + -1.8440294025e-12*T**4 + 0.683010238),
                np.where(T > 1000.0, 2.92664*np.log(T) + 0.0014879768*T + -2.84238e-07*T**2 + 3.3656793333333334e-11*T**3 + -1.68833775e-15*T**4 + 5.980528, 3.298677*np.log(T) + 0.0014082404*T + -1.981611e-06*T**2 + 1.8805050000000002e-09*T**3 + -6.112135e-13*T**4 + 3.950372),
            ])

    def get_species_gibbs_RT(self, T):
        h0_RT = self.get_species_enthalpies_RT(T)
        s0_R  = self.get_species_entropies_R(T)
        return h0_RT - s0_R

    def get_equilibrium_constants(self, T):
        RT = self.gas_constant * T
        C0 = np.log( self.one_atm / RT )

        g0_RT = self.get_species_gibbs_RT( T )
        return np.array([
                    -86*T,
                    g0_RT[2] + -1*-0.5*C0 + -1*(g0_RT[3] + 0.5*g0_RT[1]),
                    g0_RT[4] + -1*-0.5*C0 + -1*(g0_RT[5] + 0.5*g0_RT[1]),
            ])

    def get_temperature(self, H_or_E, T_guess, Y, do_energy=False):
        if do_energy == False:
            pv_fun = self.get_mixture_specific_heat_cp_mass
            he_fun = self.get_mixture_enthalpy_mass
        else:
            pv_fun = self.get_mixture_specific_heat_cv_mass
            he_fun = self.get_mixture_internal_energy_mass

        num_iter = 500
        tol = 1.0e-6
        T_i = T_guess
        dT = 1.0
        F  = H_or_E
        J  = 0.0

        for iter in range( 0, num_iter ):
            F    -= he_fun( T_i, Y )
            J    -= pv_fun( T_i, Y )
            dT    = - F / J
            T_i  += dT
            if np.abs( dT ) < tol:
                T = T_i
                break
            F = H_or_E
            J = 0.0

        T = T_i

        return T

    def get_falloff_rates(self, T, C, k_fwd):
        k_high = np.array([
        ])

        k_low = np.array([
        ])

        reduced_pressure = np.array([
        ])

        falloff_center = np.array([
        ])

        falloff_function = np.array([
        ])*reduced_pressure/(1+reduced_pressure)

        return

    def get_fwd_rate_coefficients(self, T, C):
        k_fwd = np.array([
            np.exp(26.594857854425133 + -1*(17864.293439206183 / T)),
            np.exp(12.693776816787125 + 0.7*np.log(T) + -1*(6038.634401985189 / T)),
            np.exp(18.302572655472037 + -1*(17612.683672456802 / T)),
        ])


        return k_fwd

    def get_net_rates_of_progress(self, T, C):
        k_fwd = self.get_fwd_rate_coefficients(T, C)
        log_k_eq = self.get_equilibrium_constants(T)
        k_eq = np.where(np.exp(log_k_eq) < self.big_number,
            np.exp(log_k_eq), self.big_number)
        return np.array([
                k_fwd[0]*C[0]**0.5*C[1]**0.65,
                k_fwd[1]*(C[3]*C[1]**0.5 + -1*k_eq[1]*C[2]),
                k_fwd[2]*(C[5]*C[1]**0.5 + -1*k_eq[2]*C[4]),
            ])

    def get_net_production_rates(self, rho, T, Y):
        C = self.get_concentrations(rho, Y)
        r_net = self.get_net_rates_of_progress(T, C)
        return np.array([
                -1*r_net[0],
                -1*(r_net[0] + 0.5*r_net[1] + 0.5*r_net[2]),
                r_net[1],
                2.0*r_net[0] + -1*r_net[1],
                r_net[2],
                2.0*r_net[0] + -1*r_net[2],
                0,
            ])
