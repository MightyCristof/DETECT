# scripts/luminosity_functions
from scripts.common import const, np, u
import resources.parameters as params


#########################################
#
# Quasar Luminosity Function
#
#########################################
class QuasarLuminosityFunction:
    def __init__(self, redshift, luminosity, prescription='Shen+2020', fit_type='FitA', z_ref=2.0):
        self.redshift = redshift
        self.luminosity =luminosity
        self.prescription = prescription
        self.fit_type = fit_type
        self.z_ref = z_ref
        self.fit_params = self.load_fit_parameters()
        self.log_phi = self.construct_phi_matrix()
        
    def load_fit_parameters(self):
        # Shen+2020 Table 4
        return {
            'FitA': {
                'a0': (0.8569, 0.0247, -0.0253),
                'a1': (-0.2614, 0.0162, -0.016),
                'a2': (0.0200, 0.0011, -0.0011),
                'b0': (2.5375, 0.0177, -0.0187),
                'b1': (-1.0425, 0.0164, -0.018),
                'b2': (1.1201, 0.0199, -0.0207),
                'c0': (13.0088, 0.0090, -0.009),
                'c1': (-0.5759, 0.0018, -0.002),
                'c2': (0.4554, 0.0028, -0.0027),
                'd0': (-3.5426, 0.0235, -0.020),
                'd1': (-0.3936, 0.0070, -0.007),
            },
            'FitB': {
                'a0': (0.3653, 0.0115, -0.0114),
                'a1': (-0.6006, 0.0422, -0.0417),
                'a2': (None, None, None),
                'b0': (2.4709, 0.0163, -0.0169),
                'b1': (-0.9963, 0.0167, -0.016),
                'b2': (1.0716, 0.0180, -0.0181),
                'c0': (12.9656, 0.0092, -0.008),
                'c1': (-0.5758, 0.0020, -0.001),
                'c2': (0.4698, 0.0025, -0.0026),
                'd0': (-3.6276, 0.0209, -0.020),
                'd1': (-0.3444, 0.0063, -0.006),
            },
            # ... other fits if needed
        }[self.fit_type]

    def load_model_parameters(self, z):
        # Shen+2020 Equation 14
        def T(n, x):
            if n == 0:
                return 1
            elif n == 1:
                return x
            elif n == 2:
                return 2 * x**2 - 1
            else:
                raise ValueError("Invalid value for n. It should be 0, 1, or 2.")
        if self.fit_type == 'FitA':
            gamma1 = (
                self.fit_params['a0'][0] * T(0, 1+z) +
                self.fit_params['a1'][0] * T(1, 1+z) +
                self.fit_params['a2'][0] * T(2, 1+z)
            )
        elif self.fit_type == 'FitB':
            gamma1 = self.fit_params['a0'][0] * ((1+z) / (1+self.z_ref)) ** self.fit_params['a1'][0]
        gamma2 = (
            2 * self.fit_params['b0'][0] / (
                ((1+z) / (1+self.z_ref)) ** self.fit_params['b1'][0] +
                ((1+z) / (1+self.z_ref)) ** self.fit_params['b2'][0]
            )
        )
        logL_star = (
            2 * self.fit_params['c0'][0] / (
                ((1+z) / (1+self.z_ref)) ** self.fit_params['c1'][0] +
                ((1+z) / (1+self.z_ref)) ** self.fit_params['c2'][0]
            )
        )
        logPhi_star = (
            self.fit_params['d0'][0] * T(0, 1+z) +
            self.fit_params['d1'][0] * T(1, 1+z)
        )
        return {
            'gamma1': gamma1,
            'gamma2': gamma2,
            'log_ell_star': logL_star,
            'log_phi_star': logPhi_star
        }

    def calculate_phi(self, z, luminosities):
        model_params = self.load_model_parameters(z)
        log_ell_star = model_params['log_ell_star']
        log_phi_star = model_params['log_phi_star']
        gamma1 = model_params['gamma1']
        gamma2 = model_params['gamma2']

        cgs_to_solar = const.L_sun.to(u.erg / u.s).value
        l_bol = 10**luminosities / cgs_to_solar
        phi_star = 10.0**log_phi_star
        l_star = 10.0**log_ell_star
        l_over_l_star = l_bol/l_star        
        # Shen+2020 Equation 11
        log_phi = np.log10(phi_star / (l_over_l_star**gamma1 + l_over_l_star**gamma2))
        return log_phi
    
    def construct_phi_matrix(self):
        phi = np.zeros((len(self.luminosity), len(self.redshift)))
        for i in range(len(self.redshift)):
            red = self.redshift[i]
            phi[:,i] = self.calculate_phi(red, self.luminosity)
        return phi

    def get_phi(self, z, lum=None):
        # Check if z is within the bounds
        if z < self.redshift.min() or z > self.redshift.max():
            raise ValueError(f"Redshift {z} is out of bounds [{self.redshift.min():.3f}, {self.redshift.max():.3f}]")

        if lum is None:
            # Return phi slice closest to specified z
            index_z = np.abs(self.redshift - z).argmin()
            return self.log_phi[:, index_z]
        else:
            # Check if luminosity is within the bounds
            if np.any(lum < self.luminosity.min()) or np.any(lum > self.luminosity.max()):
                raise ValueError(f"Luminosity {lum} is out of bounds [{self.luminosity.min():.3f}, {self.luminosity.max():.3f}]")
            # Find closest indices
            #index_lum = np.abs(self.luminosity - lum).argmin()
            index_lum = np.array([np.abs(self.luminosity - li).argmin() for li in lum])
            #index_z = np.abs(self.redshift - z).argmin()
            #index_z = np.array([np.abs(self.redshift - zi).argmin() for zi in z])
            index_z = np.array([np.abs(self.redshift - zi).argmin() for zi in np.atleast_1d(z)])
            return self.log_phi[index_lum, index_z]
    
    def sample_luminosity(self, z, ndraw=1):
        phi = self.get_phi(z)
        # Shift Phi to positive values
        phi_shift = phi - np.min(phi)
        phi_pdf = phi_shift / np.sum(phi_shift)
        phi_cdf = np.cumsum(phi_pdf)
        random_draw = np.random.rand(ndraw)
        samp = np.interp(random_draw, phi_cdf, self.luminosity)
        return samp
    

#########################################
#
# Galaxy Luminosity Function
#
#########################################
#class SchechterLuminosityFunction:
#    def __init__(self, phi_star=0.012, l_star=1e10, alpha=-1.2):
#        self.phi_star = phi_star*(params.H0/100)**(-3)
#        self.l_star = l_star
#        self.alpha = alpha
#    def compute_luminosity_function(self, lum):
#        return self.phi_star * (lum / self.l_star)**self.alpha * np.exp(-lum / self.l_star) / self.l_star
