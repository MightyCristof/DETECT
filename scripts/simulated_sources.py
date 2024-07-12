# scripts/simulated_source.py
from scripts.common import const, cp, np, u
from scripts.global_variables import cosmology, sed, bp
import pdb


#########################################
#
# Active Galactic Nuclei
#
#########################################
class ActiveGalacticNuclei:
    def __init__(self, rest_frequency, dluminosity, luminosity=None, ebv=0, galaxy=None):
        self.luminosity = self.set_bolometric_luminosity(luminosity)
        self.ebv = ebv
        self.eddington_ratio = self.assign_eddington_ratio(method="Suh+2015")
        self.mass = self.calculate_bh_mass()
        self.eddington_limit = self.calculate_eddington_limit()
        self.norm = self.normalize_sed(rest_frequency, dluminosity)
        self.sed_unabsorbed = self.construct_sed()
        self.sed = self.attenuate_sed()

    def set_bolometric_luminosity(self, lum):
        if lum is not None:
            return 10**lum
        else:
            raise NotImplementedError("Galaxy-to-AGN properties not yet implemented. Please define AGN bolometric luminosity.")

    def assign_eddington_ratio(self, method="Suh+2015"):
        # Suh+2015
        if method == "Suh+2015":
            # Generate bounded normal distribution
            def generate_bounded_normal(size, mean=-0.6, stddev=0.8, bounds=[-3.0,0.0]):
                data = np.random.normal(loc=mean, scale=stddev, size=size)

                if bounds is not None:
                    lower_bound, upper_bound = bounds
                    iout = np.where((data < lower_bound) | (data > upper_bound))[0]
                    if len(iout) > 0:
                        data[iout] = generate_bounded_normal(len(iout), mean=mean, stddev=stddev, bounds=bounds)
                    return data
                
            edd_ratio = (generate_bounded_normal(1))[0]
            return 10**edd_ratio
        else:
        # Implement Bernhard+2018?
            raise NotImplementedError(f"Alternative Eddington ratio method '{method}' not yet implemented. Use method='Suh+2015'.")
    
    def calculate_bh_mass(self):
        mass = self.luminosity / (1.26e38 * self.eddington_ratio)
        return mass

    def calculate_eddington_limit(self):
        edd_limit = 1.26e38 * self.mass
        return edd_limit

    def normalize_sed_components(self, coefficients):
        if len(coefficients) == 1:
            key = next(iter(coefficients))
            coefficients[key] = 1.0
        elif len(coefficients) > 1:
            raise NotImplementedError("Multiple component AGN SED not yet implemented.")
        else:
            raise ValueError("Coefficients dictionary is empty.")
        return coefficients

    def normalize_sed(self, rest_freq, dlum):
        coefficients = {component: None for component in sed.flux['agn'].keys()}
        coefficients = self.normalize_sed_components(coefficients)
        flux = sum(sed.flux['agn'][key] * coefficients[key] for key in coefficients)
        integrated_flux = np.trapz(flux[::-1], x=rest_freq[::-1])
        integrated_lum = 4 * np.pi * (dlum**2) * integrated_flux
        sed_normalization = self.luminosity / integrated_lum
        sed_normalizations = {key: coefficients[key] * sed_normalization for key in coefficients}
        return sed_normalizations
        
    def construct_sed(self):
        flux = {key: self.norm[key] * sed.flux['agn'][key] for key in self.norm}
        return flux

    def attenuate_sed(self):
        flux = cp.deepcopy(self.sed_unabsorbed)
        if 'agn' in sed.ext:
            for key in sed.ext['agn']:
                flux[key] *= 10.0**(-0.4 * self.ebv * sed.ext['agn'][key])
        return flux
    

#########################################
#
# Galaxy
#
#########################################
class Galaxy:
    def __init__(self, redshift, rest_frequency, dluminosity, luminosity=None, mass=None, ebv=0, agn=None):
        self.host_galaxy = self.is_host_galaxy(agn)
        self.ebv = ebv
        self.mass = self.calculate_stellar_mass(redshift, mass, agn)
        self.mass_to_light = self.calculate_ml_ratio()
        self.luminosity = self.set_bolometric_luminosity(luminosity)
        self.norm = self.normalize_sed(rest_frequency, dluminosity)
        self.sed_unabsorbed = self.construct_sed()
        self.sed = self.attenuate_sed()

    def is_host_galaxy(self, agn_instance=None):
        if agn_instance is not None:
            return True
        else:
            return False

    def scale_mass_ratio(self, z):
        # Correct for AGN mass overdensity 
        # Pucci+2024 Equation 2
        # ...which follows Barkana & Loeb (2001)
        # ...with definitions from Bryan & Norman 1998)
        # Ωzm === Critial matter density at redshift. Implemented in astropy.cosmology.FlatLambdaCDM
        # ∆c  === Ratio mean density / critical density at redshift
        # d   === Critial matter density at collapse redshift. Defined as "x" in Bryan & Norman (1998)
        # E   === BH mass overdensity redshift evolution 
        def calculate_delta_c(omega_m):
            d = omega_m - 1
            return 18*np.pi**2 + 82*d - 39*d**2
        
        def calculate_xi(omega_m, delta_c):
            return (cosmology.Om0/omega_m) * (delta_c/(18*np.pi**2))
        
        xi = []
        for z_i in [0, z]:
            omega_m = cosmology.Om(z_i)
            delta_c = calculate_delta_c(omega_m)
            xi.append(calculate_xi(omega_m, delta_c))

        big_e = (xi[1]**(5/6)) * ((1+z)**(5/2)) / (xi[0]**(5/6))
        return np.log10(big_e)

    def calculate_stellar_mass(self, z, mass=None, agn_instance=None):
        if self.host_galaxy is True and mass is None:
            # Define M_* through M_BH
            # Local Relation: Reines & Volonteri 2015, Equation 4
            alpha = np.random.normal(7.45, 0.08)
            beta = np.random.normal(1.05, 0.11)
            # Correct for redshift evolution on BH mass
            z_evolv = np.log10(self.scale_mass_ratio(z))
            m_stellar = 10**((np.log10(agn_instance.mass)-z_evolv-alpha)/beta + 11)
            return m_stellar
        else: 
            raise NotImplementedError("Galaxy-to-BH mass not yet implemented. Please define BH mass.")

    def calculate_ml_ratio(self):
        # Kauffmann+2003, Figure 13
        # Bounded mass-to-light ratio [2,10]
        def generate_bounded_normal(size, mean=6.0, stddev=3.0, bounds=[2.0,10.0]):
            data = np.random.normal(loc=mean, scale=stddev, size=size)

            if bounds is not None:
                lower_bound, upper_bound = bounds
                iout = np.where((data < lower_bound) | (data > upper_bound))[0]
                if len(iout) > 0:
                    data[iout] = generate_bounded_normal(len(iout), mean=mean, stddev=stddev, bounds=bounds)
                return data
        
        ml_ratio = (generate_bounded_normal(1))[0]
        return ml_ratio
    
    def set_bolometric_luminosity(self, lum=None):
        if lum is not None:
            return 10**lum
        elif self.mass is not None:
            return (self.mass / self.mass_to_light) * const.L_sun.to(u.erg/u.s).value
        else:
            # Implement López-Sanjuan+2018? 
            raise NotImplementedError("Galaxy-to-AGN properties not yet implemented. Please define AGN bolometric luminosity.")

    def relative_component_normalization(self, coefficients, method="random"):
        if method == "random":
            coeff = np.random.rand(len(coefficients))
            relative_coeff = coeff / np.sum(coeff)
        elif method == "ml_scaled":
            # Method assumes mass-to-light ratio bounded on [2,10]
            ml_ratio = self.mass_to_light
            coeff = np.zeros(len(coefficients))
            coeff[0] = (ml_ratio - 2.0) / 8.0
            coeff[1] = (1.0 - coeff[0]) * (((ml_ratio - 2.0)/8.0)* 0.25 + 0.75)
            coeff[2] = 1.0 - coeff[0] - coeff[1]
            relative_coeff = coeff / np.sum(coeff)  # Sanity check
        else:
            raise NotImplementedError(f"Alternative relative normalization method '{method}' not yet implemented. Use method='random'.")
        return {key: relative_coeff[i] for i, key in enumerate(coefficients)}

    def normalize_sed(self, rest_freq, dlum):
        coefficients = {component: None for component in sed.flux['galaxy'].keys()}
        coefficients = self.relative_component_normalization(coefficients, method='ml_scaled')
        flux = sum(sed.flux['galaxy'][key] * coefficients[key] for key in coefficients)
        integrated_flux = np.trapz(flux[::-1], x=rest_freq[::-1])
        integrated_lum = 4 * np.pi * (dlum**2) * integrated_flux
        sed_normalization = self.luminosity / integrated_lum
        sed_normalizations = {key: coefficients[key] * sed_normalization for key in coefficients}
        return sed_normalizations
    
    def construct_sed(self):
        flux = {key: self.norm[key] * sed.flux['galaxy'][key] for key in self.norm}
        return flux

    def attenuate_sed(self):
        flux = cp.deepcopy(self.sed_unabsorbed)
        if 'galaxy' in sed.ext:
            for key in sed.ext['galaxy']:
                flux[key] *= 10.0**(-0.4 * self.ebv * sed.ext['galaxy'][key])
        return flux


#########################################
#
# Source
#
#########################################
class Source:
    def __init__(self, redshift=0, agn_luminosity=0, ebv=0):
        self.z = redshift
        self.dluminosity = self.calculate_luminosity_distance()
        self.type = self.determine_type()
        self.wave = self.set_frame_data('wavelength')
        self.freq = self.set_frame_data('frequency')
        # Assign dynamically in the future
        self.agn = ActiveGalacticNuclei(self.freq['rest'],
                                        self.dluminosity,
                                        luminosity=agn_luminosity,
                                        ebv=ebv)
        self.galaxy = Galaxy(self.z,
                             self.freq['rest'],
                             self.dluminosity,
                             agn=self.agn)
        self.mass = self.agn.mass + self.galaxy.mass
        self.luminosity = self.agn.luminosity + self.galaxy.luminosity
        self.sed_unabsorbed = self.construct_sed()
        self.sed = self.attenuate_sed()
        self.photometry = self.initialize_photometry()
        self.calculate_photometry()
        self.detect = self.initialize_detections()
        self.trigger_detections()
        
    def calculate_luminosity_distance(self):
        return cosmology.luminosity_distance(self.z).to(u.cm).value
    
    def determine_type(self):
        # Placeholder method to determine if source is a star, galaxy, or galaxy+AGN
        return 'galaxy'
    
    def set_frame_data(self, frame_key):
        frame_data = {}
        if frame_key == 'wavelength':
            if not hasattr(sed, 'frame') or 'wavelength' not in sed.frame:
                raise ValueError("SED instance does not have 'frame' attribute with 'wavelength' key.")
            rest_wavelength = sed.frame['wavelength']['wave']
            frame_data['rest'] = rest_wavelength
            frame_data['obsv'] = rest_wavelength * (1 + self.z)
        elif frame_key == 'frequency':
            if not hasattr(self, 'wave'):
                raise ValueError("Source attribue 'wave' must be set before computing frequency data.")
            frame_data['rest'] = const.c.value / self.wave['rest']
            frame_data['obsv'] = const.c.value / self.wave['obsv']
        else:
            raise ValueError(f"Unknown frame_key: {frame_key}")
        return frame_data

    def construct_sed(self):
        flux = np.zeros(len(self.wave['rest']))
        if hasattr(self, 'agn'):
            flux += np.sum(np.stack([self.agn.sed_unabsorbed[key] for key in self.agn.sed_unabsorbed]), axis=0)
        if hasattr(self, 'galaxy'):
            flux += np.sum(np.stack([self.galaxy.sed_unabsorbed[key] for key in self.galaxy.sed_unabsorbed]), axis=0)
        if hasattr(self, 'star'):
            raise NotImplementedError(f"Stars not yet implemented.")
        if hasattr(self, 'neo'):
            raise NotImplementedError(f"NEOs not yet implemented.")
        return flux

    def attenuate_sed(self):
        flux = np.zeros(len(self.wave['rest']))
        if hasattr(self, 'agn'):
            flux += np.sum(np.stack([self.agn.sed[key] for key in self.agn.sed]), axis=0)
        if hasattr(self, 'galaxy'):
            flux += np.sum(np.stack([self.galaxy.sed[key] for key in self.galaxy.sed]), axis=0)
        if hasattr(self, 'star'):
            raise NotImplementedError(f"Stars not yet implemented.")
        if hasattr(self, 'neo'):
            raise NotImplementedError(f"NEOs not yet implemented.")
        return flux
    
    def initialize_photometry(self):
        return {
            'flux': {key: None for key in bp.get_filters()},
            'mag': {key: None for key in bp.get_filters()},
        }

    def calculate_photometry(self):
        bp.convolve_with_bandpass(self)
        bp.convert_flx_to_mag(self)
    
    def initialize_detections(self):
        return {key: False for key in bp.get_filters()}

    def trigger_detections(self):
        # Perform calculations and update self.photometry
        bp.compare_survey_depths(self)




#########################################
#
# Star
#
#########################################

#########################################
#
# Near-Earth Asteroid
#
#########################################


