# scripts/simulated_source.py
from scripts.common import const, cp, np, pd, u
from scripts.global_variables import cosmo, sed, select, bp
import pdb


#########################################
#
# Active Galactic Nuclei
#
#########################################
class ActiveGalacticNuclei:
    def __init__(self, source, luminosity=None, ebv=0, galaxy=None):
        self.luminosity = self.set_bolometric_luminosity(luminosity)
        self.ebv = ebv
        self.eddington_ratio = self.assign_eddington_ratio()
        self.mass = self.calculate_bh_mass()
        self.eddington_limit = self.calculate_eddington_limit()
        self.norm = self.template_normalizations(source.z, source.dluminosity, source.freq)
        self.sed_unabsorbed = self.construct_sed(source.z)
        self.sed = self.attenuate_sed()
        self.flux = self.compute_flux(source.wave)
        self.rest_flux = self.compute_flux(source.wave, z=source.z)
        self.kcorr = self.compute_kcorr(source.z)
        self.mag = self.compute_magnitude()
        self.abs_mag = self.compute_absolute_magnitude(source.dluminosity)
        
    def set_bolometric_luminosity(self, lum):
        if lum is not None:
            return 10**lum
        else:
            raise NotImplementedError("Galaxy-to-AGN properties not yet implemented. Please define AGN bolometric luminosity.")

    def assign_eddington_ratio(self, method="Weigel+2017"):
        if method == "Weigel+2017":
            # https://ui.adsabs.harvard.edu/abs/2017ApJ...845..134W/abstract
            # Radiatively efficient (x == blue/green galaxies) vs. innefficient (r == red galaxies)
            rad_x = {"log_lambda_star":-1.84, "delta1":0.47, "delta2":2.53, "log_xi_star":-1.65, "epsilon":2.22}
            rad_r = {"log_lambda_star":-2.81, "delta1":0.41, "delta2":1.22, "log_xi_star":-2.13, "epsilon":0.82}
            # Load parameters. DETERMINE WHICH CHOICE LATER
            param = rad_x
            # Build ERDF 
            edd_ratio = 10**np.arange(-8.0, 1.0, 0.01)
            erdf = 10**param['log_xi_star'] * ((edd_ratio/10**param['log_lambda_star'])**param['delta1'] + (edd_ratio/10**param['log_lambda_star'])**param['delta2'])**(-1)
            # Sample from ERDF
            dedd_ratio = np.diff(edd_ratio)
            dedd_ratio = np.append(dedd_ratio,dedd_ratio[-1])
            erdf_pdf = erdf / np.sum(erdf * dedd_ratio)
            erdf_cdf = np.cumsum(erdf_pdf * dedd_ratio)
            random_draw = np.random.rand()
            return np.interp(random_draw, erdf_cdf, edd_ratio)
        #lif method == "Wang+2017":
        #   # https://ui.adsabs.harvard.edu/abs/2017A%26A...601A..63W/abstract
        #   lambda_edd = 10**np.arange(-3.0, 0.5, 0.01)
        #   # Draw Eddington Ratio from PDF, Equation 9, parameters from Table 3
        #   def calculate_log_p_edd(z):
        #       big_a = 10**np.random.normal(loc=-2.5, scale=0.2)
        #       alpha = np.random.normal(loc=0.38, scale=0.05)
        #       gamma = np.random.normal(loc=1.8, scale=0.3)
        #       return big_a * (lambda_edd)**(-1*alpha) * (1 + z)**gamma
        #   log_p_edd = calculate_log_p_edd(z)
        #   edd_shift = log_p_edd - np.min(log_p_edd)
        #   edd_pdf = edd_shift / np.sum(edd_shift)
        #   edd_cdf = np.cumsum(edd_pdf)
        #   random_draw = np.random.rand()
        #   samp = np.interp(random_draw, edd_cdf, lambda_edd)
        #   return samp
        elif method == "Suh+2015":
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
            raise NotImplementedError(f"Alternative Eddington ratio method '{method}' not yet implemented. Use method='Weigel+2017'.")
    
    def calculate_bh_mass(self):
        mass = self.luminosity / (1.26e38 * self.eddington_ratio)
        return mass

    def calculate_eddington_limit(self):
        edd_limit = 1.26e38 * self.mass
        return edd_limit

    def relative_template_normalizations(self, coefficients):
        if len(coefficients) == 1:
            key = next(iter(coefficients))
            coefficients[key] = 1.0
        elif len(coefficients) > 1:
            raise NotImplementedError("Multiple component AGN SED not yet implemented.")
        else:
            raise ValueError("Coefficients dictionary is empty.")
        return coefficients

    def template_normalizations(self, z, dlum, freq):
        coefficients = {template: None for template in sed.components['agn'].keys()}
        coefficients = self.relative_template_normalizations(coefficients)
        flux = sum(sed.components['agn'][template] * coefficients[template] for template in coefficients)
        integrated_flux = np.trapz(flux[::-1]/(1+z), x=freq[::-1])
        integrated_lum = 4*np.pi*(dlum**2)*integrated_flux
        template_norm = self.luminosity / integrated_lum
        template_norms = {template: coefficients[template] * template_norm for template in coefficients}
        return template_norms
        
    def construct_sed(self, z):
        component_sed = {}
        for template in self.norm:
            component_sed[template] = self.norm[template]*sed.components['agn'][template]/(1+z)
        return component_sed

    def attenuate_sed(self):
        component_sed = cp.deepcopy(self.sed_unabsorbed)
        if 'agn' in sed.ext:
            for template in sed.ext['agn']:
                component_sed[template] *= 10.0**(-0.4 * self.ebv * sed.ext['agn'][template])
        return component_sed
    
    def compute_flux(self, wave, z=0):
        component_flux = {}
        for template, template_sed in self.sed.items():
            component_flux[template] = bp.convolve_with_bandpass(wave/(1+z), template_sed*(1+z))
        return component_flux
    
    def compute_kcorr(self, z):
        component_kcorr = {}
        for template in self.flux:
            component_kcorr[template] = bp.calculate_kcorrection(self.flux[template],self.rest_flux[template],z)
        return component_kcorr
    
    def compute_magnitude(self):
        component_mag = {}
        for template in self.flux:
            component_mag[template] = bp.convert_flx_to_mag(self.flux[template])
        return component_mag
    
    def compute_absolute_magnitude(self, dlum):
        component_mag = {}
        for template in self.mag:
            component_mag[template] = bp.apparent_to_absolute_magnitude(self.mag[template], dlum, self.kcorr[template])
        return component_mag
    

#########################################
#
# Galaxy
#
#########################################
class Galaxy:
    def __init__(self, source, luminosity=None, mass=None, ebv=0, agn=None):
        self.host_galaxy = self.is_host_galaxy(agn)
        self.ebv = ebv
        self.mass = self.calculate_stellar_mass(source.z, mass, agn)
        self.mass_to_light = self.calculate_ml_ratio()
        self.luminosity = self.set_bolometric_luminosity(luminosity)
        self.norm = self.template_normalizations(source.z, source.dluminosity, source.freq)
        self.sed_unabsorbed = self.construct_sed(source.z)
        self.sed = self.attenuate_sed()
        self.flux = self.compute_flux(source.wave)
        self.rest_flux = self.compute_flux(source.wave, z_rest=source.z)
        self.kcorr = self.compute_kcorr(source.z)
        self.mag = self.compute_magnitude()
        self.abs_mag = self.compute_absolute_magnitude(source.dluminosity)

    def is_host_galaxy(self, agn_instance=None):
        if agn_instance is not None:
            return True
        else:
            return False

    def scale_mass_ratio(self, z):
        # Correct for AGN mass overdensity 
        # Pacucci & Loeb 2024 Equation 2
        # https://ui.adsabs.harvard.edu/abs/2024ApJ...964..154P
        # ...which follows Barkana & Loeb (2001)
        # ...with definitions from Bryan & Norman 1998)
        # Ωzm === Critial matter density at redshift. Implemented in astropy.cosmo.FlatLambdaCDM
        # ∆c  === Ratio mean density / critical density at redshift
        # d   === Critial matter density at collapse redshift. Defined as "x" in Bryan & Norman (1998)
        # E   === BH mass overdensity redshift evolution 
        def calculate_delta_c(omega_m):
            d = omega_m - 1
            return 18*np.pi**2 + 82*d - 39*d**2
        
        def calculate_xi(omega_m, delta_c):
            return (cosmo.Om0/omega_m) * (delta_c/(18*np.pi**2))
        
        xi = []
        for z_i in [0, z]:
            omega_m = cosmo.Om(z_i)
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

    def relative_template_normalizations(self, coefficients, method="random"):
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

    def template_normalizations(self, z, dlum, freq):
        coefficients = {template: None for template in sed.components['galaxy'].keys()}
        coefficients = self.relative_template_normalizations(coefficients, method='ml_scaled')
        flux = sum(sed.components['galaxy'][template] * coefficients[template] for template in coefficients)
        integrated_flux = np.trapz(flux[::-1]/(1+z), x=freq[::-1])
        integrated_lum = 4*np.pi*(dlum**2)*integrated_flux
        template_norm = self.luminosity / integrated_lum
        template_norms = {key: coefficients[key] * template_norm for key in coefficients}
        return template_norms
    
    def construct_sed(self, z):
        component_sed = {template: self.norm[template]*sed.components['galaxy'][template]/(1+z) for template in self.norm}
        return component_sed

    def attenuate_sed(self):
        component_sed = cp.deepcopy(self.sed_unabsorbed)
        if 'galaxy' in sed.ext:
            for template in sed.ext['galaxy']:
                component_sed[template] *= 10.0**(-0.4 * self.ebv * sed.ext['galaxy'][template])
        return component_sed

    def compute_flux(self, wave, z_rest=0):
        component_flux = {}
        for template, template_sed in self.sed.items():
            component_flux[template] = bp.convolve_with_bandpass(wave/(1+z_rest), template_sed*(1+z_rest))
        return component_flux

    def compute_kcorr(self, z):
        component_kcorr = {}
        for template in self.flux:
            component_kcorr[template] = bp.calculate_kcorrection(self.flux[template],self.rest_flux[template],z)
        return component_kcorr
    
    def compute_magnitude(self):
        component_mag = {}
        for template in self.flux:
            component_mag[template] = bp.convert_flx_to_mag(self.flux[template])
        return component_mag
    
    def compute_absolute_magnitude(self, dlum):
        component_mag = {}
        for template in self.mag:
            component_mag[template] = bp.apparent_to_absolute_magnitude(self.mag[template], dlum, self.kcorr[template])
        return component_mag
    
    

#########################################
#
# Source
#
#########################################
class Source:
    def __init__(self, redshift=0, agn_luminosity=0, ebv=0, components=None):
        self.z = redshift
        self.dluminosity = self.calculate_luminosity_distance()
        self.type = self.determine_type()
        self.wave = self.set_domain_data('wavelength')
        self.freq = self.set_domain_data('frequency')
        # Dictionary of object types
        self.components = {}
        # Assign dynamically in the future
        if components is None:
            components = ['galaxy','agn']
        if 'agn' in components:
            self.components['agn'] = ActiveGalacticNuclei(self, luminosity=agn_luminosity, ebv=ebv)
        if 'galaxy' in components:
            self.components['galaxy'] = Galaxy(self, agn=self.components['agn'])
        if 'star' in components:
            self.components['star'] = Star(self)
        if 'neo' in components:
            self.components['neo'] = NearEarthObject(self)
        #self.agn = self.components['agn']
        #self.galaxy = self.components['galaxy']
        self.mass = self.combine_component_attributes('mass')
        self.luminosity = self.combine_component_attributes('luminosity')
        self.sed_unabsorbed = self.combine_component_attributes('sed_unabsorbed')
        self.sed = self.combine_component_attributes('sed')
        self.flux = self.combine_component_flux('flux')
        self.rest_flux = self.combine_component_flux('rest_flux')
        self.kcorr = self.compute_kcorrection()
        self.mag = self.compute_magnitude()
        self.abs_mag = self.compute_absolute_magnitude()
        self.detect = self.trigger_detections()
        self.select = self.trigger_selections()

    def calculate_luminosity_distance(self):
        return cosmo.luminosity_distance(self.z).to(u.cm).value
    
    def determine_type(self):
        # Placeholder method to determine if source is a star, galaxy, or galaxy+AGN
        return 'galaxy'
    
    def set_domain_data(self, domain_key):
        domain_data = {}
        if domain_key == 'wavelength':
            if not hasattr(sed, 'domain') or 'wavelength' not in sed.domain:
                raise ValueError("SED instance does not have 'domain' attribute with 'wavelength' key.")
            #domain_data['rest'] = rest_wavelength
            #domain_data['obsv'] = rest_wavelength * (1+self.z)
            domain_data = sed.domain['wavelength']['wave'] * (1+self.z)
        elif domain_key == 'frequency':
            if not hasattr(self, 'wave'):
                raise ValueError("Source attribue 'wave' must be set before computing frequency data.")
            #domain_data['rest'] = const.c.value / self.wave['rest']
            #domain_data['obsv'] = const.c.value / self.wave['obsv']
            domain_data = const.c.value / self.wave
        else:
            raise ValueError(f"Unknown domain_key: {domain_key}")
        return domain_data

    def combine_component_attributes(self, attr):
        source_attr = []
        for component in self.components.values():
            if hasattr(component, attr):
                raw_attr = getattr(component, attr)
                component_attr = np.sum(np.stack([raw_attr[template] for template in raw_attr]),axis=0) if hasattr(raw_attr,'keys') else raw_attr
                source_attr.append(component_attr)
        return np.sum(source_attr, axis=0)

    def combine_component_flux(self, attr):
        source_flux = {key: 0 for key in bp.get_filters()}
        for component in self.components.values():
            if hasattr(component, attr):
                component_flux = getattr(component, attr)
                for template_flux in component_flux.values():
                    for filter in source_flux:
                        source_flux[filter] += template_flux[filter]
        return source_flux

    #def construct_sed(self):
    #    source_sed = np.zeros(len(self.wave))
    #    for component in ['galaxy','agn']:
    #        if hasattr(self, component):
    #            source_component = getattr(self, component)
    #            source_sed += np.sum(np.stack([source_component.sed_unabsorbed[template] for template in source_component.sed_unabsorbed]), axis=0)
    #    return source_sed

    #def attenuate_sed(self):
    #    source_sed = np.zeros(len(self.wave))
    #    for component in ['galaxy','agn']:
    #        if hasattr(self, component):
    #            source_component = getattr(self, component)
    #            source_sed += np.sum(np.stack([source_component.sed[template] for template in source_component.sed]), axis=0)
    #    if hasattr(self, 'star'):
    #        raise NotImplementedError(f"Stars not yet implemented.")
    #    if hasattr(self, 'neo'):
    #        raise NotImplementedError(f"NEOs not yet implemented.")
    #    return source_sed
        
    #def compute_flux_old(self, rest=False):
    #    source_flux = {key: 0 for key in bp.get_filters()}
    #    domain = 'rest_flux' if rest else 'flux'
    #    for component in ['galaxy','agn']:
    #        if hasattr(self, component): 
    #            source_component = getattr(self, component)
    #            component_flux = getattr(source_component, domain)
    #            for template_flux in component_flux.values():
    #                for filter in source_flux:
    #                    source_flux[filter] += template_flux[filter]
    #    if hasattr(self, 'star'):
    #        raise NotImplementedError(f"Stars not yet implemented.")
    #    if hasattr(self, 'neo'):
    #        raise NotImplementedError(f"NEOs not yet implemented.")
    #    return source_flux
    
    def compute_kcorrection(self):
        return bp.calculate_kcorrection(self.flux,self.rest_flux,self.z)

    def compute_magnitude(self):
        return bp.convert_flx_to_mag(self.flux)
    
    def compute_absolute_magnitude(self):
        return bp.apparent_to_absolute_magnitude(self.mag, self.dluminosity, self.kcorr)
    
    def trigger_detections(self):
        # Perform calculations and update self.photometry
        return bp.compare_survey_depths(self.mag)

    def trigger_selections(self):
        # Perform calculations and update self.photometry
        return select.is_agn(self.mag)
    
    def get_wave(self):
        return self.wave
    
    
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
np.array([-0.15191675,5.29331781,2.87214917,1.25113931])



#########################################
#
# Methods for Analysis
#
#########################################
def process_source(source):
    row = {
        'z': source.z,
        'dluminosity': source.dluminosity,
        'mass': source.mass,
        'luminosity': source.luminosity,
        # Unpack source fluxes
        **{f"flux_{band}": value for band, value in source.flux.items()},
        # Unpack source magnitudes
        **{f"mag_{band}": value for band, value in source.mag.items()},
        # Unpack source absolute magnitudes
        **{f"abs_mag_{band}": value for band, value in source.abs_mag.items()},
        # Unpack detections by filter
        **{f"detected_{strategy}_{band}": value
           for strategy, band_dict in source.detect.items()
           for band, value in band_dict.items()},
        # Source detected in any filter
        **{f"detected_{strategy}": any(band_dict.values())
           for strategy, band_dict in source.detect.items()},
        'selected': source.select,
        'components': list(source.components.keys()),
    }
    # Galaxy component
    if 'galaxy' in source.components:
        galaxy = source.components['galaxy']
        row.update({
            'galaxy_mass': galaxy.mass,
            'galaxy_luminosity': galaxy.luminosity,
            'galaxy_mass_to_light': galaxy.mass_to_light,
            'galaxy_ebv': galaxy.ebv,
            **{f"galaxy_{comp}_norm": val
               for comp, val in galaxy.norm.items()},
            **{f"galaxy_{comp}_flux_{band}": val
               for comp, band_dict in galaxy.flux.items()
               for band, val in band_dict.items()},
            **{f"galaxy_{comp}_mag_{band}": val
               for comp, band_dict in galaxy.mag.items()
               for band, val in band_dict.items()},
            **{f"galaxy_{comp}_abs_mag_{band}": val
               for comp, band_dict in galaxy.abs_mag.items()
               for band, val in band_dict.items()},
        })
    # AGN component
    if 'agn' in source.components:
        agn = source.components['agn']
        row.update({
            'agn_mass': agn.mass,
            'agn_luminosity': agn.luminosity,
            'agn_eddington_ratio': agn.eddington_ratio,
            'agn_eddington_limit': agn.eddington_limit,
            'agn_ebv': agn.ebv,
            **{f"agn_{comp}_norm": val for comp, val in agn.norm.items()},
            **{f"agn_{comp}_flux_{band}": val
               for comp, band_dict in agn.flux.items()
               for band, val in band_dict.items()},
            **{f"agn_{comp}_mag_{band}": val
               for comp, band_dict in agn.mag.items()
               for band, val in band_dict.items()},
            **{f"agn_{comp}_abs_mag_{band}": val
               for comp, band_dict in agn.abs_mag.items()
               for band, val in band_dict.items()},
        })

    return row





#def process_source(source):
#    # Define and extract data from the source instance
#    return {
#        'z': source.z,
#        'dluminosity': source.dluminosity,
#        'mass': source.mass,
#        'luminosity': source.luminosity,
#        # Unpack source fluxes
#        **{f"{band}_flux": value for band, value in source.flux.items()},
#        # Unpack source magnitudes
#        **{f"{band}_mag": value for band, value in source.mag.items()},
#        # Unpack source absolute magnitudes
#        **{f"{band}_abs_mag": value for band, value in source.abs_mag.items()},
#        # Unpack detections by filter
#        **{f"{band}_detected_{strategy}": value
#           # Unpack survey strategy (single, coadd... or WFD, DDF... etc.)
#           for strategy, band_dict in source.detect.items()
#           # Unpack each band (e.g. 'u', etc.) and it's depth value
#           for band, value in band_dict.items()},
#        # Source detected in any filter
#        **{f"detected_{strategy}": any(band_dict.values())
#           # Unpack band dictionary for any() call
#           for strategy, band_dict in source.detect.items()},
#        # Source selected with AGN selection criteria
#        'selected': source.select,
#        # Extract data from the Galaxy instance
#        'galaxy_mass': source.galaxy.mass,
#        'galaxy_luminosity': source.galaxy.luminosity,
#        'galaxy_mass_to_light': source.galaxy.mass_to_light,
#        'galaxy_ebv': source.galaxy.ebv,
#        # Unpack individual galaxy normalizations
#        **{f"norm_{component}": value
#                for component, value in source.galaxy.norm.items()},
#        # Unpack individual galaxy fluxes
#        **{f"{band}_flux_{component}": value 
#                # Unpack each component (e.g., 'ell') and band dictionaries (e.g. {'u':value,'g':value,...})
#                for component, band_dict in source.galaxy.flux.items() 
#                # Unpack each band (e.g. 'u', etc.) and it's single value
#                for band, value in band_dict.items()},
#        # Unpack individual galaxy magnitudes
#        **{f"{band}_mag_{component}": value
#                for component, band_dict in source.galaxy.mag.items()
#                for band, value in band_dict.items()},
#        # Unpack individual galaxy absolute magnitudes
#        **{f"{band}_abs_mag_{component}": value
#                for component, band_dict in source.galaxy.abs_mag.items()
#                for band, value in band_dict.items()},
#        # Extract data from the AGN instance
#        'agn_mass': source.agn.mass,
#        'agn_luminosity': source.agn.luminosity,
#        'agn_eddington_ratio': source.agn.eddington_ratio,
#        'agn_eddington_limit': source.agn.eddington_limit,
#        'agn_ebv': source.agn.ebv,
#        # Unpack individual AGN normalizations
#        **{f"norm_{component}": value 
#            for component, value in source.agn.norm.items()},
#        # Unpack individual AGN fluxes
#        **{f"{band}_flux_{component}": value 
#                for component, band_dict in source.agn.flux.items() 
#                for band, value in band_dict.items()},
#        # Unpack individual AGN magnitudes
#        **{f"{band}_mag_{component}": value
#                for component, band_dict in source.agn.mag.items()
#                for band, value in band_dict.items()},
#        # Unpack individual AGN magnitudes
#        **{f"{band}_abs_mag_{component}": value
#                for component, band_dict in source.agn.abs_mag.items()
#                for band, value in band_dict.items()},
#        }
