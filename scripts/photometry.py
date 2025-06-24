# scripts/photometry.py
from scripts.common import const, np, pd
from resources.extinction_laws import extinction_laws
from resources.sed_templates import sed_templates
from resources.surveys import surveys
import pdb


#########################################
#
# Spectral Energy Distribution
#
#########################################
class SpectralEnergyDistribution:
    def __init__(self, sed_choice, ext_choice):
        self._meta = {"templates": self.validate_sed_selection(sed_choice),
                      "extinction": self.validate_ext_selection(ext_choice)
                      }
        self.domain = self.load_sed_templates('domain')
        # this is dreafully unclean. i hate it. make this better later
        self.domain['frequency'] = {}
        self.domain['frequency']['freq'] = const.c/self.domain['wavelength']['wave']
        #self.flux = self.load_sed_templates('components')
        self.components = self.load_sed_templates('components')
        self.ext = self.load_extinction_laws()

    def validate_sed_selection(self, sed_selection):
        template_sources = {}
        duplicates = {}

        #print(f"LOOP over SED sources...")
        for sed_source, templates in sed_selection.items():
            if not hasattr(sed_templates, sed_source.lower()):
                raise ValueError(f"SED data for '{sed_source}' not found in sed_templates module.")
            
            #print(f"  LOOP over COMPONENTS: {components}")
            for template in templates: 
                if template in template_sources:
                    if template not in duplicates:
                        duplicates[template] = [template_sources[template]]
                    duplicates[template].append(sed_source)
                else:
                    template_sources[template] = sed_source
            #print(f"    COMPONENT_SOURCES: {component_sources}")
            #print(f"    CURRENT DUPLICATES: {duplicates}")

        if duplicates:
            duplicate_messages = []
            for template, sources in duplicates.items():
                duplicate_messages.append(f"{template}: {', '.join(sources)}")
            duplicate_message = "\n".join(duplicate_messages)
            raise ValueError(f"Duplicate components detected:\n{duplicate_message}")

        return sed_selection    

    def validate_ext_selection(self, ext_selection):
        template_sources = {}
        duplicates = {}

        #print(f"LOOP over extinction laws...")
        for ext_source in ext_selection:
            #print(f"  Extinction law: {ext_name}")
            if not hasattr(extinction_laws, ext_source.lower()):
                raise ValueError(f"Extinction law '{ext_source}' not found in extinction_laws module.")

            ext_data = getattr(extinction_laws, ext_source)            
            #ext_flux = ext_data['flux']
            ext_components = ext_data['components']
            
            #print(f"  LOOP over TYPE...")
            for _, templates in ext_components.items():
                #print(f"    Checking TYPE: {flux_type}")
                if templates is not None:
                    #print(f"    For COMPONENTS: {components}")
                    for template in templates:
                        if template in template_sources:
                            if template not in duplicates:
                                duplicates[template] = [template_sources[template]]
                            duplicates[template].append(f"{ext_source}")
                            
                        else:
                            template_sources[template] = f"{ext_source}"
                #print(f"COMPONENT_SOURCES: {component_sources}")
                #print(f"CURRENT DUPLICATES: {duplicates}")
        if duplicates:
            duplicate_messages = []
            for template, sources in duplicates.items():
                duplicate_messages.append(f"{template}: {', '.join(sources)}")
            duplicate_message = "\n".join(duplicate_messages)
            raise ValueError(f"Duplicate components detected in extinction selection:\n{duplicate_message}")

        return ext_selection

    def load_sed_templates(self, category):
        # Initialize a dictionary to hold the component data
        component_data = {}

        # Retrieve the list of component types for the given key
        components = sed_templates.category_keys[category]
        #print(f"Loading components for key: {component_key}")
        #print(f"Component types to check: {component_types}")
        
        # Loop over the selected SED templates
        for sed_source, desired_templates in self._meta['templates'].items():
            #print(f"Processing SED template: {sed_name}")
            #print(f"Desired components: {desired_components}")
            # Get the sed data dictionary dynamically using getattr
            sed_data = getattr(sed_templates, sed_source.lower(), None)

            # Loop over the component types
            for component in components:
                # Fetch the specific component type
                templates = sed_data.get(category, {}).get(component, None)
                if templates is None:
                    #print(f"{component_type} not available in {sed_name} under {component_key}")
                    continue

                matching_templates = [template for template in desired_templates if template in templates]
                #print(f"Matching components for {component_type}: {matching_components}")

                if matching_templates:
                    file_paths = sed_data['file_paths']
                    columns = sed_data['columns']
                    wave_conversion = sed_data.get('wave_conversion', {})
                    flux_conversions = sed_data.get('flux_conversions', {})

                    for file_key, column_names in columns.items():
                        data_file = file_paths[file_key]
                        selected_templates = [template for template in column_names if template in matching_templates]
                        #print(f"Reading file: {data_file}")
                        #print(f"Selected columns: {selected_columns}")

                        # Read the data and apply conversions
                        #data = pd.read_csv(data_file, delim_whitespace=True, comment='#', header=None)
                        data = pd.read_csv(data_file, sep='\s+', comment='#', header=None, engine='python')
                        data.columns = column_names

                        for template in selected_templates:
                            if template in wave_conversion.get(file_key, {}):
                                #print(f"Applying wavelength conversion for {col}")
                                data[template] *= wave_conversion[file_key][template]

                            if template in flux_conversions.get(file_key, {}):
                                #print(f"Applying flux conversion for {col}")
                                data[template] *= flux_conversions[file_key][template]

                        # Initialize the nested dictionary if not already
                        if component not in component_data:
                            component_data[component] = {}

                        # Store the processed data directly under component_data
                        for template in selected_templates:
                            #print(f"Storing data for {col} under {component_type}")
                            component_data[component][template] = data[template]

        return component_data
    
    def load_extinction_laws(self):
        extinction_data = {}

        # Loop over the selected extinction laws
        #print(f"LOOP over extinction laws...")
        for ext_source in self._meta['extinction']:
            #print(f"  Extinction law: {ext_name}")
            ext_data = getattr(extinction_laws, ext_source.lower(), None)

            for component, templates in ext_data['components'].items():
                #print(f"    Component TYPE: {flux_type}")
                if component in self.components:
                    #print(f"      Matched TYPE in self.flux: {self.flux.keys()}")
                    #print(f"      Checking for COMPONENTS in TYPE: {flux_type}")
                    if templates is not None:
                        #print(f"      COMPONENTS exist: {components}")
                        matching_templates = [template for template in templates if template in self.components[component]]
                        #print(f"      Matched COMPONENTS: {matching_components}")
                        if matching_templates:
                            file_paths = ext_data['file_paths']
                            for file_key, column_names in ext_data['columns'].items():
                                #print(f"        FILE KEY: {file_key}")
                                #print(f"        Extinction law columns: {col_names}")
                                #print(f"        Reading DATA...")
                                data_file = file_paths[file_key]
                                #data = pd.read_csv(data_file, delim_whitespace=True, comment='#', header=None)
                                data = pd.read_csv(data_file, sep='\s', comment='#', header=None, engine='python')
                                data.columns = column_names

                                # Initialize the nested dictionary if not already
                                #print(f"        Adding TYPE: {flux_type}")
                                if component not in extinction_data:
                                    extinction_data[component] = {}

                                # Extract the 'kappa' column as a Series
                                #print(f"        LOOP over matched COMPONENTS...")
                                for template in matching_templates:
                                    #print(f"          Transferring KAPPA for COMPONENT: {component}")
                                    if 'kappa' in data.columns:
                                        extinction_data[component][template] = data['kappa']

        return extinction_data
    
    def meta(self):
        print("SpectralEnergyDistribution Metadata")
        for key, val in self._meta.items():
            print(f"  {key}: {val}")
            

#########################################
#
# Bandpass Filters 
#
#########################################
class Filter:
    def __init__(self, wave, freq, thru, zero_pt, wave_eff, freq_eff):
        self.wave = wave
        self.freq = freq
        self.thru = thru
        self.zero_pt = zero_pt
        self.wave_eff = wave_eff
        self.freq_eff = freq_eff

class Bandpass:
    def __init__(self, instrument):
        self.instrument = self.validate_instrument(instrument)
        self.filters = self.load_bandpass_filters()
        self.depths = self.set_survey_depths()
        self.flux_limits = self.set_survey_flux_limits()

    def validate_instrument(self, instrument_name):
        if not hasattr(surveys, instrument_name.lower()):
            raise ValueError(f"Instrument '{instrument_name}' not found in instruments module.")
        return instrument_name.lower()
    
    def load_bandpass_filters(self):
        instrument_data = getattr(surveys, self.instrument)

        filters = {}
        wave_conversion = instrument_data['wave_conversion']
        for filter in instrument_data['filters']:
            file_path = instrument_data['file_paths'].get(filter)
            if file_path is None:
                raise ValueError(f"File path not found for filter '{filter}'")
            
            data = np.loadtxt(file_path, unpack=True)
            wave = data[0] * wave_conversion  # Assumes first column is wavelength
            thru = data[1]  # Assumes second column is throughput
            zero_pt = instrument_data['zero_points'].get(filter)
            wave_eff = instrument_data['wave_eff'].get(filter) * wave_conversion
            freq = const.c.value / wave
            freq_eff = const.c.value / wave_eff
            filters[filter] = Filter(wave, freq, thru, zero_pt, wave_eff, freq_eff)
        return filters
    
    def set_survey_depths(self):
        survey_data = getattr(surveys, self.instrument)
        return survey_data['depths']
    
    def set_survey_flux_limits(self):
        flux_limit = {mode: {} for mode in self.get_depths()}
        for mode, depths in self.depths.items():
            flux_limit[mode] = self.convert_mag_to_flx(depths)
        return flux_limit

    def convolve_with_bandpass(self, source_wave, source_sed, filter=None):
        # Convolve SED with the specified bandpass filter or with all filters if filter is None
        convolved_flux = {key: None for key in self.get_filters()}
        filters_to_convolve = self.filters if filter is None else {filt: self.filters[filt] for filt in filter}

        for filter, bandpass in filters_to_convolve.items():
            # Check for overlapping wavelength
            if not (np.min(source_wave) <= np.min(bandpass.wave) or np.max(bandpass.wave) <= np.max(source_wave)):
                convolved_flux[filter] = 0#None
                continue
            
            # Interpolate the SED onto filter wavelengths
            interpolated_sed = np.interp(bandpass.wave, source_wave, source_sed, left=0.0, right=0.0)
            # Normalize throughput
            thru_norm = np.trapz(bandpass.thru[::-1], x=bandpass.freq[::-1])
            convolved_flux[filter] = np.trapz((interpolated_sed * bandpass.thru)[::-1], x=bandpass.freq[::-1]) / thru_norm
        return convolved_flux

    def calculate_kcorrection(self, obsv_flux, rest_flux, source_z, filter=None):
        # Calculate k-correction from flux
        kcorrections = {filt: None for filt in self.get_filters()}
        
        for filter in kcorrections:
            kcorrections[filter] = -2.5 * np.log10((1+source_z)*obsv_flux[filter]/rest_flux[filter])
        return kcorrections

    def convert_flx_to_mag(self, source_flux, filter=None):
        # Convert flux density to magnitude with the specified bandpass filter or with all filters if filter is None
        converted_magnitudes = {key: None for key in self.get_filters()}
        filters_to_convert = self.filters if filter is None else {filt: self.filters[filt] for filt in filter}
        
        for filter, bandpass in filters_to_convert.items():
            # Convert zero point from Jy to erg/s/cm2/Hz [1 Jy == 1e-23 erg/s/cm2/Hz]
            converted_magnitudes[filter] = -2.5 * np.log10(source_flux[filter] / (bandpass.zero_pt * 1e-23))
        return converted_magnitudes

    def convert_mag_to_flx(self, source_mag, filter=None):
        # Convert magnitude to flux density with the specified bandpass filter or with all filters if filter is None
        converted_flux = {key: None for key in self.get_filters() if filter is None}
        filters_to_convert = self.filters if filter is None else {filt: self.filters[filt] for filt in filter}

        for filter, bandpass in filters_to_convert.items():
            # Convert zero point from erg/s/cm2/Hz to Jy [1 Jy == 1e-23 erg/s/cm2/Hz]
            converted_flux[filter] = (bandpass.zero_pt * 1e-23) * 10**(-0.4 * source_mag[filter])
        return converted_flux
    
    def apparent_to_absolute_magnitude(self, apparent_mag, luminosity_distance, kcorrection=None, filter=None):
        # Compute absolute magnitudes with specified bandpass or with all filters if filter is None
        abs_magnitudes = {key: None for key in self.get_filters()}
        if kcorrection is None: kcorrection = {key: 0 for key in self.get_filters()}
        for filter in abs_magnitudes:
            # Values of luminosity distance in units of cm
            abs_magnitudes[filter] = apparent_mag[filter] - 5*np.log10(luminosity_distance) + 97.447 - kcorrection[filter]
        return abs_magnitudes

    def compare_survey_depths(self, source_mag, filter=None):
        detected_bands = {mode: {filt: False for filt in self.get_filters()} for mode in self.get_depths()}
        # Check magnitude brighter than flux limit with the specified bandpass filter or with all filters if filter is None
        if filter is None:
            filters_to_compare = self.depths
        else:
            filters_to_compare = {}
            for mode in self.depths:
                filters_to_compare[mode] = {filt: self.depths[mode] for filt in filter}
        for mode, filters in filters_to_compare.items():
            for filt, mag_depth in filters.items():
                detected_bands[mode][filt] = source_mag[filt] <= mag_depth
        return detected_bands

    def get_filters(self):
        return list(self.filters.keys())

    def get_wavelengths(self):
        return np.array([bandpass.wave_eff for bandpass in self.filters.values()])
    
    def get_frequencies(self):
        return np.array([bandpass.freq_eff for bandpass in self.filters.values()])

    def get_zero_points(self):
        return np.array([bandpass.zero_pt for bandpass in self.filters.values()])
    
    def get_depths(self):
        return list(self.depths.keys())
    

#########################################
#
# Survey
#
#########################################
#class Survey:
#    def __init__(self, survey, selection=None):
#        self.survey = self.validate_survey(survey)
#        self.bp = Bandpass(instrument=survey)
#        self.survey_depths = self.set_survey_depths()
#        self.selction = self.set_selection_criteria(selection)
#
#    def validate_survey(self, survey_name):
#        if not hasattr(surveys, survey_name.lower()):
#            raise ValueError(f"Instrument '{survey_name}' not found in instruments module.")
#        return survey_name
#    
#    def set_survey_depths(self):
#        # I DON'T HAVE TIME TO IMPLEMENT THIS LIKE I WANT
#
#        return None
#
#    def set_selection_criteria(self, selection):
#        return None












