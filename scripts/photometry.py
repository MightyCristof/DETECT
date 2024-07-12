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
        self.templates = self.validate_sed_selection(sed_choice)
        self.extinction = self.validate_ext_selection(ext_choice)
        self.frame = self.load_template_components('frame')
        self.flux = self.load_template_components('flux')
        self.ext = self.load_extinction_laws()

    def validate_sed_selection(self, sed_selection):
        component_sources = {}
        duplicates = {}

        #print(f"LOOP over SED sources...")
        for sed_name, components in sed_selection.items():
            if not hasattr(sed_templates, sed_name.lower()):
                raise ValueError(f"SED data for '{sed_name}' not found in sed_templates module.")
            
            # IMPORTANT: In Python 3.7+, dictionaries track the order in which keys are assigned.
            #            This feature is used in other areas of the code (i.e., order of galaxy components).
            #print(f"  LOOP over COMPONENTS: {components}")
            for component in components:
                if component in component_sources:
                    if component not in duplicates:
                        duplicates[component] = [component_sources[component]]
                    duplicates[component].append(sed_name)
                else:
                    component_sources[component] = sed_name
            #print(f"    COMPONENT_SOURCES: {component_sources}")
            #print(f"    CURRENT DUPLICATES: {duplicates}")

        if duplicates:
            duplicate_messages = []
            for component, sources in duplicates.items():
                duplicate_messages.append(f"{component}: {', '.join(sources)}")
            duplicate_message = "\n".join(duplicate_messages)
            raise ValueError(f"Duplicate components detected:\n{duplicate_message}")

        return sed_selection    

    def validate_ext_selection(self, ext_selection):
        component_sources = {}
        duplicates = {}

        #print(f"LOOP over extinction laws...")
        for ext_name in ext_selection:
            #print(f"  Extinction law: {ext_name}")
            if not hasattr(extinction_laws, ext_name.lower()):
                raise ValueError(f"Extinction law '{ext_name}' not found in extinction_laws module.")

            ext_data = getattr(extinction_laws, ext_name)            
            ext_flux = ext_data['flux']
            
            #print(f"  LOOP over TYPE...")
            for _, components in ext_flux.items():
                #print(f"    Checking TYPE: {flux_type}")
                if components is not None:
                    #print(f"    For COMPONENTS: {components}")
                    for component in components:
                        if component in component_sources:
                            if component not in duplicates:
                                duplicates[component] = [component_sources[component]]
                            duplicates[component].append(f"{ext_name}")
                            
                        else:
                            component_sources[component] = f"{ext_name}"
                #print(f"COMPONENT_SOURCES: {component_sources}")
                #print(f"CURRENT DUPLICATES: {duplicates}")
        if duplicates:
            duplicate_messages = []
            for component, sources in duplicates.items():
                duplicate_messages.append(f"{component}: {', '.join(sources)}")
            duplicate_message = "\n".join(duplicate_messages)
            raise ValueError(f"Duplicate components detected in extinction selection:\n{duplicate_message}")

        return ext_selection

    def load_template_components(self, component_key):
        # Initialize a dictionary to hold the component data
        component_data = {}

        # Retrieve the list of component types for the given key
        component_types = sed_templates.component_keys[component_key]
        #print(f"Loading components for key: {component_key}")
        #print(f"Component types to check: {component_types}")

        # Loop over the selected SED templates
        for sed_name, desired_components in self.templates.items():
            #print(f"Processing SED template: {sed_name}")
            #print(f"Desired components: {desired_components}")
            # Get the sed data dictionary dynamically using getattr
            sed_data = getattr(sed_templates, sed_name.lower(), None)

            # Loop over the component types
            for component_type in component_types:
                # Fetch the specific component type
                available_components = sed_data.get(component_key, {}).get(component_type, None)
                if available_components is None:
                    #print(f"{component_type} not available in {sed_name} under {component_key}")
                    continue

                matching_components = [comp for comp in desired_components if comp in available_components]
                #print(f"Matching components for {component_type}: {matching_components}")

                if matching_components:
                    file_paths = sed_data['file_paths']
                    columns = sed_data['columns']
                    wave_conversion = sed_data.get('wave_conversion', {})
                    flux_conversions = sed_data.get('flux_conversions', {})

                    for file_key, col_names in columns.items():
                        data_file = file_paths[file_key]
                        selected_columns = [col for col in col_names if col in matching_components]
                        #print(f"Reading file: {data_file}")
                        #print(f"Selected columns: {selected_columns}")

                        # Read the data and apply conversions
                        data = pd.read_csv(data_file, delim_whitespace=True, comment='#', header=None)
                        data.columns = col_names

                        for col in selected_columns:
                            if col in wave_conversion.get(file_key, {}):
                                #print(f"Applying wavelength conversion for {col}")
                                data[col] *= wave_conversion[file_key][col]

                            if col in flux_conversions.get(file_key, {}):
                                #print(f"Applying flux conversion for {col}")
                                data[col] *= flux_conversions[file_key][col]

                        # Initialize the nested dictionary if not already
                        if component_type not in component_data:
                            component_data[component_type] = {}

                        # Store the processed data directly under component_data
                        for col in selected_columns:
                            #print(f"Storing data for {col} under {component_type}")
                            component_data[component_type][col] = data[col]

        return component_data
    
    def load_extinction_laws(self):
        extinction_data = {}

        # Loop over the selected extinction laws
        #print(f"LOOP over extinction laws...")
        for ext_name in self.extinction:
            #print(f"  Extinction law: {ext_name}")
            ext_data = getattr(extinction_laws, ext_name.lower(), None)

            ext_types = ext_data['flux']
            #print(f"  LOOP over component TYPES...")
            for flux_type, components in ext_types.items():
                #print(f"    Component TYPE: {flux_type}")
                if flux_type in self.flux:
                    #print(f"      Matched TYPE in self.flux: {self.flux.keys()}")
                    #print(f"      Checking for COMPONENTS in TYPE: {flux_type}")
                    if components is not None:
                        #print(f"      COMPONENTS exist: {components}")
                        matching_components = [comp for comp in components if comp in self.flux[flux_type]]
                        #print(f"      Matched COMPONENTS: {matching_components}")
                        if matching_components:
                            file_paths = ext_data['file_paths']
                            columns = ext_data['columns']
                            #print(f"      File paths: {file_paths}")
                            #print(f"      Columns: {columns}")
                            #print(f"      LOOP over FILE KEYS....")
                            for file_key, col_names in columns.items():
                                #print(f"        FILE KEY: {file_key}")
                                #print(f"        Extinction law columns: {col_names}")
                                #print(f"        Reading DATA...")
                                data_file = file_paths[file_key]
                                data = pd.read_csv(data_file, delim_whitespace=True, comment='#', header=None)
                                data.columns = col_names

                                # Initialize the nested dictionary if not already
                                #print(f"        Adding TYPE: {flux_type}")
                                if flux_type not in extinction_data:
                                    extinction_data[flux_type] = {}

                                # Extract the 'kappa' column as a Series
                                #print(f"        LOOP over matched COMPONENTS...")
                                for component in matching_components:
                                    #print(f"          Transferring KAPPA for COMPONENT: {component}")
                                    if 'kappa' in data.columns:
                                        extinction_data[flux_type][component] = data['kappa']

        return extinction_data
    

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

    def convolve_with_bandpass(self, object, filter=None):
        # Convolve SED with the specified bandpass filter or with all filters if filter is None
        convolved_flux = object.photometry['flux']
        filters_to_convolve = self.filters if filter is None else {filter: self.filters[filter]}

        for filter, bandpass in filters_to_convolve.items():
            # Check for overlapping wavelength
            if not (np.min(object.wave['obsv']) <= np.min(bandpass.wave) or np.max(bandpass.wave) <= np.max(object.wave['obsv'])):
                convolved_flux[filter] = None
                continue
            
            # Interpolate the SED onto filter wavelengths
            interpolated_sed = np.interp(bandpass.wave, object.wave['obsv'], object.sed, left=0.0, right=0.0)
            # Normalize throughput
            thru_norm = np.trapz(bandpass.thru[::-1], x=bandpass.freq[::-1])
            convolved_flux[filter] = np.trapz((interpolated_sed * bandpass.thru)[::-1], x=bandpass.freq[::-1]) / thru_norm
            
    def convert_flx_to_mag(self, object, filter=None):
        # Convert flux to magnitude with the specified bandpass filter or with all filters if filter is None
        converted_magnitudes = object.photometry['mag']
        filters_to_convert = self.filters if filter is None else {filter: self.filters[filter]}
        
        for filter, bandpass in filters_to_convert.items():
            converted_magnitudes[filter] = -2.5 * np.log10(object.photometry['flux'][filter] / (bandpass.zero_pt * 1e-23)) # Convert zero point from Jy to erg/s/cm2/Hz

    def compare_survey_depths(self, object, mode='single', filter=None):
        # Check magnitude brighter than flux limit with the specified bandpass filter or with all filters if filter is None
        filters_to_compare = self.depths[mode] if filter is None else {filter: self.depths[mode][filter]}
        detected_bands = object.detect

        for filter, mag_depth in filters_to_compare.items():
            detected_bands[filter] = object.photometry['mag'][filter] <= mag_depth

    def get_filters(self):
        return list(self.filters.keys())

    def get_wavelengths(self):
        return np.array([bandpass.wave_eff for bandpass in self.filters.values()])
    
    def get_frequencies(self):
        return np.array([bandpass.freq_eff for bandpass in self.filters.values()])

    def get_zero_points(self):
        return np.array([bandpass.zero_pt for bandpass in self.filters.values()])
    

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












