# resources/parameters.py

# Parameters for cosmology
H0 = 70  # Hubble constant
Om0 = 0.3  # Matter density parameter

# Choice of survey/instrument
#instrument = 'LSST'
instrument = 'SDSS'

# Choice of SED templates
sed_selection = {
    'Assef': ['wave', 'agn', 'ell', 'sfg', 'irr'],
    #'Kirkpatrick': ['sfg1', 'sfg2', 'sfg3']
    # Add more sources and their associated components as needed
}

# Choice of extinction laws
ext_selection = ['mod_calzetti']
#ext_selection = {
#    'modCalzetti': ['kappa'],
#    'Carroll2023': ['kappa'],
#}

# Choice of Quasar Luminosity Function
qlf_selection = {
    'prescription': 'Shen2020QLF',
    'redshift_range': [0.001, 7.0],
    'redshift_bin': 0.001,
    'luminosity_range': [43, 48],
    'luminosity_bin': 0.001,
}

# Choice of AGN selection criteria
agn_selection = 'None'






# NOTES: AGN maximum luminosity. mass = 2.7e11, Ledd = 3.402e49 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.456L.109K/abstract)