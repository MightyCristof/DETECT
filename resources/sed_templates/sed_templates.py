# resources/sed_templates.py


# Full range of SED component...
component_keys = {
    'frame': ['wavelength', 'frequency'],
    'flux': ['agn', 'galaxy', 'star', 'neo'],
}


assef = {
    # Path to SED template files
    'file_paths': {
        'file01': 'resources/sed_templates/Assef+2013/lrt_templates.dat'
    },
    # Column names per file
    'columns': {
        'file01': ['wave', 'agn', 'agn2', 'ell', 'sfg', 'irr']
    },
    # Wavelength conversion factor to meters
    'wave_conversion': {
        'file01': {'wave': 1e-6}
    },
    # Flux conversion factors to units of [erg/s/cm^2/Hz]
    'flux_conversions': {
        'file01': {
            'agn': 1e-11,
            'agn2': 1e-14,
            'ell': 1e-16,
            'sfg': 1e-17,
            'irr': 1e-15
        }
    },
    # Specify SED components by type
    'frame': {
        'wavelength': ['wave'],
        'frequency': None,
    },
    'flux': {
        'agn': ['agn', 'agn2'],
        'galaxy': ['ell', 'sfg', 'irr'],
        'star': None, 
        'neo': None,
    },
}


kirkpatrick = {
    # Path to SED template files
    'file_paths': {
        'file01': 'resources/sed_templates/Assef+2013/lrt_templates.dat',
        'file02': 'resources/sed_templates/Assef+2013/lrt_templates.dat'
    },
    # Column names per file
    'columns': {
        'file01': ['wave', 'agn', 'agn2', 'ell', 'sfg', 'irr'],
        'file02': ['wave', 'agn', 'agn2', 'ell', 'sfg', 'irr'],
    },
    # Wavelength conversion factor to meters
    'wave_conversion': {
        'file01': {'wave': 1e-6},
        'file02': {'wave': 1e-6},
    },
    # Flux conversion factors to units of [erg/s/cm^2/Hz]
    'flux_conversions': {
        'file01': {
            'agn': 1e-11,
            'agn2': 1e-14,
            'ell': 1e-16,
            'sfg': 1e-17,
            'irr': 1e-15
        },
        'file02': {
            'agn': 1e-11,
            'agn2': 1e-14,
            'ell': 1e-16,
            'sfg': 1e-17,
            'irr': 1e-15
        },
    },
    # Specify SED components by type
    'frame': {
        'wavelength': ['wave'],
        'frequency': None,
    },
    'flux': {
        'agn': ['agn', 'agn2'],
        'galaxy': ['ell', 'sfg', 'irr'],
        'star': None, 
        'neo': None,
    },
}
