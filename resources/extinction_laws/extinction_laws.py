#scripts/extinction_laws.py


mod_calzetti = {
    # Path to extinction law files
    'file_paths': {
        'file01': 'resources/extinction_laws/Modified Calzetti.dat'
    },
    # Column names per file
    'columns': {
        'file01': ['wave', 'kappa']
    },
    # Wavelength conversion factor to meters
    'wave_conversion': {
        'file01': {'wave': 1e-10}
    },
    # Specify SED components by type to host extinction
    'frame': {
        'wavelength': ['wave'],
        'frequency': None,
    },
    'flux': {
        'agn': ['agn', 'agn2'],
        'galaxy': None,
        'stellar': None, 
        'neo': None,
    },
}

carroll2023 = {
    # Path to extinction law files
    'file_paths': {
        'file01': 'resources/extinction_laws/Carroll2023.dat'
    },
    # Column names per file
    'columns': {
        'file01': ['wave', 'kappa']
    },
    # Wavelength conversion factor to meters
    'wave_conversion': {
        'file01': {'wave': 1e-10}
    },
    # Specify SED components by type to host extinction
    'frame': {
        'wavelength': ['wave'],
        'frequency': None,
    },
    'flux': {
        'agn': None,
        'galaxy': ['ell', 'sfg', 'irr'],
        'stellar': None, 
        'neo': None,
    },
}

