# resources/surveys.py


sdss = {
    # Wavelength shorthand
    'filters': ['u', 'g', 'r', 'i', 'z'],
    # Path to filter files. Assumes column 1 == wavelength, column 2 == throughput
    'file_paths': {
        'u': 'resources/surveys/sdss/01_u.dat',
        'g': 'resources/surveys/sdss/02_g.dat',
        'r': 'resources/surveys/sdss/03_r.dat',
        'i': 'resources/surveys/sdss/04_i.dat',
        'z': 'resources/surveys/sdss/05_z.dat'
    },
    # Flux zero points in Jansky
    'zero_points': {
        'u': 3767.0,
        'g': 3631.0,
        'r': 3631.0,
        'i': 3631.0,
        'z': 3565.0
    },
    # Effective wavelength per filter
    'wave_eff': {
        'u': 0.354,
        'g': 0.475,
        'r': 0.622,
        'i': 0.763,
        'z': 0.905
    },
    # Wavelength conversion factor to meters
    'wave_conversion': 1.0e-6,
    'depths' : {
        # Single visit depths
        'single': {'u': 22.15,
                   'g': 23.13,
                   'r': 22.70,
                   'i': 22.20,
                   'z': 20.71
        }
    },
}


lsst = {
    'filters': ['u', 'g', 'r', 'i', 'z', 'y'],
    'file_paths': {
        'u': '../resources/surveys/lsst/01_u.dat',
        'g': '../resources/surveys/lsst/02_g.dat',
        'r': '../resources/surveys/lsst/03_r.dat',
        'i': '../resources/surveys/lsst/04_i.dat',
        'z': '../resources/surveys/lsst/05_z.dat',
        'y': '../resources/surveys/lsst/06_y.dat'
    },
    'zero_points': {
        'u': None,
        'g': None,
        'r': None,
        'i': None,
        'z': None,
        'y': None
    },
    # Effective wavelength per filter
    'wave_eff': {
        'u': None,
        'g': None,
        'r': None,
        'i': None,
        'z': None,
        'y': None
    },
    # Wavelength conversion factor to meters
    'wave_conversion': None,
    'depths' : {
        # Single visit depths
        'single': {'u': 23.8,
                   'g': 24.5,
                   'r': 24.03,
                   'i': 23.41,
                   'z': 22.74,
                   'y': 22.96
        },
        # 10-year survey depths
        'coadd': {'u': 25.6,
                  'g': 26.9,
                  'r': 26.9,
                  'i': 26.4,
                  'z': 25.6,
                  'y': 24.8
        },
        # Deep drilling field depths
        # SINGLE VALUES HERE AS PLACEHOLDERS
        'ddf': {'u': 23.8,
                'g': 24.5,
                'r': 24.03,
                'i': 23.41,
                'z': 22.74,
                'y': 22.96
        }
    },
}
# Values from https://www.lsst.org/scientists/keynumbers
