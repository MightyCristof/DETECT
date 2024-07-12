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
        'u': '../resources/surveys/lsst/total_u.dat',
        'g': '../resources/surveys/lsst/total_g.dat',
        'r': '../resources/surveys/lsst/total_r.dat',
        'i': '../resources/surveys/lsst/total_i.dat',
        'z': '../resources/surveys/lsst/total_z.dat',
        'y': '../resources/surveys/lsst/total_y.dat'
    },
    'zero_points': {
        'u': 3631.0,
        'g': 3631.0,
        'r': 3631.0,
        'i': 3631.0,
        'z': 3631.0,
        'y': 3631.0,
    },
    # Effective wavelength per filter
    'wave_eff': {
        'u': 367.1,
        'g': 482.7,
        'r': 622.3,
        'i': 754.6,
        'z': 869.1,
        'y': 971.2
    },
    # Wavelength conversion factor to meters
    'wave_conversion': 1.0e-9,
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
# Throughput from https://github.com/lsst/throughputs/tree/main/baseline
