#scripts/tests.py
from scripts.common import np


# Check against in src.luminosity, tests normalizations and flux integrations
def check_luminosity(flux, obsv_freq, z, dlum):
    rest_freq = obsv_freq * (1 + z)
    luminosity_density = np.trapz(flux[::-1], x=rest_freq[::-1])
    l_bol = 4 * np.pi * dlum**2 * luminosity_density
    return l_bol
