# scripts/global.py
from scripts.common import FlatLambdaCDM, np
from scripts.photometry import SpectralEnergyDistribution, Bandpass
from scripts.luminosity_functions import QuasarLuminosityFunction
import resources.parameters as params


# Set global cosmology
cosmology = FlatLambdaCDM(H0=params.H0, Om0=params.Om0)

# Set choice of SED templates 
sed = SpectralEnergyDistribution(params.sed_selection, ext_choice=params.ext_selection)

# Set choice of survey
bp = Bandpass(params.instrument)

# Construct QLF
qlf = QuasarLuminosityFunction(np.arange(*params.qlf_selection['redshift_range']),
                               np.arange(*params.qlf_selection['luminosity_range']),
                               prescription=params.qlf_selection['prescription'])


