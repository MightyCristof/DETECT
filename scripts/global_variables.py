# scripts/global.py
from scripts.common import FlatLambdaCDM, np, u
from scripts.photometry import SpectralEnergyDistribution, Bandpass
import resources.parameters as params


# Set global cosmology
cosmo = FlatLambdaCDM(H0=params.H0, Om0=params.Om0)
# Set choice of SED templates 
sed = SpectralEnergyDistribution(params.sed_selection, ext_choice=params.ext_selection)
# Set choice of survey
bp = Bandpass(params.instrument)


from scripts.luminosity_functions import QuasarLuminosityFunction
from scripts.selection_criteria import AGNSelection

# Construct QLF
log_lum = np.round(np.arange(*params.qlf_selection['luminosity_range'],params.qlf_selection['luminosity_bin']),6)
z = np.round(np.arange(*params.qlf_selection['redshift_range'],params.qlf_selection['redshift_bin']),6)
dlum = cosmo.luminosity_distance(z).to(u.cm)
vol = cosmo.comoving_volume(z)
dvol = cosmo.differential_comoving_volume(z)
qlf = QuasarLuminosityFunction(log_luminosity=log_lum, 
                               redshift=z, dluminosity=dlum, volume=vol, dvolume=dvol, 
                               prescription=params.qlf_selection['prescription'])

# Set AGN selection criteria
select = AGNSelection()