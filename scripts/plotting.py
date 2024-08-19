# scripts/plotting.py
from scripts.common import np, plt
from scripts.global_variables import bp


# Plot Source SED and template components
def plot_sed(src, phot=None):
    plt.figure(figsize=(10, 6))
    
    # Plot galaxy components
    if hasattr(src, 'galaxy'):
        galaxy = src.galaxy
        color = dict(zip(galaxy.sed.keys(), ['orange', 'green', 'blue']))
        linesty = dict(zip(galaxy.sed.keys(), ['--',':','-.']))
        for component, flux in galaxy.sed.items():
            plt.plot(src.wave['obsv'], flux * src.freq['obsv'],
                     label=component.upper(), color=color[component], linestyle=linesty[component] )
    # Plot AGN components
    if hasattr(src, 'agn'):
        agn = src.agn
        color = dict(zip(agn.sed.keys(), ['purple', 'magenta']))
        linesty = dict(zip(agn.sed.keys(), ['solid',':']))
        for component, flux in agn.sed.items():
            plt.plot(src.wave['obsv'], flux * src.freq['obsv'],
                     label=component.upper(), color=color[component], linestyle=linesty[component])
    # Plot total SED
    plt.plot(src.wave['obsv'], src.sed * src.freq['obsv'], color='black', alpha=0.5)
    # Plot data points
    plt.scatter(bp.get_wavelengths(), np.array([flux for flux in src.photometry['flux'].values()]) * bp.get_frequencies(),
                s=100, c='black', alpha=0.5)

    # Formatting
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(1e-10,1e-2)
    #plt.ylim(1e-30, 1e-1)
    plt.xlabel(r'Observed Wavelength [$\mu$m]')
    plt.ylabel(r'$\nu F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$]')
    plt.legend()
    plt.grid(True)
    plt.show()


# Plot Source color-color
def plot_cc(mags, phot=None):
    color1 = mags[0,:]-mags[1,:]
    color2 = mags[1,:]-mags[2,:]
    cut1 = 0.9
    cut2 = 0.4

    plt.figure(figsize=(10, 6))
    # Plot data points
    plt.scatter(color1[0,:], color2[0,:], s=100, c='black', alpha=0.5)
    plt.axvline(x=cut1, color='red', linestyle='--', linewidth=2)
    plt.axhline(y=cut2, color='blue', linestyle='--', linewidth=2)

    # Calculate percentiles to dynamically set limits
    x_min, x_max = np.percentile(color1, [5, 95])
    y_min, y_max = np.percentile(color2, [5, 95])

    # Formatting
    #plt.xlim(x_min, x_max)
    #plt.ylim(y_min, y_max)
    plt.xlabel(r'u - g [AB mag]')
    plt.ylabel(r'g - r [AB mag]')
    plt.legend()
    plt.grid(True)
    plt.show()
