import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skewnorm

# Generate skewed normal distribution
#def generate_skewed_distribution(size, mean=-0.6, stddev=0.8, skewness=-3, bounds=None):
#    data = skewnorm.rvs(a=skewness, loc=mean, scale=stddev, size=size)
#
#    if bounds is not None:
#        lower_bound, upper_bound = bounds
#        iout = np.where((data < lower_bound) | (data > upper_bound))[0]
#        if len(iout) > 0:
#            data[iout] = generate_skewed_distribution(len(iout), bounds=bounds)
#        return data
#
#edd_ratio = np.zeros(10000)
#for i in range(10000):
#    edd_ratio[i] = (generate_skewed_distribution(1, bounds=[-3.0, 0.0]))[0]
#
#
## Generate bounded normal distribution
#def generate_bounded_normal(size, mean=-0.6, stddev=0.8, bounds=None):
#    data = np.random.normal(loc=mean, scale=stddev, size=size)
#
#    if bounds is not None:
#        lower_bound, upper_bound = bounds
#        iout = np.where((data < lower_bound) | (data > upper_bound))[0]
#        if len(iout) > 0:
#            data[iout] = generate_bounded_normal(len(iout), bounds=bounds)
#        return data
#
#edd_ratio = np.zeros(10000)
#for i in range(10000):
#    edd_ratio[i] = (generate_bounded_normal(1, bounds=[-3.0, 0.0]))[0]
#
#edd_hist, bin_edges = np.histogram(edd_ratio, bins='auto')
#print(f"Mean: {np.mean(edd_ratio)}")
#print(f"Mode: {bin_edges[np.where(edd_hist == np.max(edd_hist))]}")
#
#
#
#plt.hist(edd_ratio, bins=50, density=True, alpha=0.75)
#plt.title('Bounded Normal Distribution')
#plt.xlabel('Values')
#plt.ylabel('Density')
#plt.show()


#def plot_hist(mn,sd,sz):
#    arr = np.random.normal(mn,sd,sz)
#    print(f"Min: {np.min(arr)}, Max: {np.max(arr)}")
#    plt.hist(arr, bins='auto')
#    plt.show()
#
#plot_hist(6,1,1000)


#coefficients = {'ell':None, 'sfg': None, 'irr':None}
#def scale_coefficients(ml_ratio, coefficients):
#    # Calculate weights based on linear interpolation between M/L ratios 2 and 10
#    if ml_ratio <= 2:
#        weight_irr = 1.0
#        weight_sfg = 0.0
#        weight_ell = 0.0
#    elif ml_ratio >= 10:
#        weight_irr = 0.0
#        weight_sfg = 0.0
#        weight_ell = 1.0
#    else:
#        weight_irr = (10 - ml_ratio) / 8.0
#        weight_sfg = 1.0 - abs(ml_ratio - 6.0) / 4.0
#        weight_ell = (ml_ratio - 2) / 8.0
#    
#    # Normalize coefficients to sum up to 1
#    total = weight_irr + weight_sfg + weight_ell
#    if total != 0:
#        weight_irr /= total
#        weight_sfg /= total
#        weight_ell /= total
#    
#    return {'irr': weight_irr, 'sfg': weight_sfg, 'ell': weight_ell}
#
## Example: Calculate scaled coefficients for different M/L ratios
#ml_ratios = np.arange(2,11)
##coefficients = {'irr': None, 'sfg': None, 'ell': None}
#
#for ml_ratio in ml_ratios:
#    scaled_coefficients = scale_coefficients(ml_ratio, coefficients)
#    print(f"M/L Ratio: {ml_ratio:.1f}")
#    print(f"Scaled Coefficients: {scaled_coefficients}")






    # Example flux limit of the survey
    flux_limit = 1e-16 * u.erg / (u.cm**2 * u.s)

    # Luminosity and redshift bins
    lum_bins = np.logspace(43, 48, 20) * u.erg / u.s
    z_bins = np.linspace(0.0, 7.0, 20)

    # Catalog of sources (these would come from your simulation)
    catalog = {
        'luminosity': np.random.uniform(1e43, 1e48, 1000) * u.erg / u.s,
        'redshift': np.random.uniform(0, 7, 1000)
    }

    # Initialize grid to hold counts
    N_detected = np.zeros((len(lum_bins) - 1, len(z_bins) - 1))

    # Loop through the luminosity and redshift bins
    for L_idx in range(len(lum_bins) - 1):
        Lmin = lum_bins[L_idx]
        Lmax = lum_bins[L_idx + 1]

        # Sources within the luminosity bin
        in_lum_bin = (catalog['luminosity'] >= Lmin) & (catalog['luminosity'] < Lmax)
        
        for z_idx in range(len(z_bins) - 1):
            zmin = z_bins[z_idx]
            zmax = z_bins[z_idx + 1]

            # Select the sources in this luminosity bin
            sources = catalog['luminosity'][in_lum_bin]
            
            # Compute the max redshift for each source:
            z_max = np.array([
                z for z in np.linspace(zmin, zmax, 100)
                if (L / (4 * np.pi * cosmo.luminosity_distance(z)**2)).to(u.erg / (u.cm**2 * u.s)) >= flux_limit
            ])
            
            # Filter only those still detectable in the redshift range
            detected_sources = (catalog['redshift'][in_lum_bin] >= zmin) & (catalog['redshift'][in_lum_bin] <= zmax)
            
            # Count the sources
            N_detected[L_idx, z_idx] = np.sum(detected_sources)










colors = ['violet','blue', 'green', 'lime green', 'goldenrod', 'red'] 
shaded_regions = []

for i, filt in enumerate(bp.filters):
    start_idx = (np.where(bp.filters[filt].thru > 1e-3))[0][0]
    end_idx = (np.where(bp.filters[filt].thru > 1e-3))[0][-1]
    shaded_regions.append((start_idx, end_idx, colors[i % len(colors)]))













fig, ax = plt.subplots()
ax.scatter(data.z,data.mag_i)
plt.show()


fig, ax = plt.subplots()
ax.scatter(data.z,data.mag_i)
xlim, ylim = ax.get_xlim(), ax.get_ylim()
buffer = 0.5
xra, yra = xlim[1]-xlim[0], ylim[1]-ylim[0]
ax.set_xlim(xlim[0]-xra*buffer,xlim[1]+xra*buffer)
ax.set_ylim(ylim[0]-yra*buffer,ylim[1]+yra*buffer)
plt.show()





