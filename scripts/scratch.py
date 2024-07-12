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



