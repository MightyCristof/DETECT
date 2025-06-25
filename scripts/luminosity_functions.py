# scripts/luminosity_functions
from scripts.common import const, np, pd, u
#import resources.parameters as params
from scripts.global_variables import bp
from scripts.global_variables import cosmo


import pdb
#########################################
#
# Quasar Luminosity Function
#
#########################################
class QuasarLuminosityFunction:
    def __init__(self, log_luminosity, redshift, dluminosity, volume, dvolume, 
                                       prescription='Shen+2020', fit_type='FitB', z_ref=2.0):
        self.log_luminosity = log_luminosity
        self.z = redshift
        self.dluminosity = dluminosity
        self.volume = volume
        self.dvolume = dvolume
        self.prescription = prescription
        self.fit_type = fit_type
        self.z_ref = z_ref
        self.fit_param = self.load_fit_parameters()
        #self.phi = self.construct_phi_matrix(self.redshift)
        self.phi = self.calculate_phi(self.log_luminosity, self.z)
        self.phi_z = self.marginalize_luminosity(phi=self.phi)

    def load_fit_parameters(self):
        # Shen+2020 Table 4
        return {'FitA': {'a0': (0.8569, 0.0247, -0.0253),
                         'a1': (-0.2614, 0.0162, -0.016),
                         'a2': (0.0200, 0.0011, -0.0011),
                         'b0': (2.5375, 0.0177, -0.0187),
                         'b1': (-1.0425, 0.0164, -0.018),
                         'b2': (1.1201, 0.0199, -0.0207),
                         'c0': (13.0088, 0.0090, -0.009),
                         'c1': (-0.5759, 0.0018, -0.002),
                         'c2': (0.4554, 0.0028, -0.0027),
                         'd0': (-3.5426, 0.0235, -0.020),
                         'd1': (-0.3936, 0.0070, -0.007),
                         },
                'FitB': {'a0': (0.3653, 0.0115, -0.0114),
                         'a1': (-0.6006, 0.0422, -0.0417),
                         'a2': (None, None, None),
                         'b0': (2.4709, 0.0163, -0.0169),
                         'b1': (-0.9963, 0.0167, -0.016),
                         'b2': (1.0716, 0.0180, -0.0181),
                         'c0': (12.9656, 0.0092, -0.008),
                         'c1': (-0.5758, 0.0020, -0.001),
                         'c2': (0.4698, 0.0025, -0.0026),
                         'd0': (-3.6276, 0.0209, -0.020),
                         'd1': (-0.3444, 0.0063, -0.006),
                         }
                }[self.fit_type]

    def load_model_parameters(self, z):
        par = self.fit_param
        zref = self.z_ref
        # Shen+2020 Equation 14
        def T(n, x):
            if n == 0:
                return 1
            elif n == 1:
                return x
            elif n == 2:
                return 2 * x**2 - 1
            else:
                raise ValueError("Invalid value for n. It should be 0, 1, or 2.")
        if self.fit_type == 'FitA': # Shen+2020 Equation 14
            gamma1 = par['a0'][0]*T(0,1+z) + par['a1'][0]*T(1,1+z) + par['a2'][0]*T(2,1+z)
        elif self.fit_type == 'FitB': # Shen+2020 Equation 16
            gamma1 = par['a0'][0]*( (1+z)/(1+zref) )**par['a1'][0]
        gamma2 = (2*par['b0'][0]) / ( ((1+z)/(1+zref))**par['b1'][0] + ((1+z)/(1+zref))**par['b2'][0] )
        logl_star = (2*par['c0'][0]) / ( ((1+z)/(1+zref))**par['c1'][0] + ((1+z)/(1+zref))**par['c2'][0] )
        logp_star = par['d0'][0]*T(0,1+z) + par['d1'][0]*T(1,1+z)
        return gamma1, gamma2, logl_star, logp_star
    
    #def calculate_phi(self, log_lum, z):
    #    gamma1, gamma2, log_ell_star, log_phi_star = self.load_model_parameters(z)
    #    ell_bol = 10**log_lum
    #    # Convert to from Solar to cgs here
    #    ell_star = 10.0**log_ell_star * const.L_sun.to(u.erg / u.s).value
    #    phi_star = 10.0**log_phi_star
    #    ell_over_ell_star = ell_bol/ell_star        
    #    # Shen+2020 Equation 11
    #    phi_at_z = phi_star / (ell_over_ell_star**gamma1 + ell_over_ell_star**gamma2)
    #    return phi_at_z
    #
    #def construct_phi_matrix(self, z):        
    #    phi = np.array([self.calculate_phi(self.log_luminosity, zi) for zi in z]).T
    #    return phi
    
    def calculate_phi(self, log_lum, z):
        # Ensure model parameters are always NumPy arrays for consistency
        gamma1, gamma2, log_ell_star, log_phi_star = map(np.atleast_1d, self.load_model_parameters(z))
        # Convert to from Solar to cgs here
        ell_bol = 10**log_lum[:,np.newaxis] / (const.L_sun.to(u.erg / u.s).value)
        ell_star = 10**log_ell_star[np.newaxis,:]
        ell_over_ell_star = ell_bol / ell_star
        phi_star = 10.0**log_phi_star[np.newaxis,:]
        # Shen+2020 Equation 11
        phi_bol = phi_star / (ell_over_ell_star**gamma1[np.newaxis,:] + ell_over_ell_star**gamma2[np.newaxis,:])
        return phi_bol
    
    def marginalize_luminosity(self, phi, lum_range=None):
        # Marginalize over specified luminosity range
        if lum_range is not None:
            ilum = np.abs(self.log_luminosity[:,np.newaxis] - lum_range).argmin(axis=0)
            ilum = np.arange(ilum[0],ilum[-1]+1)
        # Otherwise grab all indices
        else:
            ilum = np.arange(len(self.log_luminosity))
        phi_z = np.trapz(phi[ilum,:], x=self.log_luminosity[ilum], axis=0)
        return phi_z

    def get_phi(self, log_lum=None, z=None, lum_range=False, z_range=False):
        if log_lum is not None:
            # Check if luminosity is within bounds
            if np.any(log_lum < self.log_luminosity.min()) or np.any(log_lum > self.log_luminosity.max()):
                raise ValueError(f"Luminosity {log_lum} is out of bounds [{self.log_luminosity.min():.3f}, {self.log_luminosity.max():.3f}]")
            # Grab indices at closest luminosity
            ilum = np.abs(self.log_luminosity[:,None] - log_lum).argmin(axis=0)
            # If keyword set, take all indices (inclusive) between end points
            if lum_range == True: 
                ilum = np.arange(ilum[0],ilum[-1]+1)
        # Otherwise grab all indices
        else: 
            ilum = np.arange(len(self.log_luminosity))
        if z is not None:
            # Check if redshift is within bounds
            if np.any(z < self.z.min()) or np.any(z > self.z.max()):
                raise ValueError(f"Redshift {z} is out of bounds [{self.z.min():.3f}, {self.z.max():.3f}]")
            # Grab indices at closes redshift
            iz = np.abs(self.z[:,None] - z).argmin(axis=0)
            # If keyword set, take all indices (inclusive) between end points
            if z_range == True:
                iz = np.arange(iz[0],iz[-1]+1)
        # Otherwise grab all indices
        else: 
            iz = np.arange(len(self.z))
        return self.phi[np.ix_(ilum,iz)].squeeze()
    
    def get_volume(self, z):
        # Check if z is within the bounds
        if z < self.z.min() or z > self.z.max():
            raise ValueError(f"Redshift {z} is out of bounds [{self.z.min():.3f}, {self.z.max():.3f}]")
        # Return volume at closest redshift
        iz = np.abs(self.z - z).argmin()
        return self.volume[iz]
    
    def sample_luminosity(self, phi, log_lum, ndraw=1):
        dlog_lum = np.diff(log_lum)
        dlog_lum = np.append(dlog_lum, dlog_lum[0])
        pdf_phi = phi / np.sum(phi * dlog_lum[(...,) + (None,)*(phi.ndim == 2)], axis=0)
        cdf_phi = np.cumsum(pdf_phi * dlog_lum[(...,) + (None,)*(phi.ndim == 2)], axis=0)
        if phi.ndim == 2:
            ndraw = phi.shape[1]
            samp = np.zeros(ndraw)
            random_draw = np.random.rand(ndraw)
            for i in range(ndraw):
                samp[i] = np.interp(random_draw[i], cdf_phi[:,i], log_lum)
        else:
            random_draw = np.random.rand(ndraw)
            samp = np.interp(random_draw, cdf_phi, log_lum)
        return samp

    def sample_redshift(self, phi, z, dvol, ndraw=1):
        dz = np.diff(z)
        dz = np.append(dz, dz[0])
        pdf_phi = (phi * dvol) / np.sum(phi * dvol * dz)
        cdf_phi = np.cumsum(pdf_phi * dz)
        random_draw = np.random.rand(ndraw)
        samp = np.interp(random_draw, cdf_phi, z)
        return samp

    def sample_lz_pairs(self, ndraw=1):
        zsamp = self.sample_redshift(self.phi_z, self.z, self.dvolume, ndraw)
        lsamp = self.sample_luminosity(self.get_phi(z=zsamp), self.log_luminosity)
        return lsamp, zsamp

    def sample_high_luminosity(self, phi, log_lum, ndraw=1):
        dlog_lum = np.diff(log_lum)
        dlog_lum = np.append(dlog_lum, dlog_lum[0])
        log_phi = np.log10(phi)
        phi_ = log_phi - np.min(log_phi)
        pdf_phi = phi_ / np.sum(phi_ * dlog_lum[(...,) + (None,)*(phi_.ndim == 2)], axis=0)
        cdf_phi = np.cumsum(pdf_phi * dlog_lum[(...,) + (None,)*(phi_.ndim == 2)], axis=0)
        if phi_.ndim == 2:
            ndraw = phi_.shape[1]
            samp = np.zeros(ndraw)
            random_draw = np.random.rand(ndraw)
            for i in range(ndraw):
                samp[i] = np.interp(random_draw[i], cdf_phi[:,i], log_lum)
        else:
            random_draw = np.random.rand(ndraw)
            samp = np.interp(random_draw, cdf_phi, log_lum)
        return samp









    # DEFUNCT METHODS... PROBABLY
#    def sample_luminosity_OLD(self, z):
#        # Indices of closest matched redshift
#        iz = np.abs(self.redshift[:,np.newaxis] - z).argmin(axis=0)
#        dlog_lum = np.diff(self.log_luminosity)
#        dlog_lum = np.append(dlog_lum, dlog_lum[0])
#        dlog_lum = dlog_lum[:, np.newaxis]
#        phi_at_z = self.phi[:, iz]
#        pdf_phi = phi_at_z / np.sum(phi_at_z * dlog_lum, axis=0)
#        cdf_phi = np.cumsum(pdf_phi * dlog_lum, axis=0)
#        random_draw = np.random.rand(len(z))
#        samp = np.array([np.interp(random_draw[i], cdf_phi[:, i], self.log_luminosity) for i in range(len(z))])
#        return samp

#    def sample_luminosity_at_redshift(self, z, ndraw=1, test_high_lum=False):
#        phi = self.get_phi(z=z)# * self.get_volume(z)
#        dlog_lum = np.diff(self.log_luminosity)
#        dlog_lum = np.append(dlog_lum, dlog_lum[0])
#        phi_pdf = phi / np.sum(phi * dlog_lum)
#        phi_cdf = np.cumsum(phi_pdf * dlog_lum)
#        random_draw = np.random.rand(ndraw)
#        samp = np.interp(random_draw, phi_cdf, self.log_luminosity)
#        return samp
    
#    def sample_redshift_at_luminosity(self, lum, ndraw=1):
#        phi = self.get_phi(lum=lum) * self.get_volume(z)
#        dlog_lum = np.diff(self.log_luminosity)
#        dlog_lum = np.append(dlog_lum, dlog_lum[0])
#        phi_pdf = phi / np.sum(phi * dlog_lum)
#        phi_cdf = np.cumsum(phi_pdf * dlog_lum)
#        random_draw = np.random.rand(ndraw)
#        samp = np.interp(random_draw, phi_cdf, self.log_luminosity)
#        return samp


#########################################
#
# Galaxy Luminosity Function
#
#########################################
#class SchechterLuminosityFunction:
#    def __init__(self, phi_star=0.012, l_star=1e10, alpha=-1.2):
#        self.phi_star = phi_star*(params.H0/100)**(-3)
#        self.l_star = l_star
#        self.alpha = alpha
#    def compute_luminosity_function(self, lum):
#        return self.phi_star * (lum / self.l_star)**self.alpha * np.exp(-lum / self.l_star) / self.l_star



class ObservedLuminosityFunction:
    def __init__(self, log_luminosity, redshift, dluminosity, volume, dvolume, flux0_limits):
        self.log_luminosity = log_luminosity
        self.z = redshift
        self.dluminosity = dluminosity
        self.volume = volume
        self.dvolume = dvolume
        self.ledges = self.set_binning(self.log_luminosity, 0.2)
        self.zedges = self.set_binning(self.z, 0.2)
        #self.lbin = 0.2
        #self.iledge = self.set_binning(self.log_luminosity, 0.2)
        #self.ledges = np.round(self.log_luminosity[self.iledge],6)
        #self.izedge = self.set_binning(self.z, 0.2)
        #self.zbin = 0.2
        #self.zedges = np.round(self.z[self.izedge],6)
        self.flux_limits = self.set_survey_flux_limits(flux0_limits)

    def set_binning(self, arr, binsz):
        nbins = np.ceil((arr.max()-arr.min())/binsz)+1
        edges = np.arange(nbins)*binsz + arr.min()
        return edges
        #start = np.floor(arr.min() / binsz) * binsz
        #end = np.ceil(arr.max() / binsz) * binsz
        #bin_edges = np.round(np.arange(start,end,binsz),6)
        ## find indices in array closest to each bin edge
        #ind = np.searchsorted(arr,bin_edges,side="left")
        #return ind

    def set_survey_flux_limits(self, flux0_limits):
        #flux_limit0 = {key: self.convert_mag_to_flx(self.depths[key]) for key in self.depths}
        flux_limits = {mode: {} for mode in flux0_limits}
        # The distance ratio might not be right...
        #dl_ratio = (self.dluminosity/3.08567758e+19).value**2             # z=0
        #dl_ratio = (self.dluminosity/1.3215183830917299e+20).value**2     # z=1e-8 minimum for some calculations
        dl_ratio = (self.dluminosity/self.dluminosity[0]).value**2

        for mode in flux0_limits:
            for filt, flux0 in flux0_limits[mode].items():
                flux_limits[mode][filt] = flux0 * dl_ratio        
        #flux_limit = {mode: {filt: flux * distance_ratio for filt, flux in flux0_limits[mode].items()}
        #              for mode in flux0_limits}
        pdb.set_trace()
        return flux_limits
    
    def differential_number_counts(self, mags, abs_mags, z, detect, area):
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        xdr3 = np.array([15.21997,15.47582,15.72543,15.96880,16.22465,16.47426,16.71763,16.96724,17.22309,17.47270,17.72231,17.97192,18.22153,18.46490,18.71451,18.97036])
        ydr3 = np.array([0.00101,0.00202,0.00101,0.00611,0.00814,0.01509,0.04329,0.07139,0.11403,0.21367,0.33768,0.55688,0.82565,1.25047,1.89387,2.66238])
        x2qz = np.array([16.08112,16.33073,16.57410,16.82995,17.07956,17.32917,17.57878,17.82215,17.92200,18.16537,18.41498,18.66459,18.92044,19.16381,19.41342,19.66927,19.91264,20.16849,20.41186])
        y2qz = np.array([0.00658,0.00658,0.02625,0.05024,0.07940,0.12285,0.24800,0.32707,0.44534,0.73443,1.01070,1.66679,2.19819,3.09019,3.86418,4.58160,5.49034,5.66848,5.79044])
        xhiz = np.array([17.07332,17.32293,17.57878,17.82215,18.07176,18.32761,18.57722,18.82683,19.07644,19.31981,19.56942,19.81903,20.07488])
        yhiz = np.array([0.00102,0.00303,0.00905,0.00905,0.01809,0.03425,0.05187,0.07773,0.12416,0.16375,0.22060,0.28179,0.38779])
        nloz = np.array([1,3,2,10,12,23,67,107,174,327,511,849,1257,1923,2870,4028,900,216,55,20])
        nhiz = np.array([1,4,12,12,26,47,72,106,162,206,275,322,410])

        detected_mags = mags[detect]
        binsz = 0.25
        obsv = {'mags':mags, 'detected_mags':detected_mags}
        mark = dict(zip(obsv.keys(), ['o','*']))
        color = dict(zip(obsv.keys(), ['C0','C1']))
        label = dict(zip(obsv.keys(), ['Simulated—All','Simulated—Detected']))
        fig, ax = plt.subplots()
        ax.scatter(xdr3,ydr3,marker='o',color='dimgrey',label=r'SDSS DR3 $0.3<z<2.2$',zorder=1)
        ax.scatter(x2qz,y2qz,marker='s',color='lightgrey',label='2QZ/6QZ',zorder=0)
        ax.scatter(xhiz,yhiz,marker='^',color='darkgrey',label=r'SDSS DR3 $3<z<5$',zorder=2)
        for obs, mag_arr in obsv.items():
            counts, bin_edges = np.histogram(mag_arr, bins=np.arange(np.min(mag_arr), np.max(mag_arr)+binsz, binsz))
            bins = bin_edges[:-1] + 0.5*binsz
            diff_counts = counts / (area*binsz)
            ax.plot(bins,diff_counts,marker=mark[obs],color=color[obs],linestyle='None',label=label[obs])
        ax.set_yscale('log')
        plt.ylabel(r'N($i$) [deg$^{-2}$ 0.25 mag$^{-1}$]')
        plt.xlabel(r'$i$')
        plt.legend()
        plt.show()
    
    def cumulative_number_counts(self, mags, abs_mags, z, detect, area, outline=False):
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        
        xdr3 = np.array([15.21997,15.46958,15.71919,15.97504,16.22465,16.47426,16.71763,16.97348,17.22309,17.47270,17.72231,17.97192,18.22153,18.47114,18.72075,18.97036])
        ydr3 = np.array([0.00101,0.00306,0.00408,0.01022,0.01833,0.03370,0.07746,0.14601,0.26187,0.47554,0.82171,1.36798,2.22156,3.47587,5.37130,8.09681])
        x2qz = np.array([16.07488,16.33073,16.57410,16.82371,17.07956,17.32917,17.57878,17.82839,17.90952,18.15913,18.41498,18.66459,18.92044,19.16381,19.41342,19.66303,19.91264,20.16225,20.41186])
        y2qz = np.array([0.00654,0.01328,0.03961,0.09106,0.16952,0.29292,0.54529,0.87461,1.31772,2.06171,3.10793,4.80271,6.97520,10.25653,14.17420,18.63946,23.61541,29.91974,35.62590])
        xhiz = np.array([17.07332,17.32293,17.57254,17.82215,18.07176,18.32137,18.57098,18.82059,19.07020,19.31981,19.57566,19.83151,20.07488])
        yhiz = np.array([0.00103,0.00410,0.01349,0.02275,0.04130,0.07594,0.13123,0.20788,0.32932,0.49642,0.72997,1.02142,1.41154])
        nloz = np.array([1,3,2,10,12,23,67,107,174,327,511,849,1257,1923,2870,4028,900,216,55,20])
        nhiz = np.array([1,4,12,12,26,47,72,106,162,206,275,322,410])

        detected_mags = mags[detect]
        binsz = 0.25
        obsv = {'mags':mags, 'detected_mags':detected_mags}
        mark = dict(zip(obsv.keys(), ['o','*']))
        color = dict(zip(obsv.keys(), ['C0','C1']))
        label = dict(zip(obsv.keys(), ['Simulated—All','Simulated—Detected']))
        fig, ax = plt.subplots()
        ax.scatter(xdr3,ydr3,marker='o',color='dimgrey',label=r'SDSS DR3 $0.3<z<2.2$',zorder=1)
        ax.scatter(x2qz,y2qz,marker='s',color='lightgrey',label='2QZ/6QZ',zorder=0)
        ax.scatter(xhiz,yhiz,marker='^',color='darkgrey',label=r'SDSS DR3 $3<z<5$',zorder=2)
        for obs, mag_arr in obsv.items():
            counts, bin_edges = np.histogram(mag_arr, bins=np.arange(np.min(mag_arr), np.max(mag_arr)+binsz, binsz))
            bins = bin_edges[:-1] + 0.5*binsz
            diff_counts = counts / (area*binsz)
            #cum_counts = np.cumsum(counts[::-1])[::-1] / area
            cum_counts = np.cumsum(counts) / area
            ax.scatter(bins,cum_counts,marker=mark[obs],color=color[obs],label=label[obs])
        ax.set_yscale('log')
        plt.ylabel(r'N($<i$) [deg$^{-2}$]')
        plt.xlabel(r'$i$')
        plt.legend()
        plt.show()    

    def rebuild_source_sed(self, sed, src):        
        # Build source observed SED
        src_sed = []
        for comp in src.components:
            for template in sed.components[comp]:
                key = f"{comp}_{template}_norm"
                norm = src[key]
                sed_template = sed.components[comp][template]
                src_sed.append(norm * sed_template/(1+src.z))
        # OK, from here src_sed is the same as src.sed in the observed frame
        src_sed = np.sum(np.stack(src_sed),axis=0)

        # Build wave/freq domain from observed to rest
        # Build redshift array from 0 to source redshift
        # THIS MIGHT NOT WORK, THINK ABOUT "IS IT IN THIS REDSHIFT BIN"
        # NEEDS TO GO FARTHER. HOW TO DO CLEVER CUTTOFF?]
        # NOTE: you are adding z==0 into the mix
        #zarr = np.concatenate((np.array([0]),self.z[:iz+1]))
        zarr = np.concatenate((np.array([0]),self.z))
        iz = np.abs(zarr-src.z).argmin()
        # This creates a wavelength array that is build restframe OUT
        # Problem is, the src_sed is in observed frame at z=0
        # So this is essentially backwards to how we're doing it
        # I need to flip the indices
        # in this scheme, rest_wave[:,0] =close enough= obsv_wave[:,iz]
        # for now we'll go with observed wave building out to rest
        #rest_wave = np.outer(sed.domain['wavelength']['wave'],1+zarr)
        #obsv_wave = np.outer(sed.domain['wavelength']['wave']*(1+src.z),(1+zarr)**-1)
        zsed = np.outer(src_sed,1+zarr)
        wave = np.outer(sed.domain['wavelength']['wave']*(1+src.z),(1+zarr)**-1)
        # create the sed array that exists at each redshift from observed to rest

        zflx = []
        for i in range(len(zarr)):
            zflx.append(bp.convolve_with_bandpass(wave[:,i],zsed[:,i]))
        zflx = pd.DataFrame(zflx)

        import matplotlib.pyplot as plt
        colors = ['violet','dodgerblue','green','yellow','orange','red']
        labels = zflx.columns.tolist()
        for i, filt in enumerate(self.flux_limits['single'].items()):
            plt.plot(self.z,filt[1],colors[i],linestyle='--',label=labels[i]+' flux limit')
        plt.axvline(x=src.z,color='black',label='Redshift',zorder=0)
        for i in range(zflx.columns.size):
            plt.scatter(zarr,zflx.iloc[:,i],marker='.',color=colors[i],alpha=0.5,label=labels[i])
        plt.yscale('log')
        plt.xlabel(r'$z$')
        plt.ylabel(r'$F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
        plt.legend()
        plt.show()

        # OK this doesn't match the survey flux limits, so there must be an issue with how i'm creating olf.flux_limits

        pdb.set_trace()
        # then convolved to get the fluxes
        # filts = []
        # filts.append(bp.convolve_with_bandpass(wave[:,i],dat_sed[:,i]))
        # filts = pd.DataFrame(filts)
        # now you can call them like filts.u over the entire redshift array
        # and then use that to find the bisections of the flux limits
        # but you might need to consider how the flux limits should actually be at z=0 (or z=1e-8 at the very least)

        print('Rebuild in progress!')
    #    iz = np.abs(self.z - src.z).argmin()
    #    z_arr = self.z[0:iz+1]

        flux = {}
        

    def calculate_vmax(self, source_det, source_lum, source_z, source_fluxes):
        # Luminosity and redshift binning
        source_det = data.r_detected_coadd
        source_lum = np.log10(data.agn_luminosity)
        source_z = data.z
        source_flux = data.r_flux
        flux_limit = olf.flux_limits['coadd']['r']

        #phi_est = np.zeros((len(olf.ilbin)-1, len(olf.izbin)-1))

        nzbins = len(olf.zedges)
        nlbins = len(olf.ledges)
        for i in range(nzbins-1):

            #zbin = olf.z[olf.iz_edges[i]:olf.iz_edges[i+1]]
            #iz = source_det & (source_z >= zbin[0]) & (source_z < zbin[-1])
            
            iz = (source_z >= olf.zedges[i]) & (source_z < olf.zedges[i+1])
            flux_limit_bin = flux_limit[olf.izedge[i]:olf.izedge[i+1]]

            for k in range(nlbins-1):
                lmin, lmax = olf.ledges[k], olf.ledges[k+1]
                il = (source_lum >= lmin) & (source_lum < lmax)

                ibin = iz & il
                N = ibin.sum()
                if N == 0:
                    continue
                
                fluxes = source_flux[ibin]
                zmaxL = np.full(N, zbin[-1])

                for jj, f in enumerate(fluxes):
                    if (f > flux_limit_bin[0]) and (f < flux_limit_bin[-1]):
                        mask = (flux_limit_bin >= f)
                        if np.any(mask):
                            z_detect = zbin[mask]
                            f_detect = flux_limit_bin[mask]
                            z_interp = np.interp(f, f_detect[::-1], z_detect[::-1])
                            zmaxL[jj] = min(z_interp, zmax)


    #    integrating over L and z
    #        limits of integration L bins
    #        limits of integration z_min to max z where object L can still be detected in redshift bin
    #
    #        so we need differential volumes for the bin edges, and likely have to interpolate for each source
    #            actually, likely need to integrate to each source... dV/dz
        
        return None

    def calculate_vmax(self, lum, z, flux_limit):

        z = qlf.z
        dL = qlf.dluminosity.value
        F_lim = bp.flux_limits['coadd']['r']
        F_limz = F_lim * (dL/dL[0])**2
        L_lim = 4*np.pi*dL**2*F_lim # (1+z) for F_obsv
        L_0 = 4*np.pi*dL[0]**2*F_lim

        return None
