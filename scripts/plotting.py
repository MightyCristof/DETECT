# scripts/plotting.py
from scripts.common import const, cp, np, mcolors, patches, u
from scripts.common import AnnotationBbox, mcolors, mpimg, OffsetImage, plt
from scripts.common import MultipleLocator
from scripts.global_variables import bp
import pdb
import seaborn as sns

# Plot QLF vs. Luminoisty across redshift
def plot_qlf_lum(qlf):
    # Find indices of integer redshift values
    #iz_plot = np.flatnonzero(qlf.redshift % 1 == 0)
    iz_plot = np.array([np.abs(qlf.z - (w+1)).argmin() for w in range(7)])
    z_label = [f"$\\mathit{{z}} = {qlf.z[iz]:.1f}$" for iz in iz_plot]
    num_z = len(iz_plot)
    z_col = [mcolors.LinearSegmentedColormap.from_list("roygbiv_extended",
            ["red", "orange", "yellow", "green", "blue", "indigo", "violet"])(i / (num_z - 1)) for i in range(num_z)]

    # Plot QLF for each redshift
    for i, iz in enumerate(iz_plot):
        plt.plot(qlf.log_luminosity, np.log10(qlf.phi[:, iz]), color=z_col[i], alpha=0.5, label=z_label[i])

    # Formatting
    plt.xlabel(r"log L [erg/s]")
    plt.ylabel(r"log $\phi$")
    plt.legend()
    plt.show()


# Plot QLF vs. luminosity at given redshifts (Shen+2020 Fig. 5)
def plot_shen2020_fig5(qlf):
    # Find indices of redshift  to match Shen+2020
    z = np.arange(6)+1

    # Fit A
    qlfA = cp.deepcopy(qlf)
    izv = [np.abs(qlfA.z - zi).argmin() for zi in z]
    z_label = [f"$\\mathit{{z}} \\sim {qlfA.z[iz]:.1f}$" for iz in izv]
    _,_,logLA,_ = qlfA.load_model_parameters(qlfA.z[izv])
    logLA = np.log10(10**logLA * const.L_sun.to(u.erg / u.s).value)
    logPA = np.array([])
    for i in range(len(z)): logPA = np.append(logPA,np.log10(qlfA.get_phi(log_lum=logLA[i],z=z[i])))
    
    # Fit B
    qlfB = cp.deepcopy(qlf)
    qlfB.fit_type = 'FitB'
    qlfB.fit_param = qlfB.load_fit_parameters()
    qlfB.phi = qlfB.calculate_phi(qlfB.log_luminosity, qlfB.z)
    _,_,logLB,_ = qlfA.load_model_parameters(qlfA.z[izv])
    logLB = np.log10(10**logLB * const.L_sun.to(u.erg / u.s).value)
    logPB = np.array([])
    for i in range(len(z)): logPB = np.append(logPB,np.log10(qlfB.get_phi(log_lum=logLB[i],z=z[i])))
    
    # Plot QLF for each redshift
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(8, 8))
    for i in range(3):  
        for j in range(2):
            ind = 2*i + j
            ax = axes[i, j]
            # Plot Fit A  
            ax.axvline(x=logLA[ind], color='goldenrod', linewidth=2)
            ax.axhline(y=logPA[ind], color='goldenrod', linewidth=2)
            ax.plot(qlfA.log_luminosity, np.log10(qlfA.phi[:,izv[ind]]), color='purple', linewidth=4, label='Fit A')
            # Plot Fit B
            #ax.axvline(x=logLB[ind], color='yellow', linewidth=2)
            #ax.axhline(y=logPB[ind], color='yellow', linewidth=2)
            ax.plot(qlfB.log_luminosity, np.log10(qlfB.phi[:,izv[ind]]), color='pink', linewidth=4, linestyle='--', label='Fit B')
            # Formatting
            ax.set_title(z_label[ind])
            ax.set_xticks(np.arange(42,51,1))
            ax.set_xticks(np.arange(42,51,0.2), minor=True)
            ax.set_xlim(42.5, 50.3)
            ax.set_xlabel(r"log $\mathit{L}$ [erg s$^{-1}$]")
            ax.set_yticks(np.arange(-12,-1,2))
            ax.set_yticks(np.arange(-12,-1,0.5), minor=True)
            ax.set_ylim(-11.25, -2.25)
            ax.set_ylabel(r"log $\phi$ [cMpc$^{-1}$]")
            if ind == 0: ax.legend(loc='lower left')

    plt.tight_layout()
    plt.show()


# Plot QLF vs. redshift at given luminosities (Shen+2020 Fig. 8)
def plot_shen2020_fig8(qlf):
    # Find indices of log luminosity [46,47,48] to match Shen+2020
    lum = np.array([46.0, 47.0, 48.0])
    ilumv = np.abs(qlf.log_luminosity[:,np.newaxis] - lum).argmin(axis=0)
    lum_label = [f"$\\mathit{{L}} = {qlf.log_luminosity[ilum]:.1f}$" for ilum in ilumv]

    # Plot QLF for each luminosity
    plt.figure()
    for i, ilum in enumerate(ilumv):
        plt.plot(qlf.z, np.log10(qlf.phi[ilum, :]), alpha=0.5, label=lum_label[i])

    # Formatting
    plt.xticks(np.arange(0,8,1))
    plt.xticks(np.arange(0,8,0.2), minor=True)
    plt.xlim(0,7)
    plt.xlabel(r"log $\mathit{L}$ [erg s$^{-1}$]")
    plt.yticks(np.arange(-12,-2,1))
    plt.yticks(np.arange(-12,-2,0.2), minor=True)
    plt.ylim(-10.3,-3.5)
    plt.ylabel(r"log $\phi$ [cMpc$^{-1}$]")
    plt.legend()
    plt.show()


# Plot redshift sampling at specified luminosities
def plot_zsamp_at_lum(qlf):
    # Sample redshifts at each luminosity
    #zsamp_lo = qlf.sample_redshift(qlf.get_phi(log_lum=42.0),qlf.redshift,qlf.dvolume,ndraw=10000)
    zsamp_md = qlf.sample_redshift(qlf.get_phi(log_lum=44.0),qlf.z,qlf.dvolume,ndraw=10000)
    zsamp_hi = qlf.sample_redshift(qlf.get_phi(log_lum=46.0),qlf.z,qlf.dvolume,ndraw=10000)
    zsamp_xt = qlf.sample_redshift(qlf.get_phi(log_lum=48.0),qlf.z,qlf.dvolume,ndraw=10000)

    # Plot sampling histograms
    #plt.hist(zsamp_lo,bins=28,alpha=0.5,label=r"$\log \mathit{L} \sim 42.0$")
    plt.hist(zsamp_md,bins=28,alpha=0.5,label=r"$\log \mathit{L} \sim 44.0$")
    plt.hist(zsamp_hi,bins=28,alpha=0.5,label=r"$\log \mathit{L} \sim 46.0$")
    plt.hist(zsamp_xt,bins=28,alpha=0.5,label=r"$\log \mathit{L} \sim 48.0$")

    # Formatting
    plt.title(qlf.fit_type)
    plt.legend()
    plt.xlabel(r"$\mathit{z}$")
    plt.ylabel(r"Frequency")
    plt.show()


# Plot restricted luminosity comparison of marginalized QLF (number density dn) vs. redshift 
def plot_dn_restricted_lum(qlf):
    # Prepare figure
    fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(12,4))
    # Plot full L range
    ind = [0,-1]
    axes[0].plot(qlf.z,np.log10(qlf.phi_z),label = f"$\\mathit{{L}} = $ {qlf.log_luminosity[ind].astype(int).tolist()}")
    # Formatting
    axes[0].legend()
    axes[0].set_xlabel(r"$\mathit{z}$")
    axes[0].set_ylabel(r"d$\mathit{n}$ [cMPC$^{-3}$]")
    # Plot restricted L range
    qlfR = cp.deepcopy(qlf)
    l_range = (43,48)
    qlfR.phi_z = qlfR.marginalize_luminosity(qlfR.phi,lum_range=l_range)
    axes[1].plot(qlfR.z,np.log10(qlfR.phi_z),label = f'$\\mathit{{L_{{restricted}}}}$ = {list(l_range)}')
    # Formatting
    axes[1].legend()
    axes[1].set_xlabel(r"$\mathit{z}$")
    axes[1].set_ylabel(r"d$\mathit{n}$ [cMPC$^{-3}$]")
    # Fit Type
    fig.suptitle(qlf.fit_type)
    plt.show()


# Plot restricted luminosity comparison of marginalized full QLF (number density dn) vs. redshift 
def plot_dn_post_restricted_lum(qlf):
    # Prepare figure
    fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(12,4))
    # Plot full L range
    ind = [0,-1]
    axes[0].plot(qlf.z,np.log10(qlf.phi_z),label = f"$\\mathit{{L}} = $ {qlf.log_luminosity[ind].astype(int).tolist()}")
    # Formatting
    axes[0].legend()
    axes[0].set_xlabel(r"$\mathit{z}$")
    axes[0].set_ylabel(r"d$\mathit{n}$ [cMPC$^{-3}$]")
    # Plot restricted L range
    qlfR = cp.deepcopy(qlf)
    l_range = (43,48)
    qlfR.phi_z = qlfR.marginalize_luminosity(qlfR.phi,lum_range=l_range)
    axes[1].plot(qlfR.z,np.log10(qlfR.phi_z),label = f'$\\mathit{{L_{{restricted}}}}$ = {list(l_range)}')
    # Formatting
    axes[1].legend()
    axes[1].set_xlabel(r"$\mathit{z}$")
    axes[1].set_ylabel(r"d$\mathit{n}$ [cMPC$^{-3}$]")
    # Fit Type
    fig.suptitle(qlf.fit_type)
    plt.show()


# Plot L,z heatmap
def plot_luminosity_redshift_heatmap(qlf, zsamp, lsamp):
    # Compute 2D histogram (counts in each z-L bin)
    zbins = np.arange(qlf.z[0],qlf.z[-1]+1,1)
    lbins = np.arange(qlf.log_luminosity[0],qlf.log_luminosity[-1]+1,1)#42+0.5*np.arange(6.5/.5)
    counts, zedges, ledges = np.histogram2d(zsamp, lsamp, bins=[zbins, lbins])

    # Convert bin edges to bin centers for visualization
    z_centers = 0.5 * (zedges[:-1] + zedges[1:])
    l_centers = 0.5 * (ledges[:-1] + ledges[1:])

    # Create a Seaborn heatmap
    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        counts.T,
        cmap="viridis",
        xticklabels=z_centers,
        yticklabels=l_centers,
        cbar_kws={'label': 'Number of Sources'}
    ).invert_yaxis()
    #ax = sns.heatmap(counts.T, cmap="inferno", xticklabels=np.round(z_centers, 2), yticklabels=np.round(l_centers, 2), cbar_kws={'label': 'Number of Sources'}).invert_yaxis()

    # Label the axes
    #plt.gca().invert_yaxis()
    plt.xlabel("Redshift (z)")
    plt.ylabel("Luminosity (L)")
    plt.title("Luminosity vs. Redshift Heatmap")
    plt.show()


# Plot Source SED and template components
def plot_sed(source, phot=None):
    plt.figure(figsize=(10, 6))

    # Plot galaxy components
    if 'galaxy' in source.components:
        galaxy = source.components['galaxy']
        color = dict(zip(galaxy.sed.keys(), ['orange', 'green', 'blue']))
        linesty = dict(zip(galaxy.sed.keys(), ['--',':','-.']))
        for template, template_sed in galaxy.sed.items():
            plt.plot(source.wave, template_sed*source.freq,
                     label=template.upper(), color=color[template], linestyle=linesty[template])
    # Plot AGN components
    if 'agn' in source.components:
        agn = source.components['agn']
        color = dict(zip(agn.sed.keys(), ['magenta', 'darkorchid']))
        linesty = dict(zip(agn.sed.keys(), ['solid',':']))
        for template, template_sed in agn.sed.items():
            plt.plot(source.wave, template_sed*source.freq,
                     label=template.upper(), color=color[template], linestyle=linesty[template])
    # Plot total SED
    plt.plot(source.wave, source.sed*source.freq, color='black', alpha=0.5)
    # Plot data points
    plt.scatter(bp.get_wavelengths(), np.array([flux for flux in source.flux.values()])*bp.get_frequencies(),
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


def plot_kcorr_sed(source):
    colors = ['violet','dodgerblue','green','yellow','orange','red'] 
    shaded_regions = []
    for i, filt in enumerate(bp.filters):
        start_idx = (np.where(bp.filters[filt].thru > 1e-3))[0][0]
        end_idx = (np.where(bp.filters[filt].thru > 1e-3))[0][-1]
        shaded_regions.append((bp.filters[filt].wave[start_idx], bp.filters[filt].wave[end_idx], colors[i % len(colors)]))

    linesty = {'galaxy':'-.','agn':':'}
    label = {'galaxy':'Galaxy','agn':'AGN'}
    plt.plot(source.wave,source.sed,color='red',linestyle='-',label='Observed')
    for component in source.components:
        source_component = source.components[component]
        component_sed = np.sum(np.stack([template_sed for template_sed in source_component.sed.values()]), axis=0)
        plt.plot(source.wave,component_sed,color='red',linestyle=linesty[component],alpha=0.5,label=label[component])
    plt.plot(source.wave/(1+source.z),source.sed*(1+source.z),color='blue',linestyle='--',label='Restframe')
    for component in source.components:
        source_component = source.components[component]
        component_sed = np.sum(np.stack([template_sed for template_sed in source_component.sed.values()]), axis=0)
        plt.plot(source.wave/(1+source.z),component_sed*(1+source.z),color='blue',linestyle=linesty[component],alpha=0.5)
    for x_start, x_end, color in shaded_regions: plt.axvspan(x_start, x_end, color=color, alpha=0.3)
    plt.plot(bp.get_wavelengths(),np.array(list(source.flux.values())),'o',color='black',mfc='none')
    plt.plot(bp.get_wavelengths(),np.array(list(source.rest_flux.values())),'o',color='black',mfc='none')
    plt.text(0.1, 0.9, fr'$z = {source.z:.3f}$', fontsize=12, transform=plt.gca().transAxes)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Wavelength [m]')
    plt.ylabel(r'$F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
    plt.legend()
    plt.show()


def plot_flux_limits(lf):
    plt.plot(lf.z,lf.flux_limits['single']['r'],linestyle='-',color='C0')
    plt.plot(lf.z,lf.flux_limits['coadd']['r'],linestyle='--',color='C1')
    plt.text(0.2, 0.9, 'Single visit', fontsize=12, color='C0', transform=plt.gca().transAxes)
    plt.text(0.2, 0.56, 'Coadd', fontsize=12, color='C1', transform=plt.gca().transAxes)
    plt.xlabel('$z$')
    plt.ylabel(r'$F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
    plt.yscale('log')
    plt.show()


# Apparent mag:  plot_color_color(data.filter(regex='^mag_[^_]+$'),xcol='u-z',ycol='g-i')
def plot_color_color(mags, label='mag_', xcol='u-g', ycol='i-z', det=None, outline=True):
    x, y = xcol.split('-'), ycol.split('-')
    xx, yy = [f"{label}{band}" for band in x], [f"{label}{band}" for band in y]
    # Confirm single color-color (length 2)
    if len(xx) != 2 or len(yy) != 2:
        raise ValueError("Colors must be of form 'u-g'.")
    # Data and axes
    x_arr, y_arr = mags[xx[0]] - mags[xx[1]], mags[yy[0]] - mags[yy[1]]  
    x_lab, y_lab = r'$'+xcol+r'$', r'$'+ycol+r'$'

    # Filter on detections if provided
    if det is not None:
        x_det = det[xx[0]] & det[xx[1]]
        y_det = det[yy[0]] & det[yy[1]]
        s_det = x_det & y_det
    else:
        s_det = [True] * len(mags)
        
    # Create the plot
    #plt.figure(figsize=(10, 6))
    fig, ax = plt.subplots()
    if det is not None:
        ax.scatter(x_arr[~s_det], y_arr[~s_det], c='black', marker='x', alpha=0.5, zorder=3)
    ax.scatter(x_arr[s_det], y_arr[s_det], marker='o', alpha=0.5, zorder=3)

    # Outline color-color region from Richards+2006 Fig. 5 (SDSS DR3 quasar/simulated quasars)
    if outline == True:
        if xcol == 'u-g' and ycol == 'i-z':
            x1,y1 = -0.4,-0.4
            x2,y2 = 1.2,1.0
            inset_image_path = 'data/input/R06_Fig5_u-g.png'
        elif xcol == 'u-z' and ycol == 'g-i':
            x1,y1 = -1.2,-0.8
            x2,y2 = 2.5,1.4
            inset_image_path = 'data/input/R06_Fig5_u-z.png'

        if 'x1' in locals():
            wd, ht = x2-x1, y2-y1
            rectangle = patches.Rectangle((x1, y1), wd, ht, linewidth=2, edgecolor='red', facecolor='none')
            ax.add_patch(rectangle)
            ax.text(x1, y2, 'R06 (Fig. 5)', ha='left', va='bottom',
                                            color=rectangle.get_edgecolor(), fontweight='bold', 
                                            zorder=1)

        try:
            image = mpimg.imread(inset_image_path)
            imagebox = OffsetImage(image, zoom=0.3)

            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            xra, yra = xlim[1]-xlim[0], ylim[1]-ylim[0]
            ax.set_xlim(xlim[0]-xra*0.1,xlim[1]+xra*0.6)
            ax.set_ylim(ylim[0]-yra*0.1,ylim[1]+yra*0.6)
            #xpos = xlim[1] - 0.05*xra
            #ypos = ylim[1] - 0.05*yra
            ab = AnnotationBbox(imagebox, (0.98,0.98), xycoords='axes fraction', 
                                                       box_alignment=(1,1), frameon=False, 
                                                       zorder=2)
            ax.add_artist(ab)
        except Exception as e:
            print(f"Could not load inset image: {e}")

    # Formatting
    if np.diff(ax.get_xlim())[0] < 6:
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    else:
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    if np.diff(plt.gca().get_ylim())[0] < 6:
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    else:
        ax.yaxis.set_major_locator(MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    ax.set_box_aspect(1)
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)

    plt.show()

def plot_color_mag(mags, xlabel='_abs_mag', ylabel='_mag', xmag='i', ycol='g-i', xra=None, yra=None, det=None, outline=False):
    x = xmag.split('-')
    y = ycol.split('-')
    xx = [f"{band}{xlabel}" for band in x]
    yy = [f"{band}{ylabel}" for band in y]
    # Confirm single color-mag (length 2/length 1)
    if len(xx) != 1 or len(yy) != 2:
        raise ValueError("Colors must be of form 'g-r'. Magntitudes must be of form 'r'.")

    # Assign xy arrays and axes
    x_arr = mags[xx]
    y_arr = mags[yy[0]] - mags[yy[1]] 
    if '_abs' in x_arr.columns[0]:
        x_lab = r'$M_{'+xmag+r'}$'
    else: x_lab = xmag
    y_lab = r'$'+ycol+r'$'

    # Filter on detections if provided
    if det is not None:
        x_det = det[xx[0]] & det[xx[1]]
        y_det = det[yy[0]] & det[yy[1]]
        s_det = x_det & y_det
    else:
        s_det = [True] * len(mags)
        
    # Create the plot
    plt.figure(figsize=(10, 6))
    if det is not None: plt.scatter(x_arr[~s_det], y_arr[~s_det], c='black', marker='x', alpha=0.5)
    plt.scatter(x_arr[s_det], y_arr[s_det], marker='o', alpha=0.5)

    # Outline color-color region from Gavazzi+2010 Fig. 3 (SDSS DR7 at z=0 Great Wall)
    if outline == True:
        if xmag == 'i' and ycol == 'g-i':
            x1,y1 = -16.0, 0.0
            x2,y2 = -23.5, 1.5

        if 'x1' in locals():
            wd = x2-x1
            ht = y2-y1
            rectangle = patches.Rectangle((x1, y1), wd, ht, linewidth=2, edgecolor='red', facecolor='none')
            plt.gca().add_patch(rectangle)
            plt.gca().text(0.98, 0.98, 'Gavazzi+2010 (Fig. 3)', ha='right', va='top', transform=plt.gca().transAxes, 
                                                                color=rectangle.get_edgecolor(), fontweight='bold')

    # Formatting
    if xra is not None: plt.xlim(xra[0],xra[1])
    if yra is not None: plt.ylim(yra[0],yra[1])
    if np.diff(plt.gca().get_xlim())[0] < 10:
        plt.gca().xaxis.set_major_locator(MultipleLocator(1.0))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.2))
    else:
        plt.gca().xaxis.set_major_locator(MultipleLocator(5.0))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(1.0))
    if np.diff(plt.gca().get_ylim())[0] < 6:
        plt.gca().yaxis.set_major_locator(MultipleLocator(0.5))
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.1))
    else:
        plt.gca().yaxis.set_major_locator(MultipleLocator(1.0))
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.gca().invert_xaxis()
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.show()


def plot_flux_limits(lf, band='r'):
    sty = ['-','--',':']
    col = ['C0','C1','C2']
    for i, mode in enumerate(lf.flux_limits):
        plt.plot(lf.z,lf.flux_limits[mode][band],label=mode,linestyle=sty[i],color=col[i])
    plt.legend()
    plt.xlabel('$z$')
    plt.ylabel(r'$F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
    plt.yscale('log')
    plt.show()

    # PLOT SINGLE VISIT DETECTIONS
    #plt.figure(figsize=(10, 8))
    #ax = sns.heatmap(hm_depth['single'], annot=False, cmap='viridis', cbar=True, linewidths=.5, 
    #                yticklabels=['47-48', '46-47', '45-46', '44-45', '43-44', '42-43'],
    #                xticklabels=['0', '1', '2', '3', '4', '5', '6', '7'])
    #                #yticklabels=['42-43', '43-44', '44-45', '45-46', '46-47', '47-48'],
    #                #xticklabels=['7', '6', '5', '4', '3', '2', '1', '0'])
    #ax.tick_params(axis='both', which='major', labelsize=12)
    #ax.collections[0].colorbar.ax.tick_params(labelsize=12)
    #plt.title('Detection Fraction (single visit)', fontsize=24)
    #plt.ylabel(r'log $L$ [erg s$^{-1}$]', fontsize=18)
    #plt.xlabel(r'Redshift', fontsize=18)
    #if saveplot == True:
    #    plt.savefig(outdir+'detections_single_visit.png')
    #if showplot == True:
    #    plt.show()


#    # PLOT COADD VISIT DETECTIONS
#    plt.figure(figsize=(10, 8))
#    ax = sns.heatmap(hm_depth['coadd'], annot=False, cmap='viridis', cbar=True, linewidths=.5, 
#                    xticklabels=['42-43', '43-44', '44-45', '45-46', '46-47', '47-48'],
#                    yticklabels=['7', '6', '5', '4', '3', '2', '1'])
#    ax.tick_params(axis='both', which='major', labelsize=12)
#    ax.collections[0].colorbar.ax.tick_params(labelsize=12)
#    plt.title('Detection Fraction (10-yr survey)', fontsize=24)
#    plt.xlabel(r'log $L$ [erg s$^{-1}$]', fontsize=18)
#    plt.ylabel(r'Redshift', fontsize=18)
#    if saveplot == True:
#        plt.savefig(outdir+'detections_survey_coadd.png')
#    if showplot == True:
#        plt.show()


#    # PLOT SINGLE VISIT SELECTED AGN
#    plt.figure(figsize=(10, 8))
#    ax = sns.heatmap(hm_color['single'], annot=False, cmap='viridis', cbar=True, linewidths=.5, 
#                    xticklabels=['42-43', '43-44', '44-45', '45-46', '46-47', '47-48'],
#                    yticklabels=['7', '6', '5', '4', '3', '2', '1'])
#    ax.tick_params(axis='both', which='major', labelsize=12)
#    ax.collections[0].colorbar.ax.tick_params(labelsize=12)
#    plt.title('Selected AGN Fraction (single visit)', fontsize=24)
#    plt.xlabel(r'log $L$ [erg s$^{-1}$]', fontsize=18)
#    plt.ylabel(r'Redshift', fontsize=18)
#    if saveplot == True:
#        plt.savefig(outdir+'selected_single_visit.png')
#    if showplot == True:
#        plt.show()


#    # PLOT COADD VISIT SELECTED AGN
#    plt.figure(figsize=(10, 8))
#    ax = sns.heatmap(hm_color['coadd'], annot=False, cmap='viridis', cbar=True, linewidths=.5, 
#                    xticklabels=['42-43', '43-44', '44-45', '45-46', '46-47', '47-48'],
#                    yticklabels=['7', '6', '5', '4', '3', '2', '1'])
#    ax.tick_params(axis='both', which='major', labelsize=12)
#    ax.collections[0].colorbar.ax.tick_params(labelsize=12)
#    plt.title('Selected AGN Fraction (10-yr survey)', fontsize=24)
#    plt.xlabel(r'log $L$ [erg s$^{-1}$]', fontsize=18)
#    plt.ylabel(r'Redshift', fontsize=18)
#    if saveplot == True:
#        plt.savefig(outdir+'selected_survey_coadd.png')
#    if showplot == True:
#        plt.show()


#    QLF MODEL PARAMETERS
#    plt.figure(figsize=(10, 6))
#    # Plot gamma1 and gamma2
#    plt.plot(qlf.redshift, gamma1, label='gamma1', color='red', alpha=0.7)
#    plt.plot(qlf.redshift, gamma2, label='gamma2', color='orange', alpha=0.7)
#    plt.plot(qlf.redshift, log_ell_star, label='log_ell_star', color='blue', alpha=0.7)
#    plt.plot(qlf.redshift, log_phi_star, label='log_phi_star', color='green', alpha=0.7)
#    # Add labels and legend
#    plt.xlabel('Redshift')
#    plt.ylabel('Values')
#    plt.legend()
#    # Show the plot
#    plt.show()


#    plt.plot(olf.z,olf.flux_limits['single']['r'],linestyle='-',color='C0')
#    plt.plot(olf.z,olf.flux_limits['coadd']['r'],linestyle='--',color='C1')
#    plt.text(0.2, 0.9, 'Single visit', fontsize=12, color='C0', transform=plt.gca().transAxes)
#    plt.text(0.2, 0.56, 'Coadd', fontsize=12, color='C1', transform=plt.gca().transAxes)
#    plt.xlabel('$z$')
#    plt.ylabel(r'$F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
#    plt.yscale('log')
#    plt.show()
