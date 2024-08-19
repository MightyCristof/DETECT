# main.py
#import data.input.parameters as params
from scripts.common import np, const, plt, u
from scripts.global_variables import sed, bp
from scripts.simulated_sources import Source, Galaxy, ActiveGalacticNuclei
from scripts.selection_criteria import AGNSelection
from scripts.plotting import plot_sed, plot_cc
from scripts.global_variables import cosmology, qlf

import scripts.tests as test
import copy
import pdb

def main():
    #print(f"Using cosmology with H0={params.H0} and Om0={params.Om0}")
    
    nobj = 10000
    z = np.arange(1., 8.0, 1.0)
    det = np.zeros((len(bp.filters), len(z), nobj), dtype=bool)  
    sel = np.zeros((len(z), nobj), dtype=bool)

    selection = AGNSelection(criteria='Analytic')

    for i in range(len(z)):
        lagn = qlf.sample_luminosity(z[i], ndraw=nobj)
        #lgal = np.zeros(ndraw)
        #plt.hist(lagn, bins='auto')
        #plt.show()

        for j in range(nobj):
            src = Source(z[i], agn_luminosity=lagn[j])
            #lgal[i] = np.log10(src.galaxy.luminosity)
            det[:,i,j] = np.array(list(src.detect.values()))
            sel[i,j] = selection.is_agn(src.photometry['mag'])


    # passes survey depths
    det_ugriz_z = np.sum(det, axis=2) / nobj  # n obj detected in each filter across redshifts
    det_any_z = np.sum(np.any(det,axis=0),axis=1) / nobj # n obj detected in any band across redshfit
    # passes color cuts
    sel_z = np.sum(sel, axis=1) / nobj
    # total objects selected
    tot_z = np.sum(np.any(det,axis=0) & sel, axis=1) / nobj

    print(f"N objects: {nobj}")
    for i in range(len(z)):
        print(f"z:{z[i]:.1f}, det(ugriz):{det_ugriz_z[:,i]}, det(any):{det_any_z[i]:.3f}, col:{sel_z[i]:.3f}, tot:{tot_z[i]:.3f}")
        

if __name__ == "__main__":
    main()
