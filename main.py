# main.py
#import data.input.parameters as params
from scripts.common import np, const, pd, plt, sns, u
from scripts.global_variables import cosmo, qlf, bp, sed, select
from scripts.simulated_sources import Source, process_source #Galaxy, ActiveGalacticNuclei, process_source
from scripts.luminosity_functions import ObservedLuminosityFunction
#from scripts.selection_criteria import AGNSelection
from scripts.plotting import *

import scripts.tests as test
import copy
import pdb

def main():
    #print(f"Using cosmology with H0={params.H0} and Om0={params.Om0}")

    # R06 full sample: 15343
    # R06 Table 2 low-z: 12164
    # R06 TAble 2 i>=19.475: 1191
    #     full low-z sample: 13355

    sim_data = False
    if sim_data:
        data = []
        while_loop = False
        if while_loop:
            #nbreak = 15343
            nbreak = 12164
            ngood = 0
            counter = 0
            while ngood < nbreak:
                counter += 1
                lsamp, zsamp = qlf.sample_lz_pairs(1)
                src = Source(zsamp[0],agn_luminosity=lsamp[0],components='agn')
                # If (brighter than <) –OR– (dimmer than <)
                #R06 Table 2 (i<19.475) or (i>15.475)
                if (src.mag['i'] < 19.1):# and (src.mag['i'] > 15.0):
                    ngood += 1
                else:
                    continue
                data.append(process_source(src))
                if (ngood+1) % max(nbreak//10,1) == 0: print(f"{round(100*(ngood+1)/nbreak)}% complete...")
                #if ngood % 100 == 0: print(f"No. of sources: {ngood}")
        else:
            ndraw = 1000#100000
            lsamp, zsamp = qlf.sample_lz_pairs(ndraw)
            for i in range(ndraw):
                src = Source(zsamp[i],agn_luminosity=lsamp[i],components=None)
                data.append(process_source(src))
                if (i+1) % max(ndraw//10,1) == 0: print(f"{round(100*(i+1)/ndraw)}% complete...")

        data = pd.DataFrame(data)
        #match_R06_15343
        file_name = input("Please enter save file name: ")
        if file_name != "": data.to_pickle("data/output/"+file_name+".pkl")
    else:
        file_name = input("Please enter save file name: ")
        data = pd.read_pickle("data/output/"+file_name+".pkl")

        #eff_area = 2094
        #eff_area = 1622
        #nloz = 13355
        #nloz = 12164
        #temp_data = data.sample(n=nloz)
        #temp_data = data

        #plot_color_color(temp_data.filter(regex='^mag_[^_]+$'))#,xcol='u-z',ycol='g-i')
        #plot_color_color(sub_data[sub_data.abs_mag_i < -22.5].filter(regex='^mag_[^_]+$'))#,xcol='u-z',ycol='g-i')
        # Then i can do this number counts. Also restricted in apparent mag? Not sure.
        olf = ObservedLuminosityFunction(qlf.log_luminosity,qlf.z,qlf.dluminosity,qlf.volume,qlf.dvolume,bp.flux_limits)

        live_test = False
        if live_test:
            data = []
            src = Source(1.5,agn_luminosity=45.0)
            data.append(process_source(src))
            data = pd.DataFrame(data)

            dat = data.iloc[0]
            dat_sed = []                 
            for comp in dat.components:
                for template in sed.components[comp]:
                    key = f"{comp}_{template}_norm"
                    norm = dat[key]
                    sed_template = sed.components[comp][template]
                    dat_sed.append(norm * sed_template)#/(1+dat.z))
            dat_sed = np.sum(np.stack(dat_sed),axis=0)
            dat_wave = sed.domain['wavelength']['wave']

        
        olf.rebuild_source_sed(sed,data.iloc[1])
        #olf.differential_number_counts(temp_data.mag_i,temp_data.abs_mag_i,temp_data.z,temp_data.detected_single_i, eff_area)
        #olf.cumulative_number_counts(temp_data.mag_i,temp_data.abs_mag_i,temp_data.z,temp_data.detected_single_i, eff_area)

    pdb.set_trace()

    # data.to_pickle("data/output/data_snapshot_R06_12165.pkl")
    # data = pd.read_pickle("data/output/data_snapshot_R06_12165.pkl")

    log_lum = np.log10(data.luminosity_agn)
    z = data.z

    H,xedges,yedges = np.histogram2d(log_lum,z,bins=(olf.ledges,olf.zedges))
    plt.figure()
    plt.pcolormesh(xedges, yedges, H.T, shading='auto')  # Transpose H to match x-y orientation
    plt.xlabel('log(L)')
    plt.ylabel('z')
    plt.colorbar(label='Counts')
    plt.show()
#    hm_depths = {}
#    for mode in bp.get_depths():
#        cols = [filter + '_detected_' + mode for filter in bp.get_filters()]
#        hm_depths[mode] = data[cols]
#    pdb.set_trace()
#    cols = [filter + '_detected_single' for filter in bp.get_filters()]
#    hm_single = data[cols]
#    
#    # calculate fractional detections per visit mode
#    hm_depth = {mode: np.zeros((len(z), len(lrange))) for mode in det} # heatmap of passes depths
#    for i in range(len(z)):
#        for j in range(len(lrange)):
#            # Boolean mask for the current luminosity bin
#            luminosity_mask = (lagnv[i, :] >= lrange[0] + j) & (lagnv[i, :] < lrange[0] + j + 1)
#            # Count the number of Trues in detsel for the specified range
#            det_ct = {mode: np.sum(luminosity_mask & det[mode][i, :]) for mode in det}
#            # Count the number of elements in the specified range
#            num_ct = np.sum(luminosity_mask)
#            # Avoid division by zero
#            if num_ct > 0:
#                for mode in hm_depth: hm_depth[mode][i, j] = det_ct[mode] / num_ct
#            else:
#                for mode in hm_depth: hm_depth[mode][i, j] = 0
#    for mode in hm_depth: hm_depth[mode] = np.flip(np.transpose(hm_depth[mode]), axis=0)#np.flip(hm_depth[mode])   
#
#    # detections which also pass color cuts
#    detsel = {mode: np.logical_and(det[mode], sel) for mode in det_ugrizy}
#    #detsel = np.logical_and(det, sel)
#
#    # calculate fractional detections per visit mode
#    hm_color = {mode: np.zeros((len(z), len(lrange))) for mode in det} # heatmap of passes depths
#    for i in range(len(z)):
#        for j in range(len(lrange)):
#            # Boolean mask for the current luminosity bin
#            luminosity_mask = (lagnv[i, :] >= lrange[0] + j) & (lagnv[i, :] < lrange[0] + j + 1)
#            # Count the number of Trues in detsel for the specified range
#            det_ct = {mode: np.sum(luminosity_mask & detsel[mode][i, :]) for mode in det}
#            # Count the number of elements in the specified range
#            num_ct = np.sum(luminosity_mask)
#            # Avoid division by zero
#            if num_ct > 0:
#                for mode in hm_color: hm_color[mode][i, j] = det_ct[mode] / num_ct
#            else:
#                for mode in hm_color: hm_color[mode][i, j] = 0
#    for mode in hm_color: hm_color[mode] = np.flip(hm_color[mode], axis=0)
#
#
#

if __name__ == "__main__":
    main()
