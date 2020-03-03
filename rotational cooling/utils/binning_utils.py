import math
import numpy as np
from tqdm.notebook import tqdm
from .hdf_utils import load_measurement_data

def find_nearest_idx(array,value):
    """
    Get index of nearest item in array
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def combine_arrays_irregular(indices1, data1, indices2, data2):
    """
    Combine two arrays where items in either array correspond to items of a different
    index in the other.
    """
    combined = []
    for idx1, idx2 in zip(indices1, indices2):
        try:
            combined.append((data1[idx1], data2[idx2]))
        except IndexError as e:
            print(idx1, idx2)
            raise e
    return zip(*combined)

def bin_data(data_x, data_y, bins = None, width = 1):
    """
    Bin data
    """
    # check if bins are already specified, otherwise take bins from min() to max()
    # with a spacing of 1
    if isinstance(bins, type(None)):
        bins = np.arange(np.floor(np.nanmin(data_x)), np.ceil(np.nanmax(data_x))+1, width)

    # for each item in data_x get the corresponding bin index
    indices = np.digitize(data_x, bins)-1

    # put the binned data in a dictionary with the bin as index
    binned_data = {b:[] for b in bins}
    for y, idx in zip(data_y, indices):
        binned_data[bins[idx]].append(y)

    return bins, binned_data

def average_binned_data(bins, binned_data):
    """
    Computes the mean of a bin and the standard error of the means
    """
    binned_data_averaged = {}
    for b in bins:
        # only divide standard deviation by sqrt(# entries in a bin) if larger than 1
        if len(binned_data[b]) > 1:
            binned_data_averaged[b] = (np.mean(binned_data[b]), np.std(binned_data[b]/np.sqrt(len(binned_data[b]))))
        elif len(binned_data[b]) == 1:
            binned_data_averaged[b] = (np.mean(binned_data[b]), np.mean(binned_data[b]))
        else:
            binned_data_averaged[b] = (np.nan, np.nan)
    return binned_data_averaged

def bin_datasets_laserscan(fname, dset_names, descriptions, ch, laser_scanned = 'laser 1'):
    binned_integrals_dsets = {}

    for idx_dset, (dset_name, desc) in enumerate(tqdm((zip(dset_names, descriptions)), total = len(descriptions))):
        # loading data from hdf
        pxie, pxie_time, laserlock, wavelength, mirror = load_measurement_data(fname, dset_name)

        # getting corresponding indices for pmt traces and laserlock frequency data by grabbing indices with the closest
        # timestamps
        idx_pairs = {idx_pxie: find_nearest_idx(laserlock['time'], t_pxie) for idx_pxie, t_pxie in pxie_time.items()}

        # integrating the pmt traces
        integrals = [-np.trapz(pxie[idx][250:1900, ch]-np.mean(pxie[idx][50:200, ch])) for idx in range(1, max(pxie_time.keys())+1)]

        # combining the pmt and laserlock data into two lists, with each item corresponding to the closest timestamped entry for
        # each respective dataset
        freqs, integrals = combine_arrays_irregular(idx_pairs.values(), laserlock[f'{laser_scanned} frequency'], np.array(list(idx_pairs.keys()))-1, integrals)

        # bin the integrals into frequency bins of 1 MHz
        bins, binned_integrals = bin_data(freqs, integrals, width = 1)

        # average the binned integrals
        binned_integrals_averaged = average_binned_data(bins, binned_integrals)
        binned_integrals_dsets[desc] = binned_integrals_averaged

    return binned_integrals_dsets

def bin_dataset_switching(fname, dset_name, pmt_channel, switching_channel, threshold = 500):
    # loading data from HDF file
    pxie, pxie_time, laserlock, wavelength, mirror = load_measurement_data(fname, dset_name)

    switch_state = {idx_pxie:traces[:,switching_channel].mean() > threshold for idx_pxie, traces in pxie.items()}

    bins, bin_state_data = bin_data(list(switch_state.values()), [traces[:,pmt_channel] for traces in pxie.values()])
    bin_state_data = {key: [-(value - value[50:200].mean()) for value in values] for key, values in bin_state_data.items()}

    data_tmp = np.zeros(len(pxie[1]))
    switch_open = []
    switch_closed = []
    state = switch_state[1]
    cnt = 0
    for idx in range(1, max(pxie.keys())+1):
        if switch_state[idx] == state:
            data_tmp += pxie[idx][:, pmt_channel] - pxie[idx][50:200, pmt_channel].mean()
            cnt += 1
        else:
            if state:
                switch_open.append(data_tmp/cnt)
            else:
                switch_closed.append(data_tmp/cnt)
            state = switch_state[idx]
            cnt = 0
            data_tmp = pxie[idx][:, pmt_channel] - pxie[idx][50:200, pmt_channel].mean()

    # integrated averaged data with trapz
    integrated_switch_open = [np.trapz(-data[250:1750]) for data in switch_open]
    integrated_switch_closed = [np.trapz(-data[250:1750]) for data in switch_closed]

    # calculating the ratio between switching states
    ratio = [i_open/i_close for i_open, i_close in zip(integrated_switch_open, integrated_switch_closed)]
    return switch_state, bin_state_data, ratio
