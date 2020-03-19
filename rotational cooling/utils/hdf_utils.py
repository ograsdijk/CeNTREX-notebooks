import h5py

def load_measurement_data(fname, dset_name):
    # loading data from HDF file
    with h5py.File(fname, 'r') as f:
        laserlock = f[dset_name]['readout']['Laser Lock'][()]
        wavelength = f[dset_name]['readout']['Wavelength'][()]
        mirror = f[dset_name]['readout']['ZaberTMM'][()]
        pxie_time = {}
        pxie = {}
        for dset in f[dset_name]['readout']['PXIe-5171']:
            if 'events' not in dset:
                pxie[int(dset.split('_')[-1])] = f[dset_name]['readout']['PXIe-5171'][dset][()]
                pxie_time[int(dset.split('_')[-1])] = f[dset_name]['readout']['PXIe-5171'][dset].attrs['ch0 : timestamp']

    return pxie, pxie_time, laserlock, wavelength, mirror