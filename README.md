# Spectral-Element solvers
Spectral-element solvers for solving the 1D and 2D wave-equation, and 1D and 2D Helmholtz equations.

Data is written in HDF5 format and can be loaded from Python as

```python

def loadAttrFromH5(path_data):
        """ Load attributes from simulation data
            https://www.pythonforthelab.com/blog/how-to-use-hdf5-files-in-python/
        """

        with h5py.File(path_data, 'r') as f:
            settings_dict = {}
            
            settings_dict['dt'] = f['/pressures'].attrs['dt'][0]
            settings_dict['dx'] = f['/pressures'].attrs['dx'][0]
            settings_dict['c'] = f['/pressures'].attrs['c'][0]
            settings_dict['c_phys'] = f['/pressures'].attrs['c_phys'][0]
            settings_dict['rho'] = f['/pressures'].attrs['rho'][0]
            settings_dict['tmax'] = f['/pressures'].attrs['tmax'][0]
            settings_dict['sigma0'] = f['/pressures'].attrs['sigma0'][0]
            settings_dict['fmax'] = f['/pressures'].attrs['fmax'][0]
            settings_dict['tmax'] = f['/pressures'].attrs['tmax'][0]

            settings_dict['domain_minmax'] = f['/mesh'].attrs['domain_minmax']
            settings_dict['boundary_type'] = f['/mesh'].attrs['boundary_type'][0]

        return settings_dict

def loadDataFromH5(path_data, tmax=None):
    """ input
            tmax: normalized max time 
        output
            grid: is a x X y x t dimensional array
    """

    with h5py.File(path_data, 'r') as f:        

        mesh = np.asarray(f['mesh'][()])
        umesh = np.asarray(f['umesh'][()])
        ushape = f['umesh'].attrs['umesh_shape'][0]

        p = np.asarray(f['pressures'][()])
        up = np.asarray(f['upressures'][()])

        conn = np.asarray(f['conn'][()]).astype(int)
        x0_srcs = np.asarray(f['x0_srcs'][()]) if 'x0_srcs' in f else np.asarray([])

        dt = f['pressures'].attrs['dt'][0]
        t = np.asarray(f['t'][()])

        if tmax == None:
            ilast = len(t) - 1
        else:
            ilist = [i for i, n in enumerate(t) if abs(n - tmax) <= dt/2]
            if not ilist:
                raise Exception('tmax exceeds simulation data running time')
            ilast = ilist[0]            
            # crop w.r.t. time
            t = np.array(t[:ilast+1])

        if len(p.shape) == 2: # data is written differently in Matlab and C++/Python
            p = np.array([p[:ilast+1,:]])
        else:
            p = np.array(p[:,:ilast+1,:])

        assert p.shape[2] == mesh.shape[0]

        accumulators = []
        data = SimulationData(mesh, umesh, ushape, p, up, t, conn, x0_srcs, accumulators)

    return data
```
