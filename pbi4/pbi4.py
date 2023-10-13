import argparse
import json
from json import JSONEncoder
import numpy as np
import pandas as pd
import os
from math import floor
from scipy.fft import fftn,ifftn

import Utils.readgadget as readgadget
import Utils.readfof as readfof
from powerI4 import fcomb,reduce_array,assign_double_grid,assign_single_grid,\
                                compute_q2k,compute_q4k,measure_pk,measure_pk_multipoles
from bispectrumI4 import measure_bk,measure_bk_multipoles,measure_bk_cross,count_triangles

class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)
    
class IO:
    '''
    This class organises input parameters, output and reading files.
    '''
    def __init__(self, config_file = None, verbose = False, scriptmode = False):

        self.verbose = verbose
        if scriptmode: self.args = vars(self.parse_arguments())
        if config_file is None:
            self.config_file = self.args['config_file']
        else:
            self.config_file = config_file
        self.f = open(self.config_file)
        self.params = json.load(self.f)
        if scriptmode: self.params.update(dict((k,self.args[k]) for k in self.args.keys() if self.args[k] is not None))
        self.params.update(self.derived_params(self.params))
        with open(self.params['config_file'],'w') as config_file:
            json.dump(self.params,config_file,indent = 4, cls=NumpyArrayEncoder)

    def parse_arguments(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--config_file', type = str, help = 'JSON configuration file')
        parser.add_argument('--ifile', type = str, help = 'Input catalog file')
        parser.add_argument('--dfile', type = str, help = 'Input/output density file')
        parser.add_argument('--ps_outfile', type = str, help = 'Power spectrum output file')
        parser.add_argument('--bs_outfile', type = str, help = 'Bispectrum output file')
        parser.add_argument('--cat', type = str, help = 'Type of catalog.\
                                                                    gadget: Gadget particles.\
                                                                    rockstar: Rockstar Halos.\
                                                                    minerva: Minerva halos.\
                                                                    flagship: Flagship halos.\
                                                                    pinocchio: Pinocchio halos.\
                                                                    molino: Molino galaxies.\
                                                                    sancho: Sancho galaxies.')
        parser.add_argument('--box', type = float, help = 'box size of the simulations in mpc/h')
        parser.add_argument('--grid', type = float, help = 'size of the fft grid')
        parser.add_argument('--interpolation', type = int, help = 'interpolation order. \
                                                                        0: direct summation, \
                                                                        1-4: first to fourth order with interlacing')
        parser.add_argument('--interlacing', type = bool, help = 'Whether to apply interlacing')
        parser.add_argument('--kbin', type = float, help = 'bin size in units of the fundamental frequency (linear binning)')
        parser.add_argument('--kcenter', type = float, help = 'centre of the first bin in units of the fundamental frequency')
        parser.add_argument('--iopen', type = bool, help = 'Whether to compute open triangles')
        parser.add_argument('--om', type = float, help = 'matter density')
        parser.add_argument('--ol', type = float, help = 'dark energy density')
        parser.add_argument('--z', type = float, help = 'redshift of the catalog')
        parser.add_argument('--multipoles', type = bool, help = 'Whether to compute multipoles up to L=4 included')
        parser.add_argument('--los', nargs='+', type = float, help = 'Line Of Sight vector. Usage: --los X Y Z, where LOS=[X,Y,Z].')
        parser.add_argument('--cross', type = bool, help = 'Whether to compute cross power spectra/bispectra')
        parser.add_argument('--cfile', type = str, help = 'Stem of the triangle counts file')
        parser.add_argument('--npool', type = int, help = 'Number of threads for parallelized computation of FFTs')
        parser.add_argument('--savedensity', type = bool, help = 'Whether to save density')             

        args = parser.parse_args()

        return args

    def read_params(self):
        if self.verbose:
            print('Input Parameters:')
            for key,value in self.params.items():
                print(key, ' : ', value)
        return self.params
    
    def derived_params(self,params = None):
        '''
        this function computes parameters derived from the ones indicated in the 'config'
        file, such as the fundamentaly frequency.
        ###
        input: params               dict
        output: der_params          dict
        '''
        if params is None:
            params = self.params
        params['los']  = np.array(params['los'])
        der_pars = {'dk' : params['kbin'] * 2. * np.pi / params['box'],
                    'nbins_ps' : int(params['grid']/2./params['kbin']),
                    'nbins_bisp' : floor((params['grid']/3+params['kbin']-params['kcenter'])/params['kbin']),
                    'kf' : 2.*np.pi/params['box'],
                    'kN' : 2*np.pi/params['box']*params['grid']/2.,
                    'Hz' : 100.0*np.sqrt(params['om']*(1.0+params['z'])**3+params['ol']), #THIS IS REALLY H(z)
                    'los' : params['los']/(params['los'] @ params['los'])}
        kvec = np.array([(i-1)*der_pars['dk']+params['kcenter']*2.*np.pi/params['box'] for i in range(1,der_pars['nbins_ps']+1)])
        itot = 0
        if params['iopen']:
            ishift = int(params['kcenter']/params['kbin']) + 1
        else:
            ishift = int(params['kcenter']/params['kbin'])
        for l in range(1,der_pars['nbins_bisp']+1):
            for j in range(1,l+1):
                if params['cross']:
                    imax = min(der_pars['nbins_bisp'],j+l-1+ishift)
                else:
                    imax = j
                for i in range(max(1,l-j+1-ishift),imax+1):
                    itot += 1
        der_pars.update({'itot' : itot, 'kvec' : kvec})
        return der_pars

    def read_catalogues(self,params = None):
        '''
        input: params - parameters [dictionary]
                params['ifile']
                params['los']
        This function reads particles or halos or galaxy catalogs for a list of
        file structures. It also generates RSD displacements along params['los'],
        if required. After storing positions in an array, it checks for particles 
        outside the box and normalizes positions to interval [0,1].
        Formats currently supported:
            - Gadget2 particles (single type)
            - Rockstar halos
            - Minerva halos
            - Flagship halos (no RSD)
            - Pinocchio halos (no RSD)
            - Molino galaxies (no real space)
            - Sancho galaxies (no real space)
        output: npart: number of particles
                pos: array of positions with shape [3,npart]
        '''
        if params is None:
            params = self.params
        if params['cat'] == 'gadget':
            #for each file
            print('Reading Gadget snapshot.\n')
            print('FILE: %s\n' %params['ifile'])
            ptype=1
            npart = readgadget.header(params['ifile']).npart[ptype]
            pos = readgadget.read_field(params['ifile'], "POS ", ptype)     #positions in Mpc/h
            #pos /= 1e3      #for Arepo and Gadget3
            if params['multipoles']:
                print('Displacing along the LOS: (%.1f,%.1f,%.1f)'%(params['los'][0],params['los'][1],params['los'][2]))
                Vel = readgadget.read_field(params['ifile'], "VEL ", ptype)     #peculiar velocities in km/s
                pos[:] += Vel[:]*params['los']/params['Hz']*(1.+params['z'])    #note that readgadget already multiplies by the infamous sqrt(a)  
        if params['cat'] == 'rockstar':         
            print('Reading EOS Rockstar halos.\n')
            print('FILE: %s\n' %params['ifile'])
            mycat = np.load(params['ifile'])
            vec = np.array([mycat['mass'],mycat['pos'][0,:],mycat['pos'][1,:],mycat['pos'][2,:]]).T
            npart = len(vec[:,0])
            if self.verbose: print('Found %d halos.' %npart)
            if self.verbose: print('Number density: %.2e (h/Mpc)^3'%(npart/params['box']**3))
            pos = vec[:,1:4]
            if params['multipoles']:
                print('Displacing along the LOS: (%.1f,%.1f,%.1f)'%(params['los'][0],params['los'][1],params['los'][2]))
                Vel = np.array([mycat['vel'][0,:],mycat['vel'][1,:],mycat['vel'][2,:]]).T
                pos[:] += Vel[:]*params['los']/params['Hz']*(1.+params['z'])    #positions are comoving in Rockstar
        if params['cat'] == 'minerva':                 
            print('Reading Minerva halos.\n')
            print('FILE: %s\n' %params['ifile'])
            dt = np.dtype([('dummy', 'i4'),('pos', 'f4',3),('Vel', 'f4',3),('Mass', 'f4'),('dummy2', 'i4')])
            mycat = np.fromfile(params['ifile'], dtype=dt)
            vec = np.array([mycat['Mass'][:],mycat['pos'][:,0], mycat['pos'][:,1], mycat['pos'][:,2]]).T
            npart = len(vec[:,0])
            if self.verbose: print('Found %d halos.' %npart)
            if self.verbose: print('Number density: %.2e (h/Mpc)^3'%(npart/params['box']**3))
            pos = vec[:,1:4]
            if params['multipoles']:
                print('Displacing along the LOS: (%.1f,%.1f,%.1f)'%(params['los'][0],params['los'][1],params['los'][2]))
                Vel = np.array(mycat['Vel'][:,0], mycat['Vel'][:,1], mycat['Vel'][:,2])
                pos[:] += Vel[:]*params['los']   
        if params['cat'] == 'flagship':
            print('Reading Euclid Flagship galaxies.\n')
            print('FILE: %s\n' %params['ifile'])
            if params['multipoles']:
                print('redshift not implemented yet on %s'%params['cat'])
                raise SystemExit
            vec = pd.read_parquet(params['ifile'], engine='pyarrow', columns=['halo_lm','x','y','z'])
            vec = np.asarray(vec)
            vec[:,0] = 10**(vec[:,0])
            npart = len(vec[:,0])
            if self.verbose: print('Found %d halos.' %npart)
            if self.verbose: print('Number density: %.2e (h/Mpc)^3'%(npart/params['box']**3))
            pos = vec[:,1:4]
        if params['cat'] == 'quijote':
            z_dict  = {0.0:4, 0.5:3, 1.0:2, 2.0:1, 3.0:0} #{redshift:snapnum}
            snapnum = z_dict[params['z']]
            print('Reading Quijote halos.\n')
            print('FILE: %s\n' %params['ifile'])
            print('reading: ', params['ifile'],', redshift: ',params['z'],', snapnum: ',snapnum)
            # read the halo catalogue
            FoF = readfof.FoF_catalog(params['ifile'], snapnum, long_ids=False,
                            swap=False, SFR=False, read_IDs=False)
            # get the properties of the halos
            pos_h = FoF.GroupPos/1e3            #Halo positions in Mpc/h
            mass  = FoF.GroupMass*1e10          #Halo masses in Msun/h
            vec = np.array([mass[:],pos_h[:,0],pos_h[:,1],pos_h[:,2]],dtype=object).T
            npart = len(vec[:,0])
            if self.verbose: print('Found %d halos.' %npart)
            if self.verbose: print('Number density: %.2e (h/Mpc)^3'%(npart/params['box']**3))
            pos = vec[:,1:4]
            if params['multipoles']:
                print('Displacing along the LOS: (%.1f,%.1f,%.1f)'%(params['los'][0],params['los'][1],params['los'][2]))
                Vel = FoF.GroupVel*(1.0+params['z']) #Halo peculiar velocities in km/s
                pos[:] += Vel[:]*params['los']/params['Hz']*(1.+params['z'])    #note that readgadget already multiplies by the infamous sqrt(a)  
        if params['cat'] == 'molino':
            print('Reading Molino galaxies.\n')
            print('FILE: %s\n' %params['ifile'])
            if not params['multipoles']:
                print('real space not implemented yet on %s'%params['cat'])
                raise SystemExit
            #mycat = np.loadtxt(fname,unpack=True, usecols=(0,1,2,3))
            mycat = np.genfromtxt(params['ifile'])[:,[0,1,2,13]]
            vec  = np.array([mycat[:,3],mycat[:,0],mycat[:,1],mycat[:,2]]).T
            npart = len(vec[:,0])
            if self.verbose: print('Found %d galaxies.' %npart)
            if self.verbose: print('Number density: %.2e (h/Mpc)^3'%(npart/params['box']**3))
            pos = vec[:,1:4]
        if params['cat'] == 'sancho':  
            print('Reading Sancho galaxies.\n')
            print('FILE: %s\n' %params['ifile'])
            if not params['multipoles']:
                print('real space not implemented yet on %s'%params['cat'])
                raise SystemExit                      
            mycat = np.load(params['ifile'])
            vec = np.array([mycat['pos'][:,0],mycat['pos'][:,1],mycat['pos'][:,2]]).T
            npart = len(vec[:,0])
            if self.verbose: print('Found %d galaxies.' %npart)
            if self.verbose: print('Number density: %.2e (h/Mpc)^3'%(npart/params['box']**3))
            pos = vec[:,0:3]
        pos = pos.T
        pos = (pos+params['box'])%params['box']                                  # make sure particle are in the params['box']
        pos = np.minimum(1.-np.nextafter(0,1),pos/params['box'])       # normalize to interval [0,1)
        print('Done.')
        return pos
    

class DENSITIES:
    '''
    This class interpolates positions on a grid to compute a density field.
    '''
    def __init__(self, params,verbose = False):

        self.params = params
        self.verbose = verbose

    def read_density(self, params = None):
        '''
        input: params - parameters [dictionary]
                params['dfile']
                params['grid']
        This function reads a previously generated density, checks whether
        params['grid'] is lower or equal to the grid at which the density
        was computed, and if it's lower, it reduces the density array to 
        the desired grid.
        output: npart: number of particles
                dcl: density with shape [(grid, grid, grid)]
        '''
        if params is None:
            params = self.params
        print('FILE: %s\n' %params['dfile'])
        print('Reading density...\n')
        dcl_read = np.load(params['dfile'])
        header = dcl_read['header']
        print('Grid: %.0f' %header[0])
        print('Number of elements: %.0f\n' %header[1])
        if header[0] < params['grid']:
            print('\nERROR: Density file has smaller grid size than requested. \n')
            raise SystemExit
        else:
            dcl = dcl_read['data']
            if header[0] > params['grid']:
                print('Current grid size: %d' %header[0])
                print('Requested grid size: %d' %params['grid'])
                print('reducing density array to internal grid...')
                dcl = reduce_array(header[0],dcl,params['grid']) #np.resize also good but slower
            npart=header[1]
        return npart,dcl
    
    def assign_grid(self,pos, params = None, weight = None):
        if params is None:
            params = self.params
        npart = len(pos.T)
        print('Assigning %d particles on the grid...' %npart)
        if weight is None:
            weight = np.ones(npart)
        if params['interlacing']:
            dcl = assign_double_grid(params['grid'], params['interpolation'], npart, pos, weight)
        else:
            dcl = assign_single_grid(params['grid'], params['interpolation'], npart, pos, weight) 
        print('Done.')
        return npart, dcl
    
    def compute_density(self, npart, dcl, params = None):
        '''
        input: params - parameters [dictionary]
                params['interlacing']
                params['grid']
                params['kf']
                params['kN']
                params['interpolation']
                params['box']
                params['savedensity']
            npart  - number of particles (DM particles, halos, galaxies) [float]
            pos    - position of the particles in Mpc/h                  [3,npart]
            weight - density weight                                      [npart]
            FFT    - whether to Fast Fourtier Transform                  [bool]
        This function interpolates particles on a grid (or double grid if 
        interlacing is true), applying a density weight if provided.
        The interpolation is executed on an external Fortran subroutine
        included in the powerI4 library (and adapted from the powerI4 code 
        by Emiliano Sefusatti)
        If required, it computes the FFT and the fcomb routine and then
        saves the density on python binary
        output: npart: number of particles
                dcl: density with shape [(grid, grid, grid)]
        '''
        if params is None:
            params = self.params
        print('Computing final density.')
        if self.verbose:
            print('FFTs...\n')
            print('Grid = %d'%params['grid'])
            print('Fundamental frequency = %f h/Mpc' %params['kf'])
            print('Nyquist frequency = %f h/Mpc' %params['kN'])
        dcl /= npart        # FFT normalization
        dcl = fftn(dcl)
        if self.verbose: print('Running Fcomb...')
        dcl = fcomb(params['grid'], params['interpolation'], params['interlacing'], params['box'], dcl)   
        print('Done.\n')
        return dcl

    def compute_density_multipoles(self, dcl, params = None):
        '''
        input: params - parameters [dictionary]
                params['grid']
                params['box']
                params['los']
            dcl  - density from interpolation on a grid and FFT with shape [(grid,grid,grid)]
        This function computes the density multipoles 
        as explained in https://arxiv.org/pdf/1506.02729.pdf.
        Note that only a fixed LOS vector is implemented for the moment.
                dcl2, dcl4: density with shape [(grid, grid, grid)]
        '''
        if params is None:
            params = self.params
        if not params['multipoles']:
            print('WARNING: the parameter -multipoles- is set to False, so you should not compute multipoles!')
        deltax = ifftn(dcl)    
        dcl2 = -0.5*dcl
        print('Computing density multipoles...\n')
        for i in range(0,3):
            for j in range(0,i+1):
                fac = 1.
                if i != j: fac=2.   
                Q2x=deltax*params['los'][i]*params['los'][j]
                deltak=fftn(Q2x)
                Q2k=compute_q2k(params['grid'],params['box'],deltak,i+1,j+1)
                dcl2=dcl2+fac*1.5*Q2k
        print('Done with delta_2.')
        dcl4 = -7./8.*dcl-2.5*dcl2
        for i in range(0,3):
            for j in range(0,i+1):                
                for l in range(0,j+1):
                    for k in range(0,l+1): 
                        if i == j and j == l and l == k: fac = 1
                        if i == j and j == l and l > k: fac = 4
                        if i == j and j > l and l == k: fac = 6
                        if i > j and j == l and l == k: fac = 4
                        if i == j and j > l and l > k: fac = 12
                        if i > j and j == l and l > k: fac = 12
                        if i > j and j > l and l == k: fac = 12
                        Q4x=deltax*params['los'][i]*params['los'][j]*params['los'][l]*params['los'][k]
                        deltak=fftn(Q4x)
                        Q4k=compute_q4k(params['grid'],params['box'],deltak,i+1,j+1,l+1,k+1)
                        dcl4=dcl4+fac*(35./8.)*Q4k
        print('Done with delta_4.')
        return dcl2, dcl4

    def save_density(self, dcl, npart, params = None):
        if params is None:
            params = self.params
        print('FILE: %s' %params['dfile'])
        header = np.array([params['grid'], npart, params['kf'], params['kN']])
        data = dcl
        np.savez_compressed(params['dfile'], header = header, data = data)

class MEASUREMENTS:

    def __init__(self,params, verbose = False):

        self.params = params
        self.verbose = verbose

    def powerspectrum(self, npart, dcl, dcl2 = None, dcl4 = None, params = None):
        '''
        input: params - parameters [dictionary]
                params['box']
                params['dk']
                params['kcenter']
                params['nbins_ps']
                params['grid']
                params['kbin']
                params['multipoles']
        but now the measure_pk only works with no RSD although it could work with RSD. And cross PS is implemented.
        and the avgP4 is not implemented correctly in measure_pk.
        '''
        if params is None:
            params = self.params
        output_type = [('k',np.float32),
                       ('avgk',np.float32),
                       ('P0',np.float32),
                       ('P2',np.float32),
                       ('P4',np.float32),
                       ('Nmodes',np.float32),
                       ('PSN',np.float32)]
        size_output = len(params['kvec'])
        output = np.ndarray((size_output),dtype=output_type)
        psn = 1./(npart*(2.*np.pi/params['box'])**3)
        print('kmin = %f h/Mpc'%params['kvec'][0])
        print('kmax = %f h/Mpc'%params['kvec'][-1])
        print('Number of bins = %d'%params['nbins_ps'])
        if params['multipoles']:
            if params['cross']:
                print('ERROR: Cross power spectrum in redshift space not implemented yet!')
                raise SystemExit
            print('Computing redshift space power spectrum')
            pout = measure_pk_multipoles(params['grid'], params['box'], params['kbin'], \
                                         params['kcenter'],params['nbins_ps'],dcl, dcl2,dcl4)
            output['P0'] = pout[1]
            output['P2'] = pout[2]
            output['P4'] = pout[3]
            output['Nmodes'] = pout[4]
        else:
            if params['cross']:
                print('Computing reals space cross power spectrum.')
                if dcl2 is None:
                    print('ERROR: you did not include the second density field')
                    raise SystemExit
            else:
                print('Computing real space power spectrum')
                dcl2 = dcl
            pout = measure_pk(params['grid'], params['box'], params['kbin'], \
                              params['kcenter'], params['nbins_ps'],dcl, dcl2)
            output['P0'] = pout[1]
            output['Nmodes'] = pout[2]
        output['k'] = params['kvec']
        output['avgk'] = pout[0]
        output['PSN'] = psn
        print('Done.\n')
        return output
    
    def triangle_counts(self,params = None):
        if params is None:
            params = self.params
        if os.path.exists(params['cfile']) and os.path.getsize(params['cfile']) > 0:
            print('Found triangle counts file.')
            print('FILE: %s'%params['cfile'])
            counts_read = np.load(params['cfile'])
            counts = counts_read['counts']
        else:
            if params['cross']:
                print('Computing triangle counts for cross bispectra')
            counts = count_triangles(params['nbins_bisp'],params['grid'],params['kbin'],params['kcenter'],\
                                params['iopen'],params['cross'],params['npool'],params['itot']) 
        return counts
    
    def bispectrum(self,npart, counts, dcl, dcl2 = None, dcl4 = None, params = None):
        if params is None:
            params = self.params
        if params['cross']:
            output_type = [('k1',np.float32),
                           ('k2',np.float32),
                           ('k3',np.float32),
                           ('B0',np.float32),
                           ('Ntr',np.float32)]
            print('Computing cross bispectra')
            bout = measure_bk_cross(params['nbins_bisp'],params['grid'],params['kbin'],\
                    params['kcenter'], params['iopen'],params['kf'], \
                    params['npool'],params['itot'], npart,counts,dcl,dcl2)
            size_output = len(bout[0])
            output = np.ndarray((size_output),dtype=output_type)
            output['k1'] = bout[0]
            output['k2'] = bout[1]
            output['k3'] = bout[2]
            output['B0'] = bout[3]
            output['Ntr'] = bout[4]
        elif not params['multipoles']:
            output_type = [('k1',np.float32),
                           ('k2',np.float32),
                           ('k3',np.float32),
                           ('Pk1',np.float32),
                           ('Pk2',np.float32),
                           ('Pk3',np.float32),
                           ('B0+BSN',np.float32),
                           ('BSN',np.float32),
                           ('Ntr',np.float32)]
            print('Computing real space bispectrum')
            bout = measure_bk(params['nbins_bisp'],params['grid'],params['kbin'],\
                    params['kcenter'], params['iopen'],params['kf'], \
                    params['npool'],params['itot'], npart,counts,dcl)
            size_output = len(bout[0])
            output = np.ndarray((size_output),dtype=output_type)
            output['k1'] = bout[0]
            output['k2'] = bout[1]
            output['k3'] = bout[2]
            output['Pk1'] = bout[3]
            output['Pk2'] = bout[4]
            output['Pk3'] = bout[5]
            output['B0+BSN'] = bout[6]
            output['BSN'] = bout[7]
            output['Ntr'] = bout[8] 
        else:
            output_type = [('k1',np.float32),
                           ('k2',np.float32),
                           ('k3',np.float32),
                           ('Pk1',np.float32),
                           ('Pk2',np.float32),
                           ('Pk3',np.float32),
                           ('B0+BSN',np.float32),
                           ('BSN',np.float32),
                           ('Ntr',np.float32),
                           ('B2',np.float32),
                           ('B4',np.float32)]
            print('Computing redshift space bispectrum')
            bout = measure_bk_multipoles(params['nbins_bisp'],params['grid'],params['kbin'],\
                    params['kcenter'], params['iopen'],params['kf'], \
                    params['npool'],params['itot'], npart,counts,dcl,dcl2,dcl4)     
            size_output = len(bout[0])
            output = np.ndarray((size_output),dtype=output_type)
            output['k1'] = bout[0]
            output['k2'] = bout[1]
            output['k3'] = bout[2]
            output['Pk1'] = bout[3]
            output['Pk2'] = bout[4]
            output['Pk3'] = bout[5]
            output['B0+BSN'] = bout[6]
            output['BSN'] = bout[7]
            output['Ntr'] = bout[8] 
            output['B2'] = bout[9]
            output['B4'] = bout[10]  
        print('Done.\n')
        return output 

    def save_powerspectrum(self,pout, params = None):
        if params is None:
            params = self.params
        print('Saving power spectrum...\n')
        print('FILE: %s' %params['ps_outfile'])
        header = 'k     avgk     P0      P2      P4     Nmodes    PSN\n'
        np.savetxt(params['ps_outfile'], np.array([ pout['k'], pout['avgk'], pout['P0'], pout['P2'],pout['P4'],pout['Nmodes'],pout['PSN']]).T, fmt=('%0.8e','%0.8e','%0.8e','%0.8e','%0.8e','%d','%0.8e'), header = header)
        print('Done.')
        return

    def save_triangle_counts(self,counts, params = None):
        if params is None:
            params = self.params
        print('\nSaving triangle counts...')
        print('FILE %s'%params['cfile'])
        np.savez_compressed(params['cfile'], counts = counts)
        print('Done.\n') 
        return

    def save_bispectrum(self,bout, params = None):
        if params is None:
            params = self.params
        print('Saving bispectrum...\n')
        print('FILE: %s' %params['bs_outfile'])
        if params['cross']:
            header = 'k1/kF k2/kF k3/kF      B0          N_tr\n'
            np.savetxt(params['bs_outfile'], np.array([ bout['k1'], bout['k2'],bout['k3'], bout['B0'], bout['Ntr']]).T, \
                       fmt=('%d','%d','%d','%0.8e','%0.8e'), header = header)
        elif params['multipoles'] == 0:
            header = 'k1/kF k2/kF k3/kF P(k1) P(k2)       P(k3)       B0+BSN        BSN    N_tr\n'
            np.savetxt(params['bs_outfile'], np.array([ bout['k1'], bout['k2'],bout['k3'], bout['Pk1'], bout['Pk2'], bout['Pk3'],\
                                            bout['B0+BSN'],bout['BSN'], bout['Ntr']]).T, \
                       fmt=('%d','%d','%d','%0.8e','%0.8e','%0.8e','%0.8e','%0.8e','%0.8e'), header = header)
        else:
            header = 'k1/kF k2/kF k3/kF P(k1) P(k2)       P(k3)       B0+BSN        BSN    N_tr    B2      B4\n'
            np.savetxt(params['bs_outfile'], np.array([ bout['k1'], bout['k2'],bout['k3'], bout['Pk1'], bout['Pk2'], bout['Pk3'],\
                                            bout['B0+BSN'],bout['BSN'], bout['Ntr'], bout['B2'],bout['B4']]).T, \
                       fmt=('%d','%d','%d','%0.8e','%0.8e','%0.8e','%0.8e','%0.8e','%0.8e','%0.8e','%0.8e'), header = header)
        print('Done.\n')
        return
