import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#import requests
import requests
#import get()
from simulation_data import get

from .galaxy import timeaverage_stellar_formation_rate, mean_stellar_age, median_stellar_age, total_stellar_mass, halfmass_rad_stars, halflight_rad_stars, max_merger_ratio

class GalaxyPopulation():
    
    
    def __init__(self):
        self.ids = []
        self.mass_min = 0
        self.mass_max = 0
        self.redshift = 0
        self.snap = 0
        
        
    #select ids for a given redshift and mass-cut
    def select_galaxies(self, redshift, mass_min, mass_max=12):
        '''
        input params: 
            id: the simulation id of target galaxy: integer (specific to simulation, pre-check) 
            redshift: redshift of target galaxy: numerical value (default==2, specific to simulation, pre-check at https://www.tng-project.org/data/)
            mass_min: lower end of mass cut
                    [units: log10 solar mass]
            mass_max: upper end of mass cut, default 12
                    [units: log10 solar mass]
        preconditions: 
            requires get() imported from simulation_data.__init__
        output params:
            ids: array of simulation ids of galaxies in selection
                    [units: none] 
    '''
        if self.ids == [] or (self.mass_min != mass_min or self.mass_max != mass_max or self.redshift != redshift):
            h = 0.6774
            mass_minimum = 10**mass_min / 1e10 * h
            mass_maximum = 10**mass_max / 1e10 * h
            # form the search_query string by hand for once
            search_query = "?mass_stars__gt=" + str(mass_minimum) + "&mass_stars__lt=" + str(mass_maximum)
            url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + search_query
            subhalos = get(url, {'limit':5000})
            self.mass_min = mass_min
            self.mass_max = mass_max
            self.redshift = redshift
            self.ids = [ subhalos['results'][i]['id'] for i in range(subhalos['count'])]
            self.ids = np.array(self.ids, dtype=np.int32)
        return self.ids
 


    def get_galaxy_population_data(self):
        '''
        input params: 
            [none]
        preconditions: 
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
            requires self.get_median_stellar_age()
            requires self.get_halfmass_rad_stars()
            requires self.get_total_stellar_mass()
            requires self.get_halflight_rad_stars(band)
            requires self.get_timeaverage_stellar_formation_rate(timescale, binwidth)
            requires self.get_max_merger_ratio(scale)
        output params:
            checks if galaxy population file exists at target redshift. 
            if the galaxy population file does not exist, processes population data and saves relevant halo properties
            returns a dictionary with halo properties
                    output dictionary keys: 'ids': array of simulation ids of galaxies in selection
                                                    [units: none] 
                                            'median_age' : median stellar ages of galaxies in selection 
                                                    [units: Lookback time in Gyr]
                                            'halfmass_radius' : half-mass radii of galaxies in selection
                                                    [units: physical kpc]
                                            'total_mass': total stellar masses of galaxies in selection
                                                    [units: solar mass] 
                                            'halflight_radius_U' : rest-frame U-band half-light radii of galaxies in selection
                                                    [units: physical kpc]
                                            'halflight_radius_V' : rest-frame V-band half-light radii of galaxies in selection
                                                    [units: physical kpc]
                                            'halflight_radius_I' : rest-frame I-band half-light radii of galaxies in selection
                                                    [units: physical kpc]
                                            'newbin_current_SFR': current SFR averaged over last 0.01 Gyr of galaxies in selection
                                                    [units: solar mass/year]
                                            'maximum_merger_ratio_30kpc_current_fraction': greatest fraction of current stellar mass within 30 kpc that can be traced back to a single merger event of galaxies in selection
                                                    [units: none, fraction of current stellar mass within 30 kpc from galaxy CM]
    '''
        redshift = self.redshift
        galaxy_population_data = {}
        import h5py
        from pathlib import Path
        if Path('galaxy_population_data_'+str(self.redshift)+'.hdf5').is_file():
            pass
        else:
            with h5py.File('galaxy_population_data_'+str(self.redshift)+'.hdf5', 'a') as f:
                #writing data
                d1 = f.create_dataset('ids', data = self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12))
                d2 = f.create_dataset('mean_age', data = self.get_mean_stellar_age(weights=False))
                d3 = f.create_dataset('halfmass_radius', data = self.get_halfmass_rad_stars())
                d4 = f.create_dataset('total_mass', data = self.get_total_stellar_mass())
                d5 = f.create_dataset('halflight_radius_U', data = self.get_halflight_rad_stars(band='U', bound=0.5))
                d6 = f.create_dataset('halflight_radius_V', data = self.get_halflight_rad_stars(band='V', bound=0.5))
                d7 = f.create_dataset('halflight_radius_I', data = self.get_halflight_rad_stars(band='I', bound=0.5))
                d8 = f.create_dataset('newbin_current_SFR', data = self.get_timeaverage_stellar_formation_rate(timescale=0, binwidth=0.01))
                d9 = f.create_dataset('maximum_merger_ratio_30kpc_current_fraction', data = self.get_max_merger_ratio(scale=30))
                d10 = f.create_dataset('mean_age_massweighted', data = self.get_mean_stellar_age(weights=True))
                
        with h5py.File('galaxy_population_data_'+str(self.redshift)+'.hdf5', 'r') as f:
            ids = f['ids'][:]
            mean_age = f['mean_age'][:]
            mean_age_massweighted = f['mean_age_massweighted'][:]
            halfmass_radius = f['halfmass_radius'][:]
            total_mass = f['total_mass'][:]
            halflight_radius_U = f['halflight_radius_U'][:]
            halflight_radius_V = f['halflight_radius_V'][:]
            halflight_radius_I = f['halflight_radius_I'][:]
            newbin_current_SFR = f['newbin_current_SFR'][:]
            maximum_merger_ratio_30kpc_current_fraction = f['maximum_merger_ratio_30kpc_current_fraction'][:]

        galaxy_population_data = {
                                    'ids': ids,
                                    'mean_age': mean_age,
                                    'mean_age_massweighted': mean_age_massweighted,
                                    'halfmass_radius': halfmass_radius,
                                    'total_mass': total_mass,
                                    'halflight_radius_U': halflight_radius_U,
                                    'halflight_radius_V': halflight_radius_V,
                                    'halflight_radius_I': halflight_radius_I,
                                    'newbin_current_SFR': newbin_current_SFR,
                                    'maximum_merger_ratio_30kpc_current_fraction': maximum_merger_ratio_30kpc_current_fraction,
                                 }
        return galaxy_population_data

    
    #time avg SFR
    def calc_timeaverage_stellar_formation_rate(self, calc_timescale, calc_start=0, calc_binwidth=0.05):
        '''
        input params: 
            calc_timescale: length of time window for over which average SFR is calculated
                    [units: Gyr]
            calc_start: minimum lookback to which timescale is is added to get time window for calculating average SFR (default 0)
                    [units: Lookback time in Gyr]
            calc_binwidth: width of linear age bin for computing SFR (default 0.05 Gyr)
                    [units: Gyr]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires galaxy.timeaverage_stellar_formation_rate(id, redshift, timescale, start=0, binwidth=0.05)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
        output params:
            average SFR: an array of average star formation rates of galaxies in selection over a specified timescale 
                    [units: solar mass/year] 
        '''
        ids = self.ids
        time_averages = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            time_averages[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = calc_timescale, start=calc_start, binwidth=calc_binwidth)
        np.savetxt( 'z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', time_averages)
        time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', dtype=float)
        return time_avg_SFT
    
        
    def get_timeaverage_stellar_formation_rate(self, timescale, start = 0, binwidth=0.05):
        '''
        input params: 
            timescale: length of time window for over which average SFR is calculated
                    [units: Gyr]
            start: minimum lookback to which timescale is is added to get time window for calculating average SFR (default 0)
                    [units: Lookback time in Gyr]
            binwidth: width of linear age bin for computing SFR (default 0.05 Gyr)
                    [units: Gyr]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires self.calc_timeaverage_stellar_formation_rate(calc_timescale, calc_start=0, calc_binwidth=0.05)
        output params:
            checks if array of average SFR exists in temporary text file.
                if temporary text file does not exist, calculates average SFR using self.calc_timeaverage_stellar_formation_rate(calc_timescale, calc_start, calc_binwidth)
                if temporary file exists, reads array from temporary file
            returns average SFR: an array of average star formation rates of galaxies in selection over a specified timescale 
                    [units: solar mass/year] 
        '''
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr')
        if file.exists ():
            time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr', dtype=float) 
            return time_avg_SFT
        else:
            return self.calc_timeaverage_stellar_formation_rate(calc_timescale=timescale, calc_start=start, calc_binwidth=binwidth)
            
    
    
    #mean stellar age
    def calc_mean_stellar_age(self, calc_weights=False):
        '''
        input params: 
            calc_weights: When True: returns mass-weighted mean stellar age
                          When False: returns mean stellar age (default) 
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires galaxy.mean_stellar_age(id, redshift)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
        output params:
            mean stellar age: an array of mean stellar ages of galaxies in selection (mass-weighted mean ages if weights=True) 
                    [units: Lookback time in Gyr] 
        '''
        ids = self.ids
        means = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            means[i] = mean_stellar_age(redshift = self.redshift, id = id, weights=True)
        if calc_weights==False:
            filename='z='+str(self.redshift)+'_Mean_SFT'
        else:
            filename='z='+str(self.redshift)+'_Mean_massweighted_SFT'
        #save file
        np.savetxt(filename, means)
        mean_SFT = np.loadtxt(filename, dtype=float)
        return mean_SFT
    
    
    def get_mean_stellar_age(self, weights=False): 
        import pathlib
        if weights==False:
            filename='z='+str(self.redshift)+'_Mean_SFT'
        else:
            filename='z='+str(self.redshift)+'_Mean_massweighted_SFT'
        file = pathlib.Path(filename)
        if file.exists ():
            MeanSFT = np.loadtxt(filename, dtype=float) #rename pre-existing files before parameterizing further
            return MeanSFT
        else:
            return self.calc_mean_stellar_age(calc_weights=weights)
        
        

    #median stellar age
    def calc_median_stellar_age(self):
        '''
        input params: 
            [none]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires galaxy.median_stellar_age(id, redshift)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
        output params:
            median stellar age: an array of median stellar ages of galaxies in selection 
                    [units: Lookback time in Gyr] 
        '''
        ids = self.ids
        MedianSFT = np.zeros(len(ids))
        for i, id in enumerate(ids):
            MedianSFT[i] = median_stellar_age(redshift = self.redshift, id = id)
        np.savetxt('z='+ str(self.redshift) +'_Median_SFT', MedianSFT)
        median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Median_SFT', dtype=float)
        return median_SFT
    
    
    def get_median_stellar_age(self):
        '''
        input params: 
            [none]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires self.calc_median_stellar_age()
        output params:
            checks if array of median stellar age exists in temporary text file.
                if temporary text file does not exist, calculates median stellar age using self.calc_median_stellar_age()
                if temporary file exists, reads array from temporary file
            returns median stellar age: an array of median stellar ages of galaxies in selection 
                    [units: Lookback time in Gyr] 
        '''
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_Median_SFT')
        if file.exists ():
            median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Median_SFT', dtype=float) 
            return median_SFT
        else:
            return self.calc_median_stellar_age()
    

        
    #total stellar mass
    def calc_total_stellar_mass(self):
        '''
        input params: 
            [none]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires galaxy.total_stellar_mass(id, redshift)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
        output params:
            total stellar mass: an array of total stellar masses of galaxies in selection 
                    [units: log10 solar masses]
        '''
        ids = self.ids
        total_mass = np.zeros(len(ids))
        for i, id in enumerate(ids):
            total_mass[i] = total_stellar_mass(id=id, redshift=self.redshift)
        np.savetxt('z='+ str(self.redshift) +'_total_mass', total_mass)
        total_mass = np.loadtxt('z='+ str(self.redshift) +'_total_mass', dtype=float)
        return total_mass
    
    
    def get_total_stellar_mass(self):
        '''
        input params: 
            [none]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires self.calc_total_stellar_mass()
        output params:
            checks if array of total stellar mass exists in temporary text file.
                if temporary text file does not exist, calculates total stellar mass using self.calc_total_stellar_mass()
                if temporary file exists, reads array from temporary file
            returns total stellar mass: an array of total stellar masses of galaxies in selection 
                    [units: log10 solar masses]
        '''
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_total_mass')
        if file.exists ():
            total_mass = np.loadtxt('z='+ str(self.redshift) +'_total_mass', dtype=float)
            return total_mass
        else:
            return self.calc_total_stellar_mass()
        
        
        
    #half mass radius
    def calc_halfmass_rad_stars(self):
        '''
        input params: 
            [none]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires galaxy.halfmass_rad_stars(id, redshift)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
        output params:
            half mass radius: an array of stellar half-mass radii of galaxies in selection 
                    [units: physical kpc]
        '''
        ids = self.ids
        halfmass_rad = np.zeros(len(ids))
        for i, id in enumerate(ids):
            halfmass_rad[i] = halfmass_rad_stars(id=id, redshift=self.redshift)
        np.savetxt('z='+ str(self.redshift) +'_halfmass_rad', halfmass_rad)
        halfmass_rad = np.loadtxt('z='+ str(self.redshift) +'_halfmass_rad', dtype=float)
        return halfmass_rad
    
    
    def get_halfmass_rad_stars(self):
        '''
        input params: 
            [none]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires self.calc_halfmass_rad_stars()
        output params:
            checks if array of half-mass radii exists in temporary text file.
                if temporary text file does not exist, calculates half-mass radius using self.calc_halfmass_rad_stars()
                if temporary file exists, reads array from temporary file
            returns half-mass radius: an array of stellar half-mass radii of galaxies in selection 
                    [units: physical kpc]
        '''
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_halfmass_rad')
        if file.exists ():
            halfmass_rad = np.loadtxt('z='+ str(self.redshift) +'_halfmass_rad', dtype=float) 
            return halfmass_rad
        else:
            return self.calc_halfmass_rad_stars()
        
        

         #half light radius
    def calc_halflight_rad_stars(self, calc_band, calc_bound=0.5):
        '''
        input params: 
            calc_band: choice of photometric band or mass to calculate effective size in: string
                    'U': (Vega magnitude) 
                    'V': (Vega magnitude) 
                    'I': (AB magnitude)
                    'M': (solar masses)
            calc_bound: target fraction of quantity (light intensity, mass) enclosed to calculate radius (default 0.5: for half-light radius)
                    [range (0, 1]]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires galaxy.halflight_rad_stars(id, redshift, band, bound=0.5)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
        output params:
            half light radius: an array of half-light radii of galaxies in selection 
                    [units: physical kpc]
        '''
        ids = self.ids
        halflight_rad = np.zeros(len(ids))
        for i, id in enumerate(ids):
            halflight_rad[i] = halflight_rad_stars(id=id, redshift=self.redshift, band=calc_band, bound=calc_bound)
        np.savetxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad'+str(calc_bound), halflight_rad)
        halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad'+str(calc_bound), dtype=float)
        return halflight_rad
    
    
    
    def get_halflight_rad_stars(self, band, bound=0.5):
        '''
        input params: 
            band: choice of photometric band or mass to calculate effective size in: string
                    'U': (Vega magnitude) 
                    'V': (Vega magnitude) 
                    'I': (AB magnitude)
                    'M': (solar masses)
            bound: target fraction of quantity (light intensity, mass) enclosed to calculate radius (default 0.5: for half-light radius)
                    [range (0, 1]]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires self.calc_halflight_rad_stars()
        output params:
            checks if array of half-light radii exists in temporary text file.
                if temporary text file does not exist, calculates half-light radius using self.calc_halfmass_rad_stars()
                if temporary file exists, reads array from temporary file
            returns half-light radius: an array of stellar half-light radii of galaxies in selection 
                    [units: physical kpc]
        '''
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +str(band)+'_halflight_rad'+str(bound))
        if file.exists ():
            halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(band)+'_halflight_rad'+str(bound), dtype=float) 
            return halflight_rad
        else:
            return self.calc_halflight_rad_stars(calc_band=band, calc_bound=bound)
        
    
        
        #maximum merger ratio of current fraction 
    def calc_max_merger_ratio(self, calc_scale=30):
        '''
        input params: 
            calc_scale: the radial distance up to which star partickes should be considered as part of the target galaxy (default 30 kpc)
                    [units: physical kpc]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires galaxy.max_merger_ratio(id, redshift=2, scale=30)
            requires output from galaxy.get_galaxy_particle_data(id, redshift, populate_dict=True): halo file must exist
            requires output from galaxy.get_stellar_assembly_data(id, redshift, populate_dict=True): assembly data must exist
        output params:
            greatest merger contribution: the largest fraction of current stellar mass within $scale [kpc] that can be traced back to a single merger event
                    [units: none, fraction of total stellar mass within $scale [kpc] from galaxy CM]
        '''
        ids = self.ids
        mass_ratio = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            mass_ratio[i] = max_merger_ratio(id=id, redshift=self.redshift, scale=calc_scale)
        np.savetxt('z='+str(self.redshift)+'_max_merger_ratio_'+str(calc_scale), mass_ratio)
        mass_ratio = np.loadtxt('z='+str(self.redshift)+'_max_merger_ratio_'+str(calc_scale), dtype=float)
        return mass_ratio
    
    
    def get_max_merger_ratio(self, scale=30): 
        '''
                input params: 
            scale: the radial distance up to which star partickes should be considered as part of the target galaxy (default 30 kpc)
                    [units: physical kpc]
        preconditions:
            requires initialization with self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12)
            requires self.calc_max_merger_ratio()
        output params:
            checks if array of maximum merger fractions exists in temporary text file.
                if temporary text file does not exist, calculates maximum merger current fraction using self.calc_max_merger_ratio()
                if temporary file exists, reads array from temporary file
            returns greatest merger contribution: an array of the largest fraction of current stellar mass within $scale [kpc] that can be traced back to a single merger event
                    [units: none, fraction of total stellar mass within $scale [kpc] from galaxy CM]
        '''
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_max_merger_ratio_'+str(scale))
        if file.exists ():
            mass_ratio = np.loadtxt('z='+str(self.redshift)+'_max_merger_ratio_'+str(scale), dtype=float)
            return mass_ratio
        else:
            return self.calc_max_merger_ratio(calc_scale=scale)
        