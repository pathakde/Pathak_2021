import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#import requests
import requests
#import get()
from simulation_data import get

from .galaxy import timeaverage_stellar_formation_rate, median_stellar_age, total_stellar_mass, halfmass_rad_stars, halflight_rad_stars, max_merger_ratio

class GalaxyPopulation():
    
    
    def __init__(self):
        self.ids = []
        self.mass_min = 0
        self.mass_max = 0
        self.redshift = 0
        self.snap = 0
        
        
    #select ids for a given redshift and mass-cut
    def select_galaxies(self, redshift, mass_min, mass_max=12):
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
                d2 = f.create_dataset('median_age', data = self.get_median_stellar_age())
                d3 = f.create_dataset('halfmass_radius', data = self.get_halfmass_rad_stars())
                d4 = f.create_dataset('total_mass', data = self.get_total_stellar_mass())
                d5 = f.create_dataset('halflight_radius_U', data = self.get_halflight_rad_stars(band='U'))
                d6 = f.create_dataset('halflight_radius_V', data = self.get_halflight_rad_stars(band='V'))
                d7 = f.create_dataset('halflight_radius_I', data = self.get_halflight_rad_stars(band='I'))
                d8 = f.create_dataset('newbin_current_SFR', data = self.get_timeaverage_stellar_formation_rate(timescale=0, binwidth=0.01))
                d9 = f.create_dataset('maximum_merger_ratio_30kpc_current_fraction', data = self.get_max_merger_ratio(scale=30))

                
        with h5py.File('galaxy_population_data_'+str(self.redshift)+'.hdf5', 'r') as f:
            ids = f['ids'][:]
            median_age = f['median_age'][:]
            halfmass_radius = f['halfmass_radius'][:]
            total_mass = f['total_mass'][:]
            halflight_radius_U = f['halflight_radius_U'][:]
            halflight_radius_V = f['halflight_radius_V'][:]
            halflight_radius_I = f['halflight_radius_I'][:]
            newbin_current_SFR = f['newbin_current_SFR'][:]
            maximum_merger_ratio_30kpc_current_fraction = f['maximum_merger_ratio_30kpc_current_fraction'][:]

            
        galaxy_population_data = {
                                    'ids': ids,
                                    'median_age': median_age,
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
        ids = self.ids
        time_averages = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            time_averages[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = calc_timescale, start=calc_start, binwidth=calc_binwidth)
        #save file
        np.savetxt( 'z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', time_averages)
        time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', dtype=float)
        return time_avg_SFT
    
        
    def get_timeaverage_stellar_formation_rate(self, timescale, start = 0, binwidth=0.05):
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr')
        if file.exists ():
            time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr', dtype=float) 
            return time_avg_SFT
        else:
            return self.calc_timeaverage_stellar_formation_rate(calc_timescale=timescale, calc_start=start, calc_binwidth=binwidth)
            
    

    #median stellar age
    def calc_median_stellar_age(self):
        #create and populate array for mean SFT
        ids = self.ids
        MedianSFT = np.zeros(len(ids))
        for i, id in enumerate(ids):
            MedianSFT[i] = median_stellar_age(redshift = self.redshift, id = id)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_Median_SFT', MedianSFT)
        median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Median_SFT', dtype=float)
        return median_SFT
    
    
    def get_median_stellar_age(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_Median_SFT')
        if file.exists ():
            median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Median_SFT', dtype=float) 
            return median_SFT
        else:
            return self.calc_median_stellar_age()
    

        
        #total stellar mass
    def calc_total_stellar_mass(self):
        ids = self.ids
        total_mass = np.zeros(len(ids))
        for i, id in enumerate(ids):
            total_mass[i] = total_stellar_mass(id=id, redshift=self.redshift)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_total_mass', total_mass)
        total_mass = np.loadtxt('z='+ str(self.redshift) +'_total_mass', dtype=float)
        return total_mass
    
    
    def get_total_stellar_mass(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_total_mass')
        if file.exists ():
            total_mass = np.loadtxt('z='+ str(self.redshift) +'_total_mass', dtype=float)
            return total_mass
        else:
            return self.calc_total_stellar_mass()
        
        
        
        #stellar half mass radius
    def calc_halfmass_rad_stars(self):
        ids = self.ids
        halfmass_rad = np.zeros(len(ids))
        for i, id in enumerate(ids):
            halfmass_rad[i] = halfmass_rad_stars(id=id, redshift=self.redshift)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_halfmass_rad', halfmass_rad)
        halfmass_rad = np.loadtxt('z='+ str(self.redshift) +'_halfmass_rad', dtype=float)
        return halfmass_rad
    
    
    def get_halfmass_rad_stars(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_halfmass_rad')
        if file.exists ():
            halfmass_rad = np.loadtxt('z='+ str(self.redshift) +'_halfmass_rad', dtype=float) 
            return halfmass_rad
        else:
            return self.calc_halfmass_rad_stars()
        
        

         #halflight_rad_stars
    def calc_halflight_rad_stars(self, calc_band, calc_bound):
        ids = self.ids
        halflight_rad = np.zeros(len(ids))
        for i, id in enumerate(ids):
            halflight_rad[i] = halflight_rad_stars(id=id, redshift=self.redshift, band=calc_band, bound=calc_bound)
        #save file
        np.savetxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad'+str(calc_bound), halflight_rad)
        halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad'+str(calc_bound), dtype=float)
        return halflight_rad
    
    
    
    def get_halflight_rad_stars(self, band, bound):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +str(band)+'_halflight_rad'+str(bound))
        if file.exists ():
            halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(band)+'_halflight_rad'+str(bound), dtype=float) 
            return halflight_rad
        else:
            return self.calc_halflight_rad_stars(calc_band = band, calc_bound = bound)
        
    
        
        #maximum merger ration of current fraction 
    def calc_max_merger_ratio(self, calc_scale):
        ids = self.ids
        mass_ratio = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            mass_ratio[i] = max_merger_ratio(id=id, redshift=self.redshift, scale=calc_scale)
        #save file
        np.savetxt('z='+str(self.redshift)+'_max_merger_ratio_'+str(calc_scale), mass_ratio)
        mass_ratio = np.loadtxt('z='+str(self.redshift)+'_max_merger_ratio_'+str(calc_scale), dtype=float)
        return mass_ratio
    
    
    def get_max_merger_ratio(self, scale=30): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_max_merger_ratio_'+str(scale))
        if file.exists ():
            mass_ratio = np.loadtxt('z='+str(self.redshift)+'_max_merger_ratio_'+str(scale), dtype=float)
            return mass_ratio
        else:
            return self.calc_max_merger_ratio(calc_scale=scale)
        