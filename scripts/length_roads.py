"""
@author: Elco Koks
@date: Jan, 2018
"""

import pandas as pd
import os
from functions import country_road_length,create_poly_files
from multiprocess import Pool , cpu_count 
import time
import numpy as np
import urllib.request

def road_length(base_path,multiprocess=True,overwrite=True):
    
    """
    Main function to estimate the length of all the roads we are interested in. 
    
    Args:
        
        base_path: base path to location of all files and directories in this project.
        
        multiprocess: set to True by default. Set to False in the case of limited processing power.
        
        overwrite: set to True by default. This relates to all input data (i.e. .poly files, .osm.pbf files and shapefiles).
        
    Returns:
        
        an Excel file with the length of all primary, secondary, tertiary, track and other roads for each country.
    
    """
    
    print ('The calculation of road lenghts has started!')

# =============================================================================
#     """ Set path to dirs"""
# =============================================================================
    dir_out = os.path.join(base_path,'output_data')
    poly_dir = os.path.join(base_path,'poly_files')
    osm_path_in = os.path.join(base_path,'osm_continent')

# =============================================================================
#     """ create directories if they are not created yet """
# =============================================================================
    if not os.path.exists(dir_out):
         os.makedirs(dir_out)

    if not os.path.exists(poly_dir):
         os.makedirs(poly_dir)           

    if not os.path.exists(osm_path_in):
         os.makedirs(osm_path_in)     
# =============================================================================
#     """Set path to files we use """
# =============================================================================
    wb_country_in = os.path.join(base_path,'input_data','wbccodes2014.csv')
    global_shape = os.path.join(base_path,'input_data','2015_GAUL_Dataset_Mod.gdb')

# =============================================================================
#     """Load country shapes and list and only save the required countries"""
# =============================================================================
    wb_country = pd.read_csv(wb_country_in,header=0,index_col=0)
    
    #filter high income countries from country file
    country_list = wb_country[['country','continent']].loc[wb_country['wbregion']!='YHI']

    # add column to country list so we can easily look up the required continental
    # osm file for that continent    
    map_continent = {'MA': 'central-america','SA': 'south-america','EU': 'europe','AS': 'asia',
                     'AU': 'australia-oceania','AF':'africa','AM':'north-america'}
   
    country_list['osm-cont'] = country_list['continent'].map(lambda x: (map_continent[x])) 
    
# =============================================================================
#     """ create .poly files to clip countries from osm.pbf files """
# =============================================================================
    if not os.listdir(poly_dir):
        create_poly_files(base_path,global_shape,save_shapefile=overwrite)


# =============================================================================
# """ check if we have actually downloaded the openstreetmap input files. If not,
# lets download them. Note: this will take a while! """ 
# =============================================================================
    continent_list = ['central-america','south-america','europe','asia','australia-oceania','africa','north-america']

    for continent in continent_list:        
        url = 'http://download.geofabrik.de/%s-latest.osm.pbf' % continent
        if '%s-latest.osm.pbf' % (continent) not in  os.listdir(osm_path_in):
            urllib.request.urlretrieve(url, osm_path_in)
            
# =============================================================================
#     """ create extracted osm files for each country per continent """
# =============================================================================
    out = []
    countries = []
    continent_osms = []
    base_paths = []
    for country in country_list.iterrows():
        country = country[1]
        continent_osm = os.path.join(osm_path_in,'%s-latest.osm.pbf' % (country['osm-cont']))
        countries.append(country['country'])
        continent_osms.append(continent_osm)
        base_paths.append(base_path)

    # multiprocessing will start if set to True. Set to False with limited processing capacities    
    if multiprocess==True:
        pool = Pool(cpu_count()-1)
        out = pool.starmap(country_road_length, zip(countries,continent_osms,base_paths,overwrite=overwrite)) 
    
    # when multiprocessing set to False, we will just loop over the countries.
    else:
        out = []
        for country in country_list.iterrows():
            country = country[1]
            continent_osm = os.path.join(osm_path_in,'%s-latest.osm.pbf' % (country['osm-cont']))
            out.append(country_road_length(country['country'],continent_osm,base_path,overwrite=overwrite))
            
    df = pd.concat(out,axis=1).T
    
    map_country = dict(zip(wb_country['country'],wb_country['country_name']))
    df['Country'] = df.index.to_series().map(map_country)

    df.set_index('Country',inplace=True,drop=True)
    
    writer = pd.ExcelWriter(os.path.join(dir_out,'dist_roads.xlsx'))
    df.to_excel(writer,'output')
    writer.save()
    
    end = time.time()

    print('It took ' + str(np.float16((end - start))) + " seconds to finish!")         

if __name__ == "__main__":
    start = time.time()

    base_path =   os.path.join(os.path.dirname(__file__),'..')

    road_length(base_path)