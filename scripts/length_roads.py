"""
@author: Elco Koks
@date: Jan, 2018
"""

import geopandas as gpd
import pandas as pd
import os
from functions import country_road_length,get_country,create_poly_files
from multiprocess import Pool , cpu_count 
import time
import numpy as np

def road_length(base_path):
    print ('The calculation of road lenghts has started!')
# =============================================================================
#     """ Set path to dirs"""
# =============================================================================

    dir_out = os.path.join(base_path,'output_data')
    poly_dir = os.path.join(base_path,'poly_files')

# =============================================================================
#     """Set path to files we use """
# =============================================================================
    wb_country_in = os.path.join(base_path,'input_data','wbccodes2014.csv')
    global_shape = os.path.join(base_path,'WB_SHP_WGS84','2015_GAUL_Dataset_Mod.gdb')
    osm_path_in = os.path.join(base_path,'osm_continent')
    
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
        create_poly_files(base_path,global_shape,save_shapefile=True)

    """ create extracted osm files for each country per continent """
    
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

    pool = Pool(cpu_count()-1)
    out = pool.starmap(country_road_length, zip(countries,continent_osms,base_paths)) 

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
