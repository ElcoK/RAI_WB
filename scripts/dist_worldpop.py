# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 14:54:28 2017

@author: elcok
"""

import geopandas as gpd
import pandas as pd
import os
from functions import get_RAI,country_road_length#, shape_to_poly
from pathos.multiprocessing import Pool
import fiona
from shapely.geometry import mapping,shape
from rtree import index
from multiprocess import Pool , cpu_count 
import time
import numpy as np

if __name__ == "__main__":
    start = time.time()

    print ('The calculation of RAI has started!')
    
    
    """Load file paths"""
    base_path =   os.path.join(os.path.dirname(__file__),'..')

    # dir osm.pbf
    osm_path_in = os.path.join(base_path,'osm_continent')
    osm_path_out = os.path.join(base_path,'osm_country')

    # output dir
    dir_out = os.path.join(base_path,'output_data')
    ctry_out = os.path.join(base_path,'country_data')

    #input dir
    input_dir =  os.path.join(base_path,'input_data')
    # calc dir
    calc_dir = os.path.join(base_path,'calc')
    
    # WB country shapes:
    wb_country_in = os.path.join(base_path,'input_data','wbccodes2014.csv')
    grump = os.path.join(base_path,'Grump','global_urban_extent_cleaned.shp')   
    
    """Load country shapes and list and only save the required countries"""
    wb_country = pd.read_csv(wb_country_in,header=0,index_col=0)
    
    country_list = wb_country[['country','continent']].loc[wb_country['wbregion']!='YHI']
    
    map_continent = {'MA': 'central-america','SA': 'south-america','EU': 'europe','AS': 'asia',
                     'AU': 'australia-oceania','AF':'africa','AM':'north-america'}
   
    country_list['osm-cont'] = country_list['continent'].map(lambda x: (map_continent[x])) 

# =============================================================================
#     """ create .poly files to clip countries from osm.pbf files """
# =============================================================================

#   shape_to_poly(wb_poly_out,calc_dir)

            
# =============================================================================
#     """ create .poly files to clip countries from osm.pbf files """
# =============================================================================
#    if not os.listdir(poly_dir):
#        create_poly_files(base_path,global_shape,save_shapefile=True)

# =============================================================================
#       Estimate RAI
# =============================================================================
    
    # only primary and secondary
    out = []
    country_list = country_list.loc[country_list['continent'] != 'EU']
    countries = []
    continent_osms = []
    base_paths = []
    idx_paths = []
    for country in country_list.iterrows():
        country = country[1]
        continent_osm = os.path.join(osm_path_in,'%s-latest.osm.pbf' % (country['osm-cont']))
        countries.append(country['country'])
        continent_osms.append(continent_osm)
        idx_paths.append(grump)
        base_paths.append(base_path)
#        out.append(get_RAI(country['country'],continent_osm,base_path,idx_urban))

    pool = Pool(cpu_count()-1)
    out = pool.starmap(get_RAI, zip(countries,continent_osms,base_paths,idx_paths)) 

    df = pd.concat(out,axis=1).T
    
    map_dict = {0: 'Pop Rural Total', 1: 'Pop Rural < 2km', 2: 'RAI'}
    df = pd.concat([pd.DataFrame(l) for l in out],axis=1)
    df.index = df.index.map(mapper=(lambda x: map_dict[x]))

    writer = pd.ExcelWriter(os.path.join(dir_out,'RAI.xlsx'))
    df.T.to_excel(writer,'RAI_PS')
    writer.save()    
    
    end = time.time()

    print('It took ' + str(np.float16((end - start))) + " seconds to finish!")         



