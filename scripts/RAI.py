# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 16:04:29 2018

@author: cenv0574
"""

import os
from functions import unzip_worldpop,spatial_overlays,get_country,create_poly_files
import geopandas as gpd
import shutil
from rasterstats import zonal_stats
import numpy as np
import fiona
from shapely.geometry import mapping
import time
from multiprocess import Pool , cpu_count 
import urllib.request
import pandas as pd

def single_country(country,continent_osm,base_path,grump,overwrite=False,tertiary=False,track=False):

    """
    Estimation of the Rural Accesibility Index (RAI) for the specified country.
    
    The RAI is calculated by using the methodology specified in the report below:
        
        World Bank. 2016. Measuring rural access: using new technologies. 
        Washington, D.C. : World Bank Group.  http://bit.ly/2p5asME 
    
    Args:
        country: to country for which we calculate the RAI.
        
        continent_osm: the continent that the country 'belongs' to. This is required in the osm extraction.
        
        base_path: base path to location of all files.
        
        idx_urban: rTree index of the urban areas. By using an rtree index, we can quickly find the intersecting areas.
        
        tertiary: default option is 'False', but if we want to estimate the RAI for tertiary roads as well, we can set this option to 'True'.

        track: default option is 'False', but if we want to estimate the RAI for both tracks and tertiary roads as well, we can set this option to 'True'.
        
    Returns:
        A dictionary with the country name as the key and the total rural population, total population with 2km of selected roads and the RAI as values.
    """ 
    try:
    
        print('%s started!' % country)
        
    # =============================================================================
    #     """ First set all paths for output dirs"""
    # =============================================================================
        calc_dir = os.path.join(base_path,'calc')
        
    # =============================================================================
    #     """ Set the paths for the files we are going to use """
    # =============================================================================
    
        #set global file paths for worldpop
        if ('central-america' in continent_osm) or ('south-america' in continent_osm) :
            world_pop = os.path.join(base_path,'Worldpop','LAC_PPP_2015_adj_v2.tif') 
        elif ('africa' in continent_osm):
            world_pop = os.path.join(base_path,'Worldpop','AFR_PPP_2015_adj_v2.tif') 
        elif ('asia' in continent_osm):
            world_pop = os.path.join(base_path,'Worldpop','Asia_PPP_2015_adj_v2.tif') 
        elif ('europe' in continent_osm):
            world_pop = os.path.join(base_path,'Worldpop','EUROPOP_WGS84.tif') 

        # due to a few islands not included in the global worldpop data, we need to 
        # extract their data individually
        
        islands = ['FJI','KIR','MHL','FSM','PLW','PNG','WSM','SLB','TON','VUT']
        if country in islands:
            temp_path = os.path.join(base_path,'Worldpop','temp_%s' % country)
    
            if not os.path.exists(temp_path):
                os.makedirs(temp_path)
    
            unzip_worldpop(country,base_path,temp_path)
            world_pop = os.path.join(temp_path,'popmap15adj.tif')
            if country == 'KIR':
                world_pop = os.path.join(temp_path,'popmap15adj_lzw.tif')
            elif country == 'TON':
                world_pop = os.path.join(temp_path,'popmap15.tif')
            world_pop_out = os.path.join(temp_path,'popmap15adj_wgs84.tif')
            os.system('gdalwarp -t_srs EPSG:4326 -tr 0.018 0.018 -overwrite '+world_pop+' '+world_pop_out)
            world_pop = world_pop_out

        if country == 'MEX':
            world_pop = os.path.join(base_path,'Worldpop','MEX_ppp_v2c_2015_UNadj.tif')   
        elif country == 'UKR':
            world_pop = os.path.join(base_path,'Worldpop','UKR_ppp_v2b_2015_UNadj.tif')
     
        # set path for global country shapes
        shp_world = os.path.join(base_path,'input_data','country_shapes.shp')
    
        # set country shapefile paths
        shp_country =  os.path.join(calc_dir,'%s.shp' % country)
        buffer_file = os.path.join(calc_dir,'%s_buffer.shp' % country)
    
    # =============================================================================
    #     """ Get country boundary from world boundaries file and save to shp"""
    # =============================================================================
        world_boundaries = gpd.read_file(shp_world)
        country_boundary = world_boundaries.loc[world_boundaries['ISO3166_1_'] == country]
        country_boundary.crs = {'init' :'epsg:4326'}
    #    country_boundary = world_boundaries.loc[world_boundaries['ISO_Codes'] == country]
        country_boundary.to_file(shp_country)
        
        urban_geoms = gpd.read_file(grump)
        urban_geoms.crs = {'init' :'epsg:4326'}
 
        grump_country_in = os.path.join(base_path,'calc','grump_%s.shp' % country)
        grump_country_out = os.path.join(base_path,'calc','exact_grump_%s.shp' % country)

        if len(urban_geoms.loc[urban_geoms['ISO3']==country]) != 0:        
            urban_geoms = urban_geoms.loc[urban_geoms['ISO3']==country]
            urban_geoms.to_file(grump_country_in)
            os.system('ogr2ogr -skipfailures -nlt geometry -clipsrc '+shp_country+' '+grump_country_out+' '+grump_country_in)
        elif country == 'TUV':
            None
        else:
            urban_geoms = spatial_overlays(country_boundary,urban_geoms)
            urban_geoms.to_file(grump_country_in)
            os.system('ogr2ogr -skipfailures -nlt geometry -clipsrc '+shp_country+' '+grump_country_out+' '+grump_country_in)
      
    # =============================================================================
    #      """ Subtract the people living in urban areas from the total population  
    # =============================================================================

        """ Estimate the total population living in rural areas """
        country_tot = sum(item['sum'] for item in zonal_stats(country_boundary,world_pop,
            stats="sum") if item['sum'] is not None)
        if country == 'TUV':
            stats_ctry  = country_tot
        else:
            stats_ctry = country_tot- sum(item['sum'] for item in zonal_stats(gpd.read_file(grump_country_out),world_pop,
                    stats="sum") if item['sum'] is not None) 
    
    # =============================================================================
    #     """ Next step is to load roads of the country and create a buffer """
    # =============================================================================
       
        # load data, similar as when we estimate the length of the roads
        load_country = get_country(country,continent_osm,base_path,overwrite=False,RAI=True)
       
        # if we include tertiary as well, we do a different filter
        if tertiary == True:
            load_country = load_country.loc[load_country['roads'].isin(['primary','secondary','tertiary'])]
        elif track == True:
            load_country = load_country.loc[load_country['roads'].isin(['primary','secondary','tertiary','track'])]
        else:
            load_country = load_country.loc[load_country['roads'].isin(['primary','secondary'])]
            
        # to find the buffer around a road, we convert it quickly to a utm
        # coordinate system to be able to do an exact 2km around the road (it 
        # is a bit tricky to find this based on a WGS84 coordinate sytem)
        
        if len(load_country) == 0:
            if country in islands:
                try:
                    shutil.rmtree(temp_path, ignore_errors=True)
                except:
                    True
            return {country: [stats_ctry,0,0]}
        
        country_centroid = country_boundary.centroid.values[0]
        lat,lon = country_centroid.bounds[1],country_centroid.bounds[0]
        
        # formula below based on :https://gis.stackexchange.com/a/190209/80697
        epsg=int(32700-np.round((45+lat)/90,0)*100+np.round((183+lon)/6,0))
      
        load_country = load_country.to_crs(epsg=epsg)
        
        # and do the actual buffer of 2km
        load_country['geometry'] = load_country.buffer(2000)
        
        # and change the geodataframe back to WGS84 again
        load_country = load_country.to_crs(epsg=4326)
        
        # union this to one multipolygon and save to a shapefile
        try:
            poly = load_country['geometry'].unary_union
        except:
            poly = load_country['geometry'].buffer(0).unary_union
        
        # define a polygon feature geometry with one attribute
        schema = {
            'geometry': 'Polygon',
            'properties': {'id': 'int'},
        }
        
        # write a new Shapefile
        with fiona.open(buffer_file, 'w', 'ESRI Shapefile', schema) as c:
            ## If there are multiple geometries, put the "for" loop here
            c.write({
                'geometry': mapping(poly),
                'properties': {'id': 0},
            }) 
        
    # =============================================================================
    #     """ Here we exlude urban areas from the road buffer file, similar as how
    #     we did it in the country boundary file """
    # =============================================================================
        buffer_inb = os.path.join(base_path,'calc','%s_buffer.shp' % country)
        buffer_in = os.path.join(base_path,'calc','%s_buffer_exact.shp' % country)

        buffer_out = os.path.join(base_path,'calc','%s_buffer_urban.shp' % country)
        os.system('ogr2ogr -skipfailures -nlt geometry -clipsrc '+buffer_inb+' '+buffer_in+' '+shp_country)
 
        total_buffer = gpd.read_file(buffer_in)  
        if country == 'TUV':
             None   
        else:
            urban_geoms_exact = gpd.read_file(grump_country_out)
            urban_buffer = spatial_overlays(total_buffer,urban_geoms_exact,how='intersection')
            urban_buffer.to_file(buffer_out)    
   
        """ Estimate the total rural population living within 2km of selected roads"""
        tot_buf_country = sum(item['sum'] for item in zonal_stats(total_buffer,world_pop,
                stats="sum") if item['sum'] is not None) 
        if country == 'TUV':
            stats_roads = tot_buf_country
        else:
            stats_roads = tot_buf_country- sum(item['sum'] for item in zonal_stats(urban_buffer,world_pop,
                stats="sum") if item['sum'] is not None)
        
        if country in islands:
            try:
                shutil.rmtree(temp_path, ignore_errors=True)
            except:
                True
        
        """ Return the total rural population, total population with 2km of 
        selected roads and the RAI"""
        
        if stats_roads == 0:
            return {country: [stats_ctry,0,0]}
        
        stats =(stats_roads/stats_ctry)*100
        if stats > 100:
            stats = 100
            
        return {country: [stats_ctry,stats_roads,stats]}

    except Exception as e: 
        print(str(e)+' for %s' % country)
        return {country: [0,0,0]}


def all_countries(base_path,multiprocess=True,overwrite=True,tertiary=False):
    
    """
    Main function to estimate the **RAI** for all the countries we are interested in. 

    Args:
        *base_path* : Base path to the location of all files and directories in this project.
        
        *multiprocess* : Set to True by default. Set to False in the case of limited processing power.
        
        *overwrite* : Set to True by default. This relates to all input data (i.e. .poly files, .osm.pbf files and shapefiles).

        *tertiary* : Set to False by default. When set to True, the calculation of the **RAI** will run a second time, now including tertiary roads.

    Returns:
        An Excel file the **total rural population**, the **total rural population within 2km of the selected roads** and the **RAI** for each country. if *tertiary* is set to True, the Excel file will return an additional sheet with the results of that second calculation.
    
    """
    
    print ('The calculation of road lenghts has started!')
    start = time.time()

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
    overwrites = []
    tertaries_true = []
    tertaries_false = []
    for country in country_list.iterrows():
        country = country[1]
        continent_osm = os.path.join(osm_path_in,'%s-latest.osm.pbf' % (country['osm-cont']))
        countries.append(country['country'])
        continent_osms.append(continent_osm)
        base_paths.append(base_path)
        overwrites.append(overwrite)
        tertaries_true.append(True)
        tertaries_false.append(False)

    # multiprocessing will start if set to True. Set to False with limited processing capacities    
    if multiprocess==True:
        pool = Pool(cpu_count()-1)
        out = pool.starmap(single_country, zip(countries,continent_osms,base_paths,overwrites,tertaries_false)) 
    
    # when multiprocessing set to False, we will just loop over the countries.
    else:
        out = []
        i = 0
        for country in country_list.iterrows():
            country = country[1]
            continent_osm = os.path.join(osm_path_in,'%s-latest.osm.pbf' % (country['osm-cont']))
            out.append(single_country(country['country'],continent_osm,base_path,overwrites[i],tertaries_false[i]))
            i += 1
            
    map_dict = {0: 'Pop Rural Total', 1: 'Pop Rural < 2km', 2: 'RAI'}
    df = pd.concat([pd.DataFrame(l) for l in out],axis=1)
    df.index = df.index.map(mapper=(lambda x: map_dict[x]))

    map_country = dict(zip(wb_country['country'],wb_country['country_name']))
    df = df.T
    df['Country'] = df.index.to_series().map(map_country)
    df.set_index('Country',inplace=True,drop=True)

    writer = pd.ExcelWriter(os.path.join(dir_out,'RAI.xlsx'))
    df.to_excel(writer,'RAI_PS')

    if tertiary ==True:
        if multiprocess==True:
            pool = Pool(cpu_count()-1)
            out = pool.starmap(single_country, zip(countries,continent_osms,base_paths,overwrites,tertaries_false)) 
        
        # when multiprocessing set to False, we will just loop over the countries.
        else:
            out = []
            i = 0
            for country in country_list.iterrows():
                country = country[1]
                continent_osm = os.path.join(osm_path_in,'%s-latest.osm.pbf' % (country['osm-cont']))
                out.append(single_country(country['country'],continent_osm,base_path,overwrites[i],tertaries_false[i]))
                i += 1        

        df = pd.concat([pd.DataFrame(l) for l in out],axis=1)
        df.index = df.index.map(mapper=(lambda x: map_dict[x]))
    
        map_country = dict(zip(wb_country['country'],wb_country['country_name']))
        df = df.T
        df['Country'] = df.index.to_series().map(map_country)
        df.set_index('Country',inplace=True,drop=True)
    
        df.to_excel(writer,'RAI_PST')

    writer.save()    
    
    end = time.time()

    print('It took ' + str(np.float16((end - start))) + " seconds to finish!")         