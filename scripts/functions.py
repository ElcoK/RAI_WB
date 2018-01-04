
import os
from geopy.distance import vincenty
from boltons.iterutils import pairwise
from rasterstats import zonal_stats
import geopandas as gpd
from shapely.geometry import mapping,shape
import numpy as np
import fiona
from time import sleep
import pandas as pd
from rtree import index

def create_poly_files(base_path,global_shape,save_shapefile=False):

    """
    This function will create the .poly files from the world shapefile. These
    .poly files are used to extract data from the openstreetmap files.
    
    Args:
        base_path: base path to location of all files.
        
        global_shape: exact path to the global shapefile used to create the poly files.
        
        save_shape_file: when set True, the new shapefile with the countries that we include in this analysis will be saved.       
    
    Returns:
        .poly file for each country in a new dir in the working directory.
    """     
    
# =============================================================================
#     """ Create output dir for .poly files if it is doesnt exist yet"""
# =============================================================================
    poly_dir = os.path.join(base_path,'poly_files')
    
    if not os.path.exists(poly_dir):
        os.makedirs(poly_dir)
# =============================================================================
#     """ Set the paths for the files we are going to use """
# =============================================================================

    wb_poly_out = os.path.join(base_path,'input_data','country_shapes.shp')
    wb_country_in = os.path.join(base_path,'input_data','wbccodes2014.csv')
    
# =============================================================================
#   """Load country shapes and country list and only keep the required countries"""
# =============================================================================
    if global_shape.endswith('gdb'):
        wb_poly = gpd.read_file(global_shape,layer='g2015_2014_0')
    elif global_shape.endswith('.shp'):
        wb_poly = gpd.read_file(global_shape)
        
    wb_country = pd.read_csv(wb_country_in,header=0,index_col=0)
    
    # filter high income countries from country file
    country_list = wb_country['country'].loc[wb_country['wbregion']!='YHI']
    
    country_list = ['KSV']
    # filter polygon file
    wb_poly = wb_poly.loc[wb_poly['ISO3166_1_Alpha_3'].isin(country_list)]
    
    # we need to simplify the country shapes a bit. If the polygon is too diffcult,
    # osmconvert cannot handle it.
    wb_poly['geometry'] = wb_poly.simplify(tolerance = 0.005, preserve_topology=False)
    wb_poly.crs = {'init' :'epsg:4326'}

    # and save the new country shapefile if requested
    if save_shapefile == True:
        wb_poly.to_file(wb_poly_out)

# =============================================================================
#   """ The important part of this function: create .poly files to clip the country 
#   data from the openstreetmap file """    
# =============================================================================

    num = 0
    # iterate over the counties (rows) in the world shapefile
    for f in wb_poly.iterrows():
        f = f[1]
        num = num + 1
        geom=f.geometry
        
        # this will create a list of the different subpolygons
        if geom.geom_type == 'MultiPolygon':
            polygons = geom
        
        # the list will be lenght 1 if it is just one polygon
        elif geom.geom_type == 'Polygon':
            polygons = [geom]

        # define the name of the output file, based on the ISO3 code
        attr=f['ISO3166_1_Alpha_3']
        
        # start writing the .poly file
        f = open(poly_dir + "/" + attr +'.poly', 'w')
        f.write(attr + "\n")

        i = 0
        
        # loop over the different polygons, get their exterior and write the 
        # coordinates of the ring to the .poly file
        for polygon in polygons:
            polygon = np.array(polygon.exterior)

            j = 0
            f.write(str(i) + "\n")

            for ring in polygon:
                j = j + 1
                f.write("    " + str(ring[0]) + "     " + str(ring[1]) +"\n")

            i = i + 1
            # close the ring of one subpolygon if done
            f.write("END" +"\n")

        # close the file when done
        f.write("END" +"\n")
        f.close()

def get_country(country,continent_osm,base_path,overwrite=False):  

    """
    Extraction of the road data for the specified country from OSM. We use the continental OSM file, downloaded from http://download.geofabrik.de.
    

    Args:
        country: to country for which we calculate the RAI.
        
        continent_osm: the continent that the country 'belongs' to. This is 
                        required in the osm extraction.
                        
        base_path: base path to location of all files.
        
    Returns:
        load_country: a geodataframe containing all the roads and their length.
    """    

# =============================================================================
#     """ First set all paths for output dirs"""
# =============================================================================
    ctry_out = os.path.join(base_path,'country_data')
    osm_path_out = os.path.join(base_path,'osm_country')
    poly_dir = os.path.join(base_path,'poly_files')

    
# =============================================================================
#     """ Set the paths for the files we are going to use and write"""
# =============================================================================
    country_poly = os.path.join(poly_dir,country+'.poly') 
    country_shp =  os.path.join(ctry_out,country+'.shp') 
    country_pbf = os.path.join(osm_path_out,country+'.osm.pbf') 
    
# =============================================================================
#     # extract osm file for the country and write to shapefile
# =============================================================================

    if os.path.exists(country_pbf) is not True:
        clip_osm(continent_osm,country_poly,country_pbf)
        
    # get shapefile as output
#    if   (os.path.getsize(country_shp) == 143):
    if (os.path.exists(country_shp) is not True) or overwrite == True:
        extract_osm(country_shp,country_pbf)

# =============================================================================
#     Load the shapefile, remove road tags which only occur less than 15 times 
#     and estimate the length of the remaining roads    
# =============================================================================
    try:
        load_country = gpd.read_file(country_shp)
    except:
        sleep(30)
        try:
            load_country = gpd.read_file(country_shp)
        except:
            sleep(60)
            load_country = gpd.read_file(country_shp)

    uniq = load_country['highway'].value_counts()
    uniq = list(uniq[uniq > 20].index)
    load_country = load_country[load_country['highway'].isin(uniq)]
    load_country['distance'] = load_country['geometry'].apply(line_length)
    load_country = load_country[load_country['distance'] < 500]

# =============================================================================
#   Add a new column to the dataframe with the aggregated road classification
# =============================================================================
    load_country = map_roads(load_country)
    
# =============================================================================
#     And return the geodataframe
# =============================================================================
    return load_country


def get_RAI(country,continent_osm,base_path,grump):

    """
    Estimation of the Rural Accesibility Index (RAI) for the specified country.
    
    The RAI is calculated by using the methodology specified in the report 
    below:
        
        World Bank. 2016. Measuring rural access: using new technologies. 
        Washington, D.C. : World Bank Group.  http://bit.ly/2p5asME 
    
    Args:
        country: to country for which we calculate the RAI.
        
        continent_osm: the continent that the country 'belongs' to. This is required in the osm extraction.
        
        base_path: base path to location of all files.
        
        idx_urban: rTree index of the urban areas. By using an rtree index, we can quickly find the intersecting areas.
        
        tertiary: default option is 'False', but if we want to estimate the RAI for tertiary roads as well, we can set this option to 'True'.
        
    Returns:
        A dictionary with the country name as the key and the total rural population, total population with 2km of selected roads and the RAI as values.
    """ 
    try:
    
        print('%s started!' % country)
        tertiary=False
        
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
        if ('africa' in continent_osm):
            world_pop = os.path.join(base_path,'Worldpop','AFR_PPP_2015_adj_v2.tif') 
        if ('asia' in continent_osm):
            world_pop = os.path.join(base_path,'Worldpop','Asia_PPP_2015_adj_v2.tif') 

        # set path for global country shapes
        shp_world = os.path.join(base_path,'input_data','country_shapes.shp')
    
        # set country shapefile paths
        shp_country =  os.path.join(calc_dir,'%s.shp' % country)
        rural_country_shp =  os.path.join(calc_dir,'rural_%s.shp' % country)
        buffer_file = os.path.join(calc_dir,'%s_buffer.shp' % country)
        buffer_rural = os.path.join(calc_dir,'rural_roads_%s.shp' % country)
    
    # =============================================================================
    #     """ Get country boundary from world boundaries file and save to shp"""
    # =============================================================================
        world_boundaries = gpd.read_file(shp_world)
        country_boundary = world_boundaries.loc[world_boundaries['ISO3166_1_'] == country]
    #    country_boundary = world_boundaries.loc[world_boundaries['ISO_Codes'] == country]
        country_boundary.to_file(shp_country)

# =============================================================================
#     """ create spatial index of urban areas for quick intersection """
# =============================================================================
        idx_urban = index.Index()
        
        with fiona.open(grump) as shape_input:
            i = 0
    
            for poly in shape_input:
                idx_urban.insert(i, shape(poly['geometry']).bounds, poly)
                i += 1


    # =============================================================================
    #     """ Remove urban areas from the country boundary file and save to shp """
    # =============================================================================
        
        # get bounds of the country for the intersection with the rtree file
        bounds = country_boundary['geometry'].bounds
        bd_tuple = float(bounds['minx']),float(bounds['miny']),float(bounds['maxx']),float(bounds['maxy'])
        
        # we use the rtree file for this, to speed things up a bit
        hits = [n.object for n in idx_urban.intersection(bd_tuple, objects=True)]
    
        # write the country urban areas to a geopandas dataframe
    
        urban_geoms = gpd.GeoDataFrame([['urban',shape(x['geometry'])] for x in hits])
    
        urban_geoms.columns = ['urban_area','geometry']
        
        # now get the rural areas of the country and save this to a shp
        rural_country = gpd.overlay(country_boundary,urban_geoms,how='difference')
        rural_country.to_file(rural_country_shp)
    
        """ Estimate the total population living in rural areas """
        stats_ctry = sum(item['sum'] for item in zonal_stats(rural_country,world_pop,
                stats="sum") if item['sum'] is not None)        
    
    # =============================================================================
    #     """ Next step is to load roads of the country and create a buffer """
    # =============================================================================
       
        # load data, similar as when we estimate the length of the roads
        load_country = get_country(country,continent_osm,base_path)
       
        # if we include tertiary as well, we do a different filter
        if tertiary == True:
            load_country = load_country.loc[load_country['roads'].isin(['primary','secondary','tertiary'])]
        else:
            load_country = load_country.loc[load_country['roads'].isin(['primary','secondary'])]
            
        # to find the buffer around a road, we convert it quickly to a utm
        # coordinate system to be able to do an exact 2km around the road (it 
        # is a bit tricky to find this based on a coordinate sytem)
        
        country_centroid = country_boundary.centroid.values[0]
        lat,lon = country_centroid.bounds[1],country_centroid.bounds[0]
        
        # formula below based on :https://gis.stackexchange.com/a/190209/80697
        epsg=int(32700-np.round((45+lat)/90,0)*100+np.round((183+lon)/6,0))
  
        load_country = load_country.to_crs(epsg=epsg)
        
        # and do the actual buffer of 2km
        load_country['geometry'] = load_country.buffer(2000)
        
        load_country = load_country.to_crs(epsg=4326)
        
        # union this to one multipolygon and save to a shapefile
        poly = load_country['geometry'].unary_union
    
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
        
        # load files and overlay
        roads_buffer = gpd.read_file(buffer_file)
        rural_roads = gpd.overlay(roads_buffer,urban_geoms,how='difference')
        
        # and write to shapefile
        rural_roads.to_file(buffer_rural)
    
        """ Estimate the total rural population living within 2km of selected roads"""
        stats_roads = sum(item['sum'] for item in zonal_stats(buffer_rural,world_pop,
                stats="sum") if item['sum'] is not None)   
        
        """ Return the total rural population, total population with 2km of 
        selected roads and the RAI"""
        return {country: [stats_ctry,stats_roads,(stats_roads/stats_ctry)*100]}

    except:
        print('%s did not finish!' % country)

def country_road_length(country,continent_osm,base_path):
    
    """
    Function to estimate the road length for each road segment.

    Args:
        country: to country for which we calculate the RAI.
        
        continent_osm: the continent that the country 'belongs' to. This is required in the osm extraction.
        
        base_path: base path to location of all files.
        
    Returns:
        distance per road type (Primary, Secondary, Tertiary, Track, Other)
        for the specified country as pandas Series
    """
    try:
        print('%s started!' % country)
    
    # =============================================================================
    #     Load country road data
    # =============================================================================
        try:
            load_country = get_country(country,continent_osm,base_path,overwrite=False)
        except:
            sleep(30)
            load_country = get_country(country,continent_osm,base_path,overwrite=True)
            
    # =============================================================================
    #     Groupby aggregated road classifcation
    # =============================================================================
        dist_per_roadtype = load_country.groupby('roads').sum()
        dist_per_roadtype.columns = [country]
     
    # =============================================================================
    #     Return pandas Series with total road length in kilometers per road type
    # =============================================================================
        return dist_per_roadtype
    except Exception as e: print(str(e)+' for %s' % country)

def clip_osm(continent_osm,country_poly,country_pbf):
    """ Clip the country osm file from the larger continent (or planet) file and save to a new osm.pbf file. 
    This is much faster compared to clipping the osm.pbf file while extracting through ogr2ogr.
    
    This function uses the osmconvert tool, which can be found at http://wiki.openstreetmap.org/wiki/Osmconvert. 
    
    Either add the directory where this executable is located to your environmental variables or just put it in the 'scripts' directory.
    
    Args:
        continent_osm: path string to the osm.pbf file of the continent associated with the country.
        
        country_poly: path string to the .poly file, made through the 'create_poly_files' function.
        
        country_pbf: path string indicating the final output dir and output name of the new .osm.pbf file.
        
    Returns:
        a clipped .osm.pbf file.
    """ 

    os.system('osmconvert64  '+continent_osm+' -B='+country_poly+' --complete-ways -o='+country_pbf)


def extract_osm(country_shp,country_pbf):
    """Extract a shapefile with all the road information from the openstreetmap file.
    
    Args:
        country_shp: path_string indicating the final output directory and output name of the new shapefile.
        
        country_pbf: path string indicating the directory and name of the .osm.pbf file.
        
    Returns:
        A shapefile with all the roads of the clipped country. The shapefile will be in WGS84 (epsg:4326). This is the same coordinate system as Openstreetmap.
    """
    
    os.system("ogr2ogr -progress -overwrite -f \"ESRI Shapefile\" -sql \
               \"SELECT highway FROM lines WHERE highway IS NOT NULL\" \
               -lco ENCODING=UTF-8 -skipfailures -nlt geometry "+country_shp+" "+country_pbf)    
        


def line_length(line, ellipsoid='WGS-84'):
    """Length of a line in meters, given in geographic coordinates.

    Adapted from https://gis.stackexchange.com/questions/4022/looking-for-a-pythonic-way-to-calculate-the-length-of-a-wkt-linestring#answer-115285

    Args:
        line: a shapely LineString object with WGS-84 coordinates.
        
        ellipsoid: string name of an ellipsoid that `geopy` understands (see
            http://geopy.readthedocs.io/en/latest/#module-geopy.distance)

    Returns:
        Length of line in meters.
    """
    if line.geometryType() == 'MultiLineString':
        return sum(line_length(segment) for segment in line)

    return sum(
        vincenty(a, b, ellipsoid=ellipsoid).kilometers
        for a, b in pairwise(line.coords)
    )

def map_roads(load_country):

    """ 
    To create a new column with an aggregated list of road types. 
    
    Args:
        load_country: a geodataframe containing all the roads of a country.
        
    Returns:
        load_country: same geodataframe but with an additional 'roads' column containing the aggregated road types.
    """

    dict_map = {
"disused" : "other",
"dummy" : "other",
"planned" : "other",
"platform" : "other",
"unsurfaced" : "track",
"traffic_island" : "other",
"razed" : "other",
"abandoned" : "other",
"services" : "track",
"proposed" : "other",
"corridor" : "track",
"bus_guideway" : "other",
"bus_stop" : "other",
"rest_area" : "other",
"yes" : "other",
"trail" : "other",
"escape" : "track",
"raceway" : "other",
"emergency_access_point" : "track",
"emergency_bay" : "track",
"construction" : "track",
"bridleway" : "track",
"cycleway" : "other",
"footway" : "other",
"living_street" : "tertiary",
"path" : "track",
"pedestrian" : "other",
"primary" : "primary",
"primary_link" : "primary",
"residential" : "tertiary",
"road" : "secondary",
"secondary" : "secondary",
"secondary_link" : "secondary",
"service" : "tertiary",
"steps" : "other",
"tertiary" : "tertiary",
"tertiary_link" : "tertiary",
"track" : "track",
"unclassified" : "tertiary",
"trunk" : "primary",
"motorway" : "primary",
"trunk_link" : "primary",
"motorway_link" : "primary"
}
    
    load_country['roads'] = load_country['highway'].map(lambda x: (dict_map[x])) 
    
    return load_country
