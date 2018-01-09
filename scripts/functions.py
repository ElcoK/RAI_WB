
import os
from geopy.distance import vincenty
from boltons.iterutils import pairwise
from rasterstats import zonal_stats
import geopandas as gpd
from shapely.geometry import mapping
import numpy as np
import fiona
from time import sleep
import pandas as pd

import subprocess
import shutil
import warnings
from functools import reduce

import cartopy
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore", category=FutureWarning) 

def spatial_overlays(df1, df2, how='intersection'):
    """
    Compute overlay intersection of two GeoPandasDataFrames df1 and df2. This
    function is much faster compared to the original geopandas overlay method
    
    Adapted from https://github.com/geopandas/geopandas/issues/400 
    
    Args:
        df1: dataframe to be either intersected or 'erased', using the difference function.
            
        df2: dataframe to intersect or to erase, using the difference function.
        
    Return:
        df1: an either intersected or (partly) erased geopandas dataframe.
    """

    df1 = df1.copy()
    df2 = df2.copy()
    df1['geometry'] = df1.geometry.buffer(0)
    df2['geometry'] = df2.geometry.buffer(0)

# =============================================================================
#     """ the part of the function which does the intersection analysis """
# =============================================================================
    if how=='intersection':
        # Spatial Index to create intersections
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['histreg']=df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        pairs = df1['histreg'].to_dict()
        nei = []
        for i,j in pairs.items():
            for k in j:
                nei.append([i,k])
        
        pairs = gpd.GeoDataFrame(nei, columns=['idx1','idx2'], crs=df1.crs)
        pairs = pairs.merge(df1, left_on='idx1', right_index=True)
        pairs = pairs.merge(df2, left_on='idx2', right_index=True, suffixes=['_1','_2'])
        pairs['Intersection'] = pairs.apply(lambda x: (x['geometry_1'].intersection(x['geometry_2'])).buffer(0), axis=1)
        pairs = gpd.GeoDataFrame(pairs, columns=pairs.columns, crs=df1.crs)
        cols = pairs.columns.tolist()
        cols.remove('geometry_1')
        cols.remove('geometry_2')
        cols.remove('histreg')
        cols.remove('bbox')
        cols.remove('Intersection')
        dfinter = pairs[cols+['Intersection']].copy()
        dfinter.rename(columns={'Intersection':'geometry'}, inplace=True)
        dfinter = gpd.GeoDataFrame(dfinter, columns=dfinter.columns, crs=pairs.crs)
        dfinter = dfinter.loc[dfinter.geometry.is_empty==False]
        return dfinter
# =============================================================================
#     """ the part of the function which does the difference/erase analysis """
# =============================================================================
    elif how=='difference':
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['histreg']=df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        df1['new_g'] = df1.apply(lambda x: reduce(lambda x, y: x.difference(y).buffer(0), [x.geometry]+list(df2.iloc[x.histreg].geometry)) , axis=1)
        df1.geometry = df1.new_g
        df1 = df1.loc[df1.geometry.is_empty==False].copy()
        df1.drop(['bbox', 'histreg', 'new_g'], axis=1, inplace=True)
        return df1

def create_poly_files(base_path,global_shape,save_shapefile=True):

    """
    This function will create the .poly files from the world shapefile. These
    .poly files are used to extract data from the openstreetmap files.
    
    This function is adapted from the OSMPoly function in QGIS.
    
    Args:
        base_path: base path to location of all files.
        
        global_shape: exact path to the global shapefile used to create the poly files.
        
        save_shape_file: when True, the new shapefile with the countries that we include in this analysis will be saved.       
    
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
    
    # filter polygon file
    wb_poly = wb_poly.loc[wb_poly['ISO3166_1_Alpha_3'].isin(country_list)]
    wb_poly.crs = {'init' :'epsg:4326'}

    # and save the new country shapefile if requested
    if save_shapefile == True:
        wb_poly.to_file(wb_poly_out)
    
    # we need to simplify the country shapes a bit. If the polygon is too diffcult,
    # osmconvert cannot handle it.
    wb_poly['geometry'] = wb_poly.simplify(tolerance = 0.005, preserve_topology=False)


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

def get_country(country,continent_osm,base_path,overwrite=False,RAI=False):  

    """
    Extraction of the road data for the specified country from OSM. We use the continental OSM file, downloaded from http://download.geofabrik.de.
    

    Args:
        country: to country for which we calculate the RAI.
        
        continent_osm: the continent that the country 'belongs' to. This is 
                        required in the osm extraction.
                        
        base_path: base path to location of all files.
        
        overwrite: this is set on True by default, set on False if you are sure the input files are correct (it will save a lot of computation time).

        RAI: this is set on False by default. set on True if this country extracting is used for the RAI analysis. It will skip the road length calculation (saving computation time).
        
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
    uniq.extend(['primary','secondary','trunk','motorway'])
    uniq = list(set(uniq))
    load_country = load_country[load_country['highway'].isin(uniq)]
    
    if RAI is False:
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

def get_RAI(country,continent_osm,base_path,grump,overwrite=False,tertiary=False):

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
        country_boundary.crs = {'init' :'epsg:4326'}
    #    country_boundary = world_boundaries.loc[world_boundaries['ISO_Codes'] == country]
        country_boundary.to_file(shp_country)
        
        urban_geoms = gpd.read_file(grump)
        urban_geoms.crs = {'init' :'epsg:4326'}
        
    # =============================================================================
    #      """ Erase urban areas from the country shape  
    # =============================================================================
        
        if country_boundary['geometry'].values[0].geom_type == 'MultiPolygon':
            ctry_boundary = gpd.GeoDataFrame([polygon for polygon in country_boundary['geometry'].values[0]],columns=['geometry'])
            ctry_boundary = ctry_boundary.loc[ctry_boundary.geometry.area > 0.001]
        else:
            ctry_boundary = country_boundary
        # now get the rural areas of the country and save this to a shp
        rural_country =  spatial_overlays(ctry_boundary, urban_geoms, how='difference')
    #        rural_country = spatial_difference(ctry_boundary, urban_geoms)
        rural_country.to_file(rural_country_shp)
        
        """ Estimate the total population living in rural areas """
        stats_ctry = sum(item['sum'] for item in zonal_stats(rural_country,world_pop,
                stats="sum") if item['sum'] is not None)
    
    # =============================================================================
    #     """ Next step is to load roads of the country and create a buffer """
    # =============================================================================
       
        # load data, similar as when we estimate the length of the roads
        load_country = get_country(country,continent_osm,base_path,overwrite,RAI=True)
       
        # if we include tertiary as well, we do a different filter
        if tertiary == True:
            load_country = load_country.loc[load_country['roads'].isin(['primary','secondary','tertiary'])]
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
        try:
            # load files and overlay
            roads_buffer = gpd.GeoDataFrame(gpd.read_file(buffer_file)['geometry'])
            roads_buffer.reset_index(inplace=True,drop=True)
            rural_country.reset_index(inplace=True,drop=True)
            rural_roads = spatial_overlays(roads_buffer, rural_country, how='intersection')
            # and write to shapefile
            rural_roads.to_file(buffer_rural)
        except:
            if country in islands:
                try:
                    shutil.rmtree(temp_path, ignore_errors=True)
                except:
                    True
            return {country: [stats_ctry,0,0]}
    
    
        """ Estimate the total rural population living within 2km of selected roads"""
        stats_roads = sum(item['sum'] for item in zonal_stats(buffer_rural,world_pop,
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


def create_figure(country,load_country,base_path):
    """
    Create a figure of the road network
    """
    # Create figure
    plt.figure()
    
    proj_lat_lon = cartopy.crs.PlateCarree()
    ax = plt.axes([0.0,0.0,1.0, 1.0] ,facecolor='#D0E3F4', projection=proj_lat_lon)

    cmap = cm.get_cmap('Set1', 4) # Colour map (there are many others)
    cmaplist = [cmap(i) for i in range(cmap.N)]
#        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

    colors_legend = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors_legend.append(mpl.colors.rgb2hex(rgb))

    color1 =  Line2D([0], [0], linestyle="none", marker="o", markersize=10, markeredgewidth=0.0,markerfacecolor=colors_legend[0])
    color2 =  Line2D([0], [0], linestyle="none", marker="o", markersize=10, markeredgewidth=0.0,markerfacecolor=colors_legend[1])
    color3 =  Line2D([0], [0], linestyle="none", marker="o", markersize=10, markeredgewidth=0.0,markerfacecolor=colors_legend[2])
    color4 =  Line2D([0], [0], linestyle="none", marker="o", markersize=10, markeredgewidth=0.0,markerfacecolor=colors_legend[3])

    legend_handles=['primary','secondary','tertiary','track']
    
    map_linewidth = dict(zip(legend_handles,[0.8,0.5,0.3,0.1]))
    
    load_country = load_country.loc[load_country['roads'] != 'other']
    lw_all = list(load_country['roads'].map(lambda x: (map_linewidth[x])) )
    tot_bounds = list(load_country.total_bounds)
    ax.set_extent([tot_bounds[0]-0.1,tot_bounds[2]+0.1,tot_bounds[1]-0.1,tot_bounds[3]+0.1] , crs=proj_lat_lon)
    world = gpd.read_file('F:\Dropbox\WorldBank\input_data\country_shapes.shp')
    world.plot(ax=ax,color='#FEF9E0',lw=0.3,edgecolor='k')
               
    load_country.plot(ax=ax,column='roads',linewidth=lw_all,categorical=True,legend=True,cmap=cmap)

    ax.legend([color1,color2,color3,color4],legend_handles,numpoints=1, loc=1) 
    ax.background_patch.set_facecolor('#D0E3F4')
    
    admin_name = world['ADM0_NAME'].loc[world['ISO3166_1_']==country].values[0]
    plt.title('Road network of %s' % admin_name,fontweight='bold')

    figure_out= os.path.join(base_path,'Figures','%s.png' % country)
    plt.savefig(figure_out,dpi=500,bbox_inches='tight')

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
        A shapefile with all the roads of the clipped country. The shapefile will be in *WGS84* (epsg:4326). This is the same coordinate system as Openstreetmap.
    """
    
    os.system("ogr2ogr -progress -overwrite -f \"ESRI Shapefile\" -sql \
               \"SELECT highway FROM lines WHERE highway IS NOT NULL\" \
               -lco ENCODING=UTF-8 -skipfailures -nlt geometry "+country_shp+" "+country_pbf)    
        

def line_length(line, ellipsoid='WGS-84'):
    """Length of a line in meters, given in geographic coordinates.

    Adapted from https://gis.stackexchange.com/questions/4022/looking-for-a-pythonic-way-to-calculate-the-length-of-a-wkt-linestring#answer-115285

    Args:
        line: a shapely LineString object with WGS-84 coordinates.
        
        ellipsoid: string name of an ellipsoid that `geopy` understands (see http://geopy.readthedocs.io/en/latest/#module-geopy.distance).

    Returns:
        Length of line in meters.
    """
    if line.geometryType() == 'MultiLineString':
        return sum(line_length(segment) for segment in line)

    return sum(
        vincenty(a, b, ellipsoid=ellipsoid).kilometers
        for a, b in pairwise(line.coords)
    )

def unzip_worldpop(country,base_path,temp_path):
 
    archiveman = r'7z' # 7z.exe comes with 7-zip

    """Load file paths"""
    poppath = os.path.join(base_path,'Worldpop')

    islands = ['FJI','KIR','MHL','FSM','PLW','PNG','WSM','SLB','TON','VUT']
    islands_full = ['Fiji','Kiribati','Marshall Islands','Micronesia (Federated States of)','Palau','Papua New Guinea','Samoa','Solomon Islands','Tonga','Vanuatu']

    map_isl_names = dict(zip(islands,islands_full))

    island_full_name = map_isl_names[country]
    archive_name = os.path.join(poppath,'%s 100m Population.7z' % island_full_name)
    if country == 'KIR':
        _ = subprocess.check_output([archiveman, 'e','-aos', archive_name, '-o{}'.format(temp_path), 'popmap15adj_lzw.tif'])
    elif country == 'TON':
        _ = subprocess.check_output([archiveman, 'e','-aos', archive_name, '-o{}'.format(temp_path), 'popmap15.tif'])
    else:
        _ = subprocess.check_output([archiveman, 'e','-aos', archive_name, '-o{}'.format(temp_path), 'popmap15adj.tif'])

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
