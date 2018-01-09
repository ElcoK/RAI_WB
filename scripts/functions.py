
import os
from geopy.distance import vincenty
from boltons.iterutils import pairwise
import geopandas as gpd
import numpy as np
from time import sleep
import pandas as pd

import subprocess
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
    function is much faster compared to the original geopandas overlay method.
    
    Adapted from https://github.com/geopandas/geopandas/issues/400 
    
    Args:
        *df1* : The dataframe to be either intersected or 'erased', using the difference function.
            
        *df2* : The dataframe to intersect or to erase, using the difference function.
        
    Return:
        *df1*: an either intersected or (partly) erased geopandas dataframe.
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
    
    This function is adapted from the Export OSM Poly function in QGIS.
    
    Args:
        *base_path* : The base path to location of all files.
        
        *global_shape* : The exact path to the global shapefile used to create the poly files.
        
        *save_shape_file* : When set to True, the new shapefile with the countries that we include in this analysis will be saved.       
    
    Returns:
        A .poly file for each country in the 'poly_files' directory in the global working directory.
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
        *country* : The country for which we calculate the RAI.
        
        *continent_osm* : The continent the country 'belongs' to. This is required for the osm extraction.
                        
        *base_path* : The base path to location of all files and scripts.
        
        *overwrite* : This is set on True by default, set on False if you are sure the input files are correct (it will save a lot of computation time).

        *RAI* : This is set on False by default. set on True if this country extracting is used for the RAI analysis. It will skip the road length calculation (saving computation time).
        
    Returns:
        *load_country* : A geodataframe containing all the roads and their length.
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

def create_figure(country,load_country,base_path):
    """
    Create a figure of the road network
    
    Args:
        *country* : The ISO-code of the country.
        
        *load_country* : The geodataframe of the country.
        
        *base_path* : The base path to the location of all files and scripts.

    Returns:
        A .png file with a classified road system of the specified country.
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
    
    Either add the directory where this executable is located to your environmental variables or just put it in the **scripts** directory.
    
    Args:
        *continent_osm* : path string to the osm.pbf file of the continent associated with the country.
        
        *country_poly* : path string to the .poly file, made through the 'create_poly_files' function.
        
        *country_pbf* : path string indicating the final output dir and output name of the new .osm.pbf file.
        
    Returns:
        A clipped .osm.pbf file.
    """ 

    os.system('osmconvert64  '+continent_osm+' -B='+country_poly+' --complete-ways -o='+country_pbf)

def extract_osm(country_shp,country_pbf):
    """Extract a shapefile with all the road information from the openstreetmap file.
    
    Args:
        *country_shp* : The path string indicating the final output directory and output name of the new shapefile.
        
        *country_pbf* : The path string indicating the directory and name of the .osm.pbf file.
        
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
        *line* : A shapely LineString object with WGS-84 coordinates.
        
        *ellipsoid* : The string name of an ellipsoid that `geopy` understands (see http://geopy.readthedocs.io/en/latest/#module-geopy.distance).

    Returns:
        The length of the line in meters.
    """
    if line.geometryType() == 'MultiLineString':
        return sum(line_length(segment) for segment in line)

    return sum(
        vincenty(a, b, ellipsoid=ellipsoid).kilometers
        for a, b in pairwise(line.coords)
    )

def unzip_worldpop(country,base_path,temp_path):

    """Function to unzip the worldpop files for the following islands (not included in the main worldpop file):

        - Fiji
        - Kiribati
        - Marshall Islands
        - Micronesia (Federated States of)
        - Palau
        - Papua New Guinea
        - Samoa
        - Solomon Islands
        - Tonga
        - Vanuatu

    This should work out of the box if *7zip* is installed.
    
    If not: the easiest way to get this to run, is to move the *7z.exe* into the **scripts** directory. The other option would be to add the *7zip* directory to your environmental variables.

    Args:
        *country* : The ISO-code of the country.

        *base_path* : The base path to the location of all files and scripts.
        
        *temp_path* : The path to which we temporarily write the unzipped file.
        
    Returns:
        An unzipped worldpop GeoTIFF file for the specified country.
        
    """
    
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
        *load_country* : A geodataframe containing all the roads of a country.
        
    Returns:
        *load_country* : The same geodataframe but with an additional 'roads' column containing the aggregated road types.
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
