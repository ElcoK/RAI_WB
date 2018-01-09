
=========================
Getting started
=========================

The recommended option is to create a virtual environment using `Miniconda <https://conda.io/miniconda.html>`_. 

The easiest way to create a virtual environment for this project is to use the *requirements.yml* file as provided on the GitHub page of this project. 

*To create the environment using the yaml file:*

   .. code::

		conda env create -f environment.yml
	
In case of no access to the GitHub page, the other option would be to copy-paste the list below and save this to an `environment.yml` file (note the indendation):
	
   .. code::

		name: RAI_WB
		   dependencies:
			- python=3.6
			- gdal
			- numpy
			- pandas
			- shapely
			- fiona
			- geopy
			- rtree
			- geopandas
			- pathos
			- rasterio
			- rasterstats
			- boltons
			- matplotlib
			- cartopy

	
