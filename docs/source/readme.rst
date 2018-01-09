
=========================
Getting started
=========================

Recommended option is to use a `Miniconda <https://conda.io/miniconda.html>`_ environment to work in for this project, relying on conda to handle some of the trickier library dependencies.

The most easiest way to do so is to use the requirements.yml file as provided in the github page. 

*To create the environment using the yaml file:*

   .. code::

		conda env create -f environment.yml
	
In case of no access to the GitHub page, the other option would be to copy-paste the list below and save this to an `environment.yml` file (note the indendation):
	
   .. code::

		name: py36
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

	