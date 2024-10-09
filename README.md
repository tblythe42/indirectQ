Installation:
  The included .yml conda environment may be used to create the slope_area environment. If it does not work, installing conda dependencies and then geopandas via pip worked:
  conda create -n slope_area python=3.9
  conda install numpy pandas rasterio scipy gdal
  pip install geopandas
