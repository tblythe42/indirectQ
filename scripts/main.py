# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from pathlib import Path
import rasterio as rio
import geopandas as gpd

from runconfig.inputpaths import inputpaths
#from indirectQ import channelbuilder
from indirectQ import utilities

topo_pth = Path(inputpaths['topography']['dir']) / Path(inputpaths['topography']['file'])
d_xsec_pth = Path(inputpaths['drawn_cross_sections']['dir']) / Path(inputpaths['drawn_cross_sections']['file'])
chan_cnt_pth = Path(inputpaths['channel_centerline']['dir']) / Path(inputpaths['channel_centerline']['file'])
ws_pth = Path(inputpaths['water_surface_profile']['dir']) / Path(inputpaths['water_surface_profile']['file'])

chn_cnt = gpd.read_file(chan_cnt_pth)
d_xsec = gpd.read_file(d_xsec_pth)

utilities.points_along_lines(d_xsec_pth, ID_field='XSID', keep_attrs=['XSID'])

with rio.open(topo_pth, 'r') as tp_dset:
    topoZ = tp_dset.read(1)
    bnds = tp_dset.bounds

def print_hi(name):


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
