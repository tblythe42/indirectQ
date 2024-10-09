# Example usage of indirectQ for probabilistic (Monte Carlo style) multiple simulation estimate of discharge
#   based on distribution of Manning's n values

from pathlib import Path
import rasterio as rio
import geopandas as gpd
import pandas as pd

from runconfig.inputpaths import inputpaths
from indirectQ import channelbuilder
from indirectQ import utilities
from indirectQ import slopearea
from scipy.stats import norm

topo_pth = Path(inputpaths['topography']['dir']) / Path(inputpaths['topography']['file'])
d_xsec_pth = Path(inputpaths['drawn_cross_sections']['dir']) / Path(inputpaths['drawn_cross_sections']['file'])
chan_cnt_pth = Path(inputpaths['channel_centerline']['dir']) / Path(inputpaths['channel_centerline']['file'])
ws_pth = Path(inputpaths['water_surface_profile']['dir']) / Path(inputpaths['water_surface_profile']['file'])
man_reg = Path(inputpaths['Mannings_regions']['dir']) / Path(inputpaths['Mannings_regions']['file'])

surv = channelbuilder.LoadSurvey(chan_cnt_pth, ws_pth, "EPSG:6515", topo_pth=topo_pth)
surv.extractXsections(d_xsec_pth)
surv.assign_Mannings(man_reg)
surv.assign_wse()

man_r = gpd.read_file(man_reg)
rnames = man_r['Name'].tolist()
nvals = [0.12, 0.048, 0.048, 0.12, 0.02, 0.02, 0.09, 0.068, 0.064, 0.09]
samp_size = 1000

uprov = norm.rvs(loc=0.11, scale=0.01, size=samp_size)
uprch = norm.rvs(loc=0.05, scale=0.01, size=samp_size)
uplch = norm.rvs(loc=0.05, scale=0.01, size=samp_size)
uplov = norm.rvs(loc=0.11, scale=0.01, size=samp_size)
dnrov = norm.rvs(loc=0.11, scale=0.01, size=samp_size)
dnrch = norm.rvs(loc=0.05, scale=0.01, size=samp_size)
dnlch = norm.rvs(loc=0.05, scale=0.01, size=samp_size)
dnlov = norm.rvs(loc=0.11, scale=0.01, size=samp_size)

def change_manningsn(names, vals, shape):
    man_n_dict = dict(zip(names, vals))
    shape['n_Value'] = shape['Name'].map(man_n_dict)
    return shape

Qs = []
for i in range(samp_size):
    nvals[0] = uprov[i]
    nvals[1] = uprch[i]
    nvals[2] = uplch[i]
    nvals[3] = uplov[i]
    nvals[6] = dnrov[i]
    nvals[7] = dnrch[i]
    nvals[8] = dnlch[i]
    nvals[9] = dnlov[i]

    xsecs = surv.extracted_xsections
    for key, item in xsecs.items():
        change_manningsn(rnames, nvals, item['Points'])

    outQ = slopearea.calculateQ(xsecs)
    Qs.append(outQ['Q'])

FS = pd.Series(Qs, name='Discharge')

FS.to_csv(r'C:\Users\CNB968\OneDrive - MT\Conferences and Meetings\2023\AWRA\Lidar_SlopeArea.csv')

import matplotlib.pyplot as plt

XS_min = []
XS_dist = []
XS_WSE = []
for key, item in surv.extracted_xsections.items():
    min_e = item['Points']['value'].min()
    dist = item['Xsec_Cntrln_Dist']
    WSE = item['WSE_Intersects'][0][1]
    XS_min.append(min_e)
    XS_dist.append(dist)
    XS_WSE.append(WSE)

Prfl_dist = surv.channel_centerline['0']['Points']['Dwnstrm_Distance']
Prfl_el = surv.channel_centerline['0']['Points']['value']

ax = plt.axes()
ax.plot(XS_dist, XS_min, 'k-')
ax.plot(XS_dist, XS_WSE, 'b--')
ax.vlines(XS_dist, XS_min, XS_WSE)
ax.plot(Prfl_dist, Prfl_el, 'k--')
ax.set_xlabel("Distance (ft)", labelpad=15, fontsize=16)
ax.set_ylabel("Elevation (ft)", labelpad=10, fontsize=16)
#ax.set_yticklabels(ax.get_yticks(), size=14)
#ax.set_xticklabels(ax.get_xticks(), size=14)
plt.tight_layout()
plt.savefig(r'C:\Users\CNB968\OneDrive - MT\Conferences and Meetings\2023\AWRA\\Figures\SlopeAreaProfile.png')

ax = plt.axes()
for key, item in surv.extracted_xsections.items():
    ax.plot(item['Points']['Dwnstrm_Distance'], item['Points']['value'], 'k-')
    ax.plot([item['WSE_Intersects'][0][0], item['WSE_Intersects'][-1][0]], [item['WSE_Intersects'][0][1], item['WSE_Intersects'][-1][1]], 'b--')
ax.set_xlabel("Distance (ft)", labelpad=15, fontsize=16)
ax.set_ylabel("Elevation (ft)", labelpad=10, fontsize=16)
#ax.set_yticklabels(ax.get_yticks(), size=14)
#ax.set_xticklabels(ax.get_xticks(), size=14)
plt.tight_layout()


man_r = change_manningsn(rnames, nvals, man_r)

xsec_dict = utilities.points_along_lines(d_xsec_pth, ID_field='XSID', keep_attrs=['XSID'])
for key, item in xsec_dict.items():
    updt_pnts = utilities.sample_raster_at_points(item['Points'], topo_pth)
    item['Points'] = updt_pnts

comb = gpd.sjoin(xsec_dict['1']['Points'], man_r.to_crs(xsec_dict['1']['Points'].crs))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    pass

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
