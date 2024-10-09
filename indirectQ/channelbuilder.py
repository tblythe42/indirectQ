import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio

from indirectQ import utilities


class LoadSurvey:
    """
    Class to hold survey data, either directly measured points/cross sections
    or 3D topography in the form of a DEM or DSM

    Takes filepath for a .csv of survey data for cross section survey data.

    Takes filepath for a .csv of water surface profile points.

    At a minimum, .csv should be formatted to have an X, Y, Z, and Code
    columns. Code should identify each point as part of a cross section or
    HWM.


    Takes geotiff DEM or elevation raster for 3D topography data

    """
    def __init__(self, chn_cnt_pth, ws_pth, ws_crs, srvy_pth=None, topo_pth=None):

        self._topo_filepath = topo_pth
        #self._survey_points_filepath = srvy_pth
        #self._wsprofile_filepath = ws_pth
        #self._chnl_centerline_filepath = chn_cnt_pth
        self.mannings_regions = None
        self.channel_centerline = self._load_channel_centerline(chn_cnt_pth)

        # if srvy_pth is None:
        #     self.xsections = None
        # else:
        #     self.xsections = self.load_Xsections(srvy_pth)

        # TODO - remove this conditional and use exceptions in the function to either return topo array and path or None
        if topo_pth is None:
            self.terrain = None
        else:
            self._topo_filepath = topo_pth
            self.terrain = self.load_terrain(topo_pth)

        self.extracted_xsections = None
        self.WSpoints = self._load_watersurfacepoints(ws_pth, ws_crs)

    # FUTURE DEVELOPMENT - implement method to load direct cross-section survey data coded with attributes for cross-
    # ID and Manning's n for each point.
    #def load_Xsections(self, srvy_xsec_pth):

    # TODO - add exceptions here, user may load survey data and not terrain
    def load_terrain(self, topo_flpth):
        """
        Loads a raster elevation dataset
        :param topo_flpth: file path to .tif DEM
        :return: multi-dimensional numpy array with X, Y, Z values
        """
        with rasterio.open(topo_flpth, 'r') as tp_dset:
            topoz = tp_dset.read(1)
            height, width = topoz.shape
            cols, rows = np.meshgrid(np.arange(width), np.arange(height))
            xs, ys = rasterio.transform.xy(tp_dset.transform, rows, cols)
            topox = np.array(xs)
            topoy = np.array(ys)

        return np.stack((topox, topoy, topoz), axis=0)


    def _load_watersurfacepoints(self, ws_path, ws_crs):
        """
        Loads a geodataframe of surveyed water surface points with X, Y, Z coordinate fields and a CODE field for
        identifying what the point is. A coordinate reference system is needed with the survey data.
        :param ws_path: Path to water surface survey point .csv (e.g., High Water Marks).
        :param chn_cnt: Line Geometry for channel centerline (provided via load centerline or calculated).
        :return: Geodataframe of the imported points from the .csv
        """
        df = pd.read_csv(ws_path)
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.X, df.Y), crs=ws_crs)
        return gdf

    def _load_channel_centerline(self, cntrln_pth):
        """

        :param cntrln_pth: file path to centerline shapefile, sets crs for rest of imports
        :return: Geodataframe of channel centerline.
        """
        # if cntrln_pth is None:
        #     print("No filepath to centerline explicitly provided to function, using default class filepath.")
        #     raise ValueError("LoadSurvey is missing a channel centerline filepath argument.")
        chn_cnt = utilities.points_along_lines(cntrln_pth)
        for key, item in chn_cnt.items():
            updt_pnts = utilities.sample_raster_at_points(item['Points'], self._topo_filepath)
            item['Points'] = updt_pnts

        return chn_cnt


    # FUTURE - reproject all input data to same crs regardless of what it's in
    # def _datareproject(self): --- probably move to utilities


    # def _create_centerline(self): --- create a centerline from flow direction/routing of DEM (3D) option or from
    # surveyed cross section low point orthogonal vectors.

    def extractXsections(self, dxsec_pth, **kwargs):
        xsec_dict = utilities.points_along_lines(dxsec_pth, )
        for key, item in xsec_dict.items():
            updt_pnts = utilities.sample_raster_at_points(item['Points'], self._topo_filepath)
            item['Points'] = updt_pnts
        final_xsec_dict = utilities.sort_cross_sections(self.channel_centerline, xsec_dict, **kwargs)
        self.extracted_xsections = final_xsec_dict

    def assign_Mannings(self, mannings_shp):
        """
        Takes path to shapefile, loads it, and updates the cross section dictionary. Manning's shapefile should
        contain an n_Value attribute.
        :param mannings_shp: shapefile of regions with unique Manning's n values needs attributes "n_Value" and "Name"
        for each region
        :return: None, updates dictionary
        """
        man_r = gpd.read_file(mannings_shp)
        self.mannings_regions = man_r
        if self.extracted_xsections is None:
            print("No cross sections created.")
        else:
            for key, item in self.extracted_xsections.items():
                comb = gpd.sjoin(item['Points'], man_r.to_crs(item['Points'].crs))
                self.extracted_xsections[key]['Points'] = comb

    def assign_wse(self):
        wse_srtd = utilities.sort_points(self.channel_centerline, self.WSpoints)
        wses = utilities.interpz_at_intersection(self.WSpoints, self.extracted_xsections)
        self.extracted_xsections = wses


    #def plot_survey(self, plot_engine='plotly'):

    #def plot_survey3D(self, plot_engine='plotly'):

    #def plot_terrain(self, plot_engine='plotly'):

    #def plot_Xsections(self, plot_engine='plotly'):

    #def plot_profile(self, plot_engine='plotly'):