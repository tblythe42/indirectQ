import pandas as pd
import geopandas as gpd
import rasterio

import utilities

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
    def __init__(self, srvy_pth = None, topo_pth = None, ws_pth = None, chn_cnt_pth):

        self._topo_filepath = topo_pth
        self._survey_points_filepath = srvy_pth
        self._wsprofile_filepath = ws_pth
        self._chnl_centerline_filepath = chn_cnt_pth

        args = [srvy_pth, topo_pth, ws_pth]
        if all(i is None for i in args):
            raise ValueError("No input data files were specified.")


        if srvy_pth is None:
            self.xsections = None
        else:
            self.xsections = self.load_Xsections(srvy_pth)


        if topo_pth is None:
            self.terrain = None
        else:
            self.terrain = self.load_terrain(topo_pth)

        self.extracted_xsections = None
        self.channel_centerline = None
        self.WS_profile =

    # FUTURE DEVELOPMENT - implement method to load direct cross-section survey data coded with attributes for cross-
    # ID and Manning's n for each point.
    #def load_Xsections(self, srvy_xsec_pth):

    def load_terrain(self, topo_pth):
        """
        Loads a raster elevation dataset
        :param topo_pth:
        :return: multi-dimensional numpy array with X, Y, Z values
        """
        with rasterio.open(topo_pth, 'r') as tp_dset:
            topoZ = tp_dset.read(1)
            topoX =
            topoY =

        return topo_array


    def load_watersurfacepoints(self, ws_path, chn_cnt, crs):
        """
        Loads a geodataframe of surveyed water surface points with X, Y, Z coordinate fields and a CODE field for
        identifying what the point is. A coordinate reference system is needed with the survey data.
        :param ws_path: Path to water surface survey point .csv (e.g., High Water Marks).
        :param chn_cnt: Line Geometry for channel centerline (provided via load centerline or calculated).
        :return: Geodataframe of the imported points from the .csv
        """
        df = pd.read_csv(ws_path)
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.X, df.Y), crs=crs)
        return gdf

    def load_channel_centerline(self, cntrln_pth = self._chnl_centerline_filepath):
        """

        :param cntrln_pth:
        :return: Geodataframe of channel centerline.
        """
        if cntrln_pth is None:
            print("No filepath to centerline explicitly provided to function, using default class filepath.")
            raise ValueError("LoadSurvey is missing a channel centerline filepath argument.")
        chn_cnt = gpd.read_file(chan_cnt_pth)
        return chn_cnt


    # FUTURE - reproject all input data to same crs regardless of what it's in
    # def _datareproject(self): --- probably move to utilities


    # def _create_centerline(self): --- create a centerline from flow direction/routing of DEM (3D) option or from
    # surveyed cross section low point orthogonal vectors.

    def extractXsections(self, chnl_cntln_pth, drwn_xsec_pth, raster_pth=topo_pth, spacing=1):

    def plot_survey(self, plot_engine='plotly'):

    def plot_survey3D(self, plot_engine='plotly'):

    def plot_terrain(self, plot_engine='plotly'):

    def plot_Xsections(self, plot_engine='plotly'):

    def plot_profile(self, plot_engine='plotly'):