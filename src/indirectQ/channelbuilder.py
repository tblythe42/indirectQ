"""Constructs a stream channel represented by centerline and cross-sections

The channelbuilder module is used to set up an instance of a stream channel represented by a
channel centerline, cross-sections, and water surface profiles (all attributed with surveyed elevation data or
elevation data extracted from raster topography data). The ChannelBuilder instance is used as input to the
other indirect discharge calculation methods.

This module contains the following classes:
- Profile
- ChannelBuilder

-  **Author(s):** Todd Blythe, MTDNRC, CNB968
-  **Date:** Created 12/10/2025
"""
from pathlib import Path

import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio
from pyproj.crs.crs import CRS
from shapely.geometry import LineString
from affine import Affine

from indirectQ import utilities
from indirectQ.profile import Profile


class ChannelBuilder:
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
    def __init__(self,
                 centerline: str | Path | LineString | pd.DataFrame | gpd.GeoDataFrame | Profile,
                 centerline_xname: str = 'X',
                 centerline_yname: str = 'Y',
                 centerline_zname: str = 'Z',
                 crs: str | int | CRS | None = None,
                 topo: str | Path | tuple[np.ndarray, dict] | None = None,
                 topo_scalar: float | None = None):

        if isinstance(crs, (str, int)):
            self.crs = CRS.from_epsg(crs)
        else:
            self.crs = crs
        self._cntrln_cols = {
            'X': centerline_xname,
            'Y': centerline_yname,
            'elevation': centerline_zname
        }
        if topo_scalar is None:
            self._topo_metadata = {'scale': 1.0}
        else:
            self._topo_metadata = {'scale': topo_scalar}
        self._topo = None
        self._transform = None
        self._centerline = None
        self._watsurf = None
        # mannings can be added in bulk for whole cross-section, using regions (polygons), and using distance ranges
        # this attribute should be pandas df created from distance column in cross-section points attribute but with a
        # second column of Manning's n for each distance.
        self._mannings_regions = None
        self._mannings = None
        self._xsections = None

        self.terrain = topo
        self.channel_centerline = centerline

    @property
    def channel_centerline(self):
        """The channel's centerline."""
        return self._centerline

    @channel_centerline.setter
    def channel_centerline(self, chan_lin: str | Path | LineString | pd.DataFrame | gpd.GeoDataFrame | Profile) -> None:
        """Set's the channel centerline geometry."""
        if isinstance(chan_lin, (str, Path)):
            cntrln = gpd.read_file(chan_lin)
            if (isinstance(cntrln, pd.DataFrame)) & (not isinstance(cntrln, gpd.GeoDataFrame)):
                cpro = Profile.load_from_survey(cntrln, self._cntrln_cols['X'], self._cntrln_cols['Y'],
                                                self._cntrln_cols['elevation'], self.crs, 'centerline')
            else:
                if (cntrln.geom_type == 'Point').all():
                    cpro = Profile.load_from_survey(cntrln, zcoord_name=self._cntrln_cols['elevation'],
                                                    crs=self.crs, profile_name='centerline')
                elif (cntrln.geom_type == 'LineString').all():
                    if self._topo is None:
                        raise AttributeError("Topography raster data has not been loaded, set a topography dataset "
                                             "before continuing with the channel centerline.")
                    if self._transform is None:
                        raise AttributeError("A transform for topography data has not been set, set the transform "
                                             "property before continuing with the channel centerline.")

                    cpro = Profile.extract_from_raster((self._topo, self._topo_metadata), cntrln,
                                                       data_scalar=self._topo_metadata['scale'],
                                                       profile_name='centerline')
                else:
                    raise ValueError("The GeoDataFrame contains non-accepted geometry types.")
        elif isinstance(chan_lin, LineString):
            if self._topo is None:
                raise AttributeError("Topography raster data has not been loaded, set a topography dataset before "
                                     "continuing with the channel centerline.")
            if self._transform is None:
                raise AttributeError("A transform for topography data has not been set, set the transform property "
                                     "before continuing with the channel centerline.")

            cpro = Profile.extract_from_raster((self._topo, self._topo_metadata), chan_lin,
                                               data_scalar=self._topo_metadata['scale'],
                                               profile_name='centerline')
        elif (isinstance(chan_lin, pd.DataFrame)) & (not isinstance(chan_lin, gpd.GeoDataFrame)):
            cpro = Profile.load_from_survey(chan_lin, self._cntrln_cols['X'], self._cntrln_cols['Y'],
                                            self._cntrln_cols['elevation'], self.crs, 'centerline')
        elif isinstance(chan_lin, gpd.GeoDataFrame):
            if (chan_lin.geom_type == 'Point').all():
                cpro = Profile.load_from_survey(chan_lin, zcoord_name=self._cntrln_cols['elevation'],
                                                crs=self.crs, profile_name='centerline')
            elif (chan_lin.geom_type == 'LineString').all():
                if self._topo is None:
                    raise AttributeError("Topography raster data has not been loaded, set a topography dataset before "
                                         "continuing with the channel centerline.")
                if self._transform is None:
                    raise AttributeError("A transform for topography data has not been set, set the transform property "
                                         " before continuing with the channel centerline.")

                cpro = Profile.extract_from_raster((self._topo, self._topo_metadata), chan_lin,
                                                   data_scalar=self._topo_metadata['scale'],
                                                   profile_name='centerline')
            else:
                raise ValueError("The GeoDataFrame contains non-accepted geometry types.")
        elif isinstance(chan_lin, Profile):
            cpro = chan_lin
        else:
            raise TypeError("The input channel centerline is not an accepted type.")

        self._centerline = cpro

    @property
    def cross_sections(self) -> dict:
        """Cross-section profiles and attributes."""
        return self._xsections

    @property
    def terrain(self):
        """The topography raster dataset."""
        return self._topo

    @terrain.setter
    def terrain(self, rast_data: str | Path | tuple[np.ndarray, dict]) -> None:
        """Set the topography data."""
        if isinstance(rast_data, (str, Path)):
            d_arr, d_met = utilities.raster_array_from_file(rast_data, data_scalar=self._topo_metadata['scale'])
            self._topo = d_arr
            self._topo_metadata.update(d_met)
            self._transform = d_met['transform']
            if self.crs is None:
                self.crs = CRS.from_wkt(d_met['crs'].to_wkt())
            else:
                if self.crs.to_epsg() != d_met['crs'].to_epsg():
                    print("WARNING: The coordinate reference system of the input raster dataset does not match"
                          "the existing CRS for the ChannelBuilding instance.")
        else:
            if rast_data[0].shape[0] != 3:
                raise ValueError("The input array does not have axis 0 dimension equal to 3. Raster arrays must have"
                                 "axis 0 = 3 where the first position is xcoords, then ycoords, then z-values.")

            self._topo_metadata.update(rast_data[1])
            self._transform = self._topo_metadata['transform']
            self._topo = rast_data[0]

    @property
    def transform(self):
        """The affine transform for topography raster data."""
        return self._transform

    @transform.setter
    def transform(self, affn: Affine):
        """Sets the affine transform for topography data."""
        if not isinstance(affn, Affine):
            raise TypeError("The transform must be an Affine object.")

        self._transform = affn

    @property
    def water_profiles(self) -> dict:
        """Water surface profiles for the channel."""
        return self._watsurf

    @property
    def mannings(self) -> dict:
        """Manning's n values for cross-sections."""
        return self._mannings

    @mannings.setter
    def mannings(self, mannings_dict: dict) -> None:
        """Set Manning's n by cross-section"""
        for xsid, val in mannings_dict.items():
            if xsid not in list(self._xsections.keys()):
                print(f"Cross section {xsid} was not found in valid cross-sections.")
                continue
            if isinstance(val, float):
                nval = np.repeat(val, len(self._xsections[xsid]['Profile'].points))
            elif isinstance(val, (list, np.ndarray)):
                if len(val) != len(self._xsections[xsid]['Profile'].points):
                    raise ValueError("The length of input Manning's values do not match the length "
                                     "of the cross-section points.")
                nval = val
            else:
                raise ValueError("The input Manning's n values are not a single float, or list-like set of values.")

            if self._mannings is None:
                self._mannings = {xsid: nval}
            else:
                self._mannings[xsid] = nval

    def add_xsection(self,
                     xs_data: str | Path | LineString | pd.DataFrame | gpd.GeoDataFrame | Profile,
                     xs_xname: str = 'X',
                     xs_yname: str = 'Y',
                     xs_zname: str = 'elevation',
                     xs_id: str | int | None = None,
                     invert_centerline: bool = False) -> None:
        """

        Assumes the input cross-section data is in the same CRS as the ChannelBuilder instance.

        Args:
            xs_data:
            xs_xname:
            xs_yname:
            xs_zname:
            xs_id:
            invert_centerline:

        Returns:

        """
        if isinstance(xs_data, (str, Path)):
            xs = gpd.read_file(xs_data)
            if (isinstance(xs, pd.DataFrame)) & (not isinstance(xs, gpd.GeoDataFrame)):
                xspro = Profile.load_from_survey(xs, xs_xname, xs_yname,
                                                 xs_zname, self.crs, 'cross_section')
            else:
                if (xs.geom_type == 'Point').all():
                    xspro = Profile.load_from_survey(xs, zcoord_name=xs_zname,
                                                     crs=self.crs, profile_name='cross_section')
                elif (xs.geom_type == 'LineString').all():
                    if self._topo is None:
                        raise AttributeError("Topography raster data has not been loaded, set a topography dataset "
                                             "before continuing with the channel centerline.")
                    if self._transform is None:
                        raise AttributeError("A transform for topography data has not been set, set the transform "
                                             "property before continuing with the channel centerline.")

                    xspro = Profile.extract_from_raster((self._topo, self._topo_metadata), xs,
                                                        profile_name='cross_section')
                else:
                    raise ValueError("The GeoDataFrame contains non-accepted geometry types.")
        elif isinstance(xs_data, LineString):
            if self._topo is None:
                raise AttributeError("Topography raster data has not been loaded, set a topography dataset before "
                                     "continuing with the channel centerline.")
            if self._transform is None:
                raise AttributeError("A transform for topography data has not been set, set the transform property "
                                     "before continuing with the channel centerline.")

            xspro = Profile.extract_from_raster((self._topo, self._topo_metadata), xs_data,
                                                profile_name='cross_section')
        elif (isinstance(xs_data, pd.DataFrame)) & (not isinstance(xs_data, gpd.GeoDataFrame)):
            xspro = Profile.load_from_survey(xs_data, xs_xname, xs_yname,
                                             xs_zname, self.crs, 'cross_section')
        elif isinstance(xs_data, gpd.GeoDataFrame):
            if (xs_data.geom_type == 'Point').all():
                xspro = Profile.load_from_survey(xs_data, zcoord_name=xs_zname,
                                                 crs=self.crs, profile_name='cross_section')
            elif (xs_data.geom_type == 'LineString').all():
                if self._topo is None:
                    raise AttributeError("Topography raster data has not been loaded, set a topography dataset before "
                                         "continuing with the channel centerline.")
                if self._transform is None:
                    raise AttributeError("A transform for topography data has not been set, set the transform property "
                                         " before continuing with the channel centerline.")

                xspro = Profile.extract_from_raster((self._topo, self._topo_metadata), xs_data,
                                                    profile_name='cross_section')
            else:
                raise ValueError("The GeoDataFrame contains non-accepted geometry types.")
        elif isinstance(xs_data, Profile):
            xspro = xs_data
        else:
            raise TypeError("The input cross-section is not an accepted type.")

        if self._xsections is None:
            if xs_id is None:
                xs_id = 0

            xsdict = {
                xs_id: {
                    'Profile': xspro,
                    'Channel Rank': None,
                    'Xsec_Cntrln_Dist': None,
                    'WSE_Intersects': None,
                    'Mean_WSE': None
                }
            }

            self._xsections = xsdict

        else:
            if xs_id is None:
                tmp_id = 0
                while tmp_id in list(self._xsections.keys()):
                    tmp_id += 1
                xs_id = tmp_id
            else:
                if xs_id in list(self._xsections.keys()):
                    raise ValueError("The cross-section ID given already exists as a cross-section. Please choose"
                                     " a unique string or integer.")

            xsdict = {
                xs_id: {
                    'Profile': xspro,
                    'Channel Rank': None,
                    'Xsec_Cntrln_Dist': None,
                    'WSE_Intersects': None,
                    'Mean_WSE': None
                }
            }

            self._xsections.update(xsdict)

        # Sort XSections
        self._xsections = utilities.sort_cross_sections(self._centerline, self._xsections,
                                                        invert_line=invert_centerline)
        # Interp WSE
        if self._watsurf is not None:
            self._xsections = utilities.interpz_at_intersection(self._watsurf, self._xsections)

    def add_watersurface(self,
                         ws_pnts: str | Path | pd.DataFrame | gpd.GeoDataFrame,
                         ws_xname: str = 'X',
                         ws_yname: str = 'Y',
                         ws_zname: str = 'elevation',
                         reverse_centerline: bool = False,
                         ws_id: str | int | None = None) -> None:
        """

        Args:
            ws_pnts:
            ws_xname:
            ws_yname:
            ws_zname:
            reverse_centerline:
            ws_id:

        Returns:

        """
        if isinstance(ws_pnts, (str, Path)):
            ws = gpd.read_file(ws_pnts)
            if (isinstance(ws, pd.DataFrame)) & (not isinstance(ws, gpd.GeoDataFrame)):
                wsgdf = gpd.GeoDataFrame(ws, geometry=gpd.points_from_xy(ws[ws_xname], ws[ws_yname]), crs=self.crs)
                srt_pnts = utilities.sort_points(self._centerline, wsgdf, invert_line=reverse_centerline)
                ws_r_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'R'].sort_values(by='Cntrline_Dist')
                ws_l_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'L'].sort_values(by='Cntrline_Dist')
                if not ws_r_pnts.empty:
                    ws_r_cdists = ws_r_pnts['Cntrline_Dist'].values
                    ws_r_pro = Profile.load_from_survey(ws_r_pnts, ws_xname, ws_yname,
                                                        ws_zname, self.crs, 'water_surface')
                else:
                    ws_r_cdists = None
                    ws_r_pro = None
                if not ws_l_pnts.empty:
                    ws_l_cdists = ws_l_pnts['Cntrline_Dist'].values
                    ws_l_pro = Profile.load_from_survey(ws_l_pnts, ws_xname, ws_yname, ws_zname, self.crs,
                                                        'water_surface')
                else:
                    ws_l_cdists = None
                    ws_l_pro = None
            else:
                if (ws.geom_type == 'Point').all():
                    srt_pnts = utilities.sort_points(self._centerline, ws, invert_line=reverse_centerline)
                    ws_r_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'R'].sort_values(by='Cntrline_Dist')
                    ws_l_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'L'].sort_values(by='Cntrline_Dist')
                    if not ws_r_pnts.empty:
                        ws_r_cdists = ws_r_pnts['Cntrline_Dist'].values
                        ws_r_pro = Profile.load_from_survey(ws_r_pnts, zcoord_name=ws_zname,
                                                            crs=self.crs, profile_name='water_surface')
                    else:
                        ws_r_cdists = None
                        ws_r_pro = None
                    if not ws_l_pnts.empty:
                        ws_l_cdists = ws_l_pnts['Cntrline_Dist'].values
                        ws_l_pro = Profile.load_from_survey(ws_l_pnts, zcoord_name=ws_zname,
                                                            crs=self.crs, profile_name='water_surface')
                    else:
                        ws_l_cdists = None
                        ws_l_pro = None
                else:
                    raise ValueError("The GeoDataFrame contains non-accepted geometry types. Water surfaces are"
                                     " only initialized from point data.")
        elif (isinstance(ws_pnts, pd.DataFrame)) & (not isinstance(ws_pnts, gpd.GeoDataFrame)):
            wsgdf = gpd.GeoDataFrame(ws_pnts, geometry=gpd.points_from_xy(ws_pnts[ws_xname], ws_pnts[ws_yname]),
                                     crs=self.crs)
            srt_pnts = utilities.sort_points(self._centerline, wsgdf, invert_line=reverse_centerline)
            ws_r_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'R'].sort_values(by='Cntrline_Dist')
            ws_l_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'L'].sort_values(by='Cntrline_Dist')
            if not ws_r_pnts.empty:
                ws_r_cdists = ws_r_pnts['Cntrline_Dist'].values
                ws_r_pro = Profile.load_from_survey(ws_r_pnts, zcoord_name=ws_zname,
                                                    crs=self.crs, profile_name='water_surface')
            else:
                ws_r_cdists = None
                ws_r_pro = None
            if not ws_l_pnts.empty:
                ws_l_cdists = ws_l_pnts['Cntrline_Dist'].values
                ws_l_pro = Profile.load_from_survey(ws_l_pnts, zcoord_name=ws_zname,
                                                    crs=self.crs, profile_name='water_surface')
            else:
                ws_l_cdists = None
                ws_l_pro = None
        elif isinstance(ws_pnts, gpd.GeoDataFrame):
            if (ws_pnts.geom_type == 'Point').all():
                srt_pnts = utilities.sort_points(self._centerline, ws_pnts, invert_line=reverse_centerline)
                ws_r_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'R'].sort_values(by='Cntrline_Dist')
                ws_l_pnts = srt_pnts.loc[srt_pnts['Waters_Edge'] == 'L'].sort_values(by='Cntrline_Dist')
                if not ws_r_pnts.empty:
                    ws_r_cdists = ws_r_pnts['Cntrline_Dist'].values
                    ws_r_pro = Profile.load_from_survey(ws_r_pnts, zcoord_name=ws_zname,
                                                        crs=self.crs, profile_name='water_surface')
                else:
                    ws_r_cdists = None
                    ws_r_pro = None
                if not ws_l_pnts.empty:
                    ws_l_cdists = ws_l_pnts['Cntrline_Dist'].values
                    ws_l_pro = Profile.load_from_survey(ws_l_pnts, zcoord_name=ws_zname,
                                                        crs=self.crs, profile_name='water_surface')
                else:
                    ws_l_cdists = None
                    ws_l_pro = None
            else:
                raise ValueError("The GeoDataFrame contains non-accepted geometry types. Water surfaces are"
                                 " only initialized from point data.")
        else:
            raise TypeError("The input water surface data is not an accepted type.")

        if self._watsurf is None:
            if ws_id is None:
                ws_id = 0

            wsdict = {
                ws_id: {
                    'L': {
                        'Profile': ws_l_pro,
                        'Centerline_Dist': ws_l_cdists
                    },
                    'R': {
                        'Profile': ws_r_pro,
                        'Centerline_Dist': ws_r_cdists
                    }
                }
            }

            self._watsurf = wsdict

        else:
            if ws_id is None:
                tmp_id = 0
                while tmp_id in list(self._watsurf.keys()):
                    tmp_id += 1
                ws_id = tmp_id
            else:
                if ws_id in list(self._watsurf.keys()):
                    raise ValueError("The cross-section ID given already exists as water profile. Please choose"
                                     " a unique string or integer.")

            wsdict = {
                ws_id: {
                    'L': {
                        'Profile': ws_l_pro,
                        'Centerline_Dist': ws_l_cdists
                    },
                    'R': {
                        'Profile': ws_r_pro,
                        'Centerline_Dist': ws_r_cdists
                    }
                }
            }

            self._watsurf.update(wsdict)

        if self._xsections is not None:
            self._xsections = utilities.interpz_at_intersection(self._watsurf, self._xsections)

    def mannings_from_regions(self,
                              mannings_shp: str | Path | gpd.GeoDataFrame,
                              mannings_name: str = 'n_value') -> None:
        """Assigns Manning's n by intersecting polygons.

        A layer of Manning's n regions (represented as polygons) is intersected with cross-sections and associated
        Manning's value is assigned by this intersection.

        Args:
            mannings_shp:
                A file path or GeoDataFrame with polygon geometries and column/attribute with Manning's n values.
            mannings_name:
                A name of the column or attribute representing the Manning's n values.
        """
        if isinstance(mannings_shp, (str, Path)):
            man_r = gpd.read_file(mannings_shp)
        elif isinstance(mannings_shp, gpd.GeoDataFrame):
            man_r = mannings_shp
        else:
            raise TypeError("The Manning's regions are not a recognized geospatial type.")

        if self.crs is None:
            print("WARNING: no CRS has been set for the ChannelBuilder instance, Manning's regions are not allowed"
                  " to set the class CRS. Intersection will continue but may be incorrect.")
        else:
            if self.crs.to_epsg() != man_r.crs.to_epsg():
                man_r = man_r.to_crs(self.crs)

        if not (man_r.geom_type == 'Polgyon').all():
            raise TypeError("Manning's regions must all be Polygons.")

        self._mannings_regions = man_r
        if self._xsections is None:
            print("No cross sections created.")
        else:
            for key, item in self._xsections.items():
                comb = gpd.sjoin(item['Profile'].points, man_r)
                self.mannings = {key, comb[mannings_name].values}

    def mannings_by_distance(self,
                             xs_id: str | int,
                             dist_breakpoints: list | np.ndarray,
                             mannings_vals: list | np.ndarray) -> None:
        """Sets Manning's n values by distance for a cross-section.

        Distance breakpoints define where there will be changes in Manning's n values, so, there always needs to be
        n+1 Manning's values for n breakpoints. This is enforced by the function.

        Args:
            xs_id:
            dist_breakpoints:
                Breakpoints where there will be a change in Manning's n value.
            mannings_vals:
                The Manning's n values for each division of the cross-section (determined by the breakpoints).
                mannings_vals must have a length of (n + 1) where n is the length of dist_breakpoints.
        """
        pnts = self._xsections[xs_id]['Profile'].points['distance'].values
        nvals = np.empty(pnts.shape)
        nb = len(dist_breakpoints)
        if nb == 1:
            idxst = np.where(pnts <= dist_breakpoints[0])
            nvals[idxst] = mannings_vals[0]
            idxen = np.where(pnts > dist_breakpoints[0])
            nvals[idxen] = mannings_vals[1]
        elif nb == 2:
            idxst = np.where(pnts <= dist_breakpoints[0])
            nvals[idxst] = mannings_vals[0]
            idxmid = np.where((pnts > dist_breakpoints[0]) & (pnts <= dist_breakpoints[1]))
            nvals[idxmid] = mannings_vals[1]
            idxen = np.where(pnts > dist_breakpoints[1])
            nvals[idxen] = mannings_vals[2]
        else:
            for i in range(nb):
                if i == 0:
                    idxs = np.where(pnts <= dist_breakpoints[i])
                    nvals[idxs] = mannings_vals[i]
                    idxnxt = np.where((pnts > dist_breakpoints[i]) & (pnts <= dist_breakpoints[i + 1]))
                    nvals[idxnxt] = mannings_vals[i+1]
                elif i == (nb - 1):
                    idxs = np.where(pnts > dist_breakpoints[i])
                    nvals[idxs] = mannings_vals[i+1]
                else:
                    idxs = np.where((pnts > dist_breakpoints[i]) & (pnts <= dist_breakpoints[i+1]))
                    nvals[idxs] = mannings_vals[i+1]

        self.mannings = {xs_id: nvals}


    #def plot_survey(self, plot_engine='plotly'):

    #def plot_survey3D(self, plot_engine='plotly'):

    #def plot_terrain(self, plot_engine='plotly'):

    #def plot_Xsections(self, plot_engine='plotly'):

    #def plot_profile(self, plot_engine='plotly'):

