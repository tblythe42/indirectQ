"""Creates profile data for a 3-dimensional line

This module contains the following classes:
- Profile

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

from indirectQ import utilities


class Profile:
    """Stores geometry and elevation profile data for a line.

    Attributes:
        name:
            A string identifier for the profile line.
        crs:
            The coordinate reference system information for the line's geometry and profile points, as a
            pyproj CRS object.
    """

    def __init__(self,
                 geom_line: str | Path | LineString | gpd.GeoSeries | gpd.GeoDataFrame,
                 points_df: str | Path | pd.DataFrame | gpd.GeoDataFrame,
                 xcoord_name: str = 'X',
                 ycoord_name: str = 'Y',
                 zcoord_name: str = 'elevation',
                 prof_dist_name: str = 'distance',
                 crs: CRS | None = None,
                 profile_name: str | None = None
                 ):

        self.name = profile_name
        self.crs = crs
        self._line = None
        self._points = None
        self._cmap = {
            'X': xcoord_name,
            'Y': ycoord_name,
            'elevation': zcoord_name,
            'distance': prof_dist_name
        }

        self.line = geom_line
        self.points = points_df

    @property
    def line(self) -> gpd.GeoSeries:
        """The profile line geometry."""
        return self._line

    @line.setter
    def line(self, in_ln: str | Path | LineString | gpd.GeoSeries | gpd.GeoDataFrame) -> None:
        """Set the line geometries of the profiles."""
        if isinstance(in_ln, (str, Path)):
            lnfl = gpd.read_file(in_ln)
            if self.crs is None:
                self.crs = lnfl.crs
            else:
                if self.crs.to_epsg() != lnfl.crs.to_epsg():
                    lnfl = lnfl.to_crs(self.crs)
            if len(lnfl) > 1:
                print("WARNING: More than one line in the input file, will default to using only the first feature.")
            ln = lnfl.geometry.iloc[0]
        elif isinstance(in_ln, LineString):
            if self.crs is None:
                print("WARNING: No CRS is set for the profile. Make sure the input LineString coordinates"
                      "match the profile point coordinates.")
            ln = in_ln
        elif isinstance(in_ln, gpd.GeoSeries):
            if self.crs is None:
                print("WARNING: No CRS is set for the profile. Make sure the input LineString coordinates"
                      "match the profile point coordinates.")
            if len(in_ln) > 1:
                print("WARNING: More than one line in the geometry Series, only using the first feature.")
            ln = in_ln.iloc[0]
        elif isinstance(in_ln, gpd.GeoDataFrame):
            if self.crs is None:
                self.crs = in_ln.crs
            else:
                if self.crs.to_epsg() != in_ln.crs.to_epsg():
                    in_ln = in_ln.to_crs(self.crs)
            if len(in_ln) > 1:
                print("WARNING: More than one line in the GeoDataFrame, only using the first feature.")
            ln = in_ln.geometry.iloc[0]
        else:
            raise TypeError("The input line is not an accepted type.")

        self._line = ln

    @property
    def points(self) -> gpd.GeoDataFrame:
        """Elevation data at coordinate points along each profile."""
        return self._points

    @points.setter
    def points(self, newpnts: str | Path | pd.DataFrame | gpd.GeoDataFrame) -> None:
        """Define the dataframe of points for each profile."""
        if isinstance(newpnts, (str, Path)):
            gdf = gpd.read_file(newpnts)
            if (isinstance(gdf, pd.DataFrame)) & (isinstance(gdf, gpd.GeoDataFrame)):
                if self.crs is None:
                    self.crs = gdf.crs
                else:
                    if self.crs.to_epsg() != gdf.crs.to_epsg():
                        gdf = gdf.to_crs(self.crs)

                default_cols = ['distance', 'elevation']
                arg_cols = [self._cmap[default_cols[i]] for i in range(len(default_cols))]

                if all([d in gdf.columns.to_list() for d in default_cols]):
                    zname = 'elevation'
                    dist_name = 'distance'
                elif all([a in gdf.columns.to_list()] for a in arg_cols):
                    zname = self._cmap['elevation']
                    dist_name = self._cmap['distance']
                else:
                    raise KeyError("Columns for xcoords, ycoords, elevation, and distance along profile were not"
                                   "found in dataframe using default names or given names.")

                gdf = gdf.loc[:, [dist_name, zname, 'geometry']]
                gdf = gdf.rename(columns={self._cmap['distance']: 'distance', self._cmap['elevation']: 'elevation'})
            else:
                gdf = gdf.apply(pd.to_numeric, errors='coerce').fillna(gdf)
                default_cols = list(self._cmap.keys())
                arg_cols = [self._cmap[default_cols[i]] for i in range(len(default_cols))]
                # check columns here then index
                if all([d in gdf.columns.to_list() for d in default_cols]):
                    xname = 'X'
                    yname = 'Y'
                    zname = 'elevation'
                    dist_name = 'distance'
                elif all([a in gdf.columns.to_list()] for a in arg_cols):
                    xname = self._cmap['X']
                    yname = self._cmap['Y']
                    zname = self._cmap['elevation']
                    dist_name = self._cmap['distance']
                else:
                    raise KeyError("Columns for xcoords, ycoords, elevation, and distance along profile were not"
                                   "found in dataframe using default names or given names.")

                if self.crs is None:
                    print("WARNING: No coordinate reference system info is available for tabular data...")

                gdf = gpd.GeoDataFrame(gdf[[dist_name, zname]], geometry=gpd.points_from_xy(gdf[xname], gdf[yname]),
                                       crs=self.crs)
                gdf = gdf.rename(columns={self._cmap['distance']: 'distance', self._cmap['elevation']: 'elevation'})

        elif (not isinstance(newpnts, gpd.GeoDataFrame)) & (isinstance(newpnts, pd.DataFrame)):
            gdf = newpnts.apply(pd.to_numeric, errors='coerce').fillna(newpnts)
            default_cols = list(self._cmap.keys())
            arg_cols = [self._cmap[default_cols[i]] for i in range(len(default_cols))]
            if all([d in gdf.columns.to_list() for d in default_cols]):
                xname = 'X'
                yname = 'Y'
                zname = 'elevation'
                dist_name = 'distance'
            elif all([a in gdf.columns.to_list()] for a in arg_cols):
                xname = self._cmap['X']
                yname = self._cmap['Y']
                zname = self._cmap['elevation']
                dist_name = self._cmap['distance']
            else:
                raise KeyError("Columns for xcoords, ycoords, elevation, and distance along profile were not"
                               "found in dataframe using default names or given names.")

            if self.crs is None:
                print("WARNING: No coordinate reference system info is available for tabular data...")

            gdf = gpd.GeoDataFrame(gdf[[dist_name, zname]], geometry=gpd.points_from_xy(gdf[xname], gdf[yname]),
                                   crs=self.crs)
            gdf = gdf.rename(columns={self._cmap['distance']: 'distance', self._cmap['elevation']: 'elevation'})

        elif isinstance(newpnts, gpd.GeoDataFrame):
            gdf = newpnts
            if self.crs is None:
                self.crs = gdf.crs
            else:
                if self.crs.to_epsg() != gdf.crs.to_epsg():
                    gdf = gdf.to_crs(self.crs)

            default_cols = ['distance', 'elevation']
            arg_cols = [self._cmap[default_cols[i]] for i in range(len(default_cols))]

            if all([d in gdf.columns.to_list() for d in default_cols]):
                zname = 'elevation'
                dist_name = 'distance'
            elif all([a in gdf.columns.to_list()] for a in arg_cols):
                zname = self._cmap['elevation']
                dist_name = self._cmap['distance']
            else:
                raise KeyError("Columns for xcoords, ycoords, elevation, and distance along profile were not"
                               "found in dataframe using default names or given names.")

            gdf = gdf.loc[:, [dist_name, zname, 'geometry']]
            gdf = gdf.rename(columns={self._cmap['distance']: 'distance', self._cmap['elevation']: 'elevation'})

        else:
            raise TypeError("Points setter received an unexpected input type.")

        self._points = gdf

    @staticmethod
    def extract_from_raster(raster_data: str | Path | tuple[np.ndarray, dict],
                            line_geom: str | Path | LineString | gpd.GeoSeries | gpd.GeoDataFrame,
                            point_space: float = 1.0,
                            data_scalar: float = 1.0,
                            profile_name: str | None = None):
        """Extracts elevation data from a raster along a path.

        Args:
            raster_data:
                A path to geotiff or a raster (array, metadata) tuple of topography/DEM.
            line_geom:
                A line geometry that represents the path over which to create the Profile.
            point_space:
                The spacing of sampled points along the line/path.
            data_scalar:
                Conversion factor for raster data (e.g., converting meters to feet).
            profile_name:
                Name of the resulting Profile.

        Returns:
            An initialized Profile class created from raster data and a line geometry.
        """
        if isinstance(raster_data, (str, Path)):
            arr, meta = utilities.raster_array_from_file(raster_data, data_scalar=data_scalar)
        elif isinstance(raster_data, tuple):
            arr, meta = raster_data
        else:
            raise TypeError("The input raster data is not recognized as an accepted input type.")

        if isinstance(line_geom, (str, Path)):
            ln = gpd.read_file(line_geom)
        elif isinstance(line_geom, LineString):
            ln = gpd.GeoDataFrame(geometry=[line_geom])
        elif isinstance(line_geom, gpd.GeoSeries):
            ln = gpd.GeoDataFrame(geometry=line_geom)
        elif isinstance(line_geom, gpd.GeoDataFrame):
            ln = line_geom
        else:
            raise TypeError("The input line is not an accepted type.")

        if isinstance(meta['crs'], rasterio.crs.CRS):
            crs = CRS.from_wkt(meta['crs'].to_wkt())
        elif isinstance(meta['crs'], CRS):
            crs = meta['crs']
        else:
            raise AttributeError("The CRS of the raster data is not neither a rasterio.crs.CRS object nor"
                                 "a pyproj.crs.crs.CRS object.")

        if crs.to_epsg() != ln.crs.to_epsg():
            print("WARNING: raster crs and line geometry crs do not match! Reprojecting line to raster CRS.")
            ln = ln.to_crs(crs)

        pro_pnts = utilities.points_along_lines(ln, dist=point_space)
        if len(ln) > 1:
            print("WARNING: More than one line in the GeoDataFrame, only using the first feature.")
        pro_pnts = pro_pnts[0]

        if 'transform' not in list(meta.keys()):
            raise IndexError("A transform for the raster dataset was not found in the metadata dictionary.")

        updt_pnts = utilities.sample_raster_at_points(pro_pnts['Points'], arr[2, :, :], meta['transform'])
        updt_pnts = updt_pnts.rename(columns={'value': 'elevation'})
        pro_pnts['Points'] = updt_pnts

        return Profile(
            geom_line=ln,
            points_df=pro_pnts['Points'],
            crs=crs,
            profile_name=profile_name
        )

    @staticmethod
    def load_from_survey(points_df: str | Path | pd.DataFrame | gpd.GeoDataFrame,
                         xcoord_name: str = 'X',
                         ycoord_name: str = 'Y',
                         zcoord_name: str = 'elevation',
                         crs: str | int | CRS | None = None,
                         profile_name: str | None = None):
        """Loads x, y, z survey coordinates for a path or line of points.

        Args:
            points_df:
                A pandas DataFrame or geopandas GeoDataFrame that has surveyed points (as separate coordinates
                or in the case of GeoDataFrame, a geometry column of shapely Points) and elevations.
            xcoord_name:
                The column name in the dataframe that correlates with the x-coordinate data, ignored for GeoDataFrames.
            ycoord_name:
                The column name in the dataframe that correlates with the y-coordinate data, ignored for GeoDataFrames.
            zcoord_name:
                The column name in the dataframe that correlates with the elevation data.
            crs:
                A pyproj CRS object representing the coordinate reference system of the input points. If crs argument
                is passed with a GeoDataFrame, crs argument will take precedence and the GeoDataFrame will be projected
                to crs argument. If crs is None with a GeoDataFrame the GeoDataFrame's crs will be used to initialize
                the Profile instance.
            profile_name:
                The name of the Profile being created.

        Returns:
            An initialized Profile class created from the survey data.
        """
        if isinstance(crs, (str, int)):
            crs = CRS.from_epsg(crs)

        if isinstance(points_df, (str, Path)):
            gdf = gpd.read_file(points_df)
            if (isinstance(gdf, pd.DataFrame)) & (isinstance(gdf, gpd.GeoDataFrame)):
                chk = gdf.active_geometry_name
                if not (gdf.geom_type == 'Point').all():
                    raise TypeError("The GeoDataFrame does not contain all points, only points are accepted.")

                if crs is None:
                    crs = gdf.crs
                else:
                    gdf = gdf.to_crs(crs)

                ln = LineString(gdf.geometry.values)
                dists = utilities.cumdistance_from_pnts(gdf.geometry.x.values, gdf.geometry.y.values)
                gdf['distance'] = dists
                gdf = gdf.loc[:, ['distance', zcoord_name, 'geometry']]
                gdf = gdf.rename(columns={zcoord_name: 'elevation'})
            else:
                df = gdf.apply(pd.to_numeric, errors='coerce').fillna(gdf)

                if crs is None:
                    print("WARNING: No coordinate reference system info is available for tabular data...")

                pnts = gpd.points_from_xy(df[xcoord_name], df[ycoord_name])
                ln = LineString(pnts)
                dists = utilities.cumdistance_from_pnts(df[xcoord_name].values, df[ycoord_name].values)
                gdf = pd.DataFrame({'distance': dists, 'elevation': df[zcoord_name]})
                gdf = gpd.GeoDataFrame(gdf, geometry=pnts, crs=crs)

        elif (isinstance(points_df, pd.DataFrame)) & (not isinstance(points_df, gpd.GeoDataFrame)):
            df = points_df.apply(pd.to_numeric, errors='coerce').fillna(points_df)

            if crs is None:
                print("WARNING: No coordinate reference system info is available for tabular data...")

            pnts = gpd.points_from_xy(df[xcoord_name], df[ycoord_name])
            ln = LineString(pnts)
            dists = utilities.cumdistance_from_pnts(df[xcoord_name].values, df[ycoord_name].values)
            gdf = pd.DataFrame({'distance': dists, 'elevation': df[zcoord_name]})
            gdf = gpd.GeoDataFrame(gdf, geometry=pnts, crs=crs)

        elif isinstance(points_df, gpd.GeoDataFrame):
            if not (points_df.geom_type == 'Point').all():
                raise TypeError("The GeoDataFrame does not contain all points, only points are accepted.")

            if crs is None:
                crs = points_df.crs
            else:
                points_df = points_df.to_crs(crs)

            ln = LineString(points_df.geometry.values)
            dists = utilities.cumdistance_from_pnts(points_df.geometry.x.values, points_df.geometry.y.values)
            points_df['distance'] = dists
            gdf = points_df.loc[:, ['distance', zcoord_name, 'geometry']]
            gdf = gdf.rename(columns={zcoord_name: 'elevation'})

        else:
            raise TypeError("The points dataframe was an unexpected input type.")

        return Profile(
            geom_line=ln,
            points_df=gdf,
            crs=crs,
            profile_name=profile_name
        )


if __name__ == '__main__':
    pass
