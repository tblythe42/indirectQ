"""Utility classes and functions used by the indirectQ package.

This module contains the following functions:
- points_along_lines()
- xsections_along_line()
- sample_raster_at_points()
- sort_cross_sections()
- sort_points()
- interpz_at_intersection()
- cumdistance_from_pnts()
- calc_vector()
- calc_signed_angles_from_origin()
- raster_array_from_file()

-  **Author(s):** Todd Blythe, MTDNRC, CNB968
-  **Date:** Created 12/10/2025
"""
from pathlib import Path

import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
from shapely.geometry import Point, LineString
from affine import Affine

from indirectQ.profile import Profile


def points_along_lines(lines: gpd.GeoDataFrame,
                       dist: float = 1.0,
                       include_endpoint: bool = True,
                       id_field: str | None = None,
                       keep_attrs: list[str] | None = None):
    """Creates points at equal spacing along lines.

    Args:
        lines:
            A GeoDataFrame containing line geometries representing the stream channel centerline(s).
        dist:
            Distance to space points along line, interpreted in units of the input CRS.
        include_endpoint:
            Determines whether to include the endpoint of the line as a point or not, default is True.
        id_field:
            The column in the input GeoDataFrame to use as the ID or dictionary key for each line. If None is given
            the default is to assign an incremental ID as a string.
        keep_attrs:
            List of columns in the input GeoDataFrame to keep. These are included as additional keys in each line's
            dict.

    Returns:
        A dictionary with each line, attributes kept from input GeoDataFrame, the line geometry, and the points
        along the centerline as a GeoDataFrame.
    """
    if not (lines.geom_type == 'LineString').all():
        raise TypeError("Not all geometries in the input GeoDataFrame are lines.")

    lns = lines.copy()

    geodfs = {}
    for i in range(len(lns)):
        ln = lns.iloc[i, :]
        attrs = ln.index.to_list()

        if keep_attrs is not None:
            attrs_chk = [item in attrs for item in keep_attrs]
            if all(attrs_chk):
                attrs_to_keep = ln[keep_attrs].index.to_list()
            else:
                print("An attribute was listed to keep that does not exist in the available attributes. It will"
                      "be skipped.")
                attrs_to_keep = ln[attrs_chk].index.to_list()
        else:
            attrs_to_keep = None

        current_dist = dist
        pnt_coords = [Point([ln.geometry.coords[0]])]
        distances = [0.0]
        while current_dist < ln.geometry.length:
            interp_pnt = ln.geometry.interpolate(current_dist)
            pnt_coords.append(Point([interp_pnt.coords[0]]))
            distances.append(current_dist)
            current_dist += dist

        if include_endpoint:
            end_pnt = Point([ln.geometry.coords[-1]])
            end_dist = ln.geometry.length
            pnt_coords.append(end_pnt)
            distances.append(end_dist)

        pnts = gpd.GeoDataFrame({'distance': distances,
                                 'geometry': pnt_coords},
                                crs=lns.crs)

        if (attrs_to_keep is None) & (id_field is None):
            geodfs[i] = dict()
            geodfs[i]['lineID'] = i
            geodfs[i]['Line'] = ln.geometry
            geodfs[i]['Points'] = pnts
        elif (attrs_to_keep is None) & (id_field is not None):
            geodfs[i] = dict()
            geodfs[i]['lineID'] = ln[id_field]
            geodfs[i]['Line'] = ln.geometry
            geodfs[i]['Points'] = pnts
        elif (attrs_to_keep is not None) & (id_field is None):
            geodfs[i] = dict()
            geodfs[i]['lineID'] = i
            for a in attrs_to_keep:
                geodfs[i][a] = ln[a]
            geodfs[i]['Line'] = ln.geometry
            geodfs[i]['Points'] = pnts
        elif (attrs_to_keep is not None) & (id_field is not None):
            geodfs[i] = dict()
            geodfs[i]['lineID'] = ln[id_field]
            for a in attrs_to_keep:
                geodfs[i][a] = ln[a]
            geodfs[i]['Line'] = ln.geometry
            geodfs[i]['Points'] = pnts

    return geodfs


def xsections_along_line(lines: gpd.GeoDataFrame,
                         dist: float = 1.0,
                         random: bool = False,
                         id_field: str | None = None,
                         keep_attrs: list[str] | None = None):
    pass


# Does not check for consistency in coordinate systems, make sure .tif is the same crs as input points.
def sample_raster_at_points(points: gpd.GeoDataFrame, raster: np.ndarray, transform: Affine) -> gpd.GeoDataFrame:
    """Samples raster values at a point or list of points.

    Args:
        points:
            GeoDataFrame of points, must all be point geometries.
        raster:
            The data array loaded from a raster dataset.
        transform:
            The transform for the raster data to convert between grid index and spatial coordinates, must be
            an Affine object (which is what is returned by rasterio).

    Returns:
        A copy of the input points GeoDataFrame but with a new column added ['value'] which represents the
        underlying raster data for each point.
    """
    if not (points.geom_type == 'Point').all():
        raise TypeError("Not all geometries in the input GeoDataFrame are points.")

    gdf = points.copy()
    grid_ids = rio.transform.rowcol(transform, gdf.geometry.x.values, gdf.geometry.y.values)
    vals = raster[grid_ids]
    gdf['value'] = vals

    return gdf


def sort_cross_sections(profile: Profile, xsec_dict: dict, invert_line: bool = False) -> dict:
    """

    Sorts cross-section lines from upstream to downstream based on a channel centerline geometry. Channel centerline
    by default has the first vertice at the most upstream point and the last vertice at the most downstream point.

    Args:
        profile:
            An indirectQ Profile object that is used as the sorting line (perpendicular to the cross-sections).
        xsec_dict:
            A dictionary returned by the indirectQ ChannelBuilder cross_section property.
        invert_line:
            Inverts the order of the vertices of the profile (centerline), in case it is not oriented from upstream to
            downstream which is required to order the cross-sections correctly. Default is False.

    Returns:
        An updated dictionary identical to xsec_dict argument.
    """
    chn_cnt_geom = profile.line
    if invert_line:
        chn_line = chn_cnt_geom.reverse()
    else:
        chn_line = chn_cnt_geom
    wdict = xsec_dict
    dkeys = []
    xsec_dist = []
    for key in wdict:
        xsec_lin = wdict[key]['Profile'].line
        intrsct_pnt = chn_line.intersection(xsec_lin)
        xsec_prof_dist = chn_line.project(intrsct_pnt)
        dkeys.append(key)
        xsec_dist.append(xsec_prof_dist)

    fdframe = pd.DataFrame({'keys': dkeys, 'xsec_cntrln_dist': xsec_dist})
    fdframe = fdframe.sort_values(by='xsec_cntrln_dist')
    fdframe['Rank'] = np.arange(1, len(dkeys) + 1)

    for n in range(len(dkeys)):
        row = fdframe.iloc[n, :]
        wdict[row['keys']]['Channel Rank'] = row['Rank']
        wdict[row['keys']]['Xsec_Cntrln_Dist'] = row['xsec_cntrln_dist']

    return wdict


def sort_points(profile: Profile, gdframe_pnts: gpd.GeoDataFrame, invert_line: bool = False):
    """

    Args:
        profile:
        gdframe_pnts:
        invert_line:

    Returns:

    """
    """
    Sorts surveyed points from upstream to downstream based on a channel centerline geometry and assigns left and right
    edges. Channel centerline by default has the first vertice at the most upstream point and the last vertice at the
    most downstream point. Use invert_line argument if Channel centerline starts downstream and goes upstream.
    :param cntrln_pth: filepath to channel centerline shapefile
    :param xsec_dict: a cross section dictionary object formatted by the points_along_lines function.
    :param invert_line: (boolean) default False, inverts the vertices of the centerline.
    :return: an updated cross-section dictionary with longitudinal profile information added: 'Rank' - lowest rank is
    most upstream; 'Distance_to_Dwnstrm' - distance along the channel to the next downstream cross-section;
    'Distance_to_Upstrm' - distance along the channel to the next upstream cross-section.
    """
    chn_cnt_geom = profile.line
    if invert_line:
        chn_line = chn_cnt_geom.reverse()
    else:
        chn_line = chn_cnt_geom
    pnts = gdframe_pnts
    pnts['Cntrline_Dist'] = gdframe_pnts.geometry.apply(lambda row: chn_cnt_geom.project(row))

    # Use np.arctan2() to get the signed angle using coordinate points first define a function that recenters all the
    # points based on a new origin and resets the azimuth based on the centerline direction
    ws_side = []
    for i in range(len(pnts)):
        origin_pnt = chn_cnt_geom.interpolate(pnts['Cntrline_Dist'][i] - 1)
        refvec_pnt = chn_cnt_geom.interpolate(pnts['Cntrline_Dist'][i] + 1)
        refvec = calc_vector((origin_pnt.x, origin_pnt.y), (refvec_pnt.x, refvec_pnt.y))
        a = calc_signed_angles_from_origin((origin_pnt.x, origin_pnt.y), (refvec[0], refvec[1]),
                                           pnts.geometry.x[i], pnts.geometry.y[i])
        if a > 0:
            bnk = 'L'
        elif a < 0:
            bnk = 'R'
        else:
            bnk = 'UNKNOWN'

        ws_side.append(bnk)

    pnts['Waters_Edge'] = ws_side
    return pnts

    #     o_vcnt = _orth_vector(vcnt)
    #     u_vcnt = _unit_vector(vcnt)
    #     # Create Vector to WS Point
    #     vws = _calc_vector((cp.x[0], cp.y[0]), (pnts['X'][i], pnts['Y'][i]))
    #     u_vws = _unit_vector(vws)
    #     # get vector angle
    #     a = _angle_between_vectors(u_vcnt, u_vws)
    #     # check sign of dot product
    #     dot_sgn = np.dot(vws, o_vcnt)
    #     if dot_sgn > 0:
    #         bnk = 'L'
    #     elif dot_sgn < 0:
    #         bnk = 'R'
    #     else:
    #         bnk = 'UNKNOWN'
    #     ws_side.append(bnk)
    #
    # origin_pnt = chn_cnt_geom.interpolate(pnts['Cntrline_Dist'][i] - 1)
    # cpdwn = chn_cnt_geom.interpolate(pnts['Cntrline_Dist'][i] + 1)
    # vcnt = _calc_vector((cp.x[0], cp.y[0]), (cpdwn.x[0], cpdwn.y[0]))


def interpz_at_intersection(watsurf: dict, xsec_dict: dict):
    """Interpolates water surface elevation at cross-sections.

    Args:
        watsurf:
        xsec_dict:

    Returns:

    """
    for wsid, ws in watsurf.items():
        rght_ln = ws['R']['Profile'].line
        left_ln = ws['L']['Profile'].line

        for key, item in xsec_dict.items():
            if rght_ln.intersects(item['Profile'].line):
                r_pnt = rght_ln.intersection(item['Profile'].line)
                rws_dist = rght_ln.project(r_pnt)
                rpro_d = ws['R']['Profile'].points['distance'].values
                rpro_z = ws['R']['Profile'].points['elevation'].values
                rwse = np.interp(rws_dist, rpro_d, rpro_z)
            else:
                rwse = np.nan

            if left_ln.intersects(item['Profile'].line):
                l_pnt = left_ln.intersection(item['Profile'].line)
                lws_dist = left_ln.project(l_pnt)
                lpro_d = ws['L']['Profile'].points['distance'].values
                lpro_z = ws['L']['Profile'].points['elevation'].values
                lwse = np.interp(lws_dist, lpro_d, lpro_z)
            else:
                lwse = np.nan

            mn_wse = np.array([rwse, lwse])
            if np.isnan(mn_wse).all():
                print("Cross Section does not intersect available water surface profile data.")
                continue
            elif np.isnan(mn_wse).any():
                mn_wse = mn_wse[~np.isnan(mn_wse)][0]
            else:
                mn_wse = mn_wse.mean()

            xsecpntdist = item['Profile'].points['distance'].to_list()
            xsecpntz = item['Profile'].points['elevation'].to_list()
            xswse_ln = LineString([Point(xsecpntdist[0], mn_wse), Point(xsecpntdist[-1], mn_wse)])
            xspnt_ln = LineString([Point(x, y) for x, y in zip(xsecpntdist, xsecpntz)])
            wse_intsct_pnts = xspnt_ln.intersection(xswse_ln)
            wse_int_coords = [[p.x, p.y] for p in wse_intsct_pnts.geoms]
            if item['WSE_Intersects'] is None:
                item['WSE_Intersects'] = {wsid: wse_int_coords}
            else:
                item['WSE_Intersects'][wsid] = wse_int_coords
            if item['Mean_WSE'] is None:
                item['Mean_WSE'] = {wsid: mn_wse}
            else:
                item['Mean_WSE'][wsid] = mn_wse
            xsec_dict[key] = item

    return xsec_dict


#def sample_polygons_at_points(geodf_pnts, raster_pth)


#def _define_origin(point)


# def _orth_vector(vector):
#     """ Returns the vector that is orthogonal (perpendicular) to the input vector. """
#     orth_vect = np.array([vector[1], -vector[0]])
#     return orth_vect

def cumdistance_from_pnts(xcoords: list[float] | np.ndarray, ycoords: list[float] | np.ndarray) -> list[float]:
    """Calculates the sum of all distances for a set of ordered points.

    Coordinates for the set of points are assumed to be ordered (e.g., the coordinates returned from a shapely
    LineString or points sampled at equal distances along a linear path).

    Args:
        xcoords:
            A list or array of the x-coordinates.
        ycoords:
            A list or array of the y-coordinates.

    Returns:
        A list of the cumulative distance of each point from the starting point.
    """
    coords = list(zip(xcoords, ycoords))
    dists = []
    for i, pnt in enumerate(coords):
        if i == 0:
            d = 0.0
            dists.append(d)
        else:
            x1, y1 = coords[i - 1]
            x2, y2 = pnt
            d = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            dists.append(d + dists[-1])

    return dists


def calc_vector(p1: tuple[float, float], p2: tuple[float, float]) -> np.ndarray:
    """Returns the direction of a vector from a pair of points <x,y> as a numpy array.

    Args:
        p1:
            The starting or origin point.
        p2:
            The terminal point or end of the vector.
    """
    xd = p2[0] - p1[0]
    yd = p2[1] - p1[1]
    vect = np.array([xd, yd])
    return vect


# def _unit_vector(vector):
#     """ Returns the unit vector of the vector.  """
#     return vector / np.linalg.norm(vector)
#
#
# def _angle_between_vectors(v1, v2):
#     """ Returns the angle in radians between vectors 'v1' and 'v2'"""
#     v1_u = _unit_vector(v1)
#     v2_u = _unit_vector(v2)
#     return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def calc_signed_angles_from_origin(origin_pnt, ref_vector, xcoords, ycoords):
    """
    Given a new origin point and reference vector that defines the 0 azimuth, returns the signed angle between
    the reference vector and vectors formed by the origin and input x, y point coordinates.
    :param origin_pnt: (tuple) x, y coordinate of the new origin point
    :param ref_vector: (tuple) reference vector to set the 0 azimuth
    :param xcoords: (array-like) list or array/series of x coordinates
    :param ycoords: (array-like) list or array/series of y coordinates
    :return: array of signed angle to points formed by the x, y coordinate arrays
    """
    o_x, o_y = origin_pnt
    transl_x = xcoords - o_x
    transl_y = ycoords - o_y
    r_vecx, r_vecy = ref_vector
    coord_rotation = np.arctan2(r_vecy, r_vecx)
    raw_angles = np.arctan2(transl_y, transl_x)
    rotated_angles = raw_angles - coord_rotation
    return rotated_angles


def raster_array_from_file(fl_pth: str | Path, data_scalar: float = 1.0) -> tuple:
    """Creates a numpy ndarray from raster data loaded from a geotiff.

    Args:
        fl_pth:
            A pathlike object or string to a geotiff file with one band of data.
        data_scalar:
            Conversion factor for raster data (e.g., converting meters to feet for elevation data).

    Returns:
        A tuple with a multidimensional array containing x-coordinates, y-coordinates, and the raster values. The array
        shape = (X(raster_rows, raster_columns), Y(raster_rows, raster_columns), Value(raster_rows, raster_columns)).
        In the second position of the tuple is a dictionary of raster metadata including the transform and CRS.

    """
    with rio.open(fl_pth, 'r') as tp_dset:
        topoz = tp_dset.read(1)
        metdata = tp_dset.meta
        topoz = np.where(topoz == metdata['nodata'], np.nan, topoz)
        height, width = topoz.shape
        cols, rows = np.meshgrid(np.arange(width), np.arange(height))
        xs, ys = rio.transform.xy(tp_dset.transform, rows, cols)
        topox = np.array(xs)
        topox = np.reshape(topox, topoz.shape)
        topoy = np.array(ys)
        topoy = np.reshape(topoy, topoz.shape)

    return np.stack((topox, topoy, topoz * data_scalar), axis=0), metdata


def check_bounds(bounds: list | tuple | np.ndarray, compare_bounds: list | tuple | np.ndarray) -> tuple:
    """Checks if the comparison bounds are contained within another set of bounds.

    Args:
        bounds:
            A list or array-like object with geospatial bounds of the order West, South, East, North or
            (x_min, y_min, x_max, y_max)
        compare_bounds:
            A list or array-like object with geospatial bounds of the order West, South, East, North or
            (x_min, y_min, m_max, y_max). These bounds are compared to bounds to see if they fully overlap.
    Returns:
        Tuple of booleans for each bound element (e.g., if y_min of compare_bounds is within bounds, the 1 index
        of the tuple will be True.
    """
    if (compare_bounds[0] >= bounds[0]) & (compare_bounds[0] <= bounds[2]):
        xmin = True
    else:
        xmin = False

    if (compare_bounds[1] >= bounds[1]) & (compare_bounds[1] <= bounds[3]):
        ymin = True
    else:
        ymin = False

    if (compare_bounds[2] >= bounds[0]) & (compare_bounds[2] <= bounds[2]):
        xmax = True
    else:
        xmax = False

    if (compare_bounds[3] >= bounds[1]) & (compare_bounds[3] <= bounds[3]):
        ymax = True
    else:
        ymax = False

    return xmin, ymin, xmax, ymax
