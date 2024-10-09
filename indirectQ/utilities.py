# Functions and equations used in channel builder
# and Indirect Q calculations

import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
from shapely.geometry import Point, LineString


def points_along_lines(line_pth, dist=1.0, include_endpoint = True, ID_field=None, keep_attrs=None):
    """
    Creates points at equal spaced distances along lines.
    :param line_pth: path to shapefile containing line geometry(s) of the study reach channel centerline.
    :param dist: distance along line to equally space points
    :param include_endpoint: (Boolean) default is True, determines whether to include the endpoint of the line or not.
    :param ID_field: (str) attribute field to use as dictionary keys for each geometry.
    :param keep_attrs: (list of strs) attribute fields to keep with the geometry
    :return: dictionary of geodataframes containing points at specified distances along line, geometries will retain
    user specified attribute fields and be indexed by a user specified field, if none are specified it will return only
    distance and new Line_ID.
    """
    lns = gpd.read_file(line_pth)

    geodfs = {}
    for i in range(len(lns)):
        ln = lns.iloc[i, :]
        attrs = ln.index.to_list()

        if keep_attrs is not None:
            attrs_chk = [item in attrs for item in keep_attrs]
            if all(attrs_chk):
                attrs_to_keep = ln[keep_attrs].index.to_list()
            else:
                raise ValueError("An attribute was listed to keep that does not exist in the available attributes.")

        if keep_attrs is None:
            attrs_to_keep = None

        else:
            attrs_to_keep = None

        #pnts = gpd.GeoDataFrame({'Dwnstrm_Distance': 0.0, 'geometry': Point([ln.geometry.coords[0]])},
        #                        index=[0],
        #                        crs='EPSG:32100')
        current_dist = dist
        pnt_coords = [Point([ln.geometry.coords[0]])]
        distances = [0.0]
        while current_dist < ln.geometry.length:
            interp_pnt = ln.geometry.interpolate(current_dist)
            #new_pnt = gpd.GeoDataFrame({'Dwnstrm_Distance': current_dist,
            #                            'geometry': Point([interp_pnt.coords[0]])},
            #                        index=[0],
            #                        crs='EPSG:32100')
            #pnts = pd.concat([pnts, new_pnt], ignore_index=True)
            pnt_coords.append(Point([interp_pnt.coords[0]]))
            distances.append(current_dist)
            current_dist += dist

        if include_endpoint:
            #new_pnt = gpd.GeoDataFrame({'Dwnstrm_Distance': ln.geometry.length, 'geometry': Point([ln.geometry.coords[-1]])},
            #                        index=[0],
            #                        crs='EPSG:32100')
            #pnts = pd.concat([pnts, new_pnt], ignore_index=True)
            end_pnt = Point([ln.geometry.coords[-1]])
            end_dist = ln.geometry.length
            pnt_coords.append(end_pnt)
            distances.append(end_dist)

        pnts = gpd.GeoDataFrame({'Dwnstrm_Distance': distances,
                                 'geometry': pnt_coords},
                                crs=lns.crs)

        if (attrs_to_keep is None) & (ID_field is None):
            geodfs[str(i)] = dict()
            geodfs[str(i)]['lineID'] = i
            geodfs[str(i)]['Line'] = ln.geometry
            geodfs[str(i)]['Points'] = pnts
        elif (attrs_to_keep is None) & (ID_field is not None):
            geodfs[str(ln[ID_field])] = dict()
            geodfs[str(ln[ID_field])][ID_field] = ln[ID_field]
            geodfs[str(ln[ID_field])]['Line'] = ln.geometry
            geodfs[str(ln[ID_field])]['Points'] = pnts
        elif (attrs_to_keep is not None) & (ID_field is None):
            geodfs[str(i)] = dict()
            geodfs[str(i)]['lineID'] = i
            for a in attrs_to_keep:
                geodfs[str(i)][a] = ln[a]
            geodfs[str(i)]['Line'] = ln.geometry
            geodfs[str(i)]['Points'] = pnts
        elif (attrs_to_keep is not None) & (ID_field is not None):
            geodfs[str(ln[ID_field])] = dict()
            for a in attrs_to_keep:
                geodfs[str(ln[ID_field])][a] = ln[a]
            geodfs[str(ln[ID_field])]['Line'] = ln.geometry
            geodfs[str(ln[ID_field])]['Points'] = pnts

    return geodfs


# Does not check for consistency in coordinate systems, make sure .tif is the same crs as input points.
def sample_raster_at_points(geodf_pnts, raster_pth):
    """
    Takes a geodataframe with point geometries, file path to a raster, and returns the same geodataframe
    with raster value added in a "value" column.
    :param geodf_pnts: Geopandas geodataframe object that has point geometries
    :param raster_pth: file path to the raster to be sampled
    :return:
    """
    gdf = geodf_pnts
    coord_list = [(x, y) for x, y in zip(gdf["geometry"].x, gdf["geometry"].y)]

    with rio.open(raster_pth) as src:
        gdf["value"] = [x[0] / 0.3048 for x in src.sample(coord_list)]

    return gdf


def sort_cross_sections(profile_dict, xsec_dict, invert_line=False):
    """
    Sorts cross-section lines from upstream to downstream based on a channel centerline geometry. Channel centerline
    by default has the first vertice at the most upstream point and the last vertice at the most downstream point.
    Use invert_line argument if Channel centerline starts downstream and goes upstream.
    :param profile_dict: filepath to channel centerline shapefile
    :param xsec_dict: a cross section dictionary object formatted by the points_along_lines function.
    :param invert_line: (boolean) default False, inverts the vertices of the centerline.
    :return: an updated cross-section dictionary with 'Rank' - lowest rank is most upstream, and centerline
    distance.
    """
    chn_cnt_geom = profile_dict['0']['Line']
    if invert_line:
        chn_line = chn_cnt_geom.reverse()
    else:
        chn_line = chn_cnt_geom
    wdict = xsec_dict
    dkeys = []
    xsec_dist = []
    for key in wdict:
        xsec_lin = wdict[key]['Line']
        intrsct_pnt = chn_line.intersection(xsec_lin)
        xsec_prof_dist = chn_line.project(intrsct_pnt)
        dkeys.append(key)
        xsec_dist.append(xsec_prof_dist)

    Fdframe = pd.DataFrame({'keys': dkeys, 'xsec_cntrln_dist': xsec_dist})
    Fdframe.sort_values(by='xsec_cntrln_dist')
    Fdframe['Rank'] = np.arange(1, len(dkeys) + 1)

    for n in range(len(dkeys)):
        row = Fdframe.iloc[n, :]
        wdict[row['keys']]['Channel Rank'] = row['Rank']
        wdict[row['keys']]['Xsec_Cntrln_Dist'] = row['xsec_cntrln_dist']

    return wdict


def sort_points(profile_dict, gdframe_pnts, invert_line=False):
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
    chn_cnt_geom = profile_dict['0']['Line']
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
        refvec = _calc_vector((origin_pnt.x, origin_pnt.y), (refvec_pnt.x, refvec_pnt.y))
        a = _calc_signed_angles_from_origin((origin_pnt.x, origin_pnt.y), (refvec[0], refvec[1]), pnts['X'][i], pnts['Y'][i])

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


def interpz_at_intersection(prof_pnts, xsec_dict, interp_method='linear'):
    """
    Find intersection of two shapely lines and then interp the distance to the intersection point based on the
    profile line
    :param prof_line:
    :param intrsct_line:
    :return:
    """
    xsec_dict = xsec_dict
    rght_srt = prof_pnts[prof_pnts['Waters_Edge'] == 'R'].sort_values(by='Cntrline_Dist')
    rght_ln = LineString(rght_srt.geometry.to_list())
    left_srt = prof_pnts[prof_pnts['Waters_Edge'] == 'L'].sort_values(by='Cntrline_Dist')
    left_ln = LineString(left_srt.geometry.to_list())

    for key, item in xsec_dict.items():
        if rght_ln.intersects(item['Line']):
            r_pnt = rght_ln.intersection(item['Line'])
            rwse_dist = rght_ln.project(r_pnt)
            rxsec_dist = item['Line'].project(r_pnt)
            rpro_d = np.array(cumdistance_from_pnts(rght_srt['X'].to_list(), rght_srt['Y'].to_list()))
            rpro_z = rght_srt['Z'].values
            rWSE = np.interp(rwse_dist, rpro_d, rpro_z)
        else:
            rWSE = np.nan

        if left_ln.intersects(item['Line']):
            l_pnt = left_ln.intersection(item['Line'])
            lwse_dist = left_ln.project(l_pnt)
            lxsec_dist = item['Line'].project(l_pnt)
            lpro_d = np.array(cumdistance_from_pnts(left_srt['X'].to_list(), left_srt['Y'].to_list()))
            lpro_z = left_srt['Z'].values
            lWSE = np.interp(lwse_dist, lpro_d, lpro_z)
        else:
            lWSE = np.nan

        mn_WSE = np.array([rWSE, lWSE])
        if np.isnan(mn_WSE).all():
            print("Cross Section does not intersect available water surface profile data.")
            break
        elif np.isnan(mn_WSE).any():
            mn_WSE = mn_WSE[~np.isnan(mn_WSE)][0]
        else:
            mn_WSE = mn_WSE.mean()

        xsecpntdist = item['Points']['Dwnstrm_Distance'].to_list()
        xsecpntz = item['Points']['value'].to_list()
        xswse_ln = LineString([Point(xsecpntdist[0], mn_WSE), Point(xsecpntdist[-1], mn_WSE)])
        xspnt_ln = LineString([Point(x, y) for x, y in zip(xsecpntdist, xsecpntz)])
        wse_intsct_pnts = xspnt_ln.intersection(xswse_ln)
        wse_int_coords =[[p.x, p.y] for p in wse_intsct_pnts]
        item['WSE_Intersects'] = wse_int_coords
        xsec_dict[key] = item

    return xsec_dict


#def sample_polygons_at_points(geodf_pnts, raster_pth)


#def _define_origin(point)


# def _orth_vector(vector):
#     """ Returns the vector that is orthogonal (perpendicular) to the input vector. """
#     orth_vect = np.array([vector[1], -vector[0]])
#     return orth_vect

def cumdistance_from_pnts(xcoords, ycoords):
    dists = []
    for i, pnt in enumerate(zip(xcoords, ycoords)):
        if i == 0:
            d = 0.0
            dists.append(d)
        else:
            x1, y1 = list(zip(xcoords, ycoords))[i - 1]
            x2, y2 = pnt
            d = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            dists.append(d + dists[-1])

    return dists


def _calc_vector(p1, p2):
    """ Returns a vector from a pair of points (coordinate tuples). Vector direction is defined as
     from p1 -> p2"""
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


def _calc_signed_angles_from_origin(origin_pnt, ref_vector, xcoords, ycoords):
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