# USGS methods for indirect discharge calculations

import pandas as pd
import geopandas as gpd
import numpy as np
from scipy.integrate import trapezoid

from indirectQ import utilities

# class Sections:
#
#     def __init__(self, siteID, agency):
#
#         self.siteID = siteID
#         self.agency = agency
#         # FUTURE DEVELOPMENT
#         #   Define methods for looking up location information using input site ID and operating agency
#         self.siteName = None
#         self.siteLat = None
#         self.siteLon = None
#         # geometry variables
#         self.XSections = None
#         self.Centerline = None
#         self.HWMs = None
#
#     def load_Xsections(self, in_paths):
#
#
#     def load_ManningsRaster(self, man_path):
#
#
#     def load_topography(self, topo_path):
#
#     def extract_Xsections(self):
#
#     def _
#
#     def calculateQ(self):


# Currently only configure to accept 4 cross sections
def calculateQ(xsec_dict, units='US'):
    """
    Estimates Slope-Area Discharge from formatted cross section dictionary (from channelbuilder module).
    :param xsec_dict: channelbuilder.py - LoadSurvey class cross-section dictionary
    :return: float, estimated discharge
    """
    if units == 'US':
        g = 32.17405
    elif units == 'SI':
        g = 9.80665
    else:
        print("Units not recognized...assume US.")
        g = 32.17405
    # get cross section below WSE
    ranks = {}
    for key, item in xsec_dict.items():
        xsdf = item['Points'][['Dwnstrm_Distance', 'value', 'Name', 'n_Value']]
        # clip cross section to WSE
        wsclp = xsdf.copy()[xsdf['value'] <= item['WSE_Intersects'][0][1]]
        # add intersect points
        for p in item['WSE_Intersects']:
            wsclp.loc[wsclp.index[-1]+1] = [p[0], p[1], np.nan, np.nan]

        wsclp.sort_values(by='Dwnstrm_Distance', inplace=True)
        wsclp.ffill(inplace=True)
        wsclp.bfill(inplace=True)
        wsclp['depth'] = item['WSE_Intersects'][0][1] - wsclp['value']
        xsec_chnks = np.unique(wsclp['Name'])

        chnk_cols = ['a', 'wp', 'r', 'K', 'K3/a2', 'n', '1.486/n']
        chnk_vals = []
        for chnk in xsec_chnks:
            xsr = wsclp[wsclp['Name'] == chnk]
            a = trapezoid(xsr['depth'], xsr['Dwnstrm_Distance'])
            wp = utilities.cumdistance_from_pnts(xsr['Dwnstrm_Distance'].to_list(), xsr['depth'].to_list())[-1]
            if units == 'US':
                nc = 1.486 / xsr['n_Value'].iloc[2]
            elif units == 'SI':
                nc = xsr['n_Value'].iloc[2]
            else:
                print("Units not recognized...assume US.")
                nc = 1.486 / xsr['n_Value'].iloc[2]
            r = (a/wp)**(2.0/3.0)
            K = nc * a * r
            Ka = K**3 / a**2
            chnk_vals.append([a, wp, r, K, Ka, xsr['n_Value'].iloc[2], nc])

        chnkdf = pd.DataFrame(chnk_vals, columns=chnk_cols)
        chnkdf['CrossSection'] = key
        chnkdf['K/Kt'] = chnkdf['K'] / chnkdf['K'].sum()
        xs_a = chnkdf['a'].sum()
        xs_K = chnkdf['K'].sum()
        xs_alpha = (chnkdf['K3/a2'].sum()) / (xs_K**3 / xs_a**2)
        item['XS_Area'] = xs_a
        item['XS_K'] = xs_K
        item['XS_alpha'] = xs_alpha
        item['XS_Parts'] = chnkdf
        xsec_dict[key] = item
        ranks[item['Channel Rank']] = key

    rch_cols = ['Reach', 'XSs', 'Length', 'WSE_fall', 'Kw', 'hvs', 'delta_hv', 'hf', 'Sf^1/2', 'assumed_Q', 'calc_Q']
    rch_vals = []

    for rch in np.arange(1, len(xsec_dict)):
        xss = (rch, rch+1)
        rlength = xsec_dict[ranks[rch]]['Xsec_Cntrln_Dist'] - xsec_dict[ranks[rch+1]]['Xsec_Cntrln_Dist']
        K1 = xsec_dict[ranks[rch+1]]['XS_K']
        K2 = xsec_dict[ranks[rch]]['XS_K']
        A1 = xsec_dict[ranks[rch+1]]['XS_Area']
        A2 = xsec_dict[ranks[rch]]['XS_Area']
        if A1 > A2:
            k = 1.0
        elif A1 < A2:
            k = 0.5
        else:
            k = 1.0
        alpha1 = xsec_dict[ranks[rch+1]]['XS_alpha']
        alpha2 = xsec_dict[ranks[rch]]['XS_alpha']
        fall = xsec_dict[ranks[rch+1]]['WSE_Intersects'][0][1] - xsec_dict[ranks[rch]]['WSE_Intersects'][0][1]
        wght_K = np.sqrt(K1 * K2)
        assmQ = K2 * np.sqrt((fall / ((K2/K1)*rlength + ((K2**2)/(2*g*(A2**2)))*(-alpha1*k*((A2/A1)**2)+(alpha2*k)))))
        V1 = assmQ / A1
        V2 = assmQ / A2
        hv1 = alpha1 * ((V1**2) / (2 * g))
        hv2 = alpha2 * ((V2**2) / (2 * g))
        del_hv = hv1 - hv2
        if del_hv >= 0:
            hf = fall + (del_hv * 0.5)
        elif del_hv < 0:
            hf = fall + del_hv
        Sf = hf / rlength
        Sf_sqrt = np.sqrt(Sf)
        compQ = wght_K * Sf_sqrt
        rch_vals.append([rch, xss, rlength, fall, wght_K, (hv1, hv2), del_hv, hf, Sf_sqrt, assmQ, compQ])

    rchdf = pd.DataFrame(rch_vals, columns=rch_cols)

    def k_value_assgn(del_hv):
        if del_hv >= 0:
            k = 0.5
        elif del_hv < 0:
            k = 0.0
        else:
            k = 0.0

        return k

    # final equation variables
    K1 = xsec_dict[ranks[4]]['XS_K']
    K2 = xsec_dict[ranks[3]]['XS_K']
    K3 = xsec_dict[ranks[2]]['XS_K']
    K4 = xsec_dict[ranks[1]]['XS_K']
    A1 = xsec_dict[ranks[4]]['XS_Area']
    A2 = xsec_dict[ranks[3]]['XS_Area']
    A3 = xsec_dict[ranks[2]]['XS_Area']
    A4 = xsec_dict[ranks[1]]['XS_Area']
    alpha1 = xsec_dict[ranks[4]]['XS_alpha']
    alpha2 = xsec_dict[ranks[3]]['XS_alpha']
    alpha3 = xsec_dict[ranks[2]]['XS_alpha']
    alpha4 = xsec_dict[ranks[1]]['XS_alpha']
    L12 = rchdf[rchdf['Reach'] == 3].iloc[0, 2]
    L23 = rchdf[rchdf['Reach'] == 2].iloc[0, 2]
    L34 = rchdf[rchdf['Reach'] == 1].iloc[0, 2]
    fall = rchdf['WSE_fall'].sum()
    k12 = k_value_assgn(rchdf[rchdf['Reach'] == 3].iloc[0, 6])
    k23 = k_value_assgn(rchdf[rchdf['Reach'] == 2].iloc[0, 6])
    k34 = k_value_assgn(rchdf[rchdf['Reach'] == 1].iloc[0, 6])

    finalAtrm = (((K4**2)*L12)/(K1*K2)) + (((K4**2)*L23)/(K2*K3)) + (((K4**2)*L34)/(K3*K4))
    finalBtrm = ((K4**2)/((A4**2)*2*g)) * ((-alpha1*((A4/A1)**2)*(1-k12)) + (-alpha2*((A4/A2)**2)*(k23-k12)) + (-alpha3*((A4/A3)**2)*(k34-k23)) + (alpha4*(1-k34)))

    finalQ = K4 * np.sqrt(fall / (finalAtrm + finalBtrm))

    chnklst = []
    for key, item in xsec_dict.items():
        allchnks = item['XS_Parts']
        chnklst.append(allchnks)

    final_prts = pd.concat(chnklst)
    final_prts['q'] = finalQ * final_prts['K/Kt']
    final_prts['v'] = final_prts['q'] / final_prts['a']

    return {'Q': finalQ, 'Reaches': rchdf, 'XS_Parts': final_prts}