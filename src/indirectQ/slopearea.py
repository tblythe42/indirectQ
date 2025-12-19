"""Uses Slope-Area methods to estimate discharge for a channel

This module contains the following functions:
- slopearea_q()

-  **Author(s):** Todd Blythe, MTDNRC, CNB968
-  **Date:** Created 12/10/2025
"""

import pandas as pd
import numpy as np
from scipy.integrate import trapezoid

from indirectQ import utilities
from indirectQ import ChannelBuilder


def slopearea_q(channel: ChannelBuilder, units: str = 'US', ws_id: str | int | None = None) -> tuple:
    """Estimates discharge in channel using 2 or more cross-sections.

    This follows the methods for slope-area indirect discharge calculation outlined in Dalrymple and Benson (1967), from
    USGS Techniques of Water-Resources Investigations, book 3, chap. A2.

    A properly formatted indirectQ ChannelBuilder instance is required for accurate results.

    Args:
        channel:
            An indirectQ ChannelBuilder object with cross-sections added, at least one specified water profile, and
            Manning's n values assigned to each cross-section.
        units:
            Choose between 'US' (feet, cubic-feet, etc.) or 'SI' (meters, cubic-meters, etc.). This does not convert
            your survey or input data, it only specified which gravitation acceleration constant to use and whether
            to apply conversion to Manning's n values. The user must make sure their ChannelBuilder object is set up
            in the correct units.
        ws_id:
            A water surface profile ID that can be found in the ChannelBuilder object. This will run the method on
            only that water profile. Default is None, which will compute estimated discharges for all profiles found
            in the ChannelBuilder object.

    Returns:
        A tuple of discharge estimation results (DataFrame of final discharge by water profile, supplementary data used
        in the calculations)
    """
    if units == 'US':
        g = 32.17405
    elif units == 'SI':
        g = 9.80665
    else:
        print("Units not recognized...assume US.")
        g = 32.17405
    # get cross section below WSE
    if ws_id is not None:
        ws_id = [ws_id]
    else:
        ws_id = list(channel.water_profiles.keys())

    xsecs = channel.cross_sections
    finQs = []
    sup_results = {}
    for ws in ws_id:
        xsec_cols = ['XSID', 'Rank', 'Area', 'K', 'alpha', 'mean_depth']
        xsec_df = []
        chnk_dfs = []
        for key, item in xsecs.items():
            xsdf = item['Profile'].points.loc[:, ['distance', 'elevation']]
            xsdf['mannings'] = channel.mannings[key]
            # clip cross section to WSE
            wsclp = xsdf.loc[xsdf['elevation'] <= item['Mean_WSE'][ws]]
            # add intersect points
            wspnt_lst = []
            for p in item['WSE_Intersects'][ws]:
                pl = [p[0], p[1], np.nan]
                wspnt_lst.append(pl)
            wspnt_df = pd.DataFrame(wspnt_lst, columns=wsclp.columns)
            wsclp = pd.concat([wsclp, wspnt_df], ignore_index=True)

            wsclp = wsclp.sort_values(by='distance')
            wsclp = wsclp.ffill()
            wsclp = wsclp.bfill()
            wsclp['depth'] = item['Mean_WSE'][ws] - wsclp['elevation']
            bpoints = np.where(np.diff(wsclp['mannings'].values))[0] + 1
            if len(bpoints) == 0:
                wsclp['chunk'] = 'manning1'
            else:
                chnks = np.split(wsclp['mannings'].values, bpoints)
                ck_lst = []
                for i, c in enumerate(chnks):
                    cname = [f'manning{i + 1}'] * c.size
                    ck_lst.extend(cname)
                wsclp['chunk'] = ck_lst
            xsec_chnks = np.unique(wsclp['chunk'])

            chnk_cols = ['a', 'wp', 'r', 'K', 'K3/a2', 'n', '1.486/n']
            chnk_vals = []
            for chnk in xsec_chnks:
                xsr = wsclp.loc[wsclp['chunk'] == chnk]
                a = trapezoid(xsr['depth'], xsr['distance'])
                wp = utilities.cumdistance_from_pnts(xsr['distance'].to_list(), xsr['depth'].to_list())[-1]
                if units == 'US':
                    nc = 1.486 / xsr['mannings'].iloc[2]
                elif units == 'SI':
                    nc = xsr['mannings'].iloc[2]
                else:
                    print("Units not recognized...assume US.")
                    nc = 1.486 / xsr['mannings'].iloc[2]
                r = (a / wp)**(2.0 / 3.0)
                K = nc * a * r
                Ka = K**3 / a**2
                chnk_vals.append([a, wp, r, K, Ka, xsr['mannings'].iloc[2], nc])

            chnkdf = pd.DataFrame(chnk_vals, columns=chnk_cols)
            chnkdf['crosssection'] = key
            chnkdf['K/Kt'] = chnkdf['K'] / chnkdf['K'].sum()

            chnk_dfs.append(chnkdf)

            xs_a = chnkdf['a'].sum()
            xs_K = chnkdf['K'].sum()
            xs_alpha = (chnkdf['K3/a2'].sum()) / (xs_K**3 / xs_a**2)
            xsecvals = [key, item['Channel Rank'], xs_a, xs_K, xs_alpha, wsclp['depth'].mean()]
            xsec_df.append(xsecvals)

        chnk_df = pd.concat(chnk_dfs, ignore_index=True)
        xsecsdf = pd.DataFrame(xsec_df, columns=xsec_cols)

        rch_cols = ['Reach', 'Ks', 'Ke', 'As', 'Ae', 'alpha_s', 'alpha_e', 'k', 'Length', 'WSE_fall', 'Kw',
                    'hvs', 'hve', 'delta_hv', 'hf', 'Sf^1/2', 'assumed_Q', 'calc_Q']
        rch_vals = []

        for rch in np.arange(1, len(xsecsdf)):
            s = rch
            e = rch + 1
            xss = xsecsdf.loc[xsecsdf['Rank'] == s].iloc[0]
            xse = xsecsdf.loc[xsecsdf['Rank'] == e].iloc[0]
            rlength = xsecs[xse['XSID']]['Xsec_Cntrln_Dist'] - xsecs[xss['XSID']]['Xsec_Cntrln_Dist']
            #K1 = [ranks[rch+xsec_dict1]]['XS_K']
            #K2 = xsec_dict[ranks[rch]]['XS_K']
            #A1 = xsec_dict[ranks[rch+1]]['XS_Area']
            #A2 = xsec_dict[ranks[rch]]['XS_Area']
            if xse['Area'] < xss['Area']:
                ki = 0.0
            else:
                ki = 0.5
            #alpha1 = xsec_dict[ranks[rch+1]]['XS_alpha']
            #alpha2 = xsec_dict[ranks[rch]]['XS_alpha']
            fall = xsecs[xss['XSID']]['Mean_WSE'][ws] - xsecs[xse['XSID']]['Mean_WSE'][ws]
            wght_K = np.sqrt(xss['K'] * xse['K'])
            #assmQ = K2 * np.sqrt((fall / ((K2/K1)*rlength + ((K2**2)/(2*g*(A2**2)))*(-alpha1*k*((A2/A1)**2)+(alpha2*k)))))
            assmq = xse['K'] * np.sqrt((fall / ((xse['K'] / xss['K']) * rlength + ((xse['K'] ** 2) /
                                                                                   (2 * g * (xse['Area'] ** 2))) *
                                                (-xss['alpha'] * (1 - ki) * ((xse['Area'] / xss['Area']) ** 2) +
                                                 (xse['alpha'] * (1 - ki))))))
            ve = assmq / xse['Area']
            xsecsdf.loc[xsecsdf['XSID'] == xse['XSID'], 'Froude'] = ve / np.sqrt(
                g * xsecsdf.loc[xsecsdf['XSID'] == xse['XSID'], 'mean_depth'])
            vs = assmq / xss['Area']
            xsecsdf.loc[xsecsdf['XSID'] == xss['XSID'],'Froude'] = vs / np.sqrt(
                g*xsecsdf.loc[xsecsdf['XSID'] == xss['XSID'], 'mean_depth'])
            hve = xse['alpha'] * ((ve**2) / (2 * g))
            hvs = xss['alpha'] * ((vs**2) / (2 * g))
            del_hv = hvs - hve
            if del_hv >= 0:
                hf = fall + (del_hv * 0.5)
                k = 0.5
            else:
                hf = fall + del_hv
                k = 0.0
            sf = hf / rlength
            sf_sqrt = np.sqrt(sf)
            compq = wght_K * sf_sqrt
            rch_vals.append([rch, xss['K'], xse['K'], xss['Area'], xse['Area'], xss['alpha'], xse['alpha'], k,
                             rlength, fall, wght_K, hvs, hve, del_hv, hf, sf_sqrt, assmq, compq])

        rchdf = pd.DataFrame(rch_vals, columns=rch_cols)

        # final Q
        if len(rchdf) == 1:
            finq = rchdf['calc_Q']
        elif len(rchdf == 2):
            tfall = rchdf['WSE_fall'].sum()
            K3 = rchdf.loc[1, 'Ke']
            A3 = rchdf.loc[1, 'Ae']
            finq = K3 * np.sqrt((tfall / ((K3 / rchdf.loc[1, 'Ks']) * ((K3 / rchdf.loc[0, 'Ks'])
                                                                       * rchdf.loc[0, 'Length'] +
                                                                       rchdf.loc[1, 'Length']) + ((K3 ** 2) /
                                                                                                  (2 * g *
                                                                                                   (A3 ** 2)))
                                          * ((rchdf.loc[0, 'alpha_s'] * ((A3 / rchdf.loc[0, 'As']) ** 2)
                                              * (1 - rchdf.loc[0, 'k']))
                                             + (rchdf.loc[1, 'alpha_s'] * ((A3 / rchdf.loc[1, 'As']) ** 2)
                                                * (rchdf.loc[1, 'k'] - rchdf.loc[0, 'k'])) + (rchdf.loc[1, 'alpha_e']
                                                                                              * (1 - rchdf.loc[1, 'k']))
                                             ))))
        else:
            tfall = rchdf['WSE_fall'].sum()
            nrs = len(rchdf)
            Kn = rchdf.loc[nrs-1, 'Ke']
            An = rchdf.loc[nrs-1, 'Ae']
            A_lst = []
            B_lst = []
            Bp = (Kn ** 2) / ((An ** 2) * 2 * g)
            for i in range(nrs):
                Aadd = Kn * (rchdf.loc[i, 'Length'] / (rchdf.loc[i, 'Ks'] * rchdf.loc[i, 'Ke']))
                if i == 0:
                    Badd = -rchdf.loc[i, 'alpha_s'] * ((An / rchdf.loc[i, 'As']) ** 2) * (1 - rchdf.loc[i, 'k'])
                else:
                    Badd = rchdf.loc[i, 'alpha_s'] * ((An / rchdf.loc[i, 'As']) ** 2) * (rchdf.loc[i, 'k'] -
                                                                                         rchdf.loc[i-1, 'k'])
                A_lst.append(Aadd)
                B_lst.append(Badd)
            B_lst.append(rchdf.loc[nrs-1, 'alpha_e'] * (1 - rchdf.loc[nrs-1, 'k']))
            A = np.nansum(A_lst)
            Bsm = np.nansum(B_lst)
            B = Bp * Bsm
            finq = Kn * np.sqrt(tfall / (A + B))

        finQs.append([ws, finq])
        chnk_df['q'] = finq * chnk_df['K/Kt']
        chnk_df['v'] = chnk_df['q'] / chnk_df['a']
        xsecsdf['Froude'] = (finq / xsecsdf['Area']) / np.sqrt(g * xsecsdf['mean_depth'])
        sup_results[ws] = {
            'XSection_Chunks': chnk_df,
            'XSection_Calcs': xsecsdf,
            'Reach_Calcs': rchdf
        }

    return pd.DataFrame(finQs, columns=['WS_Profile', 'Q']), sup_results
