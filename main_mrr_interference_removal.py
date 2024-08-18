#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python program to perform the fuzzy-based interference line removal algorithm for MRR-2 data
Copyright (C) 2022 Kwonil Kim (KNU)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Syntax:  python <script>.py <config>.yaml <input_filename> <output_directory>

Example: python config_mrrqc.yaml mrrqc.py /foo/bar/0531.nc /path/to/save/output/without_yearmonth/

Important Notes:
 - Input file must be post-processed by Maahn and Kollias (2012) algorithm (i.e. IMProToo).
 - This altorithm assumed the rain event. Therefore, the algorithm does not process de-aliased fields, which often gives unreliable value for rain. The users should read '~~_noDA' fields.

Version history:
1.0 / Jun 2022 (kwonil): Created
1.1 / Aug 2024 (kwonil): Add a capability to update Ze MF, Fix typos, Version for the first public release

"""
__author__ = "Kwonil Kim, KNU"
__copyright__ = "Copyright (C) 2022 Kwonil Kim"
__license__ = "GPL3.0"
__version__ = "1.1"

import sys
import getopt
import copy
import kkpy
import numpy as np
import pandas as pd
import scipy.signal
import yaml
import os
import os.path
import time
import logging
import logging.config

if len(sys.argv) == 4:
    fname_cfg = sys.argv[1]
    
    # Check if extension is yaml
    _, cfg_ext = os.path.splitext(fname_cfg)
    if cfg_ext not in '.yaml':
        raise ValueError(f'Unrecognized config file: {fname_cfg}')
        sys.exit()
    
    # Read config file
    try:
        cfgs = yaml.safe_load(open(fname_cfg))
    except:
        raise ValueError(f'Failed to load {fname_cfg}')
        sys.exit()
else:
    # Syntax helper
    raise ValueError('\n\n' + \
                     '===================================================================================================\n' + \
                     'Syntax:  python <script>.py <config>.yaml <input_filename> <output_directory>\n' + \
                     'Example: python mrrqc.py config_mrrqc.yaml /foo/bar/0531.nc /path/to/save/output/without_yearmonth/\n' + \
                     '===================================================================================================\n')
    sys.exit()

# Enable logging
logging.config.dictConfig(cfgs['logging'])
logger = logging.getLogger('mrrqc')

    
def read_mrr(fnames):
    """
    Read MRR (nc).
    
    Examples
    ---------
    >>> ds_mrr = kkpy.io.read_mrr('0723.nc')
    
    >>> ds_mrr = kkpy.io.read_mrr('0723.nc', degree='')
    
    Parameters
    ----------
    fnames : str or array_like
        The MRR filename(s).
    
    Returns
    ---------
    ds : xarray dataset object
        Return MRR data.
    """
    
    return _read_mrr_mk(fnames)
    
    
def _read_mrr_mk(fnames):
    import xarray as xr
    return xr.open_mfdataset(
        fnames,
        decode_times=True,
        combine='nested'
    )


def membership_func(x, var):
    """
    Return a membership funcion value for a given x
    
    Input
    -------
    x: float or array-like
        X value
    var: str
        'Ze', 'Vr', 'PL', 'PW', 'PS', 'PSs', 'CS', 'CSs'
    
    Return
    --------
    value: float or array-like
        Membership function value
    """
    def MF_SM(x, params):
        """
        Membership function of sigmoid shape
        """
        return 1 / (1 + np.exp(params[0]*(x-params[1])))
    
    params_mf = cfgs['fuzzy']['parameter']
    if 'Ze' in var:
        logger.debug(f"params_mf: {params_mf}")
    
    return MF_SM(x, params_mf[var])

def get_membership_variables(ds):
    cfg = cfgs['mv']
    logger.debug(f"cfgs['mv']: {cfg}")
    
    try:
        PS = abs(ds['leftSlope_noDA'].values - ds['rightSlope_noDA'].values)/2.
        PSs = kkpy.util.nanconvolve2d(
            PS,
            np.zeros(cfg['PSs']['window_size'])+1) \
            / (cfg['PSs']['window_size'][0] * cfg['PSs']['window_size'][1])
        CS = scipy.signal.convolve2d(
            np.isfinite(ds['W_noDA'].values),
            np.zeros(cfg['CS']['window_size'])+1,
            boundary='symm',
            mode='same') \
            / (cfg['CS']['window_size'][0] * cfg['CS']['window_size'][1])
        CSs = kkpy.util.nanconvolve2d(
            CS,
            np.zeros(cfg['CSs']['window_size'])+1) \
            / (cfg['CSs']['window_size'][0] * cfg['CSs']['window_size'][1])
        PSs[np.isnan(ds['W_noDA'])] = np.nan
        CS[np.isnan(ds['W_noDA'])] = np.nan
        CSs[np.isnan(ds['W_noDA'])] = np.nan
    except:
        logger.error('Not enough data to calculate membership variable')

    return {
        'Ze': ds['Ze_noDA'].values,
        'Vr': ds['W_noDA'].values,
        'PL': ds['peakVelLeftBorder_noDA'].values,
        'PW': ds['peakVelRightBorder_noDA'].values - ds['peakVelLeftBorder_noDA'].values,
        'PS': PS,
        'PSs': PSs,
        'CS': CS,
        'CSs': CSs,
    }

def get_updated_mf_weight(value, data, null):
    from scipy.signal import savgol_filter, find_peaks, peak_widths
    
    def MF_GBs(x, params):
        """
        Membership function of generalized bell shape with multiple peaks
        """
        res = np.zeros(x.shape)
        for param in params:
            res += 1 / (1 + np.abs((x-param[2])/param[0])**(2*param[1]))
        res[res > 1] = 1
        return res

    # Logistic regression with sigmoid function
    null = null[np.isfinite(null)]
    data = data[np.isfinite(data)]

    if data.size > null.size:
        data = np.random.choice(data, size=null.size)
    else:
        null = np.random.choice(null, size=data.size)

    features = np.hstack([null, data])
    labels = np.hstack([np.ones(null.size), np.zeros(data.size)])

    # Histogram of reflectivity of interference lines
    hist, hist_bins = np.histogram(
        null,
        bins=np.arange(-20,50.1,1),
        weights=np.ones_like(null)/float(len(null)))
    hist_bins = np.mean([hist_bins[1:],hist_bins[:-1]], axis=0)
    # Smooth histogram
    hist = savgol_filter(hist, 11, 2)
    hist[hist < 0] = 0
    
    # Find peak locations, peak widths
    peaks, _ = find_peaks(hist, height=0.02, distance=10)
    widths = peak_widths(hist, peaks)
    
    # Generalized-Bell shape parameters 
    params_multipeak = [
        [w/1.5,3,p]
        for p, w in zip(hist_bins[peaks], widths[0])
    ]

    MF = MF_GBs(value, params_multipeak)
    weight = np.corrcoef(labels, MF_GBs(features, params_multipeak))[1][0]

    return MF, weight

def get_mask_fuzzy(dict_data, thr_A=cfgs['fuzzy']['thr_A'], weights=cfgs['fuzzy']['weight']):
    logger.debug(f"thr_A: {thr_A}")
    logger.debug(f"weights: {weights}")
    logger.debug(f"min_hgt_idx: {cfgs['fuzzy']['min_hgt_idx']}")

    MFs = [
        membership_func(dict_data[var], var) \
        for var in dict_data.keys()
    ]
    
    # Aggregate
    A = np.sum([
        MFs[i_v] * weights[var]
        for i_v, var in enumerate(dict_data.keys())
    ], axis=0) / sum(weights.values())

    # Update Ze MF, weight
    # Only if we have enough data (here we use 30)
    minhgtidx = cfgs['fuzzy']['min_hgt_idx']
    candidate_interference_line = dict_data['Ze'][:,minhgtidx:][np.where(A[:,minhgtidx:] > 0.5)]
    candidate_precipitation = dict_data['Ze'][:,:minhgtidx][np.where(A[:,:minhgtidx] < 0.5)]

    cnt_i = np.sum(np.isfinite(candidate_interference_line))
    cnt_p = np.sum(np.isfinite(candidate_precipitation))
    
    if np.logical_and(cnt_i >= 30, cnt_p >= 30):
        # Update Ze MF, weight
        new_MF_Ze, new_weight_Ze = get_updated_mf_weight(
            dict_data['Ze'],
            candidate_precipitation,
            candidate_interference_line
        )
        MFs[0] = new_MF_Ze
        
        # Update A
        A = (A * sum(weights.values()) \
            - (membership_func(dict_data['Ze'], 'Ze') * weights['Ze']) \
            + (new_MF_Ze * new_weight_Ze)) \
            / (sum(weights.values())-weights['Ze']+new_weight_Ze)

    # Height index
    hgtidx = np.zeros(dict_data['Ze'].shape)
    hgtidx[:,:] = np.arange(dict_data['Ze'].shape[1])
    
    # Masking, higher thr_MF -> weaker QC
    mask = (
        (A > thr_A) & \
        (hgtidx >= minhgtidx)
    )
    
    return mask, MFs, A

def despeckle(arr2d, kernel, thr_cnt):
    cntarr = scipy.signal.convolve2d(
        np.isfinite(arr2d),
        kernel,
        mode='same',
        boundary='symm'
    )
    arr2d[cntarr <= thr_cnt] = np.nan

    return arr2d

def process_fuzzy_qc(ds):
    logger.debug(f'Entering get_membership_variables()')
    dict_mv = get_membership_variables(ds)
    
    logger.debug(f'Entering get_mask_fuzzy()')
    mask, MFs, A = get_mask_fuzzy(dict_mv)

    # Add membership function values and AGG to output dataset
    for i_k, key in enumerate(dict_mv.keys()):
        ds[f'MF_{key}'] = (('time','range'), MFs[i_k])
    ds['AGG'] = (('time','range'), A)

    # Eliminate interference lines
    keys = [
        'Ze_noDA', 'W_noDA', 'spectralWidth_noDA',
        'skewness_noDA', 'kurtosis_noDA',
        'peakVelLeftBorder_noDA', 'peakVelRightBorder_noDA',
        'leftSlope_noDA', 'rightSlope_noDA'
    ]
    for key in keys:
        logger.debug(f'Masking {key}')
        arr = ds[key].values
        arr[mask] = np.nan
        ds[key].data = arr

    return ds

def process_despeckle(ds):
    cfg = cfgs['despeckle']
    logger.debug(f"cfgs['despeckle']: {cfg}")
    
    keys = [
        'Ze_noDA', 'W_noDA', 'spectralWidth_noDA',
        'skewness_noDA', 'kurtosis_noDA',
        'peakVelLeftBorder_noDA', 'peakVelRightBorder_noDA',
        'leftSlope_noDA', 'rightSlope_noDA'
    ]
    for i_k, key in enumerate(keys):
        logger.debug(f'Despeckling {key}')
        ds[key].data = despeckle(
            ds[key].values,
            np.zeros(cfg['window_size'])+1,
            cfg['thr_cnt']
        )
    
    return ds

def process_save_result(ds, outdir):
    dt = pd.to_datetime(ds['time'].values[0])
    logger.debug(f"dt: {dt}")
    
    outdir = f'{outdir}/{dt:%Y%m}'
    
    outfname = f'{outdir}/{dt:%m%d}.nc'
    logger.debug(f"outfname: {outfname}")
    
    os.makedirs(outdir, exist_ok=True)
    ds.to_netcdf(
        outfname,
        format='NETCDF3_CLASSIC'
    )
    
    return outfname

def main(argv):
    starttime = time.time()
    
    # Check arguments
    fname = argv[2]
    logger.info(f'INPUT: {fname}')
    outpath = argv[3]
    logger.debug(f'outpath: {outpath}')

    _, fname_ext = os.path.splitext(fname)
    if fname_ext not in '.nc':
        logger.warning(f'Recieved unrecognized file format: {fname_ext}. Make sure the input is correct.')
    if not os.path.isfile(fname):
        logger.error(f'Check if file exists: {fname}')
        return
    if not os.access(outpath, os.W_OK):
        logger.error(f'Check write permission: {outpath}')
        return
    
    # Reading mrr file
    logger.debug(f'################ read_mrr() ################')
    ds_mrr = read_mrr(fname)

    # Reading config
    cfg = cfgs['main']
    logger.debug(f"cfgs['main']: {cfg}")
    
    # Fuzzy-based removal of interference lines
    if cfg['fuzzy']:
        logger.debug(f'################ process_fuzzy_qc() ################')
        ds_mrr = process_fuzzy_qc(ds_mrr)

    # Despeckling
    if cfg['despeckle']:
        logger.debug(f"################ process_despeckle() ################")
        ds_mrr = process_despeckle(ds_mrr)

    # Save result
    if cfg['output']:
        logger.debug(f"################ process_save_result() ################")
        outfname = process_save_result(ds_mrr, outpath)
        endtime = time.time()
        logger.info(f'OUTPUT: {outfname}')
        logger.info(f'{endtime-starttime:.2f} seconds elapsed')
    else:
        endtime = time.time()
        logger.info(f'FINISHED (dry run)')
        logger.info(f'{endtime-starttime:.2f} seconds elapsed')
    
    return


if __name__ == "__main__":
    main(sys.argv)
