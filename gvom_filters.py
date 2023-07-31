#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:32:05 2023

@author: dt270490
"""


import pandas as pd

def rb_mask(pdf_cand_visibility:pd.DataFrame,rb_threshold:float):
    """
    Create a mask related to the Real/Bogus score of the optical transient
    candidates

    Parameters
    ----------
    pdf_cand_visibility : pd.DataFrame
        DESCRIPTION.
    rb_threshold : float
        threshold on the RB score to keep a transient for follow-up.

    Returns
    -------
    mask : list
        DataFrame masking the candidates not considered as real sources.

    """
    if "i:rb" in pdf_cand_visibility.columns:
        mask_rb = pdf_cand_visibility["i:rb"] > rb_threshold 
    else:
        mask_rb = pdf_cand_visibility.rb > rb_threshold
    return mask_rb

def extragal_mask(pdf_cand_visibility:pd.DataFrame,gal_lat_threshold:float):
    """
    Create a mask related to the galactic latitude of the optical transient
    candidates

    Parameters
    ----------
    pdf_cand_visibility : pd.DataFrame
        DESCRIPTION.
    gal_lat_threshold : float
        threshold on the galactic latitude to keep a transient for follow-up.

    Returns
    -------
    mask : list
        DataFrame masking the candidates too close to the galactic plane.

    """
    
    mask_north = pdf_cand_visibility.b > gal_lat_threshold 
    mask_south = pdf_cand_visibility.b< -1*gal_lat_threshold
    return mask_north | mask_south

def net_visibity_mask(pdf_cand_visibility:pd.DataFrame,duration_min):
    """
    Build the final mask to filter out the transient candidates

    Parameters
    ----------
    pdf_cand_visibility : pd.DataFrame
        DESCRIPTION.

    Returns
    -------
    mask : TYPE
        DESCRIPTION.

    """
    # The source has to visible at San Pedro Martir and Xinglong and NOT
    mask_net_visibility = (pdf_cand_visibility['OAN_SPM'] == True) & \
        (pdf_cand_visibility['Xinglong'] == True) &\
            (pdf_cand_visibility['ORM'] == True)
    mask_time_visibility= \
        (pdf_cand_visibility['OAN_SPM_max_visible_time'] > duration_min) &\
            (pdf_cand_visibility['Xinglong_max_visible_time'] > duration_min) &\
                (pdf_cand_visibility['ORM_max_visible_time'] > duration_min)
    
    return mask_net_visibility & mask_time_visibility

def mask_brightness(pdf_cand_visibility:pd.DataFrame,mag_threshold:float):
    """
    Build the mask that keeps only the brigthest candidates according to the 
    chosen mag_threshold

    Parameters
    ----------
    pdf_cand_visibility : pd.DataFrame
        DESCRIPTION.
    mag_threshold : float
        DESCRIPTION.

    Returns
    -------
    mask_brightness : TYPE
        DESCRIPTION.

    """
    if "i:magpsf" in pdf_cand_visibility.columns:
        mask_brightness = pdf_cand_visibility["i:magpsf"]<mag_threshold
    else:
        mask_brightness = pdf_cand_visibility.magpsf<mag_threshold
    
    return mask_brightness

def alert_hist_mask(pdf_cand_visibility:pd.DataFrame,jd_threshold:float):
    """
    Create a mask related to the prior detections of the optical transient
    candidates. We don't search for long transients or variable objects.

    Parameters
    ----------
    pdf_cand_visibility : pd.DataFrame
        DESCRIPTION.
    jd_threshold : float
        threshold on the prior detection delay.

    Returns
    -------
    mask : list
        DataFrame masking the candidates which have a too long detection 
        history.

    """
    if "i:jd" in pdf_cand_visibility.columns:
        mask_alert_hist = pdf_cand_visibility["i:jd"] - \
            pdf_cand_visibility["i:jdstarthist"] < jd_threshold
    else:
        mask_alert_hist = pdf_cand_visibility.jd - \
            pdf_cand_visibility.jdstarthist < jd_threshold
    return mask_alert_hist

def class_mask(pdf_cand_visibility:pd.DataFrame,class_list:list):
    """
    Create a mask related to the prior detections of the optical transient
    candidates. We don't search for long transients or variable objects.

    Parameters
    ----------
    pdf_cand_visibility : pd.DataFrame
        DESCRIPTION.
    class_list : list
        lis of transient classes to keep.

    Returns
    -------
    mask : list
        DataFrame masking the candidates which have the non desired 
        fink classification.

    """
    index_class = 0
    while index_class <len(class_list):
        if "v:classification" in pdf_cand_visibility.columns:
            mask_class = pdf_cand_visibility["v:classification"] == class_list[index_class]
        else:
            mask_class = pdf_cand_visibility.finkclass == class_list[index_class]
        if index_class == 0:
            mask_class_final = mask_class
        else:
            mask_class_final = mask_class_final | mask_class
        index_class = index_class+1
    return mask_class_final