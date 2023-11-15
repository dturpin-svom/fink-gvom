#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 19:07:22 2023

@author: dt270490
"""

import numpy as np
import pandas as pd
import utils
import astropy.units as u
from astroplan import AltitudeConstraint, AirmassConstraint, AtNightConstraint,\
MoonSeparationConstraint
from astroplan import Observer
from astropy.time import Time
from fink_too_plots import plot_common_sky, make_overlap_histo
import requests
import io
import gvom_filters
import fink_too_plots

def get_fink_alert():
    
    # Get latests alerts
    r = requests.post(
      'https://fink-portal.org/api/v1/latests',
      json={
        'class': 'Unknown',
        'n': '1000'
      }
    )
    
    # Format output in a DataFrame
    pdf = pd.read_json(io.BytesIO(r.content))
    return pdf

# Load a sample of Fink alert tests
pdf_test = get_fink_alert()
# Load the observatories info
obs_filename = "gvom_core_network.csv"
observatories = utils.load_observatories(obs_filename)
# Define the observational constraints
constraints = [
    AltitudeConstraint(30 * u.deg, 90 * u.deg),
    AirmassConstraint(2),
    MoonSeparationConstraint(min=20* u.deg)]#,AtNightConstraint.twilight_astronomical()]


all_sky_coords = utils.make_sky_grid_dict(pdf_test['i:ra'].values,
                                          pdf_test['i:dec'].values,
                                          Time(pdf_test['i:jd'].values,
                                               format='jd'))

pdf_cand_visibility = utils.obs_visibility(all_sky_coords['ra'],
                                      all_sky_coords['dec'],
                                      all_sky_coords['date_jd'],
                                      observatories,
                                      constraints,
                                      True)

# =====================================================================
# Apply the observational strategy
# =====================================================================
# Step 1: Only keep the real sources (RB score > 0.9)
mask_rb = gvom_filters.rb_mask(pdf_test,0.9)
# Step 2: Only keep the extragalactic sources

# mask the source at low galactic latitudes
mask_extragal = gvom_filters.extragal_mask(pdf_cand_visibility,
                                           15)

# Step 3: Only keep the sources both visible for more than 2 hours at SPM,
# Xinglong and Canaria Island
mask_visibility = gvom_filters.net_visibity_mask(pdf_cand_visibility,
                                                 2/24)

# =====================================================================
# Apply the science strategy
# ===================================================================== 
#Step 1 Select only the most interesting class of transients
mask_class = gvom_filters.class_mask(pdf_test,['Unknown',
                                               'SN candidate',
                                               'Ambiguous'])

#Step 2 Select only the brightest transients
mask_brightness = gvom_filters.mask_brightness(pdf_test,17.0)

# Step 3 Select only the transients "short" living transients
mask_alert_history = gvom_filters.alert_hist_mask(pdf_test,5)


mask_final_obs = mask_rb & mask_extragal & mask_visibility
mask_final_science = mask_alert_history & mask_class & mask_brightness
mask_final = mask_final_obs & mask_final_science
obs_filter_cand =\
[
    len(pdf_cand_visibility),
    len(pdf_cand_visibility[mask_rb]),
    len(pdf_cand_visibility[mask_rb&mask_extragal]),
    len(pdf_cand_visibility[mask_final_obs]),
    len(pdf_cand_visibility[mask_final_obs&mask_class]),
    len(pdf_cand_visibility[mask_final_obs&mask_class&mask_brightness]),
    len(pdf_cand_visibility[mask_final]),
]
obs_filter_cand_pdf = pd.DataFrame(obs_filter_cand)
index_mask = 2
list_masks = [[],mask_rb,mask_rb&mask_extragal,mask_final_obs,
              mask_final_obs&mask_class,
              mask_final_obs&mask_class&mask_brightness,
              mask_final]
#List the mask type names
mask_names = ['Original total', 
              'After RB filter',
              'After RB+Extragal filter', 
              'After RB+Extragal+Visibility filter (Obs.)',
              'After Obs.+Class filter',
              'After Obs.+Class+Brightness filter',
              'Obs. + Class + Brightness + Det. history filters']

list_objectId = pdf_test[list_masks[index_mask]]["i:objectId"].unique().tolist()
list_ras = pdf_test[list_masks[index_mask]]["i:ra"].values
list_decs = pdf_test[list_masks[index_mask]]["i:dec"].values

# Plot the rejection efficiency and the skymap of the remaining candidates
fink_too_plots.plot_filter_cand_hist(obs_filter_cand_pdf,False)
fink_too_plots.plot_cand_sky(obs_filter_cand_pdf.iloc[[0,index_mask]].reset_index(),
                              mask_names[index_mask],
                              [list_ras,
                              list_decs,
                              list_objectId],
                              'gal',
                              False)
