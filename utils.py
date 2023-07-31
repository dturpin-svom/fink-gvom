#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:37:49 2022

@author: dt270490
"""

import os
from astroplan import Observer, FixedTarget
from astroplan import is_observable, is_always_observable, observability_table
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
import numpy as np
import pandas as pd
import requests
import gvom_filters

def equtorad_coord(ras:np.ndarray,decs:np.ndarray):
    """
    Convert RA dec into theta,phi radian angles following the rule below
    0°<RA<+180° -> RA is positive
    RA>=+180° -> RA is negative
    DEC is ranging from 0:180° (DEC=-90° -> 180° / DEC=+90° -> 0°)

    Parameters
    ----------
    ras : np.ndarray
        array of Right ascencion.
    decs : np.ndarray
        array of Declination.

    Returns
    -------
    theta : np.ndarray
        DEC radian angles projetction.
    phi : np.ndarray
        RA radian angles projetction.

    """
    theta = np.deg2rad((decs * -1) + 90)
    phi1 = np.deg2rad(np.array(ras)[(np.array(ras)<=180.)])
    phi2 = np.deg2rad(np.array(ras)[(np.array(ras)>180.)]-360)
    phi = np.concatenate((phi1,phi2))
    
    return theta, phi


def load_observatories(obs_path_file: str):
    """
    Load the observatories for which we want to crossmacth the observable 
    skies

    Parameters
    ----------
    obs_path_file : str
        Path to the files where the observatory settings are stored.

    Returns
    -------
    observatories : Dataframes
        pandas Dataframe of the observatory configurations 
        to feed the astroplan.Observer package.

    """
    observatories = pd.read_csv(obs_path_file, delimiter=";")
    return observatories    
    
def simu_night_time_interval(ref_obs, ref_date: str, n_days: int,
                             day_bin: int):
    """
    Compute the list of the date interval during which the observable skies
    will be crossmatched (TB updated)

    Parameters
    ----------
    ref_obs : astroplan.observer.Observer
        astroplan.observer.Observer object for the Observatory chosen as the 
        reference observatory.
    ref_date : str
        Date of reference to start the simulation.
    n_days : int
        Number of simulated days.
    day_bin : int
        Day interval between two simulated nights.

    Returns
    -------
    None.

    """  
    # compute the start time of the nearest (compared to the detection time)
    # night
    date_start_night = ref_obs.twilight_evening_astronomical(ref_date,
                                                                which='next')
    
    # compute the morning time of the nearest (compared to the detection time)
    # day
    date_end_night = ref_obs.twilight_morning_astronomical(ref_date,
                                                                which='next')
    
    if isinstance(date_start_night.value,float) and\
        isinstance(date_end_night.value,float):
        # If the detection time is within the current night, use the detection
        # time as a starting date and the end of the night as an ending date
        if date_start_night < ref_date and date_end_night>ref_date:
            if date_end_night.jd - date_start_night.jd >0.5:
                date_start_night = ref_obs.twilight_evening_astronomical(ref_date,
                                                                        which='next')
                time_ranges = Time([date_start_night.iso, date_end_night.iso])
                print('Case 1a')
            else:
                time_ranges = Time([ref_date.iso, date_end_night.iso])
                print('Case 1b')
        elif date_start_night < ref_date and date_end_night<ref_date:
            # If the detection time is after the nearest night, use the next
            # night  starting date as a starting date and the end of the night
            # as an ending date
            print('Case 2')
            date_start_night = ref_obs.twilight_evening_astronomical(ref_date,
                                                                     which='next')
            date_end_night = ref_obs.twilight_morning_astronomical(date_start_night,
                                                                   which='next')
            time_ranges = Time([date_start_night.iso, date_end_night.iso])
        # If the detection time is during day time, use the night starting date
        # as a starting date and the end of the night as an ending date
        elif date_start_night > ref_date and date_end_night>ref_date:
            print('Case 3')
            date_start_night = ref_obs.twilight_evening_astronomical(ref_date,
                                                                     which='next')
            date_end_night = ref_obs.twilight_morning_astronomical(date_start_night,
                                                                   which='next')
            time_ranges = Time([date_start_night.iso, date_end_night.iso])
    # elif isinstance(ref_obs.name,np.ndarray):
    #     # If the detection time is within the current night, use the detection
    #     # time as a starting date and the end of the night as an ending date
    #     mask_night_start_before = date_start_night < ref_date
    #     mask_night_start_after = date_start_night > ref_date
    #     mask_night_end_before = date_end_night < ref_date
    #     mask_night_end_after = date_end_night > ref_date
    #     if mask_night_start_before.any() and mask_night_end_after.any():
    #         common_mask = mask_night_start_before & mask_night_end_after
    #         time_ranges_1 = Time([ref_date.iso, date_end_night[common_mask].iso])
            
    #     elif date_start_night < ref_date and date_end_night<ref_date:
    #         # If the detection time is after the nearest night, use the next
    #         # night  starting date as a starting date and the end of the night
    #         # as an ending date
    #         date_start_night = ref_obs.twilight_evening_astronomical(ref_date,
    #                                                                  which='next')
    #         date_end_night = ref_obs.twilight_morning_astronomical(date_start_night,
    #                                                                which='next')
    #         time_ranges = Time([date_start_night.iso, date_end_night.iso])
    #     # If the detection time is during day time, use the night starting date
    #     # as a starting date and the end of the night as an ending date
    #     elif date_start_night > ref_date and date_end_night>ref_date:
    #         time_ranges = Time([date_start_night.iso, date_end_night.iso])
    
    else:
        time_ranges = []
    
   
    # time_ranges = Time([ref_date.iso, (ref_date+1*u.second).iso])
    # print(time_ranges)
    return time_ranges


def build_targets(ras:np.array,decs:np.array):
    """
    Make the target objects to be ingested by the astroplan is_observable
    functions

    Parameters
    ----------
    ras : np.array
        Numpy arrays of right ascencion in units of degrees.
    decs : np.array
        Numpy arrays of declinations in units of degrees.

    Returns
    -------
    targets : 
        List of targets for which we want to estimate the observability

    """
    ras_grid, decs_grid = np.meshgrid(ras, decs)
    
    target_table = Table()
    target_table['ra'] = np.reshape(ras_grid,
                                    ras_grid.shape[0]*ras_grid.shape[1])
    target_table['dec'] = np.reshape(decs_grid,
                                      decs_grid.shape[0]*decs_grid.shape[1])
    
    targets = FixedTarget(coord=SkyCoord(ra=target_table['ra']*u.deg,
                                          dec=target_table['dec']*u.deg))
              
    
    return targets, target_table


def make_visibility_masks(constraints, observatory, targets, time_range: list):
    """
    Build the skymap mask regarding the obsverational constraints:  
    True = the position is visible
    False = the position is not visible

    Parameters
    ----------
    constraints : TYPE
        DESCRIPTION.
    observatory : TYPE
        DESCRIPTION.
    targets : TYPE
        DESCRIPTION.
    time_range : list
        time range during which the observability is estimated

    Returns
    -------
    mask_visibility : np.array
        Numpy arrays of booleans masking the invisible targets

    """
    mask_visibility = is_observable(
        constraints, observatory, targets, time_range=time_range,
        time_grid_resolution=1*u.hour)
    

    return mask_visibility

def obs_visibility(ras:list,decs:list,det_time:list,
                   observatories:pd.DataFrame,
                   constraints:list,instant_sky:bool):
    """
    This function builds a pandas DataFrame that collects all the visibility
    flags of the observatories in the network related to the positions of 
    optical transients and based on pre-defined observational constraints

    Parameters
    ----------
    ras : list
        list of right ascencion given in degrees.
    decs : list
        list of declination given in degrees.
    det_time : list
        list of first detection times expressed in Julian Date.
    observatories : pd.DataFrame
        Dataframe gathering the observatories informations.
    constraints : list
        list of astronomical constraints to be satisfied by each transient
        candidates.
    instant_sky : bool
        Flag if one wants the instantaneous visibility (True) or the 
        visibility during the entire night.

    Returns
    -------
    df : pd.DataFrame
        The Dataframe related to the observatories visibility for each candid.

    """
    mask_visibilities = []
    fraction_of_time_visibilities = []
    #---- prepare the targets
    targets = FixedTarget(coord=SkyCoord(ra=ras * u.deg,dec=decs* u.deg))
    #---- only if instantaneous sky is aksed
    if instant_sky:
        time_ranges = Time([det_time.iso, (det_time+1*u.second).iso])
    # Initialize the DataFrame     
    col_name = observatories.obs_name.values.tolist()
    df = pd.DataFrame(columns=col_name)
    for i in range(len(observatories)):
        try:
            observatory = Observer.at_site(observatories.obs_name[i])
        except:
            observatory = Observer(
                longitude=observatories.longitude.values[i] * u.deg,
                latitude=observatories.latitude.values[i] * u.deg,
                elevation=observatories.elevation.values[i] * u.m,
                name=observatories.obs_name.values[i],
            )
        
        if not instant_sky:
            time_ranges = simu_night_time_interval(
                    observatory, Time(det_time,format='jd'), 1, 1)
        print(observatory.name,time_ranges)
        if not time_ranges:
            mask_visibilities.append(False)
            fraction_of_time_visibilities.append(0.0)
        else:
            mask_visibility = make_visibility_masks(
                        constraints, observatory, targets, time_ranges)
            obs_table = observability_table(constraints, observatory,
                                            targets.coord,
                                            time_range=time_ranges,
                                            time_grid_resolution=0.5*u.hour)
            fraction_of_time_visibilities.append(obs_table['fraction of time observable'])
            mask_visibilities.append(mask_visibility)
   
        df[observatories.obs_name.values[i]] = mask_visibilities[i]
        df[observatories.obs_name.values[i]+'_max_visible_time'] = fraction_of_time_visibilities[i]
    # add the date to the data frame
    if not time_ranges:
        df['date_start'] = 'n/a'
        df['date_end'] = 'n/a'
    else:
        df['date_start'] = time_ranges[0].isot
        df['date_end'] = time_ranges[1].isot
    # add the equatorial coordinates to the data frame
    df['ra'] = ras
    df['dec'] = decs
    # add the galactic coordinates to the data frame
    coord = SkyCoord(ra=ras * u.deg, dec=decs * u.deg)
    df['l'] = coord.galactic.l.degree
    df['b'] = coord.galactic.b.degree
    mask_extragal = gvom_filters.extragal_mask(df, 15.0)
    # add the result of the extragalactic mask
    df['mask_extragal'] = mask_extragal
    return df

def network_stats(pdf_cand_visibility_ref:pd.DataFrame,
                  pdf_cand_visibility:pd.DataFrame,observatories:list,
                  save:bool):
    """
    Write into a file the results of the applied selection cuts on a bunch of
    targets or a sky fraction

    Parameters
    ----------
    pdf_cand_visibility_ref : pd.DataFrame
        Visibility of the targets by night observed at the reference Observatory.
    pdf_cand_visibility : pd.DataFrame
        Visibility of the targets by night observed at the test Observatories.
    observatories : list
        list of observatory names.
    save : bool
        Save the file. True or False.

    Returns
    -------
    d : pd.DataFrame
        Results formated into a pd.DataFrame.

    """
    d = {}
    net_stats_write = {}
    for obs_name in observatories.obs_name:
        # build the individual sky visibility mask
        mask_sky = pdf_cand_visibility[obs_name] == True
        # build the common sky visibility mask between VRO and a given obs.
        mask_common_sky = pdf_cand_visibility_ref['VRO']&\
                          pdf_cand_visibility[obs_name]
        # build the extragalactic mask
        mask_extragal = gvom_filters.extragal_mask(pdf_cand_visibility, 15.0)
        #build the final mask
        mask_final = mask_common_sky & mask_extragal
        # Compute the fraction of the sky commonly visible at VRO and SPM
        sky_frac_visibility = mask_final.sum()/pdf_cand_visibility_ref['VRO'].sum()
        mask_contents = {}
        mask_contents['mask_sky'] = mask_sky 
        mask_contents['mask_extragal'] = mask_extragal
        mask_contents['mask_common_sky'] = mask_common_sky 
        mask_contents['sky_frac_visibilities'] =  sky_frac_visibility
        
        d[obs_name] = mask_contents
        net_stats_write = {'sky_frac_visibility': [sky_frac_visibility],
              'date_start': [pdf_cand_visibility.date_start.values[0]],
              'date_end': [pdf_cand_visibility.date_end.values[0]]}

        df = pd.DataFrame(data=net_stats_write)
        if save:
            outdir = '/media/dt270490/Transcend/Workspace/SVOM/LSST-fink/FINK_ToO/'+\
                'network_sky_visibility/'+obs_name+'/'
            if not os.path.isdir(outdir):
                os.mkdir(outdir, mode = 0o777)
            if os.path.isfile(outdir+"night_skyfrac_VRO.csv"):
                df.to_csv(outdir+"night_skyfrac_VRO.csv",index=False,
                          mode='a',header=False)
            else:
                df.to_csv(outdir+"/night_skyfrac_VRO.csv",index=False,
                          mode='w')
                    
    return d

def make_sky_grid(coord_bin:int):
    """
    Define the sky RA and Dec grid with an angular bin given by coord_bin

    Parameters
    ----------
    coord_bin : integer
        bining of the RA and DEC grids.

    Returns
    -------
    ras_vec: np.array
        An array containing all the RA of the grid
    decs_vec: np.array
        An array containing all the Dec of the grid

    """
    ras = np.linspace(0, 360, coord_bin)
    decs = np.linspace(-90, 90, coord_bin)
    xx, yy = np.meshgrid(ras, decs)
    ras_vec = np.reshape(xx,xx.shape[0]*xx.shape[1])
    decs_vec = np.reshape(yy,yy.shape[0]*yy.shape[1])
    return ras_vec, decs_vec

def make_sky_grid_dict(ras_vec:list,decs_vec:list,ref_date:Time):
    """
    Build the dictionnary of all the sky coordinates to be ingested by the 
    sky observability routine

    Parameters
    ----------
    ras_vec : list
        List of RAs.
    decs_vec : list
        List of DECs.
    ref_date : Time
        Astropy.time.Time  datetime at which the grid is computed.

    Returns
    -------
    d : Dict
        List all the RA/DEC/Datetime combination of the sky grid.

    """
    # Compute the time at which the observatory visibility of 
    # each sky coordinates has to be computed
    times = np.zeros(len(decs_vec)).tolist()+ref_date
    d={
    'ra' : ras_vec,
    'dec' : decs_vec,
    'date_jd' : times[0]
        }
    return d
    