#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 18:49:13 2023

@author: dt270490
"""

import os
from healpy.newvisufunc import projview, newprojplot
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif Times New Roman')
alpha_common_sky = 0.1
import matplotlib.patches as mpatches
from mhealpy import HealpixMap
import mhealpy as hmap
import numpy as np
import healpy as hp
import random
from astropy.coordinates import get_moon, get_sun
import matplotlib as mpl
from astropy.coordinates import SkyCoord
import astropy.units as u
import utils

def get_color():
    r = random.random()
    
    b = random.random()
    
    g = random.random()
    
    color = (r, g, b)
    return color

def get_cmap_gradient_color(cmap_name):
    
    cmap = mpl.colormaps[cmap_name].resampled(64)
    new_cmap = cmap(np.linspace(0, 1, 64))
    
    return new_cmap


def plot_common_sky(pdf_cand_visibility:pd.DataFrame,
                    pdf_cand_visibility_ref:pd.DataFrame,
                    obs_name_ref,
                    obs_name_test,
                    net_stats,
                    date,
                    coord_syst:str,
                    show_indiv_sky:bool,
                    save:bool):
    """
    Plot the overlapping (and individual) sky regions between the sky visible
    by a reference Observatory at night and
    a test observatory during the next night

    Parameters
    ----------
    pdf_cand_visibility : pd.DataFrame
        DESCRIPTION.
    pdf_cand_visibility_ref : pd.DataFrame
        DESCRIPTION.
    obs_name_ref : TYPE
        DESCRIPTION.
    obs_name_test : TYPE
        DESCRIPTION.
    net_stats : TYPE
        DESCRIPTION.
    date : TYPE
        DESCRIPTION.
    coord_syst : str
        DESCRIPTION.
    show_indiv_sky : bool
        DESCRIPTION.
    save : bool
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    # plot settings
    c_common_sky = 'red'
    alpha_indiv_sky = 0.2
    fontsize_tick = 20
    fontsize_label = 24
    figsize = 19
    delta_date = (date[1] - date[0]).jd*24
    title_sky = (date[0].iso+"+ "+str(delta_date)+" hr [UTC]\n"+\
                 str("%.2f" % (net_stats[obs_name_test]\
                    ['sky_frac_visibilities']*100))+ " \% of the "+\
                obs_name_ref+" sky visible at "+obs_name_test)
                
    sdf_map = hp.read_map('/media/dt270490/Transcend/Workspace/GRB_data/'+\
                          'GRBase/catalogs/dust/EBV_SFD98_1_512.fits')
    if coord_syst =='gal':
        coord = 'G'
    elif coord_syst == 'equ':
        coord = 'C'
    else:
        coord = 'C'
    projview(sdf_map,
             xsize=1000,
             coord=["C",coord],
             graticule=True,
             graticule_labels=True,
             projection_type="mollweide",
             min=0,
             max=2,
             longitude_grid_spacing=30,
             latitude_grid_spacing=15,
             xlabel=r"$\alpha$ [deg]",
             ylabel=r"$\delta$ [deg]",
             xtick_label_color='black',
             ytick_label_color='black',
             cbar=True,
             cb_orientation='vertical',
             unit="E(B-V)",
             cmap='PuBu',
             fontsize={'title':fontsize_tick,
                       "xlabel":fontsize_label,
                       "ylabel":fontsize_label,
                       "title_label_pad":1,
                       "xtick_label": fontsize_tick,
                       "ytick_label": fontsize_tick,
                       "cbar_label":fontsize_tick,
                       "cbar_tick_label":fontsize_tick},
             override_plot_properties={"figure_width": figsize,
                                       "figure_size_ratio": 0.47},
             title=title_sky)
    
    # show the common skies
    if (net_stats[obs_name_test]['mask_common_sky']).sum()>0:
        
        mask_skies = net_stats[obs_name_test]['mask_common_sky']&\
            net_stats[obs_name_test]['mask_extragal']
        mask_gal_southern =  pdf_cand_visibility[mask_skies]["b"]<=-15
        mask_gal_northern =  pdf_cand_visibility[mask_skies]["b"]>15
        if (pdf_cand_visibility[mask_skies]["ra"]>180).sum()>0 and\
            (pdf_cand_visibility[mask_skies]["ra"]<=180).sum()>0:
            #first part of the sky map
            mask = pdf_cand_visibility[mask_skies]["ra"]>180
            mask_final_southern = mask&mask_skies&mask_gal_southern
            mask_final_northern = mask&mask_skies&mask_gal_northern
            if mask_final_southern.any():
                #plot Southern galactic sky
                ras = pdf_cand_visibility[mask_final_southern]["ra"].values
                decs = pdf_cand_visibility[mask_final_southern]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c_common_sky,
                            fmt='o',alpha = alpha_common_sky);
            if mask_final_northern.any():
                #plot Northern galactic sky
                ras = pdf_cand_visibility[mask_final_northern]["ra"].values
                decs = pdf_cand_visibility[mask_final_northern]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c_common_sky,
                            fmt='o',alpha = alpha_common_sky);
            # else:
            #     mask_final = mask&mask_skies
            #     ras = pdf_cand_visibility[mask_final]["ra"].values
            #     decs = pdf_cand_visibility[mask_final]["dec"].values
            #     theta, phi = equtorad_coord(ras,decs)
            #     newprojplot(theta=theta, phi=phi,color=c_common_sky,ls=None,
            #                 alpha=alpha_common_sky);
                
            #second part of the sky map
            mask = pdf_cand_visibility[mask_skies]["ra"]<=180
            mask_final_southern = mask&mask_skies&mask_gal_southern
            mask_final_northern = mask&mask_skies&mask_gal_northern
            if mask_final_southern.any():
                #plot Southern galactic sky
                
                ras = pdf_cand_visibility[mask_final_southern]["ra"].values
                decs = pdf_cand_visibility[mask_final_southern]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c_common_sky,
                            fmt='o',alpha = alpha_common_sky);
            if mask_final_northern.any():
                #plot Northern galactic sky
                ras = pdf_cand_visibility[mask_final_northern]["ra"].values
                decs = pdf_cand_visibility[mask_final_northern]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c_common_sky,
                            fmt='o',alpha = alpha_common_sky);
        else:
            # to smoothly plot the cropped radec due to the galactic latitude constraint
            if mask_gal_southern.any() and mask_gal_northern.any():
                #plot Southern galactic sky
                mask_final = mask_skies&mask_gal_southern
                ras = pdf_cand_visibility[mask_final]["ra"].values
                decs = pdf_cand_visibility[mask_final]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c_common_sky,
                            fmt='o',alpha = alpha_common_sky);
                #plot Northern galactic sky
                mask_final = mask_skies&mask_gal_northern
                ras = pdf_cand_visibility[mask_final]["ra"].values
                decs = pdf_cand_visibility[mask_final]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c_common_sky,
                            fmt='o',alpha = alpha_common_sky);
            else:
                mask_final = mask_skies
                ras = pdf_cand_visibility[mask_final]["ra"].values
                decs = pdf_cand_visibility[mask_final]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c_common_sky,
                            fmt='o',alpha = alpha_common_sky);
                
    if show_indiv_sky:
        obs_names = [obs_name_ref,obs_name_test]
        for obs in obs_names:
            c = get_color()
            if obs == obs_name_ref:
                mask_skies = pdf_cand_visibility_ref[obs]
            else:
                mask_skies = net_stats[obs]['mask_sky']
            if (pdf_cand_visibility[mask_skies]["ra"]>180).sum()>0 and\
                (pdf_cand_visibility[mask_skies]["ra"]<=180).sum()>0:
                #first part of the sky map
                mask = pdf_cand_visibility[mask_skies]["ra"]>180
                mask_final = mask&mask_skies
                ras = pdf_cand_visibility[mask_final]["ra"].values
                decs = pdf_cand_visibility[mask_final]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c,
                            alpha=alpha_indiv_sky );
                #second part of the sky map
                mask = pdf_cand_visibility[mask_skies]["ra"]<=180
                mask_final = mask&mask_skies
                ras = pdf_cand_visibility[mask_final]["ra"].values
                decs = pdf_cand_visibility[mask_final]["dec"].values
                theta, phi = utils.qutorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c,
                            alpha=alpha_indiv_sky );
            elif (pdf_cand_visibility[mask_skies]["ra"]>180).sum()==0:
                mask_final = mask_skies
                ras = pdf_cand_visibility[mask_final]["ra"].values
                decs = pdf_cand_visibility[mask_final]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c,lw=2,
                            alpha=alpha_indiv_sky );
            else:
                
                mask_final = mask_skies
                ras = pdf_cand_visibility[mask_final]["ra"].values
                decs = pdf_cand_visibility[mask_final]["dec"].values
                theta, phi = utils.equtorad_coord(ras,decs)
                newprojplot(theta=theta, phi=phi,color=c,lw=2,
                            alpha=alpha_indiv_sky );
    else:
        print((net_stats[obs_name_test]['mask_common_sky']).sum())
        theta = []
        phi = []
    #show the moon
    if get_moon(date[0]).ra.degree>180:
        theta_moon = np.deg2rad((get_moon(date[0]).dec.degree* -1) + 90)
        phi_moon = np.deg2rad(get_moon(date[0]).ra.degree-360)
    else:
        theta_moon = np.deg2rad((get_moon(date[0]).dec.degree* -1) + 90)
        phi_moon = np.deg2rad(get_moon(date[0]).ra.degree)
    newprojplot(theta=theta_moon, phi=phi_moon,marker='o',
                color="darkblue");
    #show the Sun
    if get_sun(date[0]).ra.degree>180:
        theta_sun = np.deg2rad((get_sun(date[0]).dec.degree* -1) + 90)
        phi_sun = np.deg2rad(get_sun(date[0]).ra.degree-360)
    else:
        theta_sun = np.deg2rad((get_sun(date[0]).dec.degree* -1) + 90)
        phi_sun = np.deg2rad(get_sun(date[0]).ra.degree)
    newprojplot(theta=theta_sun, phi=phi_sun,marker='*',ms = 10,
                color="orange");
    if save:
        outdir = '/media/dt270490/Transcend/Workspace/SVOM/LSST-fink/FINK_ToO/'+\
            'network_sky_visibility/'+obs_name_test
        if os.path.isdir(outdir): 
            plt.savefig(outdir+'/night_overlapsky_'+str(date.jd)+'_vro.png',dpi=100)
        else:
            os.mkdir(outdir, mode = 0o777)
            plt.savefig(outdir+'/night_overlapsky_'+str(date.jd)+'_vro.png',dpi=100)
    
def make_overlap_histo(sky_frac_visibilities:pd.DataFrame,
                       obs_test:str, save:bool):
    """
    Histogram of the fraction of the sky region visible by two 
    observatories
    Parameters
    ----------
    sky_frac_visibilities : pd.DataFrame
        sky_frac_visibilities : list
            Fraction of the sky visible at a given reference observatory 
            during the night at a test observatory for a given period of time.
    obs_test : str
        Observatory name.
    save : bool
        True or False.

    Returns
    -------
    None.

    """
   
    # plot settings
    fontsize_tick = 18
    fontsize_label = 22
    figsize = 14
    fig = plt.figure(figsize=(figsize,figsize))
    plt.hist(sky_frac_visibilities*100, density = True,
             bins=10, facecolor = '#2ab0ff', edgecolor='#169acf',
             linewidth=0.5)
    plt.xlabel('\% of the VRO night sky visible at '+obs_test,
               fontsize = fontsize_label)
    plt.ylabel('Probability ',
               fontsize = fontsize_label)
    plt.xticks(fontsize = fontsize_tick)
    plt.yticks(fontsize = fontsize_tick)
    plt.grid(True)
    if save:
        outdir = '/media/dt270490/Transcend/Workspace/SVOM/LSST-fink/FINK_ToO/'+\
            'network_sky_visibility/'+obs_test
        if os.path.isdir(outdir): 
            plt.savefig(outdir+'/hist_night_overlapsky_VRO_2023.png',dpi=100)
        else:
            os.mkdir(outdir, mode = 0o777)
            plt.savefig(outdir+'/hist_night_overlapsky_VRO_2023.png',dpi=100)
            
def plot_skyfrac(sky_frac_visibilities,dates, obs_names, date_ref, save):
    """
    

    Parameters
    ----------
    sky_frac_visibilities : list
        Fraction of the sky visible at a given observatory during the night at
        VRO for a given period of time

    Returns
    -------
    None.

    """
    # Compute the dates
    
    dates_day = [date-date_ref for date in dates]
    # plot settings
    fontsize_tick = 18
    fontsize_label = 22
    fontsize_legend = 18
    figsize = 14
    colors = get_cmap_gradient_color('plasma_r')
    index_col = 0
    fig = plt.figure(figsize=(figsize,figsize))
    
    for i in range(len(obs_names)-1):
        
        plt.plot(dates_day[i].jd,sky_frac_visibilities[i]*100,
                 color= colors[index_col],label=obs_names[i])
        index_col = index_col + 3
    plt.ylabel('\% of the visible VRO night sky',
               fontsize = fontsize_label)
    plt.xlabel('Day since '+date_ref.isot,
               fontsize = fontsize_label)
    plt.xticks(np.linspace(0,400,21),fontsize = fontsize_tick)
    plt.yticks(fontsize = fontsize_tick)
    plt.xlim([0,365])
    plt.grid(True)
    plt.legend(loc='best',bbox_to_anchor=(0.5, 0., 0.8, 1.0),
               fontsize=fontsize_legend)
    if save:
        outdir = '/media/dt270490/Transcend/Workspace/SVOM/LSST-fink/FINK_ToO/'+\
            'network_sky_visibility/'
        if os.path.isdir(outdir): 
            plt.savefig(outdir+'/overlapsky_VRO_2023.png',dpi=100)
        else:
            os.mkdir(outdir, mode = 0o777)
            plt.savefig(outdir+'/overlapsky_VRO_2023.png',dpi=100)

def plot_ztf_cand_hist(cand_stat:pd.DataFrame, save:bool):
    """
    Plot the distribution of the number of candidates processed by FINK during 
    the different selected month of the year

    Parameters
    ----------
    cand_stat : pd.DataFrame
        DataFrame listing the number of remaining candidates 
        passing the different selection cuts.
    save : bool
        True or False.

    Returns
    -------
    None.

    """
    
    # plot settings
    fontsize_tick = 18
    fontsize_label = 22
    fontsize_legend = 18
    fontsize_title = 26
    figsize = 14
    
    #few computation
    total_num_cand = cand_stat.iloc[0].sum()/1e6
    
    fig = plt.figure(figsize=(figsize,figsize))
    plt.bar(cand_stat.columns,cand_stat.iloc[0]/1e6, alpha=0.5)
    plt.ylabel('Number of candidates [in millions]',
               fontsize = fontsize_label)
    plt.xlabel('2022 Month',
               fontsize = fontsize_label)
    plt.xticks(fontsize = fontsize_tick)
    plt.yticks(fontsize = fontsize_tick)
    plt.grid(True)
    plt.title(str(total_num_cand)+' millions of candidates processed by FINK',
              fontsize = fontsize_title)
    if save:
        outdir = '/media/dt270490/Transcend/Workspace/SVOM/LSST-fink/FINK_ToO/'+\
            'network_sky_visibility/'
        if os.path.isdir(outdir): 
            plt.savefig(outdir+'/hist_ncand.png',dpi=100)
        else:
            os.mkdir(outdir, mode = 0o777)
            plt.savefig(outdir+'/hist_ncand.png',dpi=100)
            
def plot_filter_cand_hist(cand_stat:pd.DataFrame, save:bool):
    """
    Plot the distribution of the number of candidates processed by FINK during 
    the different selected month of the year

    Parameters
    ----------
    cand_stat : pd.DataFrame
        DataFrame listing the number of remaining candidates 
        passing the different selection cuts.
    save : bool
        True or False.

    Returns
    -------
    None.

    """
    
    # plot settings
    fontsize_tick = 18
    fontsize_label = 24
    fontsize_legend = 18
    fontsize_title = 26
    figsize = 14
    
    total_num_cand = cand_stat.iloc[0].sum()/1e6
    total_num_cand_selected = cand_stat.iloc[1].sum()
    fraction_cand_selected = cand_stat.iloc[1]/cand_stat.iloc[0]
    fig, ax = plt.subplots(figsize=(figsize,figsize))
    labels = ['Original total', 'After RB filter','After RB+Extragal filter', 'After RB+Extragal+Visibility filter (Obs.)',
              'After Obs.+Class filter','After Obs.+Class+Brightness filter','After Obs.+Class+Brightness+Det. history filter']
    for i in range(cand_stat.index.max()-4):
        if i < cand_stat.index.max()-5:
            alpha_bar = 0.2
        else:
            alpha_bar = 1.0
        ax.bar(cand_stat.columns,cand_stat.iloc[i], alpha=alpha_bar)
   
    # plt.plot(cand_stat.columns,fraction_cand_selected,':o',lw=2,markersize=10,
    #          c='black',label='_nolegend_')
    for i in range(len(fraction_cand_selected.values)):
        plt.text(ax.patches[i].get_width()+i-1.1,ax.patches[i].get_y()+5e6,
                  "{:.0f}".format(cand_stat.iloc[1][i]/1e3)+'k',
                  fontsize=fontsize_label,fontweight ='bold')
        # plt.text(ax.patches[i].get_width()+i-1.1,ax.patches[i].get_y()+1e7,
        #           "{:.4f}".format((1-fraction_cand_selected.values[i])*100),
        #           fontsize=fontsize_tick,fontweight ='bold')
    # plt.text(ax.patches[i].get_width()+i-12.5,ax.patches[i].get_y()+2e7,
    #           "Rejection fraction (\%)",
    #           fontsize=fontsize_label,fontweight ='bold')
    plt.ylabel('Number of remaining candidates',
               fontsize = fontsize_label)
    plt.xlabel('2022 Month',
               fontsize = fontsize_label)
    plt.xticks(fontsize = fontsize_tick)
    plt.yticks(ticks=np.logspace(0,8,9),fontsize = fontsize_tick)
    plt.gca().set_yticklabels(np.logspace(0,8,9))
    plt.ylim(1e0,1e7)
    plt.grid(True)
    plt.yscale('log')
    plt.legend(labels,loc='best',fontsize=fontsize_legend,
               bbox_to_anchor=(0.0, 0.5, 1.6, 0.5))
    plt.title("{:.2f}".format(total_num_cand_selected/1e6)+'/'+"{:.2f}".format(total_num_cand)+\
              ' millions of candidates \n(Rejection = '+"{:.6f}".\
                  format((1-(total_num_cand_selected*1e-6/total_num_cand))*100)+\
                  ' \%) passing the constraints',
              fontsize = fontsize_title)
    if save:
        outdir = '/media/dt270490/Transcend/Workspace/SVOM/LSST-fink/FINK_ToO/'+\
            'network_sky_visibility/'
        if os.path.isdir(outdir): 
            plt.savefig(outdir+'/hist_ncand_filter.png',dpi=100)
        else:
            os.mkdir(outdir, mode = 0o777)
            plt.savefig(outdir+'/hist_ncand_filter.png',dpi=100)

def plot_cand_sky(cand_stat:pd.DataFrame,
                  mask_type:str,
                  list_coords:list,
                  coord_syst:str,
                  save:bool):
    """
    Plot the candidates that passed the selection criteria into a custom 
    skymap

    Parameters
    ----------
    cand_stat : pd.DataFrame
        DataFrame listing the number of remaining candidates 
        passing the different selection cuts.
    mask_type : str
        Selection cut type.
    list_coords : list
        RA, dec of the remaining candidates.
    coord_syst : str
        Equatorial ("equ") or galactic ("gal").
    save : bool
        True or False.

    Returns
    -------
    None.

    """
    # Compute the rejection fraction
    total_num_cand = cand_stat.iloc[0].values[1:].sum()
    total_num_cand_selected = cand_stat.iloc[1].values[1:].sum()
    fraction_cand_selected = total_num_cand_selected/total_num_cand
    reject_frac = (1-fraction_cand_selected)*100
    # get the number of alert un objects
    n_alert = total_num_cand_selected
    n_object = len(list_coords[2])
    # plot settings
    c_common_sky = 'red'
    alpha_indiv_sky = 0.2
    fontsize_tick = 20
    fontsize_label = 24
    figsize = 19
    title_sky = ('Rejection fraction '+"{:.6f}".format(reject_frac)+'\% \n'+\
                 'Fink filters = '+mask_type+'\n '+\
                  "{:.2f}".format(n_alert/1e6)+' millions alerts related to '+\
                      "{:.2f}".format(n_object/1e6)+' millions objects')
                
    sdf_map = hp.read_map('/media/dt270490/Transcend/Workspace/GRB_data/'+\
                          'GRBase/catalogs/dust/EBV_SFD98_1_512.fits')
    if coord_syst =='gal':
        coord = 'G'
        x_axis_label = r"$l$ [deg]"
        y_axis_label = r"$b$ [deg]"
    elif coord_syst == 'equ':
        coord = 'C'
        x_axis_label = r"$\alpha$ [deg]"
        y_axis_label = r"$\alpha$ [deg]"
    else:
        coord = 'C'
        x_axis_label = r"$\alpha$ [deg]"
        y_axis_label = r"$\alpha$ [deg]"
    projview(sdf_map,
             xsize=1000,
             coord=["C",coord],
             graticule=True,
             graticule_labels=True,
             projection_type="mollweide",
             min=0,
             max=2,
             longitude_grid_spacing=30,
             latitude_grid_spacing=15,
             xlabel=x_axis_label,
             ylabel=y_axis_label,
             xtick_label_color='black',
             ytick_label_color='black',
             cbar=True,
             cb_orientation='vertical',
             unit="E(B-V)",
             cmap='PuBu',
             fontsize={'title':fontsize_tick,
                       "xlabel":fontsize_label,
                       "ylabel":fontsize_label,
                       "title_label_pad":1,
                       "xtick_label": fontsize_tick,
                       "ytick_label": fontsize_tick,
                       "cbar_label":fontsize_tick,
                       "cbar_tick_label":fontsize_tick},
             override_plot_properties={"figure_width": figsize,
                                       "figure_size_ratio": 0.47},
             title=title_sky)
    
    if coord_syst == "gal":
        ras = SkyCoord(list_coords[0]*u.degree,list_coords[1]*u.degree).galactic.l.value
        decs = SkyCoord(list_coords[0]*u.degree,list_coords[1]*u.degree).galactic.b.value
    else:
        ras = list_coords[0]
        decs = list_coords[1]
    # Plot the Northern sky only
    mask = np.array(ras) >180
    theta, phi = equtorad_coord(ras[mask],decs[mask])
    if len(theta)>0:
        newprojplot(theta=theta, phi=phi,color=c_common_sky,
                    fmt='o',alpha = 1.0);
    # Plot the Southern sky only
    mask = np.array(ras) <=180
    theta, phi = equtorad_coord(ras[mask],decs[mask])
    if len(theta)>0:
        newprojplot(theta=theta, phi=phi,color=c_common_sky,
                    fmt='o',alpha = 1.0);
    
    if save:
        outdir = '/media/dt270490/Transcend/Workspace/SVOM/LSST-fink/FINK_ToO/'+\
            'network_sky_visibility/'
        if os.path.isdir(outdir): 
            plt.savefig(outdir+'/rejection_map_vro.png',dpi=100)
        else:
            os.mkdir(outdir, mode = 0o777)
            plt.savefig(outdir+'/rejection_map_vro.png',dpi=100)