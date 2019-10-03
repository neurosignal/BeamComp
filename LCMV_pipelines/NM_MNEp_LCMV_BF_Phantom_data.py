#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 13:30:58 2019

@author: Alexandre Gramfort @ Inria, CEA/Neurospin, Universite Paris-Saclay, Paris, France
         Amit Jaiswal @ MEGIN Oy, Helsinki, Finland <amit.jaiswal@megin.fi>
                      @ School of Science, Aalto University, Espoo, Finland
         
Note: The code is used in the study:
      'Comparison of beamformer implementations for MEG source localizations' to analyze phantom data

"""

#% % Import modules
import mne
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from os.path import split, splitext
from nilearn.plotting import plot_stat_map
from nilearn.image import index_img
plt.rcParams.setdefault
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)
print(__doc__)

def my_var_cut_fn(epochs, plow, phigh, to_plot=True):
    """
    Variance base trial rejection function
    """
    trl_var, trlindx = np.empty((0,1), 'float'), np.arange(0,len(epochs))
    for trnum in range(len(epochs)):
        trl_var= np.vstack((trl_var, max(np.var(np.squeeze(epochs[trnum].get_data()), axis=1))))
    lim1 = (trl_var < np.percentile(trl_var, plow, interpolation='midpoint')).flatten()
    lim2 = (trl_var > np.percentile(trl_var, phigh, interpolation='midpoint')).flatten()
    outlr_idx = trlindx[lim1].tolist() + trlindx[lim2].tolist()
    if to_plot:
        plt.figure(), plt.scatter(trlindx, trl_var, marker='o', s=50, c='g', label='Good trials'), 
        plt.ylabel('Max. Variance accros channels')
        plt.scatter(outlr_idx, trl_var[outlr_idx], marker='o', s=50, c='r', label='Variance based bad trails'), 
        plt.xlabel('Trial number')
        plt.scatter(badtrls, trl_var[badtrls], marker='o', s=50, c='orange', label='Manually assigned bad trials')
        plt.ylim(min(trl_var)-min(trl_var)*0.01, max(trl_var)+max(trl_var)*0.01), plt.title(' Max. variance distribution') 
        plt.legend()
        plt.show()
    bad_trials = np.union1d(badtrls, outlr_idx)
    print('Removed trials: %s\n'%bad_trials)
    return bad_trials

#%% Set parameters, directories and file names
more_plots = False

data_dir     = 'BeamComp_DataRepo//MEG/Phantom//'
code_dir     = 'BeamComp_CodeRepo/LCMV_pipelines/'
subjects_dir = 'BeamComp_DataRepo//MRI/' 

for i_file in range(64):
    filenames = loadmat(code_dir + 'phantom_filenames.mat')
    filename  = filenames['phantom_filenames'][i_file][0][0]
    
    subject      = 'TriuxPhantom'
    mrifile      = subjects_dir + subject + '//mri/T1.mgz'
    surffile     = subjects_dir + subject + '//bem/watershed//' + subject + '_outer_skin_surface'
    transfile    = subjects_dir + subject + '//mri//transforms//' + subject +'-trans2.fif'
    act_dip      = loadmat(code_dir + '//triux_phantom_dipole_loc.mat') # TRIUX
    act_dip      = act_dip['biomag_phantom']
    
    #%% Read data > reject bad channels > pick only MEG channels
    fname  = data_dir + filename
    raw    = mne.io.read_raw_fif(fname, allow_maxshield=False, preload=True, verbose=True) 
    events = mne.find_events(raw, stim_channel='SYS201', min_duration=0.003)
    if i_file in range(48,64):
         events[:,1] -=3840
         events[:,2] -=3840
    
    dfname = split(splitext(raw.filenames[0])[0])[1]
    df_st  = dfname.split('_')
    dipnum = int(df_st[1][3:])
    
    bads=[]
    if not raw.info['projs']==[] or not 'sss' in raw.filenames[0]:
        bads = ['MEG1133', 'MEG1323', 'MEG0613', 'MEG1032']
    badtrls = []        
    raw.drop_channels(bads)
    raw.pick_types(meg=True)
    badtrial, badreason = [0,1], 'reset_signal'
    
    #%% Filtering for band 2-40 Hz
    raw.filter(2, 40, picks=None, filter_length='auto', l_trans_bandwidth='auto', h_trans_bandwidth='auto',
               n_jobs=1, method='fir', iir_params=None, phase='zero', fir_window='hamming', 
               fir_design='firwin', skip_by_annotation=('edge', 'bad_acq_skip'), 
               pad='reflect_limited', verbose=True)
    if more_plots:
        raw.plot()
        raw.plot_psd(average=False, spatial_colors=True, line_alpha=0.5, fmin=0.0, fmax=55.0)
    
    #%% Epoch the raw data 
    reject  = dict(mag=5000e-15, grad=300000e-15) 
    eventID = dipnum
    epochs  = mne.Epochs(raw, events, eventID, -0.5, 0.5, baseline=(None,0), picks=None, #reject=reject, 
                        preload=True, flat=None, proj=False, decim=1,reject_tmin=None, reject_tmax=None, 
                        detrend= None, on_missing='error', reject_by_annotation=True, verbose=True)
    epochs.drop(badtrial, badreason)
    #%% Find trial variance > index outliers> remove beyond plow and phigh percentile
    plow, phigh = 2.0, 98.0
    bad_trials=my_var_cut_fn(epochs, plow, phigh, to_plot=True)
    print('\n%d trial to remove from total %d trials...\nNo. of remaining trials = %d\n'%(len(bad_trials), 
                                                                                          len(epochs), 
                                                                                          len(epochs)-len(bad_trials)))
    epochs.drop(bad_trials, reason='eye_blink and high variance', verbose=True) 
    
    #%% Compute forward solution/leadfield
    bem=mne.make_sphere_model(r0=(0.0, 0.0, 0.0), head_radius=None, info=None, verbose=True)
    
    if not 'src_vol' in locals():
        src_vol = mne.setup_volume_source_space(subject=subject, pos=5.0, mri=mrifile, bem=None, surface=surffile,
                                             mindist=5.0, exclude=10.0, subjects_dir=subjects_dir, 
                                             volume_label=None, add_interpolator=True, verbose=True)
    if more_plots:
        mne.viz.plot_bem(subject=subject, subjects_dir=subjects_dir, orientation='coronal', slices=range(73,193,5), 
                     brain_surfaces=None, src=src_vol, show=True)
        mne.viz.plot_alignment(epochs.info, trans=transfile, subject=subject, subjects_dir=subjects_dir, fig=None,
                           surfaces=['head-dense', 'inner_skull'], coord_frame='head', show_axes=True,
                           meg=False, eeg='original', dig=True, ecog=True, bem=None, seeg=True,
                           src=src_vol, mri_fiducials=False,  verbose=True) 
    
    if not 'fwd' in locals():
        fwd = mne.make_forward_solution(epochs.info, trans=transfile, src=src_vol, bem=bem, 
                                        meg=True, eeg=False, mindist=2.5, n_jobs=1)
    if 'fwd' in locals() and not len(epochs.ch_names)==fwd['nchan']:
        fwd = mne.make_forward_solution(epochs.info, trans=transfile, src=src_vol, bem=bem, 
                                        meg=True, eeg=False, mindist=2.5, n_jobs=1)
    print("Leadfield size : %d sensors x %d dipoles" % fwd['sol']['data'].shape)
    
    #%% Average post-stim data and find input SNR
    evoked = epochs.average()
    evoked_pst = evoked.copy().crop(tmin=0.050, tmax=0.500)
    if more_plots:
        evoked.comment=dfname
        evoked.plot(spatial_colors=True, gfp=True, time_unit='ms')
        evoked_pst.plot(spatial_colors=True, gfp=True, time_unit='ms')
        
    noise_cov = mne.compute_covariance(epochs, tmin=-0.500, tmax=-0.050, method='empirical') 
    data_cov  = mne.compute_covariance(epochs, tmin=0.050, tmax=0.500,   method='empirical')  
        
    cov_rank = None if raw.info['proc_history']==[] else int(raw.info['proc_history'][0]['max_info']['sss_info']['nfree'])
    inverse_operator = mne.minimum_norm.make_inverse_operator(epochs.info, fwd, noise_cov, 
                                                            rank=cov_rank, loose=1.0, depth=0.199)
    snr_mnep, kk = mne.minimum_norm.estimate_snr(evoked_pst, inverse_operator, verbose=True)
    
    peak_ch, peak_time = evoked_pst.get_peak(ch_type='mag')
    tp = int((peak_time - 0.050)*evoked_pst.info['sfreq'])
    snr = snr_mnep[tp]
    snr_10log10 = 10*np.log10(snr)
    
    #%% Compute spatial filter and apply on post-stim data
    filters = mne.beamformer.make_lcmv(evoked.info, fwd, data_cov, reg=0.05, noise_cov=noise_cov, 
                                       pick_ori='max-power', rank=cov_rank, weight_norm='nai', 
                                       reduce_rank=True, verbose=True)
    stc = mne.beamformer.apply_lcmv(evoked_pst, filters, max_ori_out='signed', verbose=True)
    
    #%% Detect peak over summed power over time:: Stable over time
    stc_pow_series = np.square(stc) 
    stc_power = stc_pow_series.sum()
    
    src_peak, _ = stc_power.get_peak()
    est_loc = fwd['src'][0]['rr'][src_peak]*1000 
    
    #%% Calculate localization Error and Point spread volume
    loc_err = np.sqrt(np.sum(np.square(act_dip[dipnum-1]-est_loc))); 
    stc_power_data = stc_power.copy().data
    n_act_grid = len(stc_power_data[stc_power_data > (stc_power_data.max()*0.50)])
    PSVol      = n_act_grid*(5.0**3)
                
    print('Act_Sourceloc for %s' %dfname + '= %s' % str(act_dip[dipnum-1])) 
    print('Est_SourceLoc for %s' %dfname + '= %s' % str(np.around(est_loc,1)))
    print('Loc_error for %s' %dfname + '= %.1f mm'  %loc_err)
    print('Point Spread Volume (PSV) for %s' %dfname + '= %.1f mm^3' %PSVol)
    
    #%% Plot the activation
    if more_plots:
        plt.figure()
        ts_show = -100
        plt.plot(1e3 * stc.times, stc.data[np.argsort(stc.data.max(axis=1))[ts_show:]].T)
        _, t_peak = stc_pow_series.get_peak()
        timepoint = int(t_peak//stc_pow_series.tstep - stc_pow_series.times[0]//stc_pow_series.tstep)
        img = stc.as_volume(fwd['src'], dest='mri', mri_resolution=False, format='nifti1')
        plot_stat_map(index_img(img, timepoint), mrifile, threshold=stc.data.max()*0.50)
        plt.suptitle('%s\n'%dfname + 
                     'PeakValue= %.3f, ' % stc.data.max() + 'Est_loc= [%.1f,%.1f,%.1f], ' %tuple(est_loc) +
                     'Loc_err= %.2f mm' % loc_err, fontsize=12, color='white')
    plt.close('all')
    ################################################END####################################################