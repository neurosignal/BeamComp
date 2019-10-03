%% SPM12 LCMV beamforming pipeline for Simulated MEG data
%% Created on 23/11/2018
%% Usage: For Simulated data
%% author: Vladimir Litvak@ The Welcome Centre of Human Neuroimaging, UCL, London, UK
%%         Amit Jaiswal @ MEGIN Oy, Helsinki, Finland <amit.jaiswal@megin.fi>
%%                      @ School of Science, Aalto University, Espoo, Finland
%%         
%% Note: The code is used in the study:
%%      'Comparison of beamformer implementations for MEG source localizations'
%% 

%% Add SPM12 in path 
clear all; clc  
restoredefaultpath 
restoredefaultpath 

spm_dir  = 'BeamComp_DataRepo/Toolboxes/SPM12_24092019/spm12/';
daiss_dir= 'BeamComp_DataRepo/Toolboxes/SPM12_24092019/DAiSS-master/';
mri_dir  = 'BeamComp_DataRepo/MRI/SPM/';
data_dir = 'BeamComp_DataRepo/MEG/Simulated/';
code_dir = 'BeamComp_CodeRepo/LCMV_pipelines/';
addpath(daiss_dir)
addpath(spm_dir)
spm('defaults', 'eeg');
addpath(mri_dir)
addpath(data_dir)
addpath(code_dir)

%% Set data directory and other parameters
for i_file=1:50
    
    work_dir = '/tmp/spm_lcmv/'; mkdir(work_dir); cd(work_dir)
    unix('find . -type f -name  "*spmeeg_*"  -delete')
    
    par             = [];
    par.stimchan    = 'STI101';
    par.cov_cut     = [2 98];
    par.more_plots  = 0;
    par.gridres     = 5;
    par.trial_win   = [-0.500 0.500];
    par.bpfreq      = [2 40];
    
    load([code_dir 'simulated_filenames2.mat']);
    filename =  simulated_filenames{i_file};

    %% Select data and mri files
    mrifile     = [mri_dir 'beamcomp.nii,1'];
    segmrifname = [mri_dir 'beamcomp_1_01-brain.mat'];

    fname = [data_dir filename];
    [~, dfname, ~] = fileparts(fname);

    act_loc     = regexp(dfname,'(?<=nAm_at_).*(?=mm)','match');
    act_loc     = regexp(act_loc, '_', 'split');
    act_loc     = cellfun(@str2double,act_loc{1,1});

    load([code_dir 'channelslist.mat']);

    if isempty(strfind(dfname, 'sss'))
        badch = {'MEG2233', 'MEG2212','MEG2332', 'MEG0111'};
    else
        badch     = {};
    end
    par.trialwin = [-0.200 0.200];
    par.ctrlwin  = [-0.100 -0.010];
    par.actiwin  = [0.010 0.100];
    par.win_oi   = [-0.100, -0.010; 0.010, 0.100];
    woi          = par.win_oi*1000;    % convert into milisecond
    trial_win    = par.trialwin*1000;  % convert into milisecond

    megchan=channelslist.meg;
    megchan(ismember(megchan, badch))=[]; 

    %% Label the events
    keyset = {'Simulated'};  valueset = 5;
    evdict = containers.Map(keyset, valueset);

    %% Browse raw data
    if par.more_plots
        cfg          = [];
        cfg.channel  = [megchan par.stimchan];
        cfg.viewmode = 'vertical';
        cfg.blocksize= 15;
        cfg.preproc.demean   = 'yes';
        cfg.dataset  = fname;
        ft_databrowser(cfg);
    end

    %% Define trials
    trls = spm_make_trialdef_realdata(fname, par.stimchan, valueset, par.trial_win, 0, 0);
    if par.more_plots
        cfg        = [];
        trlss      = trls;
        cfg.trl    = trlss;
        ft_plot_events(figure, cfg, keyset, valueset)
    end  

    %% Convert fif raw data to SPM format     
    S                   = [];
    S.dataset           = fname;
    S.mode              = 'continuous';
    S.channels          = megchan; 
    S.saveorigheader    = 1;
    S.inputformat       ='neuromag_fif';
    S.conditionlabels   = char(keyset);
    S.outfile           = ['spmeeg_' dfname];
    D = spm_eeg_convert(S);

    %% Baseline correction 
    S   = [];
    S.D = D;
    S.prefix ='bc';
    D = spm_eeg_bc(S);
    %% Filter the converted raw data (band pass and bandstop)
    S       = [];
    S.D     = D;
    S.type  = 'butterworth';
    S.band  = 'bandpass';
    S.freq  = par.bpfreq;
    S.dir   = 'twopass';
    S.prefix= 'f';
    D = spm_eeg_filter(S);

    %% Epoch the filtered continuous data
    S                 = [];
    S.D               = D;           
    S.bc              = 1;
    S.trl             = trls(:,1:3); 
    S.conditionlabels = keyset;      
    S.prefix          = 'e';         
    D = spm_eeg_epochs(S);

    %% Browse converted raw data
    if par.more_plots
        cfg             = [];
        cfg.channel     = megchan;
        cfg.viewmode    = 'butterfly';
        ft_databrowser(cfg, D.ftraw);
    end

    %% Find trial variance && z-score outliers and index them
    par.badtrs    = [];
    par.bad_trials= [];
    [selecttrials, par] = NM_spm_varcut(D, par, 1);

    %% Remove bad trials
    if par.bad_trials
        D = badtrials(D, par.bad_trials, 1);
        save(D);
        S=[];
        S.D=D;
        S.prefix='r';
        D=spm_eeg_remove_bad_trials(S);
        fprintf('\nRemaining #trials = %d - %d = %d trials .........\nRemoved trials : ',...
                size(trls,1), length(par.bad_trials), D.ntrials);   disp(par.bad_trials)
    end

    D = reload(D);

    %% Model head and coregister  
    clear matlabbatch
    matlabbatch{1}.spm.meeg.source.headmodel.D = {fullfile(D)};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = 'comments';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {mrifile};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run', matlabbatch);
    close all

    %% make another temp. dir for this stim
    out_tmp_dir = [work_dir dfname];    mkdir(out_tmp_dir)
    cd(out_tmp_dir)   
    %% Prepare data    
    clear matlabbatch; tic 
    matlabbatch{1}.spm.tools.beamforming.data.dir = {pwd};
    matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(D)};
    matlabbatch{1}.spm.tools.beamforming.data.val = 1;
    matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
    matlabbatch{1}.spm.tools.beamforming.data.space = 'Head'; 
    matlabbatch{1}.spm.tools.beamforming.data.overwrite = 1;
    spm_jobman('run', matlabbatch); toc

    %% Prepare source space and compute leadfield    
    clear matlabbatch; tic 
    matlabbatch{1}.spm.tools.beamforming.sources.BF = {[out_tmp_dir '/BF.mat']};
    matlabbatch{1}.spm.tools.beamforming.sources.reduce_rank = [2 3];
    matlabbatch{1}.spm.tools.beamforming.sources.keep3d = 1;
    matlabbatch{1}.spm.tools.beamforming.sources.plugin.grid.resolution = par.gridres;
    matlabbatch{1}.spm.tools.beamforming.sources.plugin.grid.space = 'Head';
    matlabbatch{1}.spm.tools.beamforming.sources.visualise = 1;
    spm_jobman('run', matlabbatch); tocspm
    
    close all

    %% Prepare covariance matrices
    clear matlabbatch; tic 
    matlabbatch{1}.spm.tools.beamforming.features.BF = {[out_tmp_dir '/BF.mat']};
    matlabbatch{1}.spm.tools.beamforming.features.whatconditions.all = 1;
    matlabbatch{1}.spm.tools.beamforming.features.woi = woi;
    matlabbatch{1}.spm.tools.beamforming.features.modality = {
                                                              'MEG'
                                                              'MEGPLANAR'
                                                              }';
    matlabbatch{1}.spm.tools.beamforming.features.fuse = 'meg';    
    % matlabbatch{1}.spm.tools.beamforming.features.plugin.cov.foi = par.bpfreq;
    matlabbatch{1}.spm.tools.beamforming.features.plugin.cov.taper = 'hanning';
    matlabbatch{1}.spm.tools.beamforming.features.regularisation.manual.lambda = 5;
    matlabbatch{1}.spm.tools.beamforming.features.bootstrap = false;
    spm_jobman('run', matlabbatch); toc

    %% For Inverse solution (spatial filter)
    clear matlabbatch; tic
    matlabbatch{1}.spm.tools.beamforming.inverse.BF = {[out_tmp_dir '/BF.mat']};
    matlabbatch{1}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
    matlabbatch{1}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
    spm_jobman('run', matlabbatch); toc
    close all

    %% Prepare output 
    clear matlabbatch; tic
    matlabbatch{1}.spm.tools.beamforming.output.BF = {[out_tmp_dir '/BF.mat']};
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.woi = woi;
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.foi = par.bpfreq;
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.contrast = [-1 1];
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.logpower = true;
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage'; %  bycondition
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.scale = 0; %1 %2
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG';
    spm_jobman('run', matlabbatch); toc

    %% Optional for visualization
    if par.more_plots>10
        % % Write ouput map
        clear matlabbatch; tic
        matlabbatch{1}.spm.tools.beamforming.write.BF= {[out_tmp_dir '/BF.mat']};
        matlabbatch{1}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
        matlabbatch{1}.spm.tools.beamforming.write.plugin.nifti.space = 'native';%'mni'; 
        spm_jobman('run', matlabbatch); toc

        % % Plot output
        clear matlabbatch; tic
        matlabbatch{1}.spm.util.disp.data = {mrifile};
        spm_jobman('run', matlabbatch); toc
    end

    %% Load the results
    movefile('BF.mat', 'BF_LCMV.mat')
    BF = load('BF_LCMV.mat');

    %% Find hotspot/peak location  
    BF.output.image.val(isnan(BF.output.image.val))=0;
    POW = abs(BF.output.image.val);
    [~,hind] = max(POW);
    hval = BF.output.image.val(hind);
    hspot = BF.sources.pos(hind, :)*1000;
    difff = sqrt(sum((act_loc-hspot).^2));
    n_act_grid = length(POW(POW > max(POW(:))*0.50));
    PSVol = n_act_grid*(par.gridres^3);

    %% Print the results
    fprintf('Results....\n')
    fprintf('Act. Location \t\t= [%.1f, %.1f, %.1f]mm\n', act_loc)
    fprintf('Est. Location \t\t= [%.1f, %.1f, %.1f]mm\n', hspot) 
    fprintf('Localization Error \t= %.1fmm\n', difff)
    fprintf('No. of active sources \t= %d \nPoint Spread Volume(PSV)= %dmm3\n',n_act_grid, PSVol)

    %% close and delete out_tmp_dir
    close all
    cd(work_dir)
    unix(['rm -rf ' out_tmp_dir])
    %###########################END##################################################      
end