%% SPM12 LCMV beamforming pipeline for Human MEG data
%% Created on 02/12/2018
%% Usage: For Human EF MEG data
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
data_dir = 'BeamComp_DataRepo/';
code_dir = 'BeamComp_CodeRepo/LCMV_pipelines/';
addpath(daiss_dir)
addpath(spm_dir)
spm('defaults', 'eeg');
addpath(mri_dir)
addpath(data_dir)
addpath(code_dir)

work_dir = '/tmp/spm_lcmv/'; mkdir(work_dir)
cd(work_dir)
unix('find . -type f -name  "*spmeeg_*"  -delete')

%% Set data directory and other parameters
par                     = [];
par.stimchan            = 'STI 014';
par.cov_cut             = [2 98];
par.more_plots          = 0;
par.mri_seg             = 'no';
par.gridres             = 5;
par.mri_realign         = 'yes';
par.trial_win           = [-0.500 0.500];
par.bpfreq              = [2 95];

filename = 'multimodal_raw.fif'; % or 'multimodal_raw_tsss.fif'

%% Set events labels
keyset={'VEF-UR', 'VEF-LR', 'AEF-Re', 'VEF-LL', 'AEF-Le', 'VEF-UL', 'SEF-Lh', 'SEF-Rh'};
valueset=[1,2, 3, 4, 5, 8, 16, 32];
evdict=containers.Map(keyset, valueset);
%% Setting data, directories and other params
data_path       = [data_dir 'MEG/Human_EF/'];
mrifile         = [mri_dir 'beamcomp.nii,1'];
fname           = [data_path, filename];
[~, dfname,~]   = fileparts(fname);

actual_diploc = load([code_dir 'multimodal_biomag_Xfit_results.mat']);
actual_diploc = actual_diploc.multimodal_biomag_Xfit_diploc(:,4:6);
load([code_dir 'channelslist_old.mat']);

if isempty(strfind(dfname, 'sss'))
    badch       = {'MEG 0442'}; 
else
    badch       = {};
end

par.stimtype        = [par.stimchan '_up'];
par.trig_min_gap    = 0.5;
par.trial_win       = [-0.5 0.5];
par.actiwin_ph      = [0.0 0.5];
par.ctrlwin_ph      = [-0.5 0.0];
par.win_oi          = [-0.500 0.000; 0.000 0.500];
par.badtrs          = [];

megchan=channelslist.meg;
megchan(ismember(megchan, badch))=[]; 
if size(megchan, 1)>=size(megchan, 2)
    megchan=megchan';
end

%% Browse raw data
if par.more_plots
    cfg                 = [];
    cfg.channel         = [megchan par.stimchan];
    cfg.viewmode        = 'vertical';
    cfg.blocksize       = 15;
    cfg.preproc.demean  = 'yes';
    cfg.dataset         = fname;
    ft_databrowser(cfg);
end
%% Define trial definition matrix 
trls = spm_make_trialdef_realdata(fname, par.stimchan, valueset, par.trial_win, 1, 1);
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

if par.bpfreq(2)>46
    S       = [];
    S.D     = D;
    S.type  = 'butterworth';
    S.band  = 'stop';
    S.freq  = [49 51];
    S.dir   = 'twopass';
    S.prefix= 'f';
    D = spm_eeg_filter(S); 
end

%% Now run for all stim categories
for stimcat = keyset
    clear D
    D = spm_eeg_load(['ffbcspmeeg_' dfname '.mat']);
    
    stimcat=char(stimcat);
    stimval=evdict(stimcat);
    if strfind(stimcat, 'VEF')
        par.win_oi=[-0.200, -0.050; 0.050, 0.200];
    elseif strfind(stimcat, 'AEF')
        par.win_oi=[-0.150, -0.020; 0.020, 0.150];
    elseif strfind(stimcat, 'SEF')
        par.win_oi=[-0.100, -0.010; 0.010, 0.100];
    end

    woi     = par.win_oi*1000;     % convert into milisecond
    trialwin= par.trial_win*1000;  % convert into milisecond
    
    trials = trls(trls(:,4)==stimval,1:3); 
    dfname_stimcat = [dfname '_' stimcat];
    fprintf('Using: stimcat = %s\n', stimcat)
    disp(woi)

    %% Segment the filtered continuous data for this category
    S                 = [];
    S.D               = D;           
    S.bc              = 1;
    S.trl             = trials; 
    S.conditionlabels = stimcat;      
    S.prefix          = ['e' stimcat '_'];  
    D = spm_eeg_epochs(S);
    
    %% Browse converted raw data
    if par.more_plots
        cfg             = [];
        cfg.channel     = megchan;
        cfg.viewmode    = 'butterfly';
        ft_databrowser(cfg, D.ftraw);
    end
    
    %% Index bad trials
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
                size(trials,1), length(par.bad_trials), D.ntrials);   disp(par.bad_trials)
    end

    D = reload(D);
        
    %% Model head and coregister   
    clear matlabbatch     
    matlabbatch{1}.spm.meeg.source.headmodel.D = {fullfile(D)};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
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
    
    %% make another temp. dir for this stim
    out_tmp_dir = [work_dir dfname_stimcat];    mkdir(out_tmp_dir)
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
    spm_jobman('run', matlabbatch); toc
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
%     matlabbatch{1}.spm.tools.beamforming.features.plugin.cov.foi = par.bpfreq;
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
    % matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.contrast = 1;
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.contrast = [-1 1];
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.logpower = true;
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage'; %  bycondition
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.scale = 0; %1 %2
    matlabbatch{1}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG';
    spm_jobman('run', matlabbatch); toc
    
    if par.more_plots>10
        %% Write ouput map
        clear matlabbatch; tic
        matlabbatch{1}.spm.tools.beamforming.write.BF= {[out_tmp_dir '/BF.mat']};
        matlabbatch{1}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
        matlabbatch{1}.spm.tools.beamforming.write.plugin.nifti.space = 'native';%'mni'; 
        spm_jobman('run', matlabbatch); toc
    
        %% Plot output
        clear matlabbatch; tic
        matlabbatch{1}.spm.util.disp.data = {mrifile};
        %matlabbatch{1}.spm.util.disp.data = {[out_tmp_dir '/' 'uv_pow_re' stimcat '_ffbcspmeeg_multimodal_raw_tsss.nii']};
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
    difff = sqrt(sum((actual_diploc((find(valueset==evdict(stimcat))),:)-hspot).^2));
    n_act_grid = length(POW(POW > max(POW(:))*0.50));
    PSVol = n_act_grid*(par.gridres^3);

    %% Print the results
    fprintf('Results....\n')
    fprintf('Act. Location \t\t= [%.1f, %.1f, %.1f]mm\n', actual_diploc((find(valueset==evdict(stimcat))),:))
    fprintf('Est. Location \t\t= [%.1f, %.1f, %.1f]mm\n', hspot) 
    fprintf('Localization Error \t= %.1fmm\n', difff)
    fprintf('No. of active sources \t= %d \nPoint Spread Volume(PSV)= %dmm3\n',n_act_grid, PSVol)
    
    %% make some space
    close all
    clearex par out_tmp_dir trialwin woi keyset valueset evdict channelslist fname ...
            actual_diploc mrifile dfname badch trls megchan spm_dir mri_dir data_dir ...
            code_dir filename data_path stimcat work_dir
    cd(work_dir)
    unix(['rm -rf ' out_tmp_dir])
end                               
toc;
%###########################END##################################################  
        
