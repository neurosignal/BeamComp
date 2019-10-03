%% FieldTrip LCMV beamforming pipeline for Human MEG data
%% Created on Mon Jan  8 16:24:16 2018 
%% Usage: For Human MEG data
%% author: Robert Ostenveld @ Donders Institute of Brain Cognition and Behaviour, Nijmegen, The Netherlands
%%         Amit Jaiswal @ MEGIN Oy, Helsinki, Finland <amit.jaiswal@megin.fi>
%%                      @ School of Science, Aalto University, Espoo, Finland
%%         
%% Note: The code is used in the study:
%%      'Comparison of beamformer implementations for MEG source localizations'
%% 
%% Add fieldtrip in path
clear all; clc
restoredefaultpath
restoredefaultpath
ft_dir   = 'BeamComp_DataRepo/Toolboxes/fieldtrip-master_29092019/fieldtrip-master/';
mri_dir  = 'ChildBrain/BeamComp_DataRepo/MRI/FT/';
data_dir = 'BeamComp_DataRepo/';
code_dir = 'BeamComp_CodeRepo/LCMV_pipelines/';
addpath(ft_dir)
ft_defaults
addpath(mri_dir)
addpath(data_dir)
addpath(code_dir)

cd(code_dir)

%% Set data directory and other parameters
par            = [];
par.stimchan   = 'STI 014';
par.more_plots = 1;
par.mri_seg    = 'no';
par.mri_align  = 'no';
par.align_meth = 'headshape';
par.gridres    = 5.0; %mm
par.trialwin   = [-0.500 0.500];
par.bpfreq     = [2 95];
par.cov_cut    = [2, 98]; 

filename = 'multimodal_raw_tsss.fif'; % << change the file name here

actual_diploc = load([code_dir 'multimodal_biomag_Xfit_results.mat']);
actual_diploc = actual_diploc.multimodal_biomag_Xfit_diploc(:,4:6);
load([code_dir 'channelslist.mat']);

data_path   = [data_dir 'MEG/Human_EF/'];
mri_path    = [data_dir 'MRI/'];
mrifname    = [mri_path 'FT/beamcomp-amit-170204-singlecomplete.fif'];
segmrifname = [mri_path 'FT/beamcomp_1_01-brain.mat'];
fname       = [data_path filename];
[~,dfname,~]= fileparts(fname);

if isempty(strfind(dfname, 'sss'))
    badch   = {'MEG 0442'}; 
    megchan = ['MEG*',  strcat('-', badch)];
else
    badch   = {};
    megchan = 'MEG*';
end

%% Label event categories
keyset={'VEF-UR', 'VEF-LR', 'AEF-Re', 'VEF-LL', 'AEF-Le', 'VEF-UL', 'SEF-Lh', 'SEF-Rh'};
valueset=[1,2, 3, 4, 5, 8, 16, 32];
evdict=containers.Map(keyset, valueset);

%% Browse raw data
if par.more_plots
    cfg          = [];
    cfg.channel  = {megchan par.stimchan};
    cfg.viewmode = 'vertical';
    cfg.blocksize= 15;
    cfg.ylim     = [-1e-11 1e-11];%'maxmin';
    cfg.preproc.demean   = 'yes';
    cfg.dataset  = fname;
    ft_databrowser(cfg);
end
%% Define trials
cfg                     = [];       
cfg.dataset             = fname;   
cfg.channel             = megchan;
cfg.trialdef.eventtype  = par.stimchan;  
cfg.trialdef.eventvalue = valueset;      
cfg.trialdef.prestim    = 0.5;      
cfg.trialdef.poststim   = 0.5;     
cfgg = ft_definetrial(cfg);          

%% Visualize events triggers 
if par.more_plots, ft_plot_events(figure, cfgg, keyset, valueset), end

%% Preprocess continuous data
cfg                 = [];      
cfg.dataset         = fname;
cfg.channel         = megchan;
cfg.checkmaxfilter  = 0;
cfg.demean          = 'yes';
cfg.detrend         = 'yes';
cfg.bpfilter        = 'yes'; 
cfg.bpfilttype      = 'but';
cfg.bpfreq          = par.bpfreq; 
cfg.coilaccuracy    = 1;
data = ft_preprocessing(cfg);
    
cfg.bsfilter='yes';
cfg.bsfreq = [49 50];
data = ft_preprocessing(cfg, data);

%% Segment the data
cfg     = [];
cfg.trl = cfgg.trl(:,:); 
data = ft_redefinetrial(cfg, data);

%% Interactive data browser 
if par.more_plots
    cfg            = [];
    cfg.channel    = megchan;
    cfg.demean     = 'yes';
    cfg.continuous = 'no';
    cfg.viewmode   = 'butterfly';
    ft_databrowser(cfg, data);
end

%% MEG-MRI coregistration, if needed
if isequal(par.mri_align, 'yes')
    mri = ft_read_mri(mrifname);
    if isequal(par.align_meth, 'interective') % Interactive
        cfg             =[];
        cfg.method      ='interactive';
        cfg.coordsys    = 'neuromag';
        cfg.parameter   = 'anatomy';
        cfg.viewresult  =  'yes' ;
        [mri] = ft_volumerealign(cfg, mri);
    elseif isequal(par.align_meth, 'headshape') % Using ICP
        cfg                     = [];
        cfg.method              = 'headshape';
        cfg.headshape.headshape = ft_read_headshape(fname,'unit',mri.unit);
        cfg.headshape.icp       = 'yes';
        cfg.coordsys            = 'neuromag';
        cfg.parameter           = 'anatomy';
        cfg.viewresult          = 'yes';
        [mri] = ft_volumerealign(cfg, mri);
    end
end
%% Segment coregisted MRI, if needed
if isequal(par.mri_seg, 'yes')
    if ~isequal(par.mri_align, 'yes') % assuming that the given 'mrifname' is already coregistered
        mri = ft_read_mri(mrifname);
    end
    cfg          = [];  
    cfg.output   = 'brain';
    cfg.spmversion = 'spm12';
    segmri = ft_volumesegment(cfg, mri);
    segmri.transform = mri.transform;
    segmri.anatomy   = mri.anatomy;  
end

%% Reading the coregisterd and segmented mri file
if exist(segmrifname, 'file')==2 && ~isequal(par.mri_seg, 'yes')
    load(segmrifname); 
end

%% Plot segmented MRI volume
if par.more_plots
    cfg              = [];
    cfg.funparameter = 'brain';
    cfg.location     = [0,0,0];
    ft_sourceplot(cfg, segmri);
end
   
%% Compute the subject's headmodel/volume conductor model
if ~exist('headmodel', 'var')
    cfg                = [];
    cfg.method         = 'singleshell';
    % cfg.tissue         = 'brain';
    headmodel          = ft_prepare_headmodel(cfg, segmri);
end
if par.more_plots, figure, ft_plot_vol(headmodel, 'facecolor', 'brain'), rotate3d;  end

%% Create the subject specific grid (source space) > Compute forward model
if ~exist('leadfield', 'var') || length(data.label)~=length(leadfield.label)
    cfg                 = [];
    cfg.grad            = ft_convert_units(data.grad, headmodel.unit);
    cfg.headmodel       = headmodel;
    if isequal(headmodel.unit, 'mm')
        cfg.grid.resolution = par.gridres;
    elseif isequal(headmodel.unit, 'm')
        cfg.grid.resolution = par.gridres/1000;
    end
    cfg.grid.unit       = headmodel.unit;
    src_v               = ft_prepare_sourcemodel(cfg);
    
    cfg                 = [];
    cfg.grad            = ft_convert_units(data.grad, headmodel.unit);  
    cfg.headmodel       = headmodel;
    cfg.grid            = src_v;    
    cfg.channel         = data.label;
    cfg.normalize       = 'yes';    
    cfg.backproject     = 'yes';
    cfg.senstype        = 'MEG';
    leadfield           = ft_prepare_leadfield(cfg, data);
else
    disp('Using precomputed computed leadfield...')
end

%% Check whether everything is aligned and in same coordinate and unit
if par.more_plots
    ft_plot_alignment_check(fname, par, segmri, data, leadfield, headmodel)
end
clear src_v mri

%% Apply LCMV beamforming on each data category
for stimcat=keyset
    stimcat=char(stimcat);
    if strfind(stimcat, 'VEF')
        par.ctrlwin=[-0.200, -0.050];
        par.actiwin=[0.050, 0.200];
    elseif strfind(stimcat, 'AEF')
        par.ctrlwin=[-0.150, -0.020];
        par.actiwin=[0.020, 0.150];
    elseif strfind(stimcat, 'SEF')
        par.ctrlwin=[-0.100, -0.010];
        par.actiwin=[0.010, 0.100];
    end
    dfname_stim = [dfname '_' stimcat ];
    %% Extract trial for this stimulus category only && average 
    cfg = [];
    cfg.trials = find(data.trialinfo(:,1) ==  evdict(stimcat));
    data_stim=ft_preprocessing(cfg, data);
    
    %% Find trial variance outliers and index them 
    par.badtrs     = [];
    par.bad_trials = [];
    [selecttrials, par] = NM_ft_varcut3(data_stim, par, 1);

    %% Extract trial for stimulus category only > reject bad trials
    old_trs= size(data_stim.trial,2);
    cfg = [];
    cfg.trials = selecttrials;
    data_stim = ft_selectdata(cfg, data_stim);
    fprintf('\nRemaining #trials = %d - %d = %d trials .........\nRemoved trials: ',...
            old_trs, length(par.bad_trials), size(data_stim.trial,2)); disp(par.bad_trials)
    
    %% Calculate covariance for filter
    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = par.trialwin;
    cfg.vartrllength = 2;
    evoked = ft_timelockanalysis(cfg,data_stim);
        
    %% Compute noise cov.
    cfg = [];
    cfg.toilim = par.ctrlwin;
    datapre = ft_redefinetrial(cfg, data_stim);
    
    cfg = [];
    cfg.covariance='yes';
    evokedpre = ft_timelockanalysis(cfg,datapre);
    
    %% Compute data cov.
    cfg = [];
    cfg.toilim = par.actiwin;
    datapst = ft_redefinetrial(cfg, data_stim);
    
    cfg = [];
    cfg.covariance='yes';
    evokedpst = ft_timelockanalysis(cfg,datapst);

    %% Plot Evoked data and cross-check everything
    if par.more_plots>1
        ft_mixed_plots(dfname_stim, data_stim, evoked, evokedpre, evokedpst)
    end 
    
    %% Compute and apply LCMV
    cfg                  = [];
    cfg.method           = 'lcmv';
    cfg.grid             = leadfield;
    cfg.headmodel        = headmodel; 
    cfg.lcmv.keepfilter  = 'yes';
    cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
    cfg.lcmv.reducerank  = 2;
    cfg.lcmv.normalize   = 'yes';
    cfg.senstype         = 'MEG';
    cfg.lcmv.lambda      = '5%';
    source_avg           = ft_sourceanalysis(cfg, evoked); % create spatial filters

    cfg                  = [];
    cfg.method           = 'lcmv';
    cfg.senstype         = 'MEG';
    cfg.grid             = leadfield;
    cfg.grid.filter      = source_avg.avg.filter;
    cfg.headmodel        = headmodel;
    sourcepreM1=ft_sourceanalysis(cfg, evokedpre); % apply spatial filter on noise
    sourcepstM1=ft_sourceanalysis(cfg, evokedpst); % apply spatial filter on data

    %% Calculate Neural Activity Index 
    M1          = sourcepstM1;
    M1.avg.pow  = (sourcepstM1.avg.pow-sourcepreM1.avg.pow)./sourcepreM1.avg.pow;
    
    %% Find hotspot/peak location  
    M1.avg.pow(isnan(M1.avg.pow))=0;
    POW = abs(M1.avg.pow);
    [~,hind] = max(POW);
    hval = M1.avg.pow(hind);
    hspot = M1.pos(hind,:)*ft_scalingfactor(headmodel.unit,'mm');
    difff = sqrt(sum((actual_diploc((find(valueset==evdict(stimcat))),:)-hspot).^2));
    n_act_grid = length(POW(POW > max(POW(:))*0.50));
    PSVol = n_act_grid*(par.gridres^3);

    %% Print the results
    fprintf('Results....\n')
    fprintf('Act. Location \t\t= [%.1f, %.1f, %.1f]mm\n', actual_diploc((find(valueset==evdict(stimcat))),:))
    fprintf('Est. Location \t\t= [%.1f, %.1f, %.1f]mm\n', hspot) 
    fprintf('Localization Error \t= %.1fmm\n', difff)
    fprintf('No. of active sources \t= %d \nPoint Spread Volume(PSV)= %dmm3\n',n_act_grid, PSVol)
    
    clear sourcepreM1 sourcepstM1 POW hval hspot
    
end
% ##################################### END ################################