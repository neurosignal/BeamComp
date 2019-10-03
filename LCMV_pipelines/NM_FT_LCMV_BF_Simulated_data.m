%% FieldTrip LCMV beamforming pipeline for simulated data
%% Created on Mon Jan  8 16:24:16 2019 
%% Usage: For Simulated data
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
mri_dir  = 'BeamComp_DataRepo/MRI/FT/';
data_dir = 'BeamComp_DataRepo/MEG/Simulated/';
code_dir = 'BeamComp_CodeRepo/LCMV_pipelines/';
addpath(ft_dir)
ft_defaults
addpath(mri_dir)
addpath(data_dir)
addpath(code_dir)

cd(code_dir)

%% Set data directory and other parameters
par            = [];
par.mri_align  = 'no';
par.mri_seg    = 'no';
par.more_plots = 0;
par.gridres    = 5;
par.bpfreq     = [2, 40];
par.bsfreq     = [49, 51];
par.cov_cut    = [2, 98];
par.stimchan   = 'STI101';

%% Select data and mri files
for i_file = 1:50
    load([code_dir 'simulated_filenames2.mat'])    
    filename = simulated_filenames{i_file};

    mrifname    = [mri_dir 'beamcomp-amit-170204-singlecomplete.fif'];
    segmrifname = [mri_dir 'beamcomp_1_01-brain.mat'];

    fname = [data_dir filename];
    [~, dfname, ~] = fileparts(fname);

    act_loc     = regexp(dfname,'(?<=nAm_at_).*(?=mm)','match');
    act_loc     = regexp(act_loc, '_', 'split');
    act_loc     = cellfun(@str2double,act_loc{1,1});

    load('channelslist.mat');

    if isempty(strfind(dfname, 'sss'))
        badch = {'MEG2233', 'MEG2212', 'MEG2332', 'MEG0111'}; 
    else
        badch       = {};
    end

    par.trialwin = [-0.2 0.2];
    par.ctrlwin  = [-0.1 -0.01];
    par.actiwin  = [0.01 0.1];
    par.badtrs   = []; 
    megchan=channelslist.meg;
    megchan(ismember(megchan, badch))=[]; 

    % Define triggers && label them
    keyset = {'Simulated'};  valueset = 5;
    evdict=containers.Map(keyset, valueset);

    %% Browse raw data
    if par.more_plots
        cfg                 = [];
        cfg.channel         = [megchan par.stimchan];
        cfg.viewmode        = 'vertical';
        cfg.blocksize       = 15;
        cfg.ylim            = [-1e-11 1e-11]; 
        cfg.preproc.demean  = 'yes';
        cfg.dataset         = fname;
        ft_databrowser(cfg);
    end
    %% Define trials
    cfgg                     = [];                  
    cfgg.dataset             = fname;               
    cfgg.channel             = megchan;
    cfgg.trialfun            = 'ft_trialfun_general';
    cfgg.trialdef.eventtype  = par.stimchan;        
    cfgg.trialdef.eventvalue = valueset;           
    cfgg.trialdef.prestim    = abs(par.trialwin(1));
    cfgg.trialdef.poststim   = par.trialwin(2);     
    cfgg = ft_definetrial(cfgg);          

    %% Visualize events triggers
    if par.more_plots
        ft_plot_events(figure, cfgg, keyset, valueset)
    end

    %% Preprocess continuous data
    cfg               = [];
    cfg.dataset       = fname;
    cfg.channel       = megchan;
    cfg.demean        = 'yes';
    cfg.bpfilter      = 'yes'; 
    cfg.bpfilttype    = 'but';
    cfg.bpfreq        = par.bpfreq;
    cfg.coilaccuracy  = 1;
    cfg.checkmaxfilter= 0;
    data = ft_preprocessing(cfg);

    %% Epoch data
    cfg     = [];
    cfg.channel = megchan;
    cfg.trl = cfgg.trl(:,:); 
    data = ft_redefinetrial(cfg, data);

    %% Interactive data browser 
    if par.more_plots
        cfg            = [];
        cfg.channel    = megchan;
        cfg.continuous = 'no';
        cfg.viewmode   = 'butterfly';
        ft_databrowser(cfg, data);
    end

    %% Find trial variance outliers and index them 
    [selecttrials, par] = NM_ft_varcut3(data, par, 1);

    %% Extract trial for stimulus category only > reject bad trials
    old_trs= size(data.trial,2);
    cfg = [];
    cfg.trials = selecttrials;
    data = ft_selectdata(cfg, data);
    fprintf('\nRemaining #trials = %d - %d = %d trials .........\nRemoved trials: ',...
            old_trs, length(par.bad_trials), size(data.trial,2)); disp(par.bad_trials)

    %% Compute covariance for common filter
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = par.trialwin;
    cfg.vartrllength = 2;
    evoked = ft_timelockanalysis(cfg,data); 

    %% Define noise covariance
    cfg = [];
    cfg.toilim = par.ctrlwin;
    datapre = ft_redefinetrial(cfg, data);

    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = 'all';
    evokedpre = ft_timelockanalysis(cfg,datapre);

    %% Define data covariance
    cfg = [];
    cfg.toilim = par.actiwin;
    datapst = ft_redefinetrial(cfg, data);

    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = 'all';
    evokedpst = ft_timelockanalysis(cfg,datapst);

    %% Plot Evoked data and cross-check everything
    if par.more_plots
        ft_mixed_plots(dfname, data, evoked, evokedpre, evokedpst)
    end 

    %% MEG-MRI coregistration, if needed
    if isequal(par.mri_align, 'yes')
        mri = ft_read_mri(mrifname);
        if isequal(par.align_meth, 'interactive') % Interactive
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

    %% Source analysis (LCMV)
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

    %% Calculate and save Neural Activity Index
    M1          = sourcepstM1;
    M1.avg.pow  = (sourcepstM1.avg.pow-sourcepreM1.avg.pow)./sourcepreM1.avg.pow;

    %% Find hotspot/peak location  
    M1.avg.pow(isnan(M1.avg.pow))=0;
    POW = abs(M1.avg.pow);
    [~,hind] = max(POW);
    hval = M1.avg.pow(hind);
    hspot = M1.pos(hind,:)*ft_scalingfactor(headmodel.unit,'mm');
    difff = sqrt(sum((act_loc-hspot).^2));
    n_act_grid = length(POW(POW > max(POW(:))*0.50));
    PSVol = n_act_grid*(par.gridres^3);

    %% Print the results
    fprintf('Results....\n')
    fprintf('%s\n', dfname)
    fprintf('Act. Location \t\t= [%.1f, %.1f, %.1f]mm\n', act_loc)
    fprintf('Est. Location \t\t= [%.1f, %.1f, %.1f]mm\n', hspot) 
    fprintf('Localization Error \t= %.1fmm\n', difff)
    fprintf('No. of active sources \t= %d \nPoint Spread Volume(PSV)= %dmm3\n',n_act_grid, PSVol)

    clearex par ft_dir data_dir mri_dir i_file code_dir leadfield    
% ######################## END #################################
end