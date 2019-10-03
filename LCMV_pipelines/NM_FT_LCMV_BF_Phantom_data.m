%% FieldTrip LCMV beamforming pipeline for phantom data
%% Created on Mon Jan  8 16:24:16 2018 
%% Usage: For both static and moving Phantom data
%% author: Robert Ostenveld @ Donders Institute of Brain Cognition and Behaviour, Nijmegen, The Netherlands
%%         Amit Jaiswal @ MEGIN Oy, Helsinki, Finland <amit.jaiswal@megin.fi>
%%                      @ School of Science, Aalto University, Espoo, Finland
%%         
%% Note: The code is used in the study:
%%      'Comparison of beamformer implementations for MEG source localizations'

%% Add fieldtrip in path
clear all; clc
restoredefaultpath
restoredefaultpath
ft_dir   = 'BeamComp_DataRepo/Toolboxes/fieldtrip-master_29092019/fieldtrip-master/';
mri_dir  = 'BeamComp_DataRepo/MRI/FT/';
data_dir = 'BeamComp_DataRepo/MEG/Phantom/';
code_dir = 'BeamComp_CodeRepo/LCMV_pipelines/';
addpath(ft_dir)
ft_defaults
addpath(mri_dir)
addpath(data_dir)
addpath(code_dir)

cd(code_dir)

%%  Analyze all             
for i_file = 1:64
    actual_diploc = load('triux_phantom_dipole_loc.mat');
    actual_diploc = actual_diploc.biomag_phantom;
    load([code_dir 'phantom_filenames.mat']);
    fname = [data_dir phantom_filenames{i_file}];

    par=[];
    par.more_plots  = 0;
    par.stimraise   = 0;
    par.gridres     = 0.005; % m
    par.bpfreq      = [2 40];
    par.cov_cut     = [2, 98];
  
    [~, dfname,~]=fileparts(fname);
    dip     = str2double(char(regexp(dfname,'(?<=_Dip).*(?=_IASoff)','match')));
    
    phantom_ct_file= [mri_dir '../fantomi_1362_01-jne-100423-5.nii'];
       
    if isempty(strfind(dfname, 'sss'))
        bads       = {'MEG0613', 'MEG1032', 'MEG1133', 'MEG1323'}; 
    else
        bads       = {};
    end
    
    badch            = strcat('-', bads);
    megchan          = ['MEG*', badch];
    stimchan         = 'SYS201';
    par.trig_min_gap = 0.5;
    par.trialwin     = [-0.5 0.5];
    par.ctrlwin      = [-0.5 -0.05];
    par.actiwin      = [0.05 0.5];
    if ~isempty(strfind(dfname, 'movement'))
        par.stimraise = 3840;
    else
    end
    
    %% Find trigger categories && label them
    keyset={['Phantom_dip-' strcat(num2str(dip))]};
    valueset=[dip];
    evdict=containers.Map(keyset, valueset);
    %% Browse raw data
    if par.more_plots
        cfg          = [];
        cfg.dataset  = fname;
        cfg.channel  = [megchan {stimchan}];
        cfg.viewmode = 'vertical';
        cfg.blocksize= 15;
        cfg.checkmaxfilter  = 0;
        cfg.ylim     = [-1e-11 1e-11];%'maxmin';
        cfg.preproc.demean   = 'yes';
        ft_databrowser(cfg);
    end       
    %% Define trials
    hdr = ft_read_header(fname, 'checkmaxfilter', 0);
    hdr.chantype{strcmp(stimchan, hdr.label)} = 'other trigger';
    hdr.chanunit{strcmp(stimchan, hdr.label)} = 'V';
    events = ft_read_event(fname, 'header', hdr);
    cfgg = [];
    cfgg.dataset             = fname;
    cfgg.channel             = megchan;
    cfgg.checkmaxfilter      = 0;
    cfgg.trialdef.eventvalue = dip + par.stimraise;                     
    cfgg.trialdef.prestim    = abs(par.trialwin(1));
    cfgg.trialdef.poststim   = par.trialwin(2);
    cfgg.trialdef.eventtype  = stimchan;
    cfgg.trialdef.eventvalue = dip + par.stimraise; 
    cfgg.event               = events;
    cfgg = ft_definetrial(cfgg);
    
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

    %% Segment the continuous preprocessed dat
    cfg2     = []; 
    cfg2.trl = cfgg.trl(cfgg.trl(3:end,4)==dip+par.stimraise,:);% leave initial two trials
    data = ft_redefinetrial(cfg2, data); 
            
    %% Interactive data browser 
    if par.more_plots
        cfg.viewmode  = 'butterfly';
        ft_databrowser(cfg, data);
    end
    %% Find trial variance outliers and index them & remove bad trials
    par.badtrs=[];
    [selecttrials, par] = NM_ft_varcut3(data, par, 1);

    alltrials = size(data.trial,2);
    cfg = [];
    cfg.trials   = selecttrials;
    data = ft_selectdata(cfg, data);
    fprintf('\nRemaining #trials = %d - %d = %d trials .........\nRemoved trials: ',...
            alltrials, length(par.bad_trials), size(data.trial,2)); disp(par.bad_trials)

    %% Computation of covariance matrices
    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = 'all'; % confirm this 
    cfg.vartrllength = 2;
    timelock = ft_timelockanalysis(cfg,data);
    
    % Noise
    cfg = [];            
    cfg.toilim = par.ctrlwin;
    datapre = ft_redefinetrial(cfg, data);
    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = 'all';
    cfg.vartrllength = 2;
    avgpre = ft_timelockanalysis(cfg,datapre);
    
    % Data
    cfg = [];
    cfg.toilim = par.actiwin;
    datapst = ft_redefinetrial(cfg, data);
    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = 'all';
    cfg.vartrllength = 2;
    avgpst = ft_timelockanalysis(cfg,datapst);
    
    if par.more_plots
        ft_mixed_plots(dfname, data, timelock, avgpre, avgpst)
    end

    %clear datapre datapst data hdr

    %% Define headmodel
    headmodel=[];
    headmodel.o       = [0,0,0];   
    headmodel.r       = .070;      
    headmodel= ft_convert_units(headmodel, 'm');

    %%  Prepare forward model
    if ~exist('leadfield', 'var')
        cfg                 = [];
        cfg.grad            = ft_convert_units(timelock.grad, headmodel.unit);  
        cfg.headmodel       = headmodel; % volume conduction headmodel
        cfg.channel         = timelock.label;
        cfg.grid.xgrid      = -headmodel.r:par.gridres:headmodel.r;
        cfg.grid.ygrid      = -headmodel.r:par.gridres:headmodel.r;
        cfg.grid.zgrid      = -par.gridres:par.gridres:headmodel.r;
        cfg.normalize       = 'yes';    % to remove depth bias
        cfg.normalizeparam  = 0.5;
        cfg.backproject     = 'yes';
        cfg.senstype        = 'MEG';
        leadfield           = ft_prepare_leadfield(cfg, timelock);
    end
    if exist('leadfield', 'var')
        if size(leadfield.label,1)~=size(timelock.label,1)
            cfg                 = [];
            cfg.grad            = ft_convert_units(timelock.grad, headmodel.unit);  
            cfg.headmodel       = headmodel; % volume conduction headmodel
            cfg.channel         = timelock.label;
            cfg.grid.xgrid      = -headmodel.r:par.gridres:headmodel.r;
            cfg.grid.ygrid      = -headmodel.r:par.gridres:headmodel.r;
            cfg.grid.zgrid      = -par.gridres:par.gridres:headmodel.r;
            cfg.normalize       = 'yes';    % to remove depth bias
            cfg.normalizeparam  = 0.5;
            cfg.backproject     = 'yes';
            cfg.senstype        = 'MEG';
            leadfield           = ft_prepare_leadfield(cfg, timelock);
        else
            disp('Using previously computed leadfield >>>>>')
        end
    end
        
    %% Cross-check all alignments
    if par.more_plots
        figure
        ft_plot_headshape(ft_convert_units(ft_read_headshape(fname), 'm'))
        ft_plot_vol(headmodel, 'facecolor', 'skin'); alpha 0.3; camlight
        ft_plot_mesh(leadfield.pos(leadfield.inside,:), 'vertexcolor', 'skin', 'vertexsize', 20)
        ft_plot_sens(ft_convert_units(ft_read_sens(fname, 'senstype', 'meg'), 'm'), 'coilshape', 'circle', 'coilsize', 0.015, 'facecolor', [1 1 1])
        rotate3d 
    end

    %% Campute and apply LCMV 
    cfg                  = [];
    cfg.method           = 'lcmv';
    cfg.grid             = leadfield; 
    cfg.headmodel        = headmodel; 
    cfg.senstype         = 'MEG';
    cfg.lcmv.keepfilter  = 'yes';
    cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
    cfg.lcmv.reducerank  = 2;
    cfg.lcmv.normalize   = 'yes';
    cfg.lcmv.lambda      = '5%';
    source_avg           = ft_sourceanalysis(cfg, timelock); % create spatial filters

    cfg                  = [];
    cfg.method           = 'lcmv';
    cfg.senstype         = 'MEG';
    cfg.grid             = leadfield;
    cfg.grid.filter      = source_avg.avg.filter;
    cfg.headmodel        = headmodel;
    sourcepreM1 = ft_sourceanalysis(cfg, avgpre); % apply spatial filter on noise
    sourcepstM1 = ft_sourceanalysis(cfg, avgpst); % apply spatial filter on data

    %% Calculate and save NAI
    M1          = sourcepstM1;
    M1.avg.pow  = (sourcepstM1.avg.pow-sourcepreM1.avg.pow)./sourcepreM1.avg.pow;
        
    %% Find hotspot/peak location  
    M1.avg.pow(isnan(M1.avg.pow))=0;
    POW = abs(M1.avg.pow);
    [~,hind] = max(POW);
    hval = M1.avg.pow(hind);
    hspot = M1.pos(hind,:)*ft_scalingfactor(headmodel.unit,'mm');
    difff = sqrt(sum((actual_diploc(dip,:)-hspot).^2));
    n_act_grid = length(POW(POW > max(POW(:))*0.50));
    PSVol = n_act_grid*((par.gridres*ft_scalingfactor(headmodel.unit,'mm'))^3);
    fprintf('%s\t\t => [%.1f %.1f %.1f]\t\t %.1f \t%.1f \t%.1f \t%.1f \n',...
                        dfname, hspot, hval, difff, n_act_grid, PSVol)
    %% Print the results
    fprintf('Results....\n')
    fprintf('%s\n', dfname)
    fprintf('Act. Location \t\t= [%.1f, %.1f, %.1f]mm\n', actual_diploc(dip,:))
    fprintf('Est. Location \t\t= [%.1f, %.1f, %.1f]mm\n', hspot) 
    fprintf('Localization Error \t= %.1fmm\n', difff)
    fprintf('No. of active sources \t= %d \nPoint Spread Volume(PSV)= %dmm3\n',n_act_grid, PSVol)
       
    %% Plot result 
    if par.more_plots
        if ~exist('segmri', 'var');  segmri= ft_convert_units(ft_read_mri(phantom_ct_file), 'm'); end
        cfg              = [];
        cfg.parameter    = 'pow'; %'avg.pow';
        cfg.downsample   = 3; 
        cfg.interpmethod = 'nearest';
        source_int  = ft_sourceinterpolate(cfg, M1, segmri);
        source_int = ft_convert_units(source_int, 'mm');

        source_int.mask = source_int.pow > max(source_int.pow(:))*0.50; % Set threshold for plotting
        cfg                 = [];
        cfg.method          = 'ortho';
        cfg.funparameter    = 'pow'; %'avg.pow';
        cfg.maskparameter   = 'mask';
        cfg.funcolormap     = 'hot';
        cfg.colorbar        = 'yes';
        cfg.location        = hspot;
        ft_sourceplot(cfg, source_int);
    end
    
    clearex par dip actual_diploc leadfield ft_dir mri_dir data_dir code_dir i_file
    close all
end 
