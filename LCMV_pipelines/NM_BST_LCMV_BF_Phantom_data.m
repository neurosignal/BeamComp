%% Brainstorm LCMV beamforming pipeline for Phantom MEG data
%% Created on 02/10/2018
%% Usage: For Phantom data
%% author: John Mosher@ University of Texas Health Science Center, Huston, Texas, USA
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
bst_dir  = 'BeamComp_DataRepo/Toolboxes/brainstorm3_26092019/brainstorm3-master/';
mri_dir  = 'BeamComp_DataRepo/MRI/';
data_dir = 'BeamComp_DataRepo/MEG/';
code_dir = 'BeamComp_CodeRepo/LCMV_pipelines/';
addpath(bst_dir)
addpath(mri_dir)
addpath(data_dir)
addpath(code_dir)

%% Set parameters
par                 = [];
par.Baseline        = [-0.500 -0.050];
par.DataWindow      = [0.050 0.500];
par.sampling_rate   = 1000; 
par.stimchan        = 'SYS201';
par.bads            = 'MEG1133, MEG1323, MEG0613, MEG1032';
par.more_plots      = 0;
par.bpfreq          = [2 40];
par.gridres         = 5; % mm
par.cov_cut         = [2, 98];
par.badtrs          = [1,2];
par.moving_maxf     = '';

%% Start brainstorm
if ~brainstorm('status')
    brainstorm nogui
end

%% Run to analyze all data (make i_file selection for specific data)
for i_file = 1:64
    
    %% ===== CREATE PROTOCOL =====
    ProtocolName = 'ELEKTA_PHANTOM';

    % Delete existing protocol
    gui_brainstorm('DeleteProtocol', ProtocolName);
    gui_brainstorm('DeleteProtocol', ProtocolName);
    gui_brainstorm('DeleteProtocol', ProtocolName);

    % Create new protocol
    gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

    % Start a new report
    bst_report('Start');

    % Phantom anatomy
    SubjectName = ['Phantom_' date];
    DipoleFile = generate_phantom_elekta(SubjectName);

    % Load all filenames
    load([code_dir 'phantom_filenames.mat'])

    filename = phantom_filenames{i_file};

    data_path = [data_dir 'Phantom/'];

    RawFile       = [data_path filename];
    [~, dfname,~] = fileparts(RawFile);

    dipnum = str2double(char(regexp(dfname,'(?<=_Dip).*(?=_IASoff)','match')));

    Act_locc=load(file_fullpath(DipoleFile));
    act_diploc = [Act_locc.Dipole.Loc]';

    %% Process: Create link to raw file
    clear sFiles sFilesEpochs sFilesAvg sAvgLCMV
    sFiles = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
        'subjectname',    SubjectName, ...
        'datafile',       {RawFile, 'FIF'}, ...
        'channelreplace', 1, ...
        'channelalign',   1, ...
        'evtmode',        'value');
    
    %% Mark bad channels
    if isempty(strfind(dfname, 'sss'))
        sFiles = bst_process('CallProcess', 'process_channel_setbad', sFiles, [], ...
            'sensortypes', par.bads);
    end

    %% Get subject definition
    sSubject = bst_get('Subject', SubjectName);
    % Get MRI file and surface files
    MriFile    = sSubject.Anatomy(sSubject.iAnatomy).FileName;
    CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
    HeadFile   = sSubject.Surface(sSubject.iScalp).FileName;

    if par.more_plots
        % Display MRI
        hFigMri1 = view_mri(MriFile);
        hFigMri3 = view_mri_3d(MriFile, [], [], 'NewFigure');
        hFigMri2 = view_mri_slices(MriFile, 'x', 20); 
        pause(0.5);
        % Close figures
        close([hFigMri1 hFigMri2 hFigMri3]);
        % Display scalp and cortex
        hFigSurf = view_surface(HeadFile);
        hFigSurf = view_surface(CortexFile, [], [], hFigSurf);
        hFigMriSurf = view_mri(MriFile, CortexFile);

        % Figure configuration
        iTess = 2;
        panel_surface('SetShowSulci',     hFigSurf, iTess, 1);
        panel_surface('SetSurfaceColor',  hFigSurf, iTess, [1 0 0]);
        panel_surface('SetSurfaceSmooth', hFigSurf, iTess, 0.5, 0);
        panel_surface('SetSurfaceTransparency', hFigSurf, iTess, 0.8);
        figure_3d('SetStandardView', hFigSurf, 'left');
        pause(0.5);
        % Close figures
        close([hFigSurf hFigMriSurf]);
    end 
    
    %fig1 = view_headpoints(sFiles.ChannelFile, HeadFile); % check the MEG/MRI alignment
    % Process: Refine registration
    sFiles = bst_process('CallProcess', 'process_headpoints_refine', sFiles, []);
    %fig2 = view_headpoints(sFiles.ChannelFile, HeadFile); % check the MEG/MRI alignment
    %pause(2.0)
    %close([fig1 fig2])
 
    %% Process: DC offset correction: [All file]
    sFiles = bst_process('CallProcess', 'process_baseline', sFiles, [], ...
        'baseline',    [], ...
        'sensortypes', 'MEG', ...
        'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
        'read_all',    0);

    %% Process: Remove linear trend: All file
    sFiles = bst_process('CallProcess', 'process_detrend', sFiles, [], ...
        'timewindow',  [], ...
        'sensortypes', 'MEG', ...
        'read_all',    0);

    %% Process: Band-pass
    sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
        'sensortypes', 'MEG', ...
        'highpass',    par.bpfreq(1), ...
        'lowpass',     par.bpfreq(2), ...
        'attenuation', 'strict', ...  % 60dB
        'mirror',      0, ...
        'useold',      0, ...
        'read_all',    0);

    %% %% ===== READ EVENTS =====
    % Process: Read from channel
    sFiles = bst_process('CallProcess', 'process_evt_read', sFiles, [], ...
        'stimchan',  par.stimchan, ...
        'trackmode', 1, ...  % Value: detect the changes of channel value
        'zero',      0);
    % Process: Delete spurious other events unrelated to dipoles
    sFiles = bst_process('CallProcess', 'process_evt_delete', sFiles, [], ...
        'eventname', '256, 768, 1792, 3072, 3840, 3584, 4096, 6144, 7168, 7680, 7936, transient, transient_bandpass');
    % Process: Rename events to have a leading zero, for proper sorting
    sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
        'src',  num2str(dipnum), ...
        'dest', sprintf('%02d ',dipnum));
    % Delete the first event of the first category (there is always an artifact)
    LinkFile = file_fullpath(sFiles.FileName);
    LinkMat = load(LinkFile, 'F');
    if ~isempty(LinkMat.F.events) && ~isempty(LinkMat.F.events(1).times)
        LinkMat.F.events(1).times(1)   = [];
        try LinkMat.F.events(1).samples(1) = [];end
        LinkMat.F.events(1).epochs(1)  = [];
    end
    bst_save(LinkFile, LinkMat, 'v6', 1);

    eventlabel = sprintf('%02d ',dipnum);

    if ~isempty(strfind(dfname, 'movement'))
        eventlabel = sprintf('%02d ',dipnum+3840);
    end
    fprintf('\neventlabel = %s\n', eventlabel)
    
    %% Make Epochs 
    clear sFilesEpochs
    sFilesEpochs = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
        'subjectname', SubjectName, ...
        'condition',   '', ...
        'eventname',   eventlabel, ...
        'timewindow',  [], ...
        'epochtime',   [par.Baseline(1), par.DataWindow(2)], ...
        'createcond',  1, ...
        'ignoreshort', 1, ...
        'usectfcomp',  0, ...
        'usessp',      0, ...
        'freq',        par.sampling_rate, ...
        'baseline',    [par.Baseline]);

    %% Remove bad trials
    par.bad_trials = [];
    [selecttrials, par] = NM_bst_varcut(sFilesEpochs, par, dfname, par.more_plots);

    badtrials = fliplr(par.bad_trials);
    for ii=badtrials
        sFilesEpochs(ii)=[];
    end
    fprintf('Remaing number of trials = %d\n',length(sFilesEpochs))
    %%  Process: Average: By trial group (folder average)
    sFilesAvg = bst_process('CallProcess', 'process_average', sFilesEpochs, [], ...
        'avgtype',    5, ...  
        'avg_func',   1, ...  
        'weighted',   0, ...
        'keepevents', 0);
    
    %% Process: Compute Head Model
    sFilesAvg = bst_process('CallProcess', 'process_headmodel', sFilesAvg, [], ...
                    'Comment',     'phantom headmodel', ...
                    'sourcespace', 2, ...  % MRI volume
                    'meg',         2, ...  % single sphere
                    'volumegrid',  struct(...
                    'Method',        'isotropic', ...
                    'nLayers',       17, ...
                    'Reduction',     3, ...
                    'nVerticesInit', 4000, ...
                    'Resolution',    0.005, ...
                    'FileName',      []));

    [subj_datadir, ~,~ ] = fileparts(file_fullpath(sFilesAvg.FileName));
    HeadModelFile = [subj_datadir '//headmodel_vol_meg_sphere.mat'];
    
    if par.more_plots
        view_spheres(HeadModelFile, sFilesAvg.ChannelFile, sSubject)
        view_gridloc(HeadModelFile, 'V'); view(172,6), hold on
    end
    %% Process: Compute covariance (noise or data)
    bst_process('CallProcess', 'process_noisecov', sFilesEpochs, [], ...
        'baseline',       par.Baseline, ...
        'datatimewindow', par.DataWindow, ... 
        'sensortypes',    'MEG', ...
        'target',         1, ...  %noise
        'dcoffset',       1, ...  
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'replacefile',    1);  

    bst_process('CallProcess', 'process_noisecov', sFilesEpochs, [], ...
        'baseline',       par.Baseline, ...
        'datatimewindow', par.DataWindow, ... 
        'sensortypes',    'MEG', ...
        'target',         2, ...  %data
        'dcoffset',       1, ...  
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'replacefile',    1);  
    
    %% Extract averaged data
    [~, avgtime, avgdata, ~, ~] = NM_bst_extract_data(par, dfname, sFilesAvg);
    if par.more_plots
        figure(); ax1=subplot(2,2,[1 2]); ax2=subplot(2,2,3); ax3=subplot(2,2,4);
        plot(ax1, avgtime, avgdata), xlim(ax1,[avgtime(1), avgtime(end)])
        title(ax1, [dfname '- Avg.data (sensor)'])
    end
    
    %% Make LCMV
     sAvgLCMV = bst_process('CallProcess', 'process_inverse_2018', sFilesAvg, [], ...
    'output',  3, ...  % Full results: one per file
    'inverse', struct(...
         'Comment',        'PNAI: MEG ALL', ...
         'InverseMethod',  'lcmv', ...
         'InverseMeasure', 'nai', ...
         'SourceOrient',   {{'free'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       1, ...
         'WeightExp',      0.5, ...
         'WeightLimit',    10, ...
         'NoiseMethod',    'reg', ...
         'NoiseReg',       0.05, ...
         'SnrMethod',      'rms', ...
         'SnrRms',         0, ...
         'SnrFixed',       3, ...
         'ComputeKernel',  0, ...
         'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));

    %% Convert 3D moment in power and find peak 
    result_lcmv = load(file_fullpath(sAvgLCMV.FileName));
    
    tlim = [0.050 0.500];
    tcrop = intersect(find(result_lcmv.Time>tlim(1)), find(result_lcmv.Time<tlim(2)));
    stc_xyz = result_lcmv.ImageGridAmp(:,tcrop).^2;
    cnt=0; stc=[];
    for ii=1:3:length(stc_xyz)
        cnt=cnt+1;
        stc(cnt,:) = stc_xyz(ii,:) + stc_xyz(ii+1,:) + stc_xyz(ii+2,:);
    end
    
    grand_pow    = mean(stc,2); 
    [valx, indx] = sort(abs(grand_pow),'descend');
    n_act_grid   = length(grand_pow(grand_pow > max(grand_pow(:))*0.50));
    PSVol        = n_act_grid*(par.gridres^3);
    Est_loc      = result_lcmv.GridLoc(indx(1), :)*1000;
    
    Est_locNM    = [-Est_loc(2),Est_loc(1),Est_loc(3)]; % NM_loc= [-(CTF_loc_y),(CTF_loc_x),(CTF_loc_z)]
    Est_val      = valx(1);

    difff = sqrt(sum((act_diploc(dipnum,:)*1000-Est_loc).^2));
    disp(fix([Est_locNM, Est_val, difff, n_act_grid, PSVol]))    
    
    %% Plot output
    if par.more_plots
        plot(ax2, valx,'-bo', 'MarkerEdgeColor','r','MarkerSize',10), title(ax2, 'Activation order')
        set(ax2, 'XScale', 'log'), xlabel(ax2, 'Source points in decaying activation order'), ylabel(ax2, 'Activation')
        plot(ax3, grand_pow,'--bo', 'MarkerEdgeColor','r','MarkerSize',5), title(ax3, 'Activation at each grid location')
        xlabel(ax3, 'Source points'), ylabel(ax3, 'Activation')
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15)
        set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'none')
    end
    %%
    clearex phantom_filenames DipoleFile SubjectName ProtocolName par ...
            bst_dir mri_dir data_dir code_dir i_file
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%