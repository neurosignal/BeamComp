function NM_bst_view_dataplots(par, sFilesAvg, pausetime, to_plot)
    if to_plot
        hFigMeg = view_timeseries(sFilesAvg.FileName, 'MEG');
        hFigMegGrad = view_timeseries(sFilesAvg.FileName, 'MEG GRAD');
        hFigMegMag = view_timeseries(sFilesAvg.FileName, 'MEG MAG');
        % hFigEeg = view_timeseries(sFilesRaw.FileName, 'Misc');
        hFigStim = view_timeseries(sFilesAvg.FileName, 'Stim', {par.stimchan});
        % hFigSel = view_timeseries(sFilesRaw.FileName, 'Stim', {'MLT11','MLT12','MLT13'});
        % Figure configuration
        pause(0.5);
        %panel_record('SetTimeLength', 2);
        %panel_record('SetStartTime', 100);
        panel_record('SetDisplayMode', hFigMeg, 'column');
        panel_montage('SetCurrentMontage', hFigMeg, 'MEG');
        % Set filters: panel_filter('SetFilters', LowPassEnabled, LowPassValue, HighPassEnabled, HighPassValue, SinRemovalEnabled, SinRemovalValue, MirrorEnabled, FullSourcesEnabled)
        % panel_filter('SetFilters', 1, 45, 1, 1, 0, [], 0, 0);
        pause(0.5);
        panel_record('SetDisplayMode', hFigMeg, 'butterfly');
        panel_montage('SetCurrentMontage', hFigMeg, '');

        panel_record('SetDisplayMode', hFigMegGrad, 'butterfly');
        panel_record('SetDisplayMode', hFigMegMag, 'butterfly');



        fig1 = view_topography(sFilesAvg.FileName, 'MEG', '2DDisc');
        fig2 = view_topography(sFilesAvg.FileName, 'MEG', '2DSensorCap');%
        fig3 = view_topography(sFilesAvg.FileName, 'MEG', '3DSensorCap');
        fig4 = view_topography(sFilesAvg.FileName, 'MEG', '2DLayout'); % topoplot  

        pause(pausetime)

        % Close figures
        close([hFigMeg hFigMegGrad hFigMegMag hFigStim]);
        close([fig1 fig2 fig3 fig4]);
    end
    
end