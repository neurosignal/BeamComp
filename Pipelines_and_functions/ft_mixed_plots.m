function ft_mixed_plots(dfname, data, evoked, evokedpre, evokedpst)
    figure(), FS = 11;
    subplot_tight(4,4,[1,2],0.05);  
    plot(evoked.time, evoked.avg); xlim([evoked.time(1) evoked.time(end)])
    title(dfname, 'FontSize', FS, 'Interpreter', 'none')
    subplot_tight(4,4,[5,6],0.05); clear trl_var;
    for trl = 1:size(data.trial,2), trl_var(trl,:) = max(var(data.trial{1,trl}(:,:)'));  end
    scatter(1:size(data.trial,2), trl_var, 50, 'go', 'filled'); xlim([1 size(data.trial,2)]), title('Max. variance for all selected trials', 'FontSize', FS)           
    subplot_tight(2,4,3,0.05);
    imagesc(evokedpre.cov), title(['Noise Cov [' num2str(min(min(evokedpre.cov))) ' ' num2str(max(max(evokedpre.cov))) ']'], 'FontSize', FS)
    colorbar('South')
    subplot_tight(2,4,4,0.05);
    imagesc(evokedpst.cov), title(['Data Cov [' num2str(min(min(evokedpst.cov))) ' ' num2str(max(max(evokedpst.cov))) ']'], 'FontSize', FS)
    colorbar('South')
    subplot_tight(2,4,5,0.05);   
    cfg = []; cfg.layout = 'neuromag306mag.lay';
    ft_multiplotER(cfg, evoked); 
    title('Magnetometer', 'FontSize', FS)
    subplot_tight(2,4,6,0.05);   
    cfg = []; cfg.layout = 'neuromag306planar.lay';
    ft_multiplotER(cfg, evoked); 
    title('Gradiometer', 'FontSize', FS)
    cfg=[];  cfg.method = 'mtmfft'; cfg.output = 'pow';
    cfg.foi  = 1:1:25; cfg.taper = 'hanning';
    tfr= ft_freqanalysis(cfg, evokedpst);
    subplot_tight(2,4,7,0.05);
    cfg = []; cfg.layout = 'neuromag306mag.lay'; 
    ft_topoplotTFR(cfg, tfr); title('Post-stim Mags', 'FontSize', FS)
    subplot_tight(2,4,8,0.05);
    cfg = []; cfg.layout = 'neuromag306planar.lay';
    ft_topoplotTFR(cfg, tfr); title('Post-stim Grads', 'FontSize', FS), clear tfr
    set(gcf, 'Position', get(0, 'Screensize'));
end