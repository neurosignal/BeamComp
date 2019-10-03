function [selecttrials, par] = NM_bst_varcut(sFilesEpochs, par, fname, to_plot)
%% Author: Amit Jaiswal @MEGIN, Helsinki.

    fprintf('\nDetecting the trials with excessive max/min variance...\n')
    
    Channels = load(file_fullpath(sFilesEpochs(1).ChannelFile),'Channel');
    idxmeg = double(~cellfun(@isempty, strfind({Channels.Channel.Name}', 'MEG')));
    if isempty(strfind(fname, 'sss')) 
        badch = par.badch; else badch={}; end
    idxbadch  = -double(ismember({Channels.Channel.Name}', badch));
    megchindx = logical(idxbadch + idxmeg);
    clear trl_var
    trlindx = 1:length(sFilesEpochs);
    for trl = trlindx
        data = load(file_fullpath(sFilesEpochs(trl).FileName), 'F');
        data = data.F(megchindx,:)';
        trl_var(trl,:) = max(var(data));
        clear data
    end
    clear Channels idxmeg idxbadch badch megchindx
    percentiles = prctile(trl_var, par.cov_cut);
    outlr_idx = trl_var < percentiles(1) | trl_var > percentiles(2);
    bd_trl_var = trlindx(outlr_idx);
    bd_trls = bd_trl_var;
    par.bad_trials = union([par.badtrs, bd_trls], 1);
    if to_plot
        try
            figure
            scatter(trlindx, trl_var, 100, 'b*'); 
            title('Max. variance'), hold on
            scatter(par.bad_trials, trl_var(par.bad_trials), 100, 'ro'); 
            xlim([trlindx(1)-1 trlindx(end)+1]), 
            ylim([min(trl_var)-min(trl_var)*0.1 max(trl_var)+max(trl_var)*0.1]),
            title('Max. variance', 'FontSize', 15)
            legend({'all trials','bad trials'}, 'FontSize', 15)
        catch err
            if exist('err', 'var')
                scatter([trlindx', trl_var], 'b*'); 
                title('Max. variance'), hold on
                scatter([par.bad_trials', trl_var(par.bad_trials)], 'ro'); 
                xlim([trlindx(1)-1 trlindx(end)+1]), 
                ylim([min(trl_var)-min(trl_var)*0.1 max(trl_var)+max(trl_var)*0.1]),
                title('Max. variance', 'FontSize', 15)  
                legend({'all trials','bad trials'}, 'FontSize', 15)
            end
        end
    end
     fprintf('Found total %d trials to remove : > >\n',length(par.bad_trials))
     disp(par.bad_trials)
     selecttrials = setdiff(trlindx, par.bad_trials);
     
end