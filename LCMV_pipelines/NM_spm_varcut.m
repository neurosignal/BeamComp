function [selecttrials, par] = NM_spm_varcut(D, par, to_plot)
    clear trl_var
    trlindx = 1:D.ntrials;
    timelk = D.fttimelock.trial;
    for trl = 1:D.ntrials
        trl_var(trl,:)= max(var(squeeze(timelk(trl,:,:))'));
    end
    percentiles = prctile(trl_var, par.cov_cut);
    outlr_idx = trl_var < percentiles(1) | trl_var > percentiles(2);
    bd_trl_var = trlindx(outlr_idx);
    bd_trls = union(bd_trl_var, []);
    par.bad_trials = [par.badtrs, bd_trls];

    if to_plot    
        figure
        subplot(111),scatter(1:D.ntrials, trl_var, 50, 'b*'); %xlim([1 D.ntrials]), 
        title('Max. variance'), hold on
        scatter(par.bad_trials, trl_var(par.bad_trials), 70, 'ro', 'linewidth',2); 
        xlim([1 D.ntrials]), title('Max. variance')
        legend({'all trials', 'bad trials'})
    end
    selecttrials = setdiff(trlindx, par.bad_trials);
end