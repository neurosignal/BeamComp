function [trl]= spm_make_trialdef_realdata(dataset, stimchan, stimval, trialwin, rmtrialstart, rmtrialend)
%% it's for real data
    cfg                     = [];                   
    cfg.dataset             = dataset;
    cfg.trialdef.eventtype  = stimchan;
    cfg.trialdef.eventvalue = stimval;                     
    cfg.trialdef.prestim    = abs(trialwin(1));                     
    cfg.trialdef.poststim   = trialwin(2); 
    cfg.minlength = cfg.trialdef.prestim + cfg.trialdef.poststim;
    cfg = ft_definetrial(cfg);
    % Leave the first and last trigger to avoid data deficiency error
    cfg.trl(1:rmtrialstart,:)  =[];
    cfg.trl(end-rmtrialend:end,:)=[];

    trl=cfg.trl(:,1:4);

    %fprintf(['**Subtracted first ' num2str(rmtrialstart) ' and last ' num2str(rmtrialend) ' trials to remove the reset artifacts. Resulting ' num2str(size(trl,1)) ' trials in total...........\n\n'])

end
