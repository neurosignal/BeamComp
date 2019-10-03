function [Chanlabel, avgtime1, avgdata1, avgtime2, avgdata2] = NM_bst_extract_data(par, dfname, sFilesAvg)
    % Extract evoked and post-stim-evoked(active) data

    Channels =load(file_fullpath(sFilesAvg.ChannelFile), 'Channel');
    Chanlabel = {Channels.Channel.Name};
    avgtime1  = load(file_fullpath(sFilesAvg.FileName), 'Time');
    avgtime1  = avgtime1.Time; avgtime2 = avgtime1(avgtime1>par.DataWindow(1));
    idxmeg = double(~cellfun(@isempty, strfind({Channels.Channel.Name}', 'MEG')));
    if isempty(strfind(dfname, 'sss')) 
        if length(regexp(par.badch(~isspace(par.badch)), ',', 'split'))==1
            badch = regexp(par.badch, ',', 'split');
        else
            badch = regexp(par.badch(~isspace(par.badch)), ',', 'split'); 
        end
    else
        badch={};
    end
    idxbadch  = -double(ismember({Channels.Channel.Name}', badch));
    megchindx = idxbadch + idxmeg;
    avgdata2 = load(file_fullpath(sFilesAvg.FileName), 'F');
    avgdata1 = avgdata2.F(logical(megchindx),:);
    avgdata2 = avgdata2.F(logical(megchindx),avgtime1>par.DataWindow(1));
end