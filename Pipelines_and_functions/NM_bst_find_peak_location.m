function [Est1, Est2]=NM_bst_find_peak_location(result_lcmv, tlim)
% Just to make the main script much cleaner
    tcrop = foi_nrm  = nearest(gradfreq_h_base.freq, freq_nrm); result_lcmv.Time
    if ~isempty(result_lcmv.ImageGridAmp)
        stc_xyz = result_lcmv.ImageGridAmp.^2;
        cnt=0; clear stc
        for ii=1:3:length(stc_xyz)
            cnt=cnt+1;
            stc(cnt,:) = stc_xyz(ii,:) + stc_xyz(ii+1,:) + stc_xyz(ii+2,:);
        end
        [val, ind] = max(stc(:)); % Taking close to MNE (finding the max on stc series)
        [row_ind, ~] = ind2sub(size(stc),ind);
        Est_loc_lcmv = result_lcmv.GridLoc(row_ind, :)*1000;
        % Note: Neuromag_location = [-(CTF_loc_y), (CTF_loc_x), (CTF_loc_z)] *****************
        Est_loc_lcmv_neuromag = [-Est_loc_lcmv(2),Est_loc_lcmv(1),Est_loc_lcmv(3)]; 
        Est_val_lcmv = val;
        stimcat_act_loc_ctf = act_diploc(dipnum,:)*1000;
        difff_lcmv = sqrt(sum((stimcat_act_loc_ctf-Est_loc_lcmv).^2));
        disp([Est_loc_lcmv_neuromag, Est_val_lcmv, difff_lcmv])
        
    else
        % 1st way to find location    
        stc_xyz = result_lcmv.ImagingKernel*avgdata_pst;
        stc_xyz = stc_xyz.^2;
        cnt=0; clear stc
        for ii=1:3:length(stc_xyz)
            cnt=cnt+1;
            stc(cnt,:) = stc_xyz(ii,:) + stc_xyz(ii+1,:) + stc_xyz(ii+2,:);
            % disp([ii, cnt])
        end
        [val, ind] = max(stc(:)); % Taking close to MNE (finding the max on stc series)
        [row_ind, ~] = ind2sub(size(stc),ind);
        Est_loc_lcmv = result_lcmv.GridLoc(row_ind, :)*1000;
        % Note: Neuromag_location = [-(CTF_loc_y), (CTF_loc_x), (CTF_loc_z)] *****************
        Est_loc_lcmv_neuromag = [-Est_loc_lcmv(2),Est_loc_lcmv(1),Est_loc_lcmv(3)]; 
        Est_val_lcmv = val;
        stimcat_act_loc_ctf = act_diploc(dipnum,:)*1000;
        difff_lcmv = sqrt(sum((stimcat_act_loc_ctf-Est_loc_lcmv).^2));
        disp([Est_loc_lcmv_neuromag, Est_val_lcmv, difff_lcmv])

        % 2nd way to find location    
        grand_pow = sum(stc,2); %% similar to FieldTrip || SPM12 (finding the max of the summed series)
        [valx, indx] = sort(abs(grand_pow),'descend');
        n_hind = indx(1:50);
        n_act_grid2 = length(grand_pow(grand_pow > max(grand_pow(:))*0.50));
        PSVol2 = n_act_grid2*(par.gridres^3);
        % pow_ = sum(stc, 2);
        % [val_, ind_] = max(grand_pow);
        Est_loc_lcmv_ = result_lcmv.GridLoc(indx(1), :)*1000;
        % Note: Neuromag_location = [-(CTF_loc_y), (CTF_loc_x), (CTF_loc_z)] *****************
        Est_loc_lcmv_neuromag_ = [-Est_loc_lcmv_(2),Est_loc_lcmv_(1),Est_loc_lcmv_(3)]; 
        Est_val_lcmv_ = valx(1);
        % figure, plot(valx)
        difff_lcmv_ = sqrt(sum((stimcat_act_loc_ctf-Est_loc_lcmv_).^2));
        disp([Est_loc_lcmv_neuromag_, Est_val_lcmv_, difff_lcmv_])
        disp([difff_lcmv, difff_lcmv_])
    end
end