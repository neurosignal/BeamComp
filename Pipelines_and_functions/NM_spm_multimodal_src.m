function grid = NM_spm_multimodal_src(grid, resolution, BF)
    fprintf('\nNote: Using multimodal data specific constrained source space....\n')
    clear vv vv_ ww ww_ yy yy_ zz
    if strfind(char(BF.data.D.conditions(1)), 'VEF')
        [vv,~]  = find(grid.pos(:,2)<0.0-resolution*5);
        zz = vv;
    elseif strfind(char(BF.data.D.conditions(1)), 'AEF')
        [vv,~]  = find(grid.pos(:,2)>0.0-resolution*5);
        [vv_,~] = find(grid.pos(:,2)<0.0+resolution*5);
        [ww,~] = find(grid.pos(:,3)>0.0+resolution*2);  
        zz = intersect(intersect(vv,vv_),ww);
    elseif strfind(char(BF.data.D.conditions(1)), 'SEF')
        [vv,~]  = find(grid.pos(:,2)>0.0-resolution*5);
        [vv_,~] = find(grid.pos(:,2)<0.0+resolution*5);
        [ww,~] = find(grid.pos(:,3)>0.0+resolution*6); 
        zz = intersect(intersect(vv,vv_),ww);
    end
    grid.inside = zz;
    grid.pos = grid.pos(grid.inside,:);
    grid.inside = [];
    grid.inside = zz;
    
    disp(['Using ' num2str(length(grid.inside)) ' points in source space. . .'])
    
    
    
