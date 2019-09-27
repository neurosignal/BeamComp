function grid = NM_spm_simulated_src(grid, resolution, act_loc)
    XX=act_loc(1)/1000; YY=act_loc(2)/1000; ZZ=act_loc(3)/1000;
    resolution2 = resolution*6;

    clear vv vv_ ww ww_ yy yy_ zz
    [vv,~]  = find(grid.pos(:,1)>=XX-resolution2);
    [vv_,~] = find(grid.pos(:,1)<=XX+resolution2);
    [ww,~]  = find(grid.pos(:,2)>=YY-resolution2);
    [ww_,~] = find(grid.pos(:,2)<=YY+resolution2);
    [yy,~]  = find(grid.pos(:,3)>=ZZ-resolution2);        
    [yy_,~] = find(grid.pos(:,3)<=ZZ+resolution2);  
    zz=intersect(intersect(intersect(intersect(intersect(vv,vv_),ww),ww_),yy),yy_);
    grid.inside = zz;
    grid.pos = grid.pos(grid.inside,:);
    grid.inside = [];
    grid.inside = zz;
    
    disp(['Taking grid only ' num2str(resolution2*1000) ' mm around(6D) the actual source. . .'])
    disp(['Using ' num2str(length(grid.inside)) ' points in source space. . .'])
