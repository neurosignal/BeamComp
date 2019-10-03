function grid = NM_spm_phantom_src(grid, BF, resolution)
    dnum = str2double(char(regexp(char(BF.data.D.conditions(1)),'(?<=Phantom_dip-).*(?=)','match')));
    quard1 = 1:8;     quard2 = 9:16;
    quard3 = 17:24;   quard4 = 25:32;
    resolution2 = resolution*2;
    if any(dnum==quard1)
%         clear vv ww xx yy zz
%         vv = grid.pos(:,1)>=-resolution2;
%         ww = grid.pos(:,2)>=-resolution2;
%         xx = grid.pos(:,2)<=resolution2;
%         yy = grid.pos(:,3)>=-resolution2;
%         zz = and(and(and(vv,ww),xx),yy);
%         grid.pos = grid.pos(zz,:);
        clear vv ww xx yy zz
        [vv,~] = find(grid.pos(:,1)>=-resolution2);
        [ww,~] = find(grid.pos(:,2)>=-resolution2);
        [xx,~] = find(grid.pos(:,2)<=resolution2);
        [yy,~] = find(grid.pos(:,3)>=-resolution2);        
    elseif  any(dnum==quard2)
        clear vv ww xx yy zz
        [vv,~] = find(grid.pos(:,2)<=resolution2);
        [ww,~] = find(grid.pos(:,1)>=-resolution2);
        [xx,~] = find(grid.pos(:,1)<=resolution2);
        [yy,~] = find(grid.pos(:,3)>=-resolution2);
    elseif  any(dnum==quard3)
        clear vv ww xx yy zz
        [vv,~] = find(grid.pos(:,1)<=resolution2);
        [ww,~] = find(grid.pos(:,2)>=-resolution2);
        [xx,~] = find(grid.pos(:,2)<=resolution2);
        [yy,~] = find(grid.pos(:,3)>=-resolution2);
    elseif  any(dnum==quard4)
        clear vv ww xx yy zz
        [vv,~] = find(grid.pos(:,2)>=-resolution2);
        [ww,~] = find(grid.pos(:,1)>=-resolution2);
        [xx,~] = find(grid.pos(:,1)<=resolution2);
        [yy,~] = find(grid.pos(:,3)>=-resolution2);
    end
    zz=intersect(intersect(intersect(vv,ww),xx),yy);
    grid.inside = zz;
    grid.pos = grid.pos(grid.inside,:);
    grid.inside = [];
    grid.inside = zz;
    
    disp(['Using ' num2str(length(grid.inside)) ' points in source space. . .'])
