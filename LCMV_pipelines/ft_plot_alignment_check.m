function ft_plot_alignment_check(fname, par, segmri, data, leadfield, headmodel)
        fprintf('\n\nCross check the alignment of isotrack points, headmodel and leadfield inside sensors array >>>>\n\n')
        headshape=ft_read_headshape(fname,'unit',segmri.unit);
        figure
        ft_plot_headshape(headshape, 'fidlabel','no', 'vertexsize', 20)
        text(headshape.fid.pos(:,1), headshape.fid.pos(:,2), headshape.fid.pos(:,3), headshape.fid.label(:), 'fontsize',10);
        text(headshape.pos(1:4,1), headshape.pos(1:4,2), headshape.pos(1:4,3), headshape.label(1:4), 'color', 'b', 'fontsize',10);
        clear headshape
        ft_plot_vol(headmodel, 'facecolor', 'skin'); alpha 0.3; camlight
        ft_plot_sens(ft_convert_units(ft_read_sens(fname, 'senstype', 'meg'), segmri.unit), 'coilshape', 'circle', 'coilsize', 0.015, 'facecolor', [1 1 1])
        ft_plot_mesh(leadfield.pos(leadfield.inside,:));
        grad=ft_convert_units(data.grad, 'm'); ft_plot_topo3d(grad.chanpos, ones(306,1)), alpha 0.2; 
        title('Aligned plot : Scanning grid | head model | isotrack points | Sensor array')
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15)
        
        if data.grad.unit==leadfield.unit && data.grad.unit==headmodel.unit
            disp('Everything in the same unit..')
        else
            error('Unit difference between data.grad, headmodel and leadfield...')
        end
        
end