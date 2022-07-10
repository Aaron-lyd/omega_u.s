function figure_hf = scatter_plotting(data1, data2, data1_rms, data2_rms,txt_x, txt_y, title_text,OPTS_plot, OPTS_FIGS)
nx = OPTS_plot.nx;
ny = OPTS_plot.ny;
font_size = OPTS_plot.font_size;
OPTS_AXES = OPTS_plot.AXES;
n_figures = nx*ny;
nb = OPTS_plot.nbins;


for i = 1:n_figures
    
    ax = subaxis(nx, ny, i, OPTS_AXES{:});
    
    data1_2d = squeeze(data1(:,:,i));
    data2_2d = squeeze(data2(:,:,i));
    good = ~isnan(data1_2d) & ~isnan(data2_2d);
    xvec = data1_2d(good);
    yvec = data2_2d(good);
    hist2D(ax,log10(abs(xvec)), log10(abs(yvec)),'nbins',nb, 'RescaleFcn', @(x) sqrt(x))
    plot(log10(data1_rms(i)), log10(data2_rms(i)), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
    patch_x = [-12 log10(3e-9) log10(3e-9) -12];
    patch_y = [log10(1e-8) log10(1e-8) -20 -20];
    blue = [0.5843 0.8157 0.9882];
    patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
    grid on
    colorbar
    xlim([-12, -4])
    ylim([-20, -5])
    refline(1,0)
    xline(log10(3e-9),'-',{'Acceptable','Limit'});
    yline(log10(1e-8),'-',{'Acceptable','Limit'});
    xlabel(sprintf(txt_x ,data1_rms(i)) , 'fontsize',font_size,'Interpreter','latex');
    ylabel(sprintf(txt_y ,data2_rms(i)) , 'fontsize',font_size,'Interpreter','latex');
    title(title_text(i,:), 'fontsize',font_size,'Interpreter','latex')
    set(gca,'fontsize', font_size);
    
end
figure_hf = [];
end