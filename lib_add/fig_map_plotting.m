function figure_hf = fig_map_plotting(data, data_rms,title_text,OPTS_plot, OPTS_FIGS)
nx = OPTS_plot.nx;
ny = OPTS_plot.ny;
X = OPTS_plot.X;
Y = OPTS_plot.Y;
land_mask = OPTS_plot.land_mask;
n_lim = OPTS_plot.n_lim;
p_lim = OPTS_plot.p_lim;
bwr = OPTS_plot.bwr;
font_size = OPTS_plot.font_size;
OPTS_AXES = OPTS_plot.AXES;
n_figures = nx*ny;
u_log10 = OPTS_plot.upb;
l_log10 = OPTS_plot.lowb;
logscale = OPTS_plot.log;


for i = 1:n_figures
    
    ax = subaxis(nx, ny, i, OPTS_AXES{:});
    hf = fig_map(ax, X, Y, squeeze(data(:,:,i)), land_mask, OPTS_FIGS);
    caxis([n_lim(i), p_lim(i)])
    if bwr(i)
        if logscale
            pcol_logsigned(ax, X, Y, squeeze(data(:,:,i)), l_log10, u_log10);
        else
            colormap(ax, bluewhitered), colorbar;
        end

    else
        colorbar(ax)
    end
    title(sprintf(title_text(i,:) ,data_rms(i)) , 'fontsize',font_size,'Interpreter','latex');
    set(gca,'fontsize', font_size);
    if mod(i,2) == 1 && i~=n_figures - 1
        ax.XTickLabel = [];
    elseif mod(i,2) == 0 && i~=n_figures
        ax.XTickLabel = [];
        ax.YTickLabel = [];
    elseif i==n_figures
        ax.YTickLabel = [];
    end
end
figure_hf = [];
end