%% figure_excitatory_postsynaptic_potential.m
% clearvars;
% close all;
% clc;

if isempty(who('save_all_figures*'))
	save_all_figures_as_pdfs = 0;
	save_all_figures_as_figs = 0;
end

tau_m = 10;
tau_s = 2.5;
K = (4^(4/3))/3;
t_j = 10;

t = 0 : 0.01 : tau_m*5;
t_j_ind = find(t==t_j);

u = zeros(length(t),1);


for k = t_j_ind:length(t)
   u(k) = K * ( exp(-(t(k)-t_j)/tau_m) - exp(-(t(k)-t_j)/tau_s));
end

ymin = -0.2;
ymax = 2.2;

fig_w = 12;
fig_h = 5;

fig_epsp = figure;
set(fig_epsp, 'Color', 'w');
set(fig_epsp, 'Units', 'centimeters');
set(fig_epsp, 'PaperUnits', 'centimeters');
set(fig_epsp, 'Position', [2, 2, fig_w, fig_h]);
set(fig_epsp, 'PaperPosition', [0, 0, fig_w, fig_h]);
set(fig_epsp, 'PaperSize', [fig_w, fig_h]);
set(fig_epsp, 'PaperPositionMode', 'manual');

plot(t, u, 'k', 'LineWidth', 1)
line(t_j*[1,1], [ymin,ymax], 'LineStyle', ':', 'Color', 'k');

ylim([ymin, ymax]);

xlabel('Time [ms]');
ylabel('Postsynaptic Potential');

%% save figure and PDF
fig_file_name = fullfile('figures','fig_ps_epsp');

if save_all_figures_as_pdfs
	FT_PDF = '.pdf';
	saveas(fig_epsp,[fig_file_name, FT_PDF]);
	disp(['Figure saved as: ', fig_file_name, FT_PDF]);
end
if save_all_figures_as_figs
	FT_FIG = '.fig';
	saveas(fig_epsp,[fig_file_name, FT_FIG]);
	disp(['Figure saved as: ', fig_file_name, FT_FIG]);
end
