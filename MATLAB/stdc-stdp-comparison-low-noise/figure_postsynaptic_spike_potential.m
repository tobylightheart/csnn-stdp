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
theta = 500;
K1 = 2;
K2 = 4;
t_i = 10;

t = 0 : 0.01 : tau_m*5;
t_i_ind = find(t==t_i);

u = zeros(length(t),1);

u(1:t_i_ind) = 0:(500/(t_i_ind-1)):500;

for k = t_i_ind:length(t)
   u(k) = theta * (K1* exp(-(t(k)-t_i)/tau_m) ...
	   - K2*(exp(-(t(k)-t_i)/tau_m) - exp(-(t(k)-t_i)/tau_s)));
end

ymin = -600;
ymax = 1200;

fig_w = 12;
fig_h = 6;

fig_spike = figure;
set(fig_spike, 'Color', 'w');
set(fig_spike, 'Units', 'centimeters');
set(fig_spike, 'PaperUnits', 'centimeters');
set(fig_spike, 'Position', [2, 2, fig_w, fig_h]);
set(fig_spike, 'PaperPosition', [0, 0, fig_w, fig_h]);
set(fig_spike, 'PaperSize', [fig_w, fig_h]);
set(fig_spike, 'PaperPositionMode', 'manual');

plot(t, u, 'k', 'LineWidth', 1)
line(t_i*[1,1], [ymin,ymax], 'LineStyle', ':', 'Color', 'k');
line([t(1),t(end)], theta*[1,1], 'LineStyle', ':', 'Color', 'k');
line([t(1),t(end)], [0,0], 'LineStyle', '--', 'Color', 'k');

ylim([ymin, ymax]);

xlabel('Time [ms]');
ylabel('Postsynaptic Potential');

%% save figure and PDF
fig_file_name = fullfile('figures','fig_ps_postsynaptic_spike');

if save_all_figures_as_pdfs
	FT_PDF = '.pdf';
	saveas(fig_spike,[fig_file_name, FT_PDF]);
	disp(['Figure saved as: ', fig_file_name, FT_PDF]);
end
if save_all_figures_as_figs
	FT_FIG = '.fig';
	saveas(fig_spike,[fig_file_name, FT_FIG]);
	disp(['Figure saved as: ', fig_file_name, FT_FIG]);
end
