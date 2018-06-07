%% figure_presynaptic_neuron_activty.m
% clc;
% close all;
% clearvars;

if isempty(who('save_all_figures*'))
	save_all_figures_as_pdfs = 0;
	save_all_figures_as_figs = 0;
end

% choose background rate
r_b = 10;

% initialise random number generator for repeatability
load('rng_states_spike_list');
rng(rng_states{1});

% generate spike list
gen_spike_list_lif_p_all;

% select presynaptic neuron activity to plot
plot_neurons = [(1:50)'; n_pattern+(1:50)'];

% find spikes in the first segment time
plot_spikes = spike_list(plot_neurons,:);
plot_spikes = plot_spikes(:, any(plot_spikes < segment_time));
plot_spikes(plot_spikes > segment_time) = Inf;

[neuron_s, spike_num_s] = find(plot_spikes < Inf);
plot_spike_list = [plot_spikes(plot_spikes<Inf), neuron_s];

num_spikes = length(spike_num_s);

% find spikes in the first pattern time
plot_pattern = spike_repeat_list(plot_neurons,:);
plot_pattern = plot_pattern(:, any(plot_pattern < segment_time));
plot_pattern(plot_pattern > segment_time) = Inf;

[neuron_p, spike_num_p] = find(plot_pattern < Inf);
plot_pattern_list = [plot_pattern(plot_pattern<Inf), neuron_p];

num_pattern = length(spike_num_p);

fig_w = 12;
fig_h = 8;

fig_ex_activity = figure;
set(fig_ex_activity, 'Color', 'w');
set(fig_ex_activity, 'Units', 'centimeters');
set(fig_ex_activity, 'PaperUnits', 'centimeters');
set(fig_ex_activity, 'Position', [2, 2, fig_w, fig_h]);
set(fig_ex_activity, 'PaperPosition', [0, 0, fig_w, fig_h]);
set(fig_ex_activity, 'PaperSize', [fig_w, fig_h]);
set(fig_ex_activity, 'PaperPositionMode', 'manual');

hold on;
box on;

plot(plot_spike_list(:,1), plot_spike_list(:,2), 'k.');

plot(plot_pattern_list(:,1), plot_pattern_list(:,2), ...
      'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4, ...
      'MarkerEdgeColor', 'k');

xlim([-0.002, 0.202]);
ylim([-1, 102]);

set(gca, 'YTick', [1, 25:25:100]);

xlabel('Time [s]');
ylabel('Presynaptic Neurons');

%% save figure and PDF
fig_file_name = fullfile('figures','fig_ps_example_presynaptic_activity_rb10');
if save_all_figures_as_pdfs
	FT_PDF = '.pdf';
	saveas(fig_ex_activity,[fig_file_name, FT_PDF]);
	disp(['Figure saved as: ', fig_file_name, FT_PDF]);
end
if save_all_figures_as_figs
	FT_FIG = '.fig';
	saveas(fig_ex_activity,[fig_file_name, FT_FIG]);
	disp(['Figure saved as: ', fig_file_name, FT_FIG]);
end
