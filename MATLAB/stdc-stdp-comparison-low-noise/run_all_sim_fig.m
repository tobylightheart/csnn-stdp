%% run_all_sim_fig.m
close all;
clearvars;
clc;

save_all_figures_as_pdfs = 0;
save_all_figures_as_figs = 0;
if save_all_figures_as_pdfs || save_all_figures_as_figs
	if ~isdir('figures')
		mkdir('figures')
	end
end

disp('Presenting presynaptic and postsynaptic activity models.');
figure_presynaptic_neuron_activity;
figure_excitatory_postsynaptic_potential;
figure_postsynaptic_spike_potential;

disp('Running postsynaptic neuron simulations with STDP.');
run_lif_sim_batches_addSTDP;

disp('Running simulations with neuron construction estimating STDP.');
run_lif_sim_batches_addSTDC;

disp('Displaying results from STDP simulation and neuron construction.')
figure_postsynaptic_spike_times;
figure_weight_change_examples;
figure_weight_histograms;