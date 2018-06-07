%% n_post_activity_analysis.m
% 
clc;
clearvars;
close all;

%% postsynaptic activity
disp('Load postsynaptic spike-time lists and parameters.');
load('n_post_activity_results.mat', ...
  't_f_post_lists', 't_f_post_count', ...
  'sim_time', 'num_sim', ...
  'input_rate_list', 'num_ir', ...
  'synapse_weight_list', 'num_sw', ...
  'threshold_list', 'num_th' );

%% first spike time
disp('Calculating average first spike times.');
results_first_spike_times = zeros(num_th, num_sw, num_sim, num_ir);
results_first_spike_avg = zeros(num_th, num_sw, num_ir);
for irn = 1:num_ir
  for swn = 1:num_sw
    for thn = 1:num_th
      t_f_temp = [t_f_post_lists{thn, swn, :, irn}];
      results_first_spike_times(thn, swn, :, irn) = t_f_temp(1,:);
      not_inf = results_first_spike_times(thn, swn, :, irn) ~= Inf;
      results_first_spike_avg(thn, swn, irn) = ...
          sum(results_first_spike_times(thn, swn, not_inf, irn)) ...
          /sum(not_inf);
    end
  end
end

%% average postsynaptic spike-rates (taking all 1s simulations)
disp('Calculating average spike rates.');
results_spike_rates_avg = zeros(num_th, num_sw, num_ir);
for irn = 1:num_ir
  for swn = 1:num_sw
    for thn = 1:num_th
      results_spike_rates_avg(thn, swn, irn) = ...
          sum(t_f_post_count(thn, swn, :, irn)) / (num_sim * sim_time);
    end
  end
end
