%% Poisson patterns with noise
% Modified from Masquelier et al (2008)
% requires r_b in the workspace (taken from run_lif_sim_batches*)
disp('Initialising spike time generation parameters.');

% pattern repetition segment length and time-steps
pattern_time = 0.02; % pattern activity lasts 20ms
segment_time_pre = 0.1; % 100ms before the pattern
segment_time_post = 0.1 - pattern_time; % time after the pattern

segment_time = segment_time_pre + pattern_time + segment_time_post;
sim_time = 21*segment_time;   % simulation time

dt = 0.001;       % time step (1 millisecond), time units seconds

n_pre = 2000;     % number of neurons
n_pattern = 1000; % number of neurons in each repeating pattern
n_pattern_list = 1:n_pattern;
n_background = 1000;
n_background_list = n_pattern + 1:n_background;

% defining pattern segments parameters
% select pattern neuron (logical array)
pattern_neuron = [true(n_pattern,1); ...
    false(n_background,1)];

pattern_start_step = round(segment_time_pre/dt);
pattern_stop_step = round((segment_time_pre + pattern_time)/dt);
segment_stop_step = round(segment_time/dt);
sim_stop_step = round(sim_time/dt);

repeat_n = round(sim_time/segment_time);

%% Generating neuron spike times
disp('Generating neuron pattern spike-times.')

% number of pattern spikes
num_ps = 80 * n_pattern * pattern_time;  % 80Hz pattern activity  
spike_pattern_list_max = round(300*pattern_time);
spike_pattern_list_inc = 10;
spike_pattern_list = Inf(n_pre, spike_pattern_list_max);
             
count_spl = [ones(n_pattern,1); zeros(n_background,1)];  % count spike pattern list
spike_pattern_list(n_pattern_list,1) = rand(n_pattern,1)*pattern_time;

num_ps_rem = num_ps - n_pattern;
spl = [randi(n_pattern,num_ps_rem,1), rand(num_ps_rem,1)*pattern_time];

for k = 1:num_ps_rem
  count_spl(spl(k,1)) = count_spl(spl(k,1)) + 1;
  
  % increase the spike pattern list length if necessary
  if any(count_spl >= spike_pattern_list_max)
    spike_pattern_list = [spike_pattern_list, ...
      Inf(n_pre, spike_pattern_list_inc)];
    spike_pattern_list_max = spike_pattern_list_max ...
      + spike_pattern_list_inc;
  end
  
  spike_pattern_list(spl(k,1), count_spl(spl(k,1))) = spl(k,2);
end

disp('Create pattern repetition spike times.')
% add repeat patterns to spike train
% find occupied spike pattern list
spl_end = find(sum(any(spike_pattern_list<Inf),3),1,'last');
spike_repeat_list = Inf(n_pre, repeat_n*spl_end);

for k = 1:repeat_n
  spike_repeat_list(:,(1+(k-1)*spl_end):(k*spl_end)) = ...
      spike_pattern_list(:,1:spl_end) ...
      + segment_time*(k-1) ...
      + segment_time_pre ...
      + dt*randn(n_pre,spl_end);
end

%% Generate background spike times
disp('Generating background spike-times.');

%r_b maximum background firing rate (Hz) taken from run_lif_sim_batches

spike_list_max = 100*sim_time;
spike_list_inc = 100;
spike_list = Inf(n_pre, spike_list_max);
count_sl = zeros(n_pre,1);    % count spike list entries

f_ind = (1:n_pre)';  
  
%simulate neuron activity
for t_step = 0:sim_stop_step
  segment_step = mod(t_step,segment_stop_step);
  
  % background neuron spikes (non-pattern)
  f_b = r_b*dt > rand(n_pre,1);
 
  % clear background spikes for pattern neurons in pattern time
  if (segment_step >= pattern_start_step) ...
      && (segment_step < pattern_stop_step)
    
    f_b(pattern_neuron) = false;
  end
  
  f_b_neuron = find(f_b);
  f_b_num = sum(f_b);
  t_f_b = dt*(t_step + rand(f_b_num, 1));

  count_sl(f_b) = count_sl(f_b) + 1;
  if any(count_sl > spike_list_max)
    spike_list = [spike_list, Inf(n_pre, spike_list_inc)];
    spike_list_max = spike_list_max + spike_list_inc;
  end

  for k = 1:f_b_num
    spike_list(f_b_neuron(k), count_sl(f_b_neuron(k))) = ...
      t_f_b(k);
  end      
end

spike_list = [spike_list, spike_repeat_list];

disp('Sort and trim spike lists for storage.');
% spike list: sort and trim (sparse increases storage requirements)
spike_list = sort(spike_list,2);
% trim spike list and terminate with infinity column
spike_list = [spike_list(:,1:find(any(spike_list<Inf),1,'last')), Inf(n_pre,1)];

% spike repeat list: sort, 
spike_repeat_list = sort(spike_repeat_list,2);

disp('Spike generation complete.');
