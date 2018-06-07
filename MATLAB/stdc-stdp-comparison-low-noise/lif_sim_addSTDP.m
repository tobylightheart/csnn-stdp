%% lif_sim_addSTDP.m
% receives a presynaptic neuron spike times 't_f_pre_list' and simulates
% postsynaptic neurons with addSTDP

%% Network parameters
n_pre = 2000;       % number of presynaptic neurons (set by gen_poisson)
n_post = 1;   % number of additive STDP postsynaptic neurons

%% Synapse weights
w_max = 1;    % maximum synapse weight
w_min = 0;    % minimum synapse weight

w_IJ = 0.5*ones(n_pre,n_post); %0.475*ones(n_pre,n_post);

%% Spike-Timing-Dependent Plasticity parameters
long_time = 1e6;

% Additive STDP
A_p = 0.03125; 
A_n = 0.85*A_p; 
STDP_tau_p = 0.0168;
STDP_tau_n = 0.0337;

% addSTDP updates for first presynaptic and postsynaptic spike
first_f_pre = zeros(n_pre, n_post); 
first_f_post = zeros(n_pre, n_post); 

%% Neuron activity
% replicate and store list of spike times for simulation
%t_f_pre_list = sort([spike_list, spike_list+sim_time, spike_list+2*sim_time],2);

t_f_post_list_max = sim_time*20; %(repeat_n+1)*3 + 1000;
t_f_post_list_inc = 20;

% addSTDP neuron spike times (to be recorded)
t_f_post_list = Inf(n_post, t_f_post_list_max);
t_f_post_list_count = zeros(n_post,1);

% postsynaptic neuron spike delay to record synapse weights
theta_b = 2*STDP_tau_n;

w_record = false;
w_IJ_record = Inf(n_pre, t_f_post_list_max, n_post);

%% LIF neuron parameters
p_thr = 0.25*n_pre;       % neuron potential threshold
p_rest = 0;               % neuron potential resting value

p_reset_m = p_thr*(2-4);  % neuron membrane reset potential
p_reset_s = p_thr*4;      % neuron synapse reset potential

p_tau_m = 0.010;      % neuron membrane potential time constant
p_tau_s = 0.0025;     % neuron postsynaptic potential time constant

K_EPSP = 4^(4/3)/3;   % normalisation factor for synaptic transmission

refractory_period = 0.001;   % refractory period length
t_refrac = -long_time*ones(n_post,1);   % refractory time after spiking

% spiking neuron potentials
p_post_m = p_rest*ones(n_post,1);
p_post_s = p_rest*ones(n_post,1);

%% simulation
t_now = 0;  % time of current neuron spike

t_b = -long_time; % record weights time

t_f_post = -long_time*ones(n_post,1);

t_f_pre = -long_time*ones(n_pre,1);
t_f_pre_list_count = ones(n_pre,1);
t_f_pre_list_temp = t_f_pre_list(:,1);

t_out = 0;
t_out_period = 0.1;

disp(['Simulation time: ' num2str(t_out, '%3.1f') ...
  '; # Spikes: ' num2str(t_f_post_list_count)]);

% iterate through each presynaptic spike
while any(t_f_pre_list_temp(:)~=inf)
  %% 1. find next presynaptic neuron spike
  t_next = min(t_f_pre_list_temp(:));
  Dt = t_next-t_now;
  t_now = t_next;  
  f_pre_neuron = find(t_f_pre_list_temp == t_now);
  
  % 1a. display time
  if floor(t_now/t_out_period) >= round(t_out/t_out_period)+1
    t_out = t_out+t_out_period;       
    disp(['Simulation time: ' num2str(t_out, '%3.1f') ...
      '; # Spikes: ' num2str(t_f_post_list_count)]);
  end
  
  % 1b. record weights if sufficient time after a postsynaptic neuron spike
  if any(w_record) && (t_now > t_b) 
    % record new weight vector
    w_IJ_record(:,t_f_post_list_count(w_record)) = w_IJ(:,w_record);     
    w_record(w_record) = false;
  end
  
  %% 2. update the postsynaptic neuron potential and transmission
  % STDP only neurons
  %  decay factor
  p_dec_m = exp(-Dt/p_tau_m);
  p_dec_s = exp(-Dt/p_tau_s);
  
  if n_post > 0
    % update potential from presynaptic spike
    p_post_m = p_post_m * p_dec_m;
    p_post_s = p_post_s * p_dec_s;
  end
    
  %% 3. detect postsynaptic neuron spikes
  % STDP only neurons
  if n_post > 0
    % check for neurons in refractory period
    n_refrac = t_now <= t_refrac;
    % find neurons that exceed threshold
    f_post = ((p_post_m + p_post_s) >= p_thr) & ~n_refrac;
    if any(f_post)
      % reset neuron potential
      p_post_m(f_post) = p_reset_m;      
      p_post_s(f_post) = p_reset_s;

      % record postsynpatic neurons spike time
      t_f_post(f_post) = t_now;

      % record refractory period
      t_refrac(f_post) = t_now + refractory_period;
      
      % check synapse weight record status
      if w_record(f_post)
        disp('WARNING: Postsynaptic spike occurred before weight recording time.');
        w_IJ_record(:,t_f_post_list_count(w_record)) = w_IJ(:,w_record);
      else
        w_record(f_post) = true;
      end
      t_b = t_now + theta_b;

      % store postsynaptic spike time in the list
      t_f_post_list_count = t_f_post_list_count+1;
      while t_f_post_list_count > t_f_post_list_max
        disp('WARNING: t_f_post_list bounds exceeded: size increased.');
        t_f_post_list_max = t_f_post_list_max + t_f_post_list_inc;
        t_f_post_list = [t_f_post_list, Inf(n_post, t_f_post_list_inc)];
      end
      t_f_post_list(t_f_post_list_count) = t_now;
    end
  end

  %% 4. STDP: potentiate connections of active postsynaptic neurons
  % update STDP only neurons
  if n_post > 0
    if any(f_post)
      % difference in spike times
      Dt_f = t_now - t_f_pre;
      
      % weight update for each synapse
      Dw_IJ = A_p.*exp(-Dt_f/STDP_tau_p);
      
      % update selected weights
      w_IJ(:,f_post) = w_IJ(:,f_post) ...
        + bsxfun(@times, first_f_post(:,f_post), Dw_IJ);
      
      % enforce maximum weight limit
      w_IJ(:,f_post) = min(w_max, w_IJ(:,f_post));
      
      % the next presynaptic spike will be the first 
      first_f_pre(:,f_post) = 1;
      
      % the next postsynaptic spike will not be first
      first_f_post(:,f_post) = 0;
    end
  end

  %% 5. update presynaptic activity times, spike count and temp spike list
  t_f_pre(f_pre_neuron) = t_now;
  t_f_pre_list_count(f_pre_neuron) = t_f_pre_list_count(f_pre_neuron)+1;
  for n = f_pre_neuron'
    t_f_pre_list_temp(n) = t_f_pre_list(n,t_f_pre_list_count(n));
  end
  
  %% 6. STDP: depress connections of active presynaptic neurons  
  % STDP only neurons
  if n_post > 0
    % find presynaptic neurons that are spiking first after the
    % postsynpatic neuron
    first_f_pre_neuron = f_pre_neuron(any(first_f_pre(f_pre_neuron,:),2));
    
    % cancel updates if no spike is the first
    if any(first_f_pre_neuron) 
      % calculate the difference in post- and pre-synaptic spike times
      Dt_f = t_f_post - t_now;
      
      % update for each synapse
      Dw_IJ = -A_n*exp(Dt_f/STDP_tau_n);
      
      % update selected weights (first presynaptic neuron spikes only)
      w_IJ(first_f_pre_neuron,:) = w_IJ(first_f_pre_neuron,:) ...
        + bsxfun(@times, first_f_pre(first_f_pre_neuron,:), Dw_IJ');
      
      % enforce minimum weight
      w_IJ(first_f_pre_neuron,:) = max(w_min, w_IJ(first_f_pre_neuron,:));
    end
    
    % no longer first presynaptic spike
    first_f_pre(f_pre_neuron,:) = 0;
    % now the next postsynaptic spike will cause an update
    first_f_post(f_pre_neuron,:) = 1;
  end  
  
  %% 7. Start presynaptic spike transmission
  % STDP only neurons
  if n_post > 0
    if size(f_pre_neuron,1) > 1
      p_post_m = p_post_m + K_EPSP*sum(w_IJ(f_pre_neuron,:))';
      p_post_s = p_post_s - K_EPSP*sum(w_IJ(f_pre_neuron,:))';
    else
      p_post_m = p_post_m + K_EPSP*w_IJ(f_pre_neuron,:)'; 
      p_post_s = p_post_s - K_EPSP*w_IJ(f_pre_neuron,:)'; 
    end   
  end
  
end

% if simulation ends before recording
if any(w_record)
  % record new weight vector
  w_IJ_record(:,t_f_post_list_count(w_record)) = w_IJ(:,w_record);     
  w_record(w_record) = false;
end