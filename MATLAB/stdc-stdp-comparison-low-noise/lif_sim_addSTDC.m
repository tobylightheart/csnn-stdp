%% addSTDC_sim
% receives presynaptic neuron spike times 't_f_pre_list' and simulates
% postsynaptic neuron construction with addSTDP

%% Construction parameters
stdp_estimate = [1; 5; 10; 15; 20];   % # plasticity iterations estimated
num_stdp_e = length(stdp_estimate);  % number of estimate values

c_max = 100;    % maximum number of neurons allowed to be constructed
                % first constructed neuron has anomalous weights
%c_update = 1;   % acceleration of weight updates
c_status = 0;   % construction status, 0: no construction, 1: in progress 

%% Network parameters
n_pre = 2000;       % number of presynaptic neurons (must match poisson_hidden_pattern_gen)

n_post_trig = 1;    % number of constructive trigger neurons
n_post = 0;         % number of constructed neurons

%% Synapse weights
w_max = 1;    % maximum synapse weight
w_min = 0;    % minimum synapse weight

%w_mu = (w_max-w_min)/2; % initial weight mean
%w_sd = 0.1;   % synapse weight standard deviation
%w_IJ = w_sd*randn + w_mu*ones(n_pre,n_post); % initial weights of addSTDP neurons
%w_IJ = (w_max - w_min)*rand(n_pre,n_post_astdp) + w_min;

w_0 = 0.5;
w_IJ_0 = w_0*ones(n_pre,1);

w_IJ_trig = w_0*ones(n_pre,n_post_trig); % weights of trigger neurons

w_IJ_initial = mat2cell(zeros(num_stdp_e*n_pre,c_max),n_pre*ones(num_stdp_e,1),c_max);
w_IJ = mat2cell(zeros(num_stdp_e*n_pre,c_max),n_pre*ones(num_stdp_e,1),c_max);

%% Spike-Timing-Dependent Plasticity parameters
long_time = 1e6;

% Additive STDP
A_p = 0.03125; %0.03125;
A_n = 0.85*A_p; %0.85*A_p;
STDP_tau_p = 0.0168;
STDP_tau_n = 0.0337;

% addSTDP updates for first presynaptic and postsynaptic spike
first_f_pre = mat2cell(zeros(num_stdp_e*n_pre,c_max),n_pre*ones(num_stdp_e,1),c_max);
first_f_post = mat2cell(zeros(num_stdp_e*n_pre,c_max),n_pre*ones(num_stdp_e,1),c_max); % initially no neurons

% Spike-timing-dependent plasticity estimation
theta_a = 3*STDP_tau_p;
theta_b = 2*STDP_tau_n;


%% Neuron activity
% replicate and store list of spike times for simulation
%t_f_pre_list = sort([spike_list, spike_list+sim_time, spike_list+2*sim_time],2);

t_f_post_list_max = sim_time*20; %(repeat_n+1)*3 + 1000;
t_f_post_list_trig_max = sim_time*20;
t_f_post_list_inc = 20;

t_f_post_list = mat2cell(Inf(num_stdp_e*c_max,t_f_post_list_max),c_max*ones(num_stdp_e,1),t_f_post_list_max); % addSTDC neuron spike times (to be recorded)
t_f_post_list_trig = Inf(n_post_trig,t_f_post_list_trig_max);    % trigger neuron spike times (to be recorded)

t_f_post_list_count = mat2cell(zeros(num_stdp_e*c_max,1), c_max*ones(num_stdp_e,1), 1);
t_f_post_list_trig_count = zeros(n_post_trig,1);

%% LIF neuron parameters
p_thr = 500;       % neuron potential threshold
p_rest = 0;               % neuron potential resting value

p_reset_m = p_thr*(2-4);  % neuron membrane reset potential
p_reset_s = p_thr*4;      % neuron synapse reset potential

p_tau_m = 0.010;      % neuron membrane potential time constant
p_tau_s = 0.0025;     % neuron postsynaptic potential time constant

K_EPSP = 4^(4/3)/3; % normalisation factor for synaptic transmission

refractory_period = 0.001;   % refractory period length
t_refrac = mat2cell(-long_time*ones(num_stdp_e*c_max,1),c_max*ones(num_stdp_e,1),1);    % refractory time after spiking

% spiking neuron potentials
p_post_m_trig = p_rest*ones(n_post_trig,1);
p_post_s_trig = p_rest*ones(n_post_trig,1);

p_post_m = mat2cell(p_rest*ones(num_stdp_e*c_max,1),c_max*ones(num_stdp_e,1),1);
p_post_s = mat2cell(p_rest*ones(num_stdp_e*c_max,1),c_max*ones(num_stdp_e,1),1);

%% simulation
t_now = 0;  % time of current neuron spike

t_a = -long_time; % construction time window start
t_b = -long_time; % construction time window end

% neuron spike times within window for construction
Dt_f_trig = [long_time*ones(n_pre,1), -long_time*ones(n_pre,1)];
TRIG_POT = 1; % potentiating spike time matrix column
TRIG_DEP = 2; % depressing spike time matrix column

t_f_post = mat2cell(-long_time*ones(num_stdp_e*c_max,1),c_max*ones(num_stdp_e,1),1);
t_f_post_trig = -long_time*ones(n_post_trig,1);

t_f_pre = -long_time*ones(n_pre,1);
t_f_pre_list_count = ones(n_pre,1);
t_f_pre_list_temp = t_f_pre_list(:,1);

t_out = 0;
t_out_period = 0.1;

disp(['Simulation time: ' num2str(t_out, '%3.1f') ...
  '; # Spikes: ' num2str(t_f_post_list_count{1}(1))]);

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
      '; # Spikes: ' num2str(t_f_post_list_count{1}(1))]);
  end
  
  %% 2. check neuron construction state
  % construction ends if time has expired
  if c_status && (t_now > t_b)
    % reset construction status
    c_status = 0;
    
    % Does the simulation have less than the maximum number of
    % constructed neurons?
    if n_post < c_max
      % yes: construct a neuron
      n_post = n_post + 1;      
      
      % calculated synapse weight estimates 
      % addSTDP estimation
      u_p = Dt_f_trig(:,TRIG_POT)<long_time;
      u_n = Dt_f_trig(:,TRIG_DEP)>-long_time; 

      for se = 1:num_stdp_e

        w_IJ_new = w_IJ_0;

        % potentiation update only
        w_IJ_new(u_p&~u_n) = min(w_max, ...
            w_IJ_0(u_p&~u_n) + stdp_estimate(se) * A_p ...
            * exp(-Dt_f_trig(u_p&~u_n,TRIG_POT)/STDP_tau_p));  
        % depression update only
        w_IJ_new(~u_p&u_n) = max(w_min, ...
            w_IJ_0(~u_p&u_n) - stdp_estimate(se) * A_n ...
            * exp(Dt_f_trig(~u_p&u_n,TRIG_DEP)/STDP_tau_n));
        % potentiation and depression updates
        w_IJ_new(u_p&u_n) = min(w_max - ...
            A_n*exp(Dt_f_trig(u_p&u_n,TRIG_DEP)/STDP_tau_n), ...
            max(w_min, w_IJ_0(u_p&u_n) + stdp_estimate(se) * ...
            (A_p * exp(-Dt_f_trig(u_p&u_n,TRIG_POT)/STDP_tau_p) - ...
            A_n * exp(Dt_f_trig(u_p&u_n,TRIG_DEP)/STDP_tau_n))));

        % first-spike matrix values are initialised assuming a postsynaptic 
        % spike at trigger time
        first_f_pre{se}(:,n_post) = ~u_n; 
        first_f_post{se}(:,n_post) = u_n;

        % add neuron weight estimates
        w_IJ_initial{se}(:,n_post) = w_IJ_new;
        w_IJ{se}(:,n_post) = w_IJ_new;
      end
    end
  end
  
  %% 3. update the postsynaptic neuron potential
  %  decay factor
  p_dec_m = exp(-Dt/p_tau_m);
  p_dec_s = exp(-Dt/p_tau_s);
  
  % constructed neurons
  if n_post > 0
    for se = 1:num_stdp_e
      p_post_m{se} = p_post_m{se} * p_dec_m;
      p_post_s{se} = p_post_s{se} * p_dec_s;
    end
  end
  
  % increase potential
  % construction trigger neurons
  if n_post_trig > 0
    if c_status
      % construction is in progress:
      %   pause neuron potential updates
      %   record the first presynaptic spike for weight depression
      if Dt_f_trig(f_pre_neuron,TRIG_DEP) < t_f_post_trig-t_now
        Dt_f_trig(f_pre_neuron,TRIG_DEP) = t_f_post_trig-t_now;
      end
    else
      p_post_m_trig = p_post_m_trig * p_dec_m;
      p_post_s_trig = p_post_s_trig * p_dec_s;
    end
  end
 
  %% 4. detect postsynaptic neuron spikes
  % constructed neurons
  if n_post > 0
    f_post = mat2cell(zeros(num_stdp_e*n_post,1),n_post*ones(num_stdp_e,1),1);
    for se = 1:num_stdp_e
      % check for neurons in refractory period 
      n_refrac = t_now <= t_refrac{se};
      % find neurons that exceed threshold
      f_post{se} = ((p_post_m{se} + p_post_s{se}) >= p_thr) & ~n_refrac;
      if any(f_post{se})
        % reset neuron potential
        p_post_m{se}(f_post{se}) = p_reset_m;
        p_post_s{se}(f_post{se}) = p_reset_s;

        % record postsynpatic neurons spike time
        t_f_post{se}(f_post{se}) = t_now;

        % record refractory period
        t_refrac{se}(f_post{se}) = t_now + refractory_period;

        % store postsynaptic spike time in the list
        t_f_post_list_count{se}(f_post{se}) = t_f_post_list_count{se}(f_post{se})+1;
        while any(t_f_post_list_count{se} > t_f_post_list_max)
          disp('WARNING: t_f_post_list bounds exceeded: size increased.')
          %disp('Paused: Press any key to continue.');
          %pause;
          t_f_post_list_max = t_f_post_list_max + t_f_post_list_inc;
          for list_num = 1:num_stdp_e
            t_f_post_list{list_num} = [t_f_post_list{list_num}, Inf(c_max, t_f_post_list_inc)];
          end
        end

        list_update = find(f_post{se})';
        for lu = list_update
          t_f_post_list{se}(lu,t_f_post_list_count{se}(lu)) = t_now;
        end
      end
    end
  end
  
  %% 4a. detect trigger neuron spikes: start construction
  % trigger neurons
  if n_post_trig > 0
    % check if enoughtime has elapsed since last construction
    % find neurons that exceed threshold
    f_post_trig = ((p_post_m_trig + p_post_s_trig) >= p_thr) & (t_now > t_b);
    if any(f_post_trig)
      % start construction window for depressed weights
      c_status = 1;
      
      % reset spike-time difference lists
      Dt_f_trig = [long_time*ones(n_pre,1), -long_time*ones(n_pre,1)];

      % record difference in pre- and postsynaptic neuron spike times
      Dt_f_trig(:,TRIG_POT) = t_now - t_f_pre;
      Dt_f_trig(Dt_f_trig(:,TRIG_POT)>theta_a,TRIG_POT) = long_time;

      % reset neuron potential
      p_post_m_trig(f_post_trig) = p_rest;
      p_post_s_trig(f_post_trig) = p_rest;

      % construction refractory period
      t_b = t_now + theta_b;

      % record postsynpatic neurons spike time
      t_f_post_trig(f_post_trig) = t_now;
      t_f_post_list_trig_count(f_post_trig) = t_f_post_list_trig_count(f_post_trig)+1;
      while any(t_f_post_list_trig_count > t_f_post_list_trig_max)
        disp('WARNING: t_f_post_list_trig bounds exceeded: size increased.')
        t_f_post_list_trig_max = t_f_post_list_trig_max + t_f_post_list_inc;
        t_f_post_list_trig = [t_f_post_list_trig, Inf(n_post_trig, t_f_post_list_inc)];
      end
      
      list_update = find(f_post_trig)';
      for lu = list_update
        t_f_post_list_trig(lu,t_f_post_list_trig_count(lu)) = t_now;
      end
    end
  end
  
  %% 5. STDP: potentiate connections of active postsynaptic neurons
  % update constructed neurons 
  if n_post > 0
    for se = 1:num_stdp_e
      if any(f_post{se})
        % calculate difference in spike times
        Dt_f_astdc = t_now - t_f_pre;

        % weight update for each synapse
        Dw_IJ = A_p.*exp(-Dt_f_astdc/STDP_tau_p);

        % update selected weights
        w_IJ{se}(:,f_post{se}) = w_IJ{se}(:,f_post{se}) + ...
          bsxfun(@times, first_f_post{se}(:,f_post{se}), Dw_IJ);

        % enforce maximum weight limit
        w_IJ{se}(:,f_post{se}) = min(w_max, w_IJ{se}(:,f_post{se}));

        % the next presynaptic spikes will be the first
        first_f_pre{se}(:,f_post{se}) = 1;

        % the next postsynaptic spikes will not be first
        first_f_post{se}(:,f_post{se}) = 0;

  %       for n_loop_post = find(f_post_astdc)'
  %         % presynaptic neurons that have not fired since the previous
  %         % postsynaptic spike do not have synapses updated
  %         Dt_f_astdc(~first_f_post_astdc(:,n_loop_post)) = inf;
  % 
  %         % update weights
  %         w_IJ_astdc(:,n_loop_post) = min(w_max, w_IJ_astdc(:,n_loop_post) + A_p.*exp(-Dt_f_astdc/STDP_tau_p));
  %       end    
      end 
    end
  end

  %% 6. update presynaptic activity times, spike count and temp spike list
  t_f_pre(f_pre_neuron) = t_now;
  t_f_pre_list_count(f_pre_neuron) = t_f_pre_list_count(f_pre_neuron)+1;
  for n = f_pre_neuron'
    t_f_pre_list_temp(n) = t_f_pre_list(n,t_f_pre_list_count(n));
  end
  
  %% 7. STDP: depress connections of active presynaptic neurons
  % trigger neurons do not receive weight updates
  
  % constructed neurons
  if n_post > 0
    for se = 1:num_stdp_e
      % find presynaptic neurons that are spiking first after the
      % postsynpatic neuron
      first_f_pre_neuron = f_pre_neuron(any(first_f_pre{se}(f_pre_neuron,:),2));

      % cancel updates if no spike is the first    
      if any(first_f_pre_neuron)
        % difference in post- and pre-synaptic spike times
        Dt_f_astdc = t_f_post{se} - t_now;  

        % weight update for each synapse
        Dw_IJ = -A_n*exp(Dt_f_astdc/STDP_tau_n);

        % update selected weights (first presynaptic spike only)
        w_IJ{se}(first_f_pre_neuron,:) = w_IJ{se}(first_f_pre_neuron,:) ...
          + bsxfun(@times, first_f_pre{se}(first_f_pre_neuron,:), Dw_IJ');

        % enforce minimum weight
        w_IJ{se}(first_f_pre_neuron,:) = max(w_min, w_IJ{se}(first_f_pre_neuron,:));

      end
      % no longer first presynaptic spike
      first_f_pre{se}(f_pre_neuron,:) = 0;
      % now the next postsynaptic spike will cause an update
      first_f_post{se}(f_pre_neuron,:) = 1;
    end
  end
  
  %% 8. start transmission of postsynaptic potential to neurons
  if n_post > 0
    for se = 1:num_stdp_e
      if size(f_pre_neuron,1) > 1
        p_post_m{se} = p_post_m{se} + K_EPSP*sum(w_IJ{se}(f_pre_neuron,:))';
        p_post_s{se} = p_post_s{se} - K_EPSP*sum(w_IJ{se}(f_pre_neuron,:))';
      else
        p_post_m{se} = p_post_m{se} + K_EPSP*w_IJ{se}(f_pre_neuron,:)';
        p_post_s{se} = p_post_s{se} - K_EPSP*w_IJ{se}(f_pre_neuron,:)';
      end
    end
  end

  if (n_post_trig > 0) && (t_now > t_b)
    if size(f_pre_neuron,1) > 1
      p_post_m_trig = p_post_m_trig + K_EPSP*sum(w_IJ_trig(f_pre_neuron,:))';
      p_post_s_trig = p_post_s_trig - K_EPSP*sum(w_IJ_trig(f_pre_neuron,:))';
    else
      p_post_m_trig = p_post_m_trig + K_EPSP*w_IJ_trig(f_pre_neuron,:)';
      p_post_s_trig = p_post_s_trig - K_EPSP*w_IJ_trig(f_pre_neuron,:)';
    end
  end  
end

% if the simulation ends before construction completes
if c_status
  % reset construction status
  c_status = 0;

  % Does the simulation have less than the maximum number of
  % constructed neurons?
  if n_post < c_max
    % yes: construct a neuron
    n_post = n_post + 1;      

    % calculated synapse weight estimates 
    % addSTDP estimation
    u_p = Dt_f_trig(:,TRIG_POT)<long_time;
    u_n = Dt_f_trig(:,TRIG_DEP)>-long_time; 

    for se = 1:num_stdp_e

      w_IJ_new = w_IJ_0;

      % potentiation update only
      w_IJ_new(u_p&~u_n) = min(w_max, ...
          w_IJ_0(u_p&~u_n) + stdp_estimate(se) * A_p ...
          * exp(-Dt_f_trig(u_p&~u_n,TRIG_POT)/STDP_tau_p));  
      % depression update only
      w_IJ_new(~u_p&u_n) = max(w_min, ...
          w_IJ_0(~u_p&u_n) - stdp_estimate(se) * A_n ...
          * exp(Dt_f_trig(~u_p&u_n,TRIG_DEP)/STDP_tau_n));
      % potentiation and depression updates
      w_IJ_new(u_p&u_n) = min(w_max - ...
          A_n*exp(Dt_f_trig(u_p&u_n,TRIG_DEP)/STDP_tau_n), ...
          max(w_min, w_IJ_0(u_p&u_n) + stdp_estimate(se) * ...
          (A_p * exp(-Dt_f_trig(u_p&u_n,TRIG_POT)/STDP_tau_p) - ...
          A_n * exp(Dt_f_trig(u_p&u_n,TRIG_DEP)/STDP_tau_n))));

      % first-spike matrix values are initialised assuming a postsynaptic 
      % spike at trigger time
      first_f_pre{se}(:,n_post) = ~u_n; 
      first_f_post{se}(:,n_post) = u_n;

      % add neuron weight estimates
      w_IJ_initial{se}(:,n_post) = w_IJ_new;
      w_IJ{se}(:,n_post) = w_IJ_new;
    end
  end
end
