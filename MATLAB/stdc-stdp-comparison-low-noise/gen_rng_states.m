% batch random number generator state generation
num_rng_states = 100;
rng_states = cell(num_rng_states,1);

for k = 1:num_rng_states
  
  rng('shuffle');
  rng_states{k} = rng;

end

save('rng_states_spike_list.mat', 'rng_states');
