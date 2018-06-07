%% figure_postsynaptic_spike_times.m
% clc;
% clearvars;
% close all;

if isempty(who('save_all_figures*'))
	save_all_figures_as_pdfs = 0;
	save_all_figures_as_figs = 0;
end

%% set up simulation constants and data containers
SEC_TO_MSEC = 1000;

num_sim = 100;
sim_list = 1:100;
sim_list_l = length(sim_list);

r_b_list = 0:10:20;
rate_list = 1:3;
rate_list_l = length(rate_list);

pattern_time = 0.02;
segment_time_pre = 0.1;
segment_time_post = 0.08;
segment_time = segment_time_pre + pattern_time + segment_time_post;

dir_expd = 'lif_sim_data_addSTDC'; 
dir_stdp = 'lif_sim_data_addSTDP'; 

what_expd = what(dir_expd);
what_stdp = what(dir_stdp);

load(fullfile(dir_expd, ...
    what_expd.mat{1}), ...
    'num_stdp_e', ...
    'stdp_estimate', ...
    't_f_post_list_count', ...
    't_f_post_list_trig_count');

n_plot = sum(t_f_post_list_count{1}>0); 

t_f_post_expd = zeros(sim_list_l, ...
    n_plot, ...
    num_stdp_e, ...
    rate_list_l);
  
t_f_post_expd_trig = zeros(sim_list_l * t_f_post_list_trig_count, ...
    rate_list_l);    
    
load(fullfile(dir_stdp, ...
  what_stdp.mat{1}), ...
  't_f_post_list_count');

t_f_post_stdp = zeros(sim_list_l, t_f_post_list_count, rate_list_l);
post_lists = zeros(n_plot*sim_list_l, num_stdp_e, rate_list_l);

% figure parameters
x_lim_stdp = [-1 21; -1 21; -1 21];
y_lim_time = [7 21; 4 14; 4 14];

fig_h = [10; 7; 7];
fig_w = [12; 12; 12];

box_w_expd = 0.7;
box_w_stdp = 1;

MS = 'MarkerSize';
marker_s_expd = 4;
marker_s_stdp = 4;
MEC = 'MarkerEdgeColor';
color_stdp = [0.4, 0.4, 0.4];
color_expd = [0, 0, 0];

%% extract the postsynaptic spike-times from simulation results
for r = 1:rate_list_l
  for k = 1:sim_list_l
    disp(['Rate: ' num2str(r_b_list(rate_list(r))) ...
      '; Simulation: ' num2str(k)]);
    load(fullfile(dir_expd, ...
      what_expd.mat{(rate_list(r)-1)*num_sim + k}), ...
      't_f_post_list_count', ...
      't_f_post_list', ...
      't_f_post_list_trig_count', ...
      't_f_post_list_trig');

    for nse = 1:num_stdp_e
      for np = 1:n_plot
        t_f_post_expd(k,np,nse,r) = ...
          rem(t_f_post_list{nse}(np,1), segment_time) - segment_time_pre;
      end
    end
    
    expd_trig_ind = ...
      1+(k-1)*t_f_post_list_trig_count:(k)*t_f_post_list_trig_count;
    t_f_post_expd_trig(expd_trig_ind,r) = ...
      rem(t_f_post_list_trig(1:t_f_post_list_trig_count), segment_time) ...
      - segment_time_pre;
    
    load(fullfile(dir_stdp, ...
      what_stdp.mat{(rate_list(r)-1)*num_sim + k}), ...
      't_f_post_list_count', ...
      't_f_post_list');

    t_f_post_stdp(k,:,r) = t_f_post_list(1:t_f_post_list_count) ...
      - (0:t_f_post_list_count-1)*segment_time ...
      - segment_time_pre;
  end

  for nse = 1:num_stdp_e
    post_lists(:,nse,r) = reshape(t_f_post_expd(:,:,nse,r), ...
      [n_plot*sim_list_l, 1]);
  end
end

%% generate box plots for postsynaptic spike-time statistics
fig_boxplot = gobjects(3,1);
for r = 1:rate_list_l
  fig_boxplot(r) = figure;
  hold on;
  set(gcf, 'Color', 'w');
  set(gcf, 'Units', 'centimeters');
  set(gcf, 'Position', [5, 5, fig_w(r), fig_h(r)]);
  set(gcf, 'PaperPosition', [0, 0, fig_w(r), fig_h(r)]);
  set(gcf, 'PaperSize', [fig_w(r), fig_h(r)]);
  
  boxplot(t_f_post_stdp(:,:,r)*SEC_TO_MSEC, ...
    'width', box_w_stdp, ...
    'positions', (1:t_f_post_list_count)-1, ...
    'whisker', 0, ...
    'labels', (1:t_f_post_list_count)-1, ...
    'symbol', '', ...
    'colors', color_stdp, ...
    'boxstyle', 'filled');
  set(gca, 'XTickLabel', {' '});
  stdp_lists_min = min(t_f_post_stdp(:,:,r))*SEC_TO_MSEC;
  stdp_lists_max = max(t_f_post_stdp(:,:,r))*SEC_TO_MSEC;
  plot((1:t_f_post_list_count)-1, stdp_lists_min, '+', ...
    MEC, color_stdp, MS, marker_s_stdp);
  plot((1:t_f_post_list_count)-1, stdp_lists_max, '+', ...
    MEC, color_stdp, MS, marker_s_stdp);  
  
  boxplot(t_f_post_expd_trig(:,r)*SEC_TO_MSEC, ...
    'width', box_w_expd, ...
    'positions', 0, ...    
    'whisker', 0, ...
    'symbol', '', ...
    'colors', color_expd);     
  set(gca, 'XTickLabel', {' '});
  plot(0, min(t_f_post_expd_trig(:,r))*SEC_TO_MSEC, 'd', ... '^', ...
    MEC, color_expd, MS, marker_s_expd);
  plot(0, min(t_f_post_expd_trig(:,r))*SEC_TO_MSEC, '.', ...
    MEC, color_expd, MS, marker_s_expd);
  plot(0, max(t_f_post_expd_trig(:,r))*SEC_TO_MSEC, 'd', ... 'v', ...
    MEC, color_expd, MS, marker_s_expd);
  plot(0, max(t_f_post_expd_trig(:,r))*SEC_TO_MSEC, '.', ...
    MEC, color_expd, MS, marker_s_expd);
   
  boxplot(post_lists(:,:,r)*SEC_TO_MSEC, ...
    'width', box_w_expd, ...
    'positions', stdp_estimate, ...    
    'whisker', 0, ...
    'symbol', '', ...
    'colors', color_expd);    
  set(gca, 'XTickLabel', {' '});
  post_lists_min = min(post_lists(:,:,r))*SEC_TO_MSEC;
  post_lists_max = max(post_lists(:,:,r))*SEC_TO_MSEC;
  plot(stdp_estimate', post_lists_min, 'd', ... '^', ...
    MEC, color_expd, MS, marker_s_expd);
  plot(stdp_estimate', post_lists_min, '.', ...
    MEC, color_expd, MS, marker_s_expd);  
  plot(stdp_estimate', post_lists_max, 'd', ... 'v', ...
    MEC, color_expd, MS, marker_s_expd);
  plot(stdp_estimate', post_lists_max, '.', ...
    MEC, color_expd, MS, marker_s_expd);
  
  if r == 1
    record_stdp_max = stdp_lists_max;
    record_expd_max = post_lists_max;
  end
  
  if r == 1
    rect_p = [11.2,16.8,9,3.6];
    box_r = 1.2;
    box_x = rect_p(1)+1;
    box_y_stdp = rect_p(2)+2.0;
    box_y_expd = rect_p(2)+0.4;
    text_x = box_x+1.2;
  else
    if r == 2
      rect_p = [-0.2,4.6,9,3.9];
    else
      rect_p = [11.2,9.5,9,3.9];
    end
    box_r = 1.3;
    box_x = rect_p(1)+1;
    box_y_stdp = rect_p(2)+2.2;
    box_y_expd = rect_p(2)+0.4;
    text_x = box_x+1.2;
  end
  
  % start custom legend
  rectangle('Position', rect_p, 'FaceColor', 'w', 'LineWidth', 0.5);

  boxplot(box_y_stdp+box_r*[0:0.01:1], ...
      'width', box_w_stdp, ...
      'positions', box_x, ...
      'whisker', 0, ...
      'labels', box_x, ...
      'symbol', '', ...
      'colors', color_stdp, ...
      'boxstyle', 'filled'); 
  set(gca, 'XTickLabel', {' '});
  plot(box_x, box_y_stdp+box_r, '+', ...
      MEC, color_stdp, MS, marker_s_stdp);
  plot(box_x, box_y_stdp, '+', ...
      MEC, color_stdp, MS, marker_s_stdp);  

  boxplot(box_y_expd+box_r*[0:0.01:1], ...
      'width', box_w_expd, ...
      'positions', box_x, ...    
      'whisker', 0, ...
      'symbol', '', ...
      'colors', color_expd);    
  set(gca, 'XTickLabel', {' '});
  plot(box_x, box_y_expd, 'd', ...  '^', ...
      MEC, color_expd, MS, marker_s_expd);
  plot(box_x, box_y_expd, '.', ...
      MEC, color_expd, MS, marker_s_expd);  
  plot(box_x, box_y_expd+box_r, 'd', ... 'v', ...
      MEC, color_expd, MS, marker_s_expd);
  plot(box_x, box_y_expd+box_r, '.', ...
      MEC, color_expd, MS, marker_s_expd);

  text(text_x, box_y_stdp+box_r/2, 'Simulated STDP');
  text(text_x, box_y_expd+box_r/2, 'Construction');
  % custom legend complete
  
  set(gca, 'XTick', 0:5:t_f_post_list_count-1, ...
    'XTickLabel', 0:5:t_f_post_list_count-1);
  
  xlim(x_lim_stdp(r,:));
  ylim(y_lim_time(r,:));
  
  xlabel('Pattern Iteration, m');
  ylabel('Spike Latency [ms]');
  
  grid on;

end

%% save figures
fig_file_names = cell(3,1);
for r = 1:rate_list_l
  fig_file_names{r} = fullfile('figures', ...
      ['fig_ps_spike_times_rb',...
      num2str(r_b_list(rate_list(r)),'%02d')]);
end

if save_all_figures_as_pdfs
  for r = 1:rate_list_l
    FT_PDF = '.pdf';  
    saveas(fig_boxplot(r), [fig_file_names{r}, FT_PDF]);
    disp(['File saved as: ', fig_file_names{r}, FT_PDF]);
  end
end
if save_all_figures_as_figs	
  for r = 1:rate_list_l
    FT_FIG = '.fig';
    saveas(fig_boxplot(r), [fig_file_names{r}, FT_FIG]);
    disp(['File saved as: ', fig_file_names{r}, FT_FIG]);
  end
end
