%% figure_weights_lif_sim.m
% clc;
% clearvars;
% close all;

if isempty(who('save_all_figures*'))
	save_all_figures_as_pdfs = 0;
	save_all_figures_as_figs = 0;
end

CHAR_A = 'A';

num_sim = 100;
plot_sim = 1;

r_b_list = 0:10:20;
plot_rate = 2;

n_pattern = 1000;
w_0 = 0.5;

% background rate 0
% pattern neurons
% 1: STDP pos (and increasing); estimate pos
% 2: STDP first pos then neg; estimate pos and neg
% 6: STDP pos then neg; estimate pos and middle
% 11: STDP pos then neg; estimate pos
% background neurons
% none

% background rate 10
% pattern neurons
% 1: STDP pos; estimate wide pos distribution
% 2: STDP neg; estimate narrow neg distrubtion
% 22: STDP pos then neg; estimate all pos
% 23: STDP mixed then neg; estimate bimodal pos and neg
% background neurons
% 9: STDP small pos; estimate pos and middle
% 25: STDP neg; estimate pos, neg, and middle
% 1001: STDP mixed in middle; estimate wide dist over full range
% 1009: STDP slight neg; estimate all neg dist
% [excluded] 13: STDP one pos then neg; estimate bimodal pos and neg

% background rate 20
% pattern neurons
% 1: STDP start pos then neg; estimate bimodal pos and neg
% 2: STDP neg; estimate narrow neg
% 7: STDP mixed pos and neg; estimate bimodal pos and slight neg
% 8: STDP pos; estimate pos
% 10: STDP middle; estimate narrow middle
% background neurons
% 1001: STDP mixed in middle; estimate wide dist over full range
% 1016: STDP slight neg; estimate dist mostly neg, one pos outlier

%select_w = [11:20, 1,2,6];
switch plot_rate
  case 1
    select_w = [1,11,2,6];
  case 2
    select_w = [1,2,22,23,9,25,1001,1009]; %[1,2,9,13,22,1001,1009];
  case 3
    select_w = [1,2,7,8,10,1001];
end

num_w = length(select_w);

dir_stdp = 'lif_sim_data_addSTDP'; %_const_spl';
dir_expd = 'lif_sim_data_addSTDC'; %_const_spl';

what_stdp = what(dir_stdp);
what_expd = what(dir_expd);

sub_cols = 2;
sub_rows = ceil(num_w/2);

fig_h = 4 * sub_rows;
fig_w = 12;

MS = 'MarkerSize';
marker_s_expd = 5;
marker_s_stdp = 3;
MEC = 'MarkerEdgeColor';
color_stdp = [0.4, 0.4, 0.4];
color_expd = [0, 0, 0];

% Simulated STDP data
load(fullfile(dir_stdp, ...
  what_stdp.mat{(plot_rate-1)*num_sim + plot_sim}), ...
  't_f_post_list_count', ...
  'w_IJ_record');  

% Estimated STDP data
load(fullfile(dir_expd, ...
  what_expd.mat{(plot_rate-1)*num_sim + plot_sim}), ...
  'num_stdp_e', ...
  'stdp_estimate', ...
  't_f_post_list_trig_count', ...
  'w_IJ_initial');

w_IJ_expd_lists = zeros(t_f_post_list_trig_count,num_stdp_e);

fig_weights = figure;
%hold on;
set(fig_weights, 'Color', 'w');
set(fig_weights, 'Units', 'centimeters');
set(fig_weights, 'Position', [5, 5, fig_w, fig_h]);
set(fig_weights, 'PaperPosition', [0, 0, fig_w, fig_h]);
set(fig_weights, 'PaperSize', [fig_w, fig_h]); 

sp_h = cell(num_w,1);

for nw = 1:num_w %sw = select_w
   
  for nse = 1:num_stdp_e
    w_IJ_expd_lists(:,nse) = ...
      w_IJ_initial{nse}(select_w(nw),1:t_f_post_list_trig_count)';
  end
  
  sp_h{nw} = subplot(sub_rows, sub_cols, nw);
  hold on;
  
  plot(stdp_estimate, w_IJ_expd_lists, 'k.', ...
    MS, marker_s_expd, MEC, color_expd);
  
  stdp_update = 0:t_f_post_list_count-1;
  stdp_weight = [w_0, w_IJ_record(select_w(nw),1:t_f_post_list_count-1)];
  
  plot(stdp_update, stdp_weight, 'k+', ...
    MS, marker_s_stdp, MEC, color_stdp);
  
  set(sp_h{nw}, 'XTick', 0:5:t_f_post_list_count-1);
  set(sp_h{nw}, 'YTick', 0:0.25:1);
  xlim([-1 21]);
  ylim([-0.05 1.05]);
  grid on;
  set(gca, 'Box', 'on');
  
  if rem(num_w,2) == 0
    if ceil((nw-0.01)/2) >= sub_rows 
      xl_h = xlabel('STDP Update');
    else
      set(sp_h{nw}, 'XTickLabel', []);
    end
  else
    if ceil((nw+0.01)/2) >= sub_rows 
      xl_h = xlabel('STDP Update');
    else
      set(sp_h{nw}, 'XTickLabel', []);
    end
  end
  if rem(nw,2) == 1  
    yl_h = ylabel('Weight');
  else
    set(sp_h{nw}, 'YTickLabel', []);
  end
 
  text(0.75, 0.9, char(CHAR_A+nw-1));
end

switch plot_rate
  case 1
    for nw = 1:num_w
      set(sp_h{nw}, 'Position', get(sp_h{nw}, 'Position') + [0 0 0.05 0.05]);
      if nw > 2
        set(sp_h{nw}, 'Position', get(sp_h{nw}, 'Position') + [0 0.02 0 0]);
      end
      if rem(nw,2) == 0
        set(sp_h{nw}, 'Position', get(sp_h{nw}, 'Position') + [-0.01,0,0,0]);
      end      
    end
  case 2
    for nw = 1:num_w
      set(sp_h{nw}, 'Position', get(sp_h{nw}, 'Position') + [0,0,0.05,0.05]);
      set(sp_h{nw}, 'Position', ...
        get(sp_h{nw}, 'Position') + [0,-floor((nw-0.1)/2)*0.010,0,0]);
      if rem(nw,2) == 0
        set(sp_h{nw}, 'Position', get(sp_h{nw}, 'Position') + [-0.01,0,0,0]);
      end
    end
  case 3
    for nw = 1:num_w
      set(sp_h{nw}, 'Position', get(sp_h{nw}, 'Position') + [0,0,0.05,0.05]);
      if rem(nw,2) == 0
        set(sp_h{nw}, 'Position', get(sp_h{nw}, 'Position') + [-0.01,0,0,0]);
      end   
    end
end

%% save figure
fig_file_name = fullfile('figures', ...
	['fig_ps_synapse_updates_rb' ...
		num2str(r_b_list(plot_rate),'%02d')]);

if save_all_figures_as_pdfs
	FT_PDF = '.pdf';
	saveas(fig_weights,[fig_file_name, FT_PDF]);
	disp(['Figure saved as: ', fig_file_name, FT_PDF]);
end
if save_all_figures_as_figs
	FT_FIG = '.fig';
	saveas(fig_weights,[fig_file_name, FT_FIG]);
	disp(['Figure saved as: ', fig_file_name, FT_FIG]);
end
