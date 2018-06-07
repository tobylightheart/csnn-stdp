%% figure_weight_histograms_patterns.m
% clc;
% clearvars;
% close all;

if isempty(who('save_all_figures*'))
	save_all_figures_as_pdfs = 0;
	save_all_figures_as_figs = 0;
end

num_sim = 100;
sim_list = 1:100;
sim_list_l = length(sim_list);

r_b_list = 0:10:20;
rate_list = 1:3;
rate_list_l = length(rate_list);

stdp_estimate = [1, 5, 10, 15, 20];
num_stdp_e = length(stdp_estimate);

stdp_label_yf = 0.85;
stdp_label_xf = [0.11, 0.11; 0.11, 0.11; 0.74, 0.71];

hist_bin_centres = 0.02:0.04:0.98;
x_tick_pos = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1];

hist_bin_centres_diff = -1:0.1:1;
x_tick_pos_diff = -1:0.5:1;

y_lims = zeros(num_stdp_e, 2, rate_list_l);
y_lims(:,:,1) = [0, 1; 0, 0.4; 0, 0.4; 0, 0.4; 0, 0.4];
y_lims(:,:,2) = [0, 1; 0, 0.4; 0, 0.2; 0, 0.2; 0, 0.2];
y_lims(:,:,3) = [0, 1; 0, 0.4; 0, 0.4; 0, 0.2; 0, 0.2];

y_tick_pos = cell(num_stdp_e, rate_list_l);
y_tick_pos{1,1} = [0, 0.2, 0.4, 0.6, 1];
y_tick_pos{2,1} = [0, 0.1, 0.2, 0.4];
y_tick_pos{3,1} = [0, 0.1, 0.2, 0.4];
y_tick_pos{4,1} = [0, 0.1, 0.2, 0.4];
y_tick_pos{5,1} = [0, 0.1, 0.2, 0.4];
y_tick_pos{1,2} = [0, 0.2, 0.4, 0.6, 1];
y_tick_pos{2,2} = [0, 0.1, 0.2, 0.4];
y_tick_pos{3,2} = [0, 0.1, 0.2];
y_tick_pos{4,2} = [0, 0.1, 0.2];
y_tick_pos{5,2} = [0, 0.1, 0.2];
y_tick_pos{1,3} = [0, 0.2, 0.4, 0.6, 1];
y_tick_pos{2,3} = [0, 0.1, 0.2, 0.4];
y_tick_pos{3,3} = [0, 0.1, 0.2, 0.4];
y_tick_pos{4,3} = [0, 0.1, 0.2];
y_tick_pos{5,3} = [0, 0.1, 0.2];

axes_dims = cell(num_stdp_e,2);
axes_handles = cell(num_stdp_e,2);
fig_h = 16;
fig_w = 12;
bar_color = [0.4, 0.4, 0.4];

FT_FIG = '.fig';
FT_PDF = '.pdf';

dir_expd = 'lif_sim_data_addSTDC';
dir_stdp = 'lif_sim_data_addSTDP';

what_expd = what(dir_expd);
what_stdp = what(dir_stdp);

load(fullfile(dir_expd, what_expd.mat{1}), ...
	'n_pre', 'num_stdp_e', 'stdp_estimate', 'n_post');

n_p_list = 1:1000;
n_p_list_l = length(n_p_list);
n_b_list = 1001:2000;
n_b_list_l = length(n_b_list);

w_IJ_all_expd = Inf(n_pre, n_post, sim_list_l, num_stdp_e, rate_list_l);
w_IJ_all_stdp = Inf(n_pre, sim_list_l, num_stdp_e, rate_list_l);

w_IJ_diff_p = Inf(n_p_list_l, n_post, sim_list_l, num_stdp_e, rate_list_l);
w_IJ_diff_b = Inf(n_b_list_l, n_post, sim_list_l, num_stdp_e, rate_list_l);

f_expd = zeros(sim_list_l,1);
f_stdp = zeros(sim_list_l,1);

for r = 1:rate_list_l
	for k = 1:sim_list_l
		disp(['Rate=' num2str(r_b_list((rate_list(r)))) ...
			'; Simulation: ' num2str(sim_list(k))]);
		load(fullfile(dir_stdp, ...
			what_stdp.mat{sim_list(k) + (rate_list(r)-1)*sim_list_l}), ...
			'w_IJ_record', 't_f_post_list_count');
		load(fullfile(dir_expd, ...
			what_expd.mat{sim_list(k) + (rate_list(r)-1)*sim_list_l}), ...
			'w_IJ_initial', 't_f_post_list_trig_count');
		
		f_stdp(k) = t_f_post_list_count;
		f_expd(k) = t_f_post_list_trig_count;
		
		for se = 1:num_stdp_e
			w_IJ_all_stdp(:,k,se,r) = w_IJ_record(:,stdp_estimate(se));
			w_IJ_all_expd(:,:,k,se,r) = w_IJ_initial{se}(:,1:f_expd(k));
			
			w_IJ_diff_p(:,:,k,se,r) = bsxfun(@minus, ...
				w_IJ_initial{se}(n_p_list,1:f_expd(k)), ...
				w_IJ_record(n_p_list,stdp_estimate(se)));
			w_IJ_diff_b(:,:,k,se,r) = bsxfun(@minus, ...
				w_IJ_initial{se}(n_b_list,1:f_expd(k)), ...
				w_IJ_record(n_b_list,stdp_estimate(se)));
		end
	end
end

%% synapse weight histograms (pattern neurons): columns-f,M, rows-STDP-STDC
LWmean = 1;
LWmed = 1.5;
LWq = 1;

fig_w = 12;
fig_h = 8;
fig_ind = 0;
hist_edges_col_M = [-0.02, 1e-9, 0.04:0.04:0.96, 1-1e-9, 1.02];

fig_hist_rate = gobjects(6,1);
h_axes_col_M = gobjects(4,3,2);

for pb = 1:2
	for r = 1:rate_list_l
		fig_ind = fig_ind + 1;
		fig_hist_rate(fig_ind) = figure;
		set(gcf, 'Color', 'w');
		set(gcf, 'Units', 'centimeters');
		set(gcf, 'Position', [5+fig_ind, 10-fig_ind, fig_w, fig_h]);
		set(gcf, 'PaperPosition', [0, 0, fig_w, fig_h]);
		set(gcf, 'PaperSize', [fig_w, fig_h]);
		h_ind = 0;
		
		% subplots counted in rows
		se_ind = 0;
		for se = [3,5] % selecting f,M = 10 and 20
			se_ind = se_ind+1;
			
			for k = 1:2
				switch k
					case 1
						if pb == 1
							w_IJ_temp = w_IJ_all_stdp(n_p_list, :, se, r);
						else
							w_IJ_temp = w_IJ_all_stdp(n_b_list, :, se, r);
						end
					case 2
						if pb == 1
							w_IJ_temp = w_IJ_all_expd(n_p_list, :, :, se, r);
						else
							w_IJ_temp = w_IJ_all_expd(n_b_list, :, :, se, r);
						end
				end
				w_IJ_temp = sort(w_IJ_temp(:));
				
				h_axes_col_M(se_ind+(k-1)*2,r,pb) = subplot(2,2,se_ind+(k-1)*2);
				histogram(w_IJ_temp(:), hist_edges_col_M, ...
					'Normalization', 'probability', ...
					'FaceColor', bar_color);
				hold on;
				grid on;
				box on;
				
				h_count_total = length(w_IJ_temp(:));
				plot(w_IJ_temp(round(h_count_total/2))*ones(2,1), ...
					[0; 1], '--k', 'LineWidth', LWmed);
				plot(w_IJ_temp(round(h_count_total/4))*ones(2,1), ...
					[0; 1], ':k', 'LineWidth', LWq);
				plot(w_IJ_temp(round(3*h_count_total/4))*ones(2,1), ...
					[0; 1], ':k', 'LineWidth', LWq);
				plot(mean(w_IJ_temp)*ones(2,1), ...
					[0; 1], 'k', 'LineWidth', LWmean);
				
				xlim([-0.05,1.05]);
				
				if pb == 1
					if r == 1
						ylim([0,0.3]);
					elseif r == 2
						ylim([0,0.2]);
					elseif r == 3
						ylim([0,0.3]);
					end
				else
					if r == 1
						ylim([0,1.05]);
					elseif r == 2
						ylim([0,0.45]);
					elseif r == 3
						ylim([0,0.45]);
					end
				end
				switch k
					case 1
						if se_ind == 1
							ylabel({'STDP';'Frequency'});
						end
						title(['f,M=', num2str(stdp_estimate(se))], 'FontWeight', 'normal');
					case 2
						if se_ind == 1
							ylabel({'Constructed';'Frequency'});
						end
				end
				set(gca, 'XTick', 0:0.25:1);
				if k == 1
					set(gca, 'XTickLabel', {});
				else
					xlabel('Synapse Weight');
				end
				if se_ind == 2
					set(gca, 'YTickLabel', {});
				end
				
			end
		end
		
		for h = 1:length(h_axes_col_M)
			if h == 1 || h == 2
				set(h_axes_col_M(h,r,pb), 'Position', ...
					get(h_axes_col_M(h,r,pb),'Position') + [0.015,-0.005,0.05,0.015]);
			else
				set(h_axes_col_M(h,r,pb), 'Position', ...
					get(h_axes_col_M(h,r,pb),'Position') + [0.015,0.038,0.05,0.015]);
			end
		end
		if pb == 1
			fig_file_name = fullfile('figures', ...
				['fig_ps_w_hist_p_rb', num2str(r_b_list(rate_list(r)),'%02d')]);
		else
			fig_file_name = fullfile('figures', ...
				['fig_ps_w_hist_b_rb', num2str(r_b_list(rate_list(r)),'%02d')]);
		end
		if save_all_figures_as_figs
			saveas(fig_hist_rate(fig_ind), [fig_file_name, FT_FIG]);
			disp(['File saved as: ', fig_file_name, FT_FIG]);
		end
		if save_all_figures_as_pdfs
			saveas(fig_hist_rate(fig_ind), [fig_file_name, FT_PDF]);
			disp(['File saved as: ', fig_file_name, FT_PDF]);
		end
	end
end

%% synapse weight difference histograms (pattern neurons)
fig_h = 11.5;
hist_edges_col_M_D = [-1.08:0.08:1.08];

h_axes_col_diff_p = gobjects(3,2);
h_ind = 0;

fig_hist_diff_p = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5+fig_ind, 10-fig_ind, fig_w, fig_h]);
set(gcf, 'PaperPosition', [0, 0, fig_w, fig_h]);
set(gcf, 'PaperSize', [fig_w, fig_h]);

for r = 1:rate_list_l
	
	% subplots counted in rows
	se_ind = 0;
	for se = [3,5] % selecting f,M = 10 and 20
		se_ind = se_ind+1;
		w_IJ_temp = w_IJ_diff_p(:,:,:,se,r);
		w_IJ_temp = sort(w_IJ_temp(:));
		
		h_axes_col_diff_p(r,se_ind) = subplot(3,2,se_ind+(r-1)*2);
		histogram(w_IJ_temp(:), hist_edges_col_M_D, ...
			'Normalization', 'probability', ...
			'FaceColor', bar_color);
		hold on;
		grid on;
		box on;
		
		h_count_total = length(w_IJ_temp(:));
		%bar(h_centres, h_count/h_count_total, 1, 'FaceColor', bar_color);
		plot(w_IJ_temp(round(h_count_total/2))*ones(2,1), ...
			[0; 1], '--k', 'LineWidth', LWmed);
		plot(w_IJ_temp(round(h_count_total/4))*ones(2,1), ...
			[0; 1], ':k', 'LineWidth', LWq);
		plot(w_IJ_temp(round(3*h_count_total/4))*ones(2,1), ...
			[0; 1], ':k', 'LineWidth', LWq);
		plot(mean(w_IJ_temp)*ones(2,1), ...
			[0; 1], 'k', 'LineWidth', LWmean);
		
		if r == 1
			title(['f,M=', num2str(stdp_estimate(se))], 'FontWeight', 'normal');
		end
		
		xlim([-1.1,1.1]);
		ylim([0,0.55]);
		set(gca,'XTick', -1:0.5:1);
		set(gca,'YTick', [0,0.1,0.3,0.5]);
		if r ~= 3
			set(gca,'XTickLabel', {});
		else
			xlabel('Weight Error');
		end
		if se_ind == 1
			ylabel('Frequency');
		end
		if se_ind == 2
			text(1.25,0.17+((3-r)*0.01),['r=' num2str(r_b_list(r)) 'Hz'],...
				'FontSize', 11, 'Rotation', 90);
			set(gca,'YTickLabel', {});
		end
	end
end

for r = 1:rate_list_l
	for k = 1:2
		set(h_axes_col_diff_p(r,k),'Position', ...
			get(h_axes_col_diff_p(r,k),'Position') + [-0.02,-0.0,0.05,0.03]);
		
	end
end

fig_file_name = fullfile('figures','fig_ps_Dw_hist_p');
if save_all_figures_as_pdfs
	saveas(fig_hist_diff_p, [fig_file_name, FT_FIG]);
	disp(['File saved as: ', fig_file_name, FT_FIG]);
end
if save_all_figures_as_figs
	saveas(fig_hist_diff_p, [fig_file_name, FT_PDF]);
	disp(['File saved as: ', fig_file_name, FT_PDF]);
end

%% synapse weight difference histograms (background noise neurons)
fig_h = 8;
h_axes_col_diff_b = gobjects(2,2);
h_ind = 0;

fig_hist_diff_b = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5+fig_ind, 10-fig_ind, fig_w, fig_h]);
set(gcf, 'PaperPosition', [0, 0, fig_w, fig_h]);
set(gcf, 'PaperSize', [fig_w, fig_h]);

for r = 2:rate_list_l
	
	% subplots counted in rows
	se_ind = 0;
	for se = [3,5] % selecting f,M = 10 and 20
		se_ind = se_ind+1;
		w_IJ_temp = w_IJ_diff_b(:,:,:,se,r);
		w_IJ_temp = sort(w_IJ_temp(:));
		
		h_axes_col_diff_b(r-1,se_ind) = subplot(2,2,se_ind+(r-2)*2);
		histogram(w_IJ_temp(:), hist_edges_col_M_D, ...
			'Normalization', 'probability', ...
			'FaceColor', bar_color);
		hold on;
		grid on;
		box on;
		
		h_count_total = length(w_IJ_temp(:));
		plot(w_IJ_temp(round(h_count_total/2))*ones(2,1), ...
			[0; 1], '--k', 'LineWidth', LWmed);
		plot(w_IJ_temp(round(h_count_total/4))*ones(2,1), ...
			[0; 1], ':k', 'LineWidth', LWq);
		plot(w_IJ_temp(round(3*h_count_total/4))*ones(2,1), ...
			[0; 1], ':k', 'LineWidth', LWq);
		plot(mean(w_IJ_temp)*ones(2,1), ...
			[0; 1], 'k', 'LineWidth', LWmean);
		
		if r == 2
			title(['f,M=', num2str(stdp_estimate(se))], 'FontWeight', 'normal');
		end
		
		xlim([-1.1,1.1]);
		ylim([0,0.4]);
		set(gca,'XTick', -1:0.5:1);
		set(gca,'YTick', [0,0.2,0.4]);
		if r ~= 3
			set(gca,'XTickLabel', {});
		else
			xlabel('Weight Error');
		end
		if se_ind == 1
			ylabel('Frequency');
			
		end
		if se_ind == 2
			text(1.25,0.12+((3-r)*0.01),['r=' num2str(r_b_list(r)) 'Hz'],...
				'FontSize', 11, 'Rotation', 90);
			set(gca,'YTickLabel', {});
		end
	end
end

for r = 2:rate_list_l
	for k = 1:2
		if r == 2
			set(h_axes_col_diff_b(r-1,k), 'Position', ...
				get(h_axes_col_diff_b(r-1,k),'Position') + [-0.02,-0.005,0.05,0.015]);
		else
			set(h_axes_col_diff_b(r-1,k), 'Position', ...
				get(h_axes_col_diff_b(r-1,k),'Position') + [-0.02,0.038,0.05,0.015]);
		end
		
	end
end

% save figures
fig_file_name = fullfile('figures','fig_ps_Dw_hist_b');
if save_all_figures_as_pdfs
	saveas(fig_hist_diff_b, [fig_file_name, FT_PDF]);
	disp(['File saved as: ', fig_file_name, FT_PDF]);
end
if save_all_figures_as_figs
	saveas(fig_hist_diff_b, [fig_file_name, FT_FIG]);
	disp(['File saved as: ', fig_file_name, FT_FIG]);
end
