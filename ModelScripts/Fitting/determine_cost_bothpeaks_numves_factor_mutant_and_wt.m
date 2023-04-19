function [cost_value_bothpeaks, num_ves_factor_bothpeaks] = determine_cost_bothpeaks_numves_factor_mutant_and_wt(peak1s_wt, peak2s_extr_wt, peak1s_mutant, peak2s_extr_mutant)

%loading experimental data:
load('new_experiments_mutant_and_wt_peaks12_ppr.mat')

%defining cost function with finding the optimal number of vesicles:
cost_function = @(x)(sum(((peak1s_wt*x - peak1_mean_all_cells_exp_wt).^2)./abs(peak1_mean_all_cells_exp_wt))...
                    +sum(((peak2s_extr_wt*x - peak2_extr_mean_all_cells_exp_wt).^2)./abs(peak2_extr_mean_all_cells_exp_wt))...
                    +sum(((peak1s_mutant*x - peak1_mean_all_cells_exp_mutant).^2)./abs(peak1_mean_all_cells_exp_mutant))...
                    +sum(((peak2s_extr_mutant*x - peak2_extr_mean_all_cells_exp_mutant).^2)./abs(peak2_extr_mean_all_cells_exp_mutant)));

%running cost minimalization function:
[num_ves_factor_bothpeaks, cost_value_bothpeaks] = fminsearch(cost_function, 2);