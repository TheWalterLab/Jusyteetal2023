function [cost_value, num_ves_factor, peak1s_wt_corrected, peak2s_extr_wt_corrected, pprs_extr_wt, peak1s_mutant_corrected, peak2s_extr_mutant_corrected, pprs_extr_mutant] = cost_model4_mutant_and_wt(pars, prim_coop, save_data)

%%%This script determines the cost value of unpriming model.
%%%The num_ves is corrected for each choice of [Q_max, kD_wt, prim_rate, unprim_rate, m] values to minimize the
%%%cost function

%%NOTE: This optimizes the perimeter models (all SVs at equal distance).

Q_max = pars(1);
m = pars(2);
Ca_prim_type = prim_coop*m;
prim_rate_const = pars(3);
unprim_rate_const = pars(4);
kD_wt_par = pars(5); %for log file
kD_wt = ((pars(5))^2)^m;

kD_mutant = kD_wt;
kD_mutant_par = pars(5); %for log file
CaMax_rest = 190e-9;

pVr2_hack = 0;

CaExtracellular = [0.4 0.75 1.5 3 6];
save_calc_loc = 1;

%displaying parameter values for current iteration
disp('Running next simulation')
disp(['kD_wt (fitting)= (' num2str(kD_wt_par) ')^2'])
disp(['Q_max (fitting)= ' num2str(Q_max)])
disp(['m (fitting)= ' num2str(m)])
disp(['prim_rate_const (fitting)= ' num2str(prim_rate_const)])
disp(['unprim_rate_const (fitting)= ' num2str(unprim_rate_const)])
disp('Number of vesicles will be corrected later')

tic;

%Running determenistic simulations:
stoch_on_off = 0;
rand_ves_on_off = 1;
CalC_on_off = 1;

if sum(pars <= 0)
    cost_value = 1e6;
    num_ves_factor = NaN;
    peak1s = NaN;
    peak1s_corrected = NaN;

else    
    num_stim = 2;
    stim_freq = 100;
    par_free_wt = [Q_max Ca_prim_type kD_wt prim_rate_const unprim_rate_const CaMax_rest 1 0];
    par_free_mutant = [Q_max Ca_prim_type kD_mutant prim_rate_const unprim_rate_const CaMax_rest 1 1];
    model_type = 47;
    [par_init_wt, savefilename_wt] = parameter_choices(par_free_wt, model_type, 0, rand_ves_on_off);    
    [~, ~, ~, peak1s_wt, peak2s_extr_wt, pprs_extr_wt, ~, pVr1s_wt, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init_wt, CaExtracellular, save_data, savefilename_wt, save_calc_loc, pVr2_hack, stim_freq, num_stim); 
    [par_init_mutant, savefilename_mutant] = parameter_choices(par_free_mutant, model_type, 0, rand_ves_on_off);
    CalC_on_off = 0;
    [~, ~, ~, peak1s_mutant, peak2s_extr_mutant, pprs_extr_mutant, ~, pVr1s_mutant, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init_mutant, CaExtracellular, save_data, savefilename_mutant, save_calc_loc, pVr2_hack, stim_freq, num_stim); 
    
    %Calculating cost value:
    [cost_value, num_ves_factor] = determine_cost_bothpeaks_numves_factor_mutant_and_wt(peak1s_wt, peak2s_extr_wt, peak1s_mutant, peak2s_extr_mutant);

    peak1s_wt_corrected = peak1s_wt*num_ves_factor;
    peak2s_extr_wt_corrected = peak2s_extr_wt*num_ves_factor;
    peak1s_mutant_corrected = peak1s_mutant*num_ves_factor;
    peak2s_extr_mutant_corrected = peak2s_extr_mutant*num_ves_factor;
end

num_ves_factor
cost_value

%Write log file for each iteration

    filecont = sprintf(     ['kD_wt: (' num2str(kD_wt_par) ')^2' '\n' ...
                             'Q_max: ' num2str(Q_max) '\n' ...
                             'm: ' num2str(m) '\n' ...
                             'prim_rate_const: ' num2str(prim_rate_const) '\n' ...
                             'unprim_rate_const: ' num2str(unprim_rate_const) '\n' ...
                             'Num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                             'Cost: ' num2str(cost_value) '\n' ...
                            '------------------- \n', ...
                            'PVR_wt: ' num2str(pVr1s_wt) '\n' ...
                            'peak1s_wt:' num2str(peak1s_wt_corrected) '\n', ...
                            'peak2s_extr_wt:' num2str(peak2s_extr_wt_corrected) '\n', ...
                            'PPRs_extr_wt: ' num2str(pprs_extr_wt) '\n'   ...
                            'PVR_mutant: ' num2str(pVr1s_mutant) '\n' ...
                            'peak1s_mutant:' num2str(peak1s_mutant_corrected) '\n', ...
                            'peak2s_extr_mutant' num2str(peak2s_extr_mutant_corrected) '\n', ...
                            'PPRs_extr_mutant: ' num2str(pprs_extr_mutant) '\n'   ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

    fullfile = ['./Sim_data/Log_files/Log_Optimization_model4_mutant_and_wt_primcoop' num2str(prim_coop)];

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');



toc