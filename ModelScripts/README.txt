ReadMe file for script (Jusyte et al.)

About the script:

 - The model and some of the following instructions are based on the 2020, Kobbersmed at al. eLife paper "Rapid regulation of vesicle priming explains synaptic facilitation despite heterogeneous vesicle:Ca2+ channel distances". The paper and the original code can be found here: https://elifesciences.org/articles/51032/figures#content. In the original code different variations of the model are explored. For this paper model type 47 was updated with a new unpriming reaction. Since some of the old parameters were updated and new parameters were added the other model types may not work anymore.

 - All simulations were carried out in Matlab (R2020a) on a computer allowing parallel computing.

Prior to simulations:
 
 - Before simulations the calcium calculator CalC needs to be installed in the same folder as the scripts (the outermost folder). The program as well as a manual on installation and usage can be found here: https://web.njit.edu/~matveev/calc.html. After installation, remember to set permissions to allow execution of CalC. Also make sure that the path to the CalC program is correct in the bottom of the RunCalC_det.m script. Otherwise, no calcium files will be created, and the Matlab script will return an error as it cannot open the calcium file. CalC version 6.8.6 for Windows was used in the simulations presented in the paper. This version of the program (calc_686.exe, with plug-ins: cygcc_s-she.dll, cygstdc++-6.dll, cygwin1.dll) can be in the folder containing the scripts and can be used without further installation. 
 
 - As the scripts are sorted in folders (Common, Deterministic, Stochastic, etc.). Remember to add all folders to the search path.
 
 - Lines 59-61 in testing_the_system.m are responsible for deleting old CalC files when running new ones for the current simulation. In line 59 PCWIN64 must be changed to the name of the computer (this can be acquired by simply typing “computer” to the Command Window). Lines 60-61 may need other syntax on Mac OS and Linux systems.

Calling a single simulation:

 - The script test_run.m in the main folder defines all relevant parameters for (both a deterministic and a stochastic) simulation of the unpriming model presented in the paper with the best fitted parameters. As default it runs 200 repetitions of the stochastic simulation. This can be adjusted in line 10 of model_parameters_det.m. Function run_more_reps.m can also be used for running more repetitions.

Other relevant files and functions:

 - parameter_choices.m: Given a model choice and initial parameters, this script provides an input vector of parameters to be fed into model_parameters_DET. Also provides an initial result file name (modified in other scripts).  Lines 419-426 of parameter_choices.m contain the parameters of the model that are set from test_tun.m as default. KD is raised to the power of n*m (K_D^(n*m)) in the parameter initialization and it is passed down in that form for further use (prim_kM). m is multiplied by n and it is passed down as Ca_prim_type. Unprim_mutant defines the type of unpriming (0: Ca2+ dependent unpriming rate defined in equation 18 of the article, 1: Constant low Ca2+ independent unpriming rate defined in equation 19 of the article, 2: Constant high unpriming rate).

 - calculate_Caprim_rate.m: Function calculating the priming and unpriming rates.

 - fit_model4_mutant_and_wt.m, cost_model4_mutant_and_wt.m and determine_cost_bothpeaks_numves_factor_mutant_and_wt.m: functions used for fitting the model parameters to the experimental data.

 - new_experiments_mutant_and_wt_peaks12_ppr.mat and full_exp_data.mat: files that contain the experimental data used for the fitting and generating the figures.

 - figures.m: script used for generating figures 2 and s1.

 - new_result: folder containing the simulation results used for generating the figures.
