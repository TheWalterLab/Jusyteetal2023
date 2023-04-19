function test_run()

%%THIS SCRIPTS RUNS THE MODEL AT EXPERIMENTAL CALCIUM CONCENTRATIONS WITH
%%THE BEST FIT PARAMETERS OF THE UNPRIMING MODEL

model_type = 47; %The choice of model (and which parameters to define) as defined and explained in parameter_choices.m

Q_max = 14.047;
m = 2.0431;
Ca_prim_type = 2*m;
prim_rate_const = 92.1614;
unprim_rate_const = 265.0334;
kD = (4.0757e-08);
kD_par_free = kD^Ca_prim_type;
CaMax_rest = 190e-9;
num_ves_factor = 1.3356;
unprim_mutant = 0; %0 for Ca dependent unpriming, 1 for constant low unpriming, 2 for constant high unpriming

par_free = [Q_max Ca_prim_type kD_par_free prim_rate_const unprim_rate_const CaMax_rest num_ves_factor unprim_mutant];

stoch_on_off = 1; %0 for deterministic simulation, 1 for stochastic, 2 for both
rand_ves_on_off = 1; %This defines the vesicle placement as described in determ_vesicle_distances.m
CalC_on_off = 1; %0: No CalC simulation (if Ca files are already generated), >0 for CalC simulation. If == 1 all older calcium files are deleted
CaExtracellular = [0.4 0.75 1.5 3 6]; %Extracellular calcium concentrations (mM) to simulate
save_data = 1; %Defines how much data to be saved (1: A lot, including all eEJCs. 2: Data summary). See end of simulation_call_det.m and simulation_call_stoch.m
save_calc_loc = 1; %%Defines the location of calcium files (to be able to run more procedures in parallel without interfering with each other). If ==1, Calc simulation is 20.5 ms, if ==2 163.5 ms
pVr2_hack = 0; %Set this to zero. Only used for estimation of pVr2 (by putting new vesicles into the system right before second stimulation)

[par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off); %Defines the parameter vector and result filename to be inputted in the following function call. 
% savefilename = ['results_' num2str(unprim_mutant)]; %Sometimes Matlab cannot write to a file if the name of the file is too long. Uncomment this line to resolve the issue.
testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack)
