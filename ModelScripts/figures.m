%% This script generates plots that were used for Figure 2 and Figure S1 in the Jusyte et al., 2023 Cell Reports paper   

%% Load files
% 
clear all
close all

load('.\Exp_data\full_exp_data.mat');

%load simulation results
wt = load('.\Sim_data\new_results\*.mat');
mutant = load('.\Sim_data\new_results\*.mat');

save_folder =  '.\NewFigures\';
CaExtracellular = [0.4 0.75 1.5 3 6];

%% plot EPSC1 peaks (Figure 2 G)

figure()

title('Peak EPSC_1 simulation')
hold on

epsc1_wt = errorbar(CaExtracellular,wt.peak1_means,std(wt.peak1s_all,0,2)./sqrt(200));
epsc1_mutant = errorbar(CaExtracellular,mutant.peak1_means,std(mutant.peak1s_all,0,2)./sqrt(200));

epsc1_wt.Marker = '.';
epsc1_mutant.Marker = '.';
epsc1_wt.MarkerSize = 15;
epsc1_mutant.MarkerSize = 15;
epsc1_wt.MarkerFaceColor = "auto";
epsc1_mutant.MarkerFaceColor = "auto";
epsc1_wt.LineWidth = 1;
epsc1_mutant.LineWidth = 1;

epsc1_wt_exp = errorbar(CaExtracellular,-wt_exp.epsc1mean,wt_exp.sem_epsc1);
epsc1_mutant_exp = errorbar(CaExtracellular,-mutant_exp.epsc1mean,mutant_exp.sem_epsc1);

epsc1_wt_exp.Marker = '.';
epsc1_mutant_exp.Marker = '.';
epsc1_wt_exp.MarkerSize = 15;
epsc1_mutant_exp.MarkerSize = 15;
epsc1_wt_exp.MarkerFaceColor = "auto";
epsc1_mutant_exp.MarkerFaceColor = "auto";
epsc1_wt_exp.LineWidth = 1;
epsc1_mutant_exp.LineWidth = 1;

xticks(CaExtracellular)
xticklabels(CaExtracellular)
xlim([0 6.2])
ylim([-150 0])
ylabel('Peak EPSC_1 (nA)')
set(gca,'ydir','reverse')
xlabel('[Ca^{2+}]{ext} (mM)')
legend('wt sim', 'mutant sim', 'wt exp', 'mutant exp' ,'Location','best')
set(gca, 'TickDir', 'out')
box off
hold off

saveas(gcf,[save_folder ,date ,'_epsc1peaks'],'epsc')

%% plot PPR (Figure 2 H)

figure()

title('PPR simulation')
hold on

ppr_wt = errorbar(CaExtracellular,wt.ppr_extr_means,std(wt.ppr_extrs_all,0,2)./sqrt(200));
ppr_mutant = errorbar(CaExtracellular,mutant.ppr_extr_means,std(mutant.ppr_extrs_all,0,2)./sqrt(200));

ppr_wt.Marker = '.';
ppr_mutant.Marker = '.';
ppr_wt.MarkerSize = 15;
ppr_mutant.MarkerSize = 15;
ppr_wt.MarkerFaceColor = "auto";
ppr_mutant.MarkerFaceColor = "auto";
ppr_wt.LineWidth = 1;
ppr_mutant.LineWidth = 1;

PPR_wt_exp = errorbar(CaExtracellular,wt_exp.pprmean,wt_exp.sem_ppr);
PPR_mutant_exp = errorbar(CaExtracellular,mutant_exp.pprmean,mutant_exp.sem_ppr);

PPR_wt_exp.Marker = '.';
PPR_mutant_exp.Marker = '.';
PPR_wt_exp.MarkerSize = 15;
PPR_mutant_exp.MarkerSize = 15;
PPR_wt_exp.MarkerFaceColor = "auto";
PPR_mutant_exp.MarkerFaceColor = "auto";
PPR_wt_exp.LineWidth = 1;
PPR_mutant_exp.LineWidth = 1;

xticks(CaExtracellular)
xticklabels(CaExtracellular)
xlim([0 6.2])
ylim([0 2.5])
ylabel('PPR')
xlabel('[Ca^{2+}]_{ext} (mM)')
legend('wt sim','mutant sim','wt exp','mutant exp','Location','best')
set(gca, 'TickDir', 'out')
box off
hold off

saveas(gcf,[save_folder ,date ,'ppr'],'epsc')

%% Calcium transient (Figure 2 E)
% Plot the Calcium simulated via the CalC program over the simulation time

figure();
hold on

%plot(wt.Ca_time_equi_cell {1, 1}(1:100:end) , (wt.Ca_R_equi_cell{1,1}(1:100:end,25)))
plot(wt.Ca_time_equi_cell {5, 1}(1:100:end) , (wt.Ca_R_equi_cell{5, 1}(1:100:end,25)))

%legend('low 0.4','high 6','Location','best')
ylabel( '[Ca^{2+}]_i (microM)')
ylim([0.01 200])
set(gca,'YTick',[0.01 0.1 1 10 100])
xlabel('Time(ms)')
xlim([0 20])
xticks([0 5 10 15 20]);
xticklabels({'0','5','10','15', '20'})

title('high')

set(gca, 'YScale', 'log')
set(gca, 'TickDir', 'out')
hold off

saveas(gcf,[save_folder ,date ,'calcium_transient_high'],'epsc')


%% Low primed Vesicles (Figure 2 E1)

% Calculate the number of primed vesicles. All states (0 to 5) are considered
% State 100 correspond to fusion
% The time considered is the time of the simulation 

figure()

for i = 1: 200
    states = wt.ves_states_cell{1,i};
    t = wt.stoch_time_cell{1,i};
    primed = sum(states~=100);
    wt.primed_low(:,i) = interp1(t, primed, wt.time_vector);
end

for i = 1: 200
    states = mutant.ves_states_cell{1,i};
    t = mutant.stoch_time_cell{1,i};
    primed = sum(states~=100);
    mutant.primed_low(:,i) = interp1(t, primed, mutant.time_vector);
end 

%
hold on

plot(wt.time_vector(1:100:end), mean(wt.primed_low(1:100:end,:),2))
plot(mutant.time_vector(1:100:end), mean(mutant.primed_low(1:100:end,:),2))

xlabel('Time(s)')
xlim([0 0.020])
xticks([0 0.005 0.010 0.015 0.020])
xticklabels({'0','0.005','0.010','0.015', '0.020'})

ylim([0 300])
yticks([0 50 100 150 200 250 300])
yticklabels({'0','50','100','150', '200', '250', '300'})
ylabel('Number of docked vesicles')

legend('Wt','on', 'off' ,'Location','best')
title('low')
% 
set(gca, 'TickDir', 'out')
hold off

saveas(gcf,[save_folder ,date ,'primed_vesicles_low04'],'epsc')
saveas(gcf,[save_folder ,date ,'primed_vesicles_low04'],'fig')

%% High primed Vesicles (Figure 2 E2)
% Calculate the number of primed vesicles. All states (0 to 5) are considered
% State 100 correspond to fusion
% The time considered is the time of the simulation 

figure()

for i = 1: 200
    states = wt.ves_states_cell{5,i};
    t = wt.stoch_time_cell{5,i};
    primed = sum(states~=100);
    wt.primed_high(:,i) = interp1(t, primed, wt.time_vector);
end

for i = 1: 200
    states = mutant.ves_states_cell{5,i};
    t = mutant.stoch_time_cell{5,i};
    primed = sum(states~=100);
    mutant.primed_high(:,i) = interp1(t, primed, mutant.time_vector);
end 


hold on
plot(wt.time_vector(1:100:end), mean(wt.primed_high(1:100:end,:),2))
plot(mutant.time_vector(1:100:end), mean(mutant.primed_high(1:100:end,:),2))

xlabel('Time(s)')
xlim([0 0.020])
xticks([0 0.005 0.010 0.015 0.020])
xticklabels({'0','0.005','0.010','0.015', '0.020'})

ylim([0 300])
yticks([0 50 100 150 200 250 300])
yticklabels({'0','50','100','150', '200', '250', '300'})
ylabel('Number of docked vesicles')

legend('Wt','on', 'off' ,'Location','best')
title('high')
% 
set(gca, 'TickDir', 'out')
hold off
xlim([0 0.020])
saveas(gcf,[save_folder ,date ,'primed_vesicles_high6'],'epsc')

%% eEPSCs
% Plot of the stochastic eEPSC over simulation time for the 5 external
% calcium concentration
% !!! Warning !!! Important to load again simulation results 
clear all
%load simulation results

wt = load('.\Sim_data\new_results\*.mat');
mutant = load('.\Sim_data\new_results\*.mat');

save_folder =  '.\Figures\';
CaExtracellular = [0.4 0.75 1.5 3 6];

%% (Figure 2 F)
figure ();
title('EPSC')
%eEPSC 0.4 mM 
subplot(1,5,1)
hold on;
plot(wt.time_vector(1:100:end), wt.stoch_EPSC_means(1,1:100:end), 'LineWidth', 3 )
plot(mutant.time_vector(1:100:end), mutant.stoch_EPSC_means(1,1:100:end), 'LineWidth', 2)
%plot(off.time_vector(1:100:end), off.stoch_EPSC_means(1,1:100:end), 'LineWidth', 2 )
line([ 0.005 0.015],[-100 -100],'LineWidth', 3, 'color', 'k' )
line([ 0.005 0.005],[-80 -100],'LineWidth', 3, 'color', 'k' )
set(gca,'XColor', 'none','YColor','none')
ylabel('20 nA','Color','k')
xlabel('10 ms', 'Color','k')
ylim([-150 0])
title('0.4 mM')
%legend('WT','Mutant','Location','best')
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca, 'TickDir', 'out')
box off
hold off

% EPSC 0.75 mM 
subplot(1,5,2)
%figure()
hold on
plot(wt.time_vector(1:100:end), wt.stoch_EPSC_means(2,1:100:end), 'LineWidth', 3 )
plot(mutant.time_vector(1:100:end), mutant.stoch_EPSC_means(2,1:100:end), 'LineWidth', 2 )
%plot(off.time_vector(1:100:end), off.stoch_EPSC_means(2,1:100:end), 'LineWidth', 2 )
ylim([-150 0])
title('0.75 mM')
legend('WT','on','off','Location','best')
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca, 'TickDir', 'out')
set(gca,'XColor', 'none','YColor','none')
box off
hold off


% EPSC 1.5 mM
subplot(1,5,3)
%figure()
hold on
plot(wt.time_vector(1:100:end), wt.stoch_EPSC_means(3,1:100:end), 'LineWidth', 3 )
plot(mutant.time_vector(1:100:end), mutant.stoch_EPSC_means(3,1:100:end), 'LineWidth', 2 )
%plot(off.time_vector(1:100:end), off.stoch_EPSC_means(3,1:100:end), 'LineWidth', 2 )
ylim([-150 0])
title('1.5 mM')
%legend('WT','Mutant','Location','best')
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca, 'TickDir', 'out')
set(gca,'XColor', 'none','YColor','none')
box off
hold off


% EPSC 3 mM 
subplot(1,5,4)
%figure()
hold on
plot(wt.time_vector(1:100:end), wt.stoch_EPSC_means(4,1:100:end), 'LineWidth', 3 )
plot(mutant.time_vector(1:100:end), mutant.stoch_EPSC_means(4,1:100:end), 'LineWidth', 2 )
%plot(off.time_vector(1:100:end), off.stoch_EPSC_means(4,1:100:end), 'LineWidth', 2 )
ylim([-150 0])
title('3 mM')
%legend('WT','Mutant','Location','best')
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca, 'TickDir', 'out')
set(gca,'XColor', 'none','YColor','none')
box off
hold off

% EPSC 6 mM 
subplot(1,5,5)
%figure()
hold on
plot(wt.time_vector(1:100:end), wt.stoch_EPSC_means(5,1:100:end), 'LineWidth', 3 )
plot(mutant.time_vector(1:100:end), mutant.stoch_EPSC_means(5,1:100:end), 'LineWidth', 2 )
%plot(off.time_vector(1:100:end), off.stoch_EPSC_means(5,1:100:end), 'LineWidth', 2 )
ylim([-150 0])
title('6 mM')
%legend('WT','Mutant','Location','best')
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca, 'TickDir', 'out')
set(gca,'XColor', 'none','YColor','none')
box off
hold off

saveas(gcf,[save_folder ,date ,'epscall'],'epsc')

%%
% %% 
% % wt
% for n = 1:200
%     fused = wt.fused_ves_cell{5,n};
%     time_stoch = wt.stoch_time_cell{5, n}; 
%     wt.time_vector = [0:5e-4:0.025];
%     for i = 1:(length(wt.time_vector)-1)
%          v_fused_wt(i,n) = sum(fused(time_stoch > wt.time_vector(i) & time_stoch<wt.time_vector(i+1)));
%          
%     end
% end
% 
% % on
% for n = 1:200
%     fused = mutant.fused_ves_cell{5,n};
%     time_stoch = mutant.stoch_time_cell{5, n}; 
%     mutant.time_vector = [0:5e-4:0.025];
%     for i = 1:length(mutant.time_vector)-1
%          v_fused_mutant(i,n) = sum(fused(time_stoch > mutant.time_vector(i) & time_stoch<mutant.time_vector(i+1)));
%          
%     end
% end

%%
% 
% cum_v_fused_wt = cumsum(v_fused_wt);
% cum_v_fused_mutant = cumsum(v_fused_mutant);
% 
% %--plots
% hold on 
% 
% plot(wt.time_vector(1:end-1), mean(cum_v_fused_wt,2))
% plot(mutant.time_vector(1:end-1), mean(cum_v_fused_mutant,2))
% 
% xlim([0 0.022])
% xticks([0 0.005 0.010 0.015 0.020])
% xticklabels({'0','0.005','0.010','0.015', '0.020'})
% xlabel('Time(s)');
% 
% ylim([0 450])
% yticks([0 100 200 300 400])
% yticklabels({'0', '100', '200', '300', '400'})
% ylabel('Fusion rate high')
% 
% legend('WT','On','Off','Location','best')
% 
% set(gca, 'TickDir', 'out')
% hold off
% 
% 
% saveas(gcf,[save_folder ,date ,'fusion_rate_high6'],'epsc')

%% calculating fusions
for j = 1:length(CaExtracellular)
    for n = 1:200
        fused = wt.fused_ves_cell{j,n};
        time_stoch = wt.stoch_time_cell{j, n}; 
        fusion10ms(n)=sum(fused(time_stoch < 0.01));
    end

    fusedat10wt(j)=mean(fusion10ms);
    for n = 1:200
        fused = mutant.fused_ves_cell{j,n};
        time_stoch = mutant.stoch_time_cell{j, n}; 
        fusion10ms(n)=sum(fused(time_stoch < 0.01));
    end
    fusedat10on(j)=mean(fusion10ms);
end


%% loading simulation parameters
parstruct = load('.\Sim_data\par.mat');
par = parstruct.par;
m = 2.0431;
Ca_prim_type = 2*m;
prim_kM = (40.0757e-9)^Ca_prim_type;
prim_rate_const = 92.1614; % par(32)
unprim_rate_const = 265.0334; % par(33)
unprim_rate_const_0 = 0; % par(44)
CaExt_extr = [0:0.05:6] * 1e-3;
CaIntra = [0:0.5:100]*10e-09;
Km = 2.679 * 1e-3; % Km,current ; Table 1. 
Ca_max = 190e-9;
save_folder =  '.\NewFigures\';
num_ves = round(1.3356*180);

%%
% Calculate internal calcium concentration from external calcium
% concentration
for i = (1:length(CaExt_extr))
   Ca = Ca_max * (CaExt_extr(i) / (Km + CaExt_extr(i)));
   Ca_bas(i) = Ca;
end

% Calculate priming rate and unpriming rate for the calculated internal calcium concentration 
for i = (1:length(Ca_bas))
    [Caprim_rate_ext(:,i), Caunprim_rate_ext(:,i)] = calculate_Caprim_rate(Ca_bas(i), Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, 0);
end

% Calculate Steady State for the calculated internal calcium concentrations
for i = (1:length(Caunprim_rate_ext))
     [SS_normal_ext]= calculate_steady_state(par, Ca_bas(i), 1);
      SS_normal_ext(end+1) = SS_normal_ext(1)*(Caunprim_rate_ext(i)/Caprim_rate_ext(i));
      SS_total_ext{i} = SS_normal_ext/sum(SS_normal_ext);
end

% Calculate % occupied release site for the calculated internal calcium concentrations
for i = (1:length(SS_total_ext))
   percent_wt_ext(:,i) = (sum(SS_total_ext{1,i})-SS_total_ext{1, i}(end))*100;
end


% Calculate priming rate and unpriming rate for the calculated internal calcium concentration 
for i = (1:length(Ca_bas))
    [Caprim_rate_ext_on(:,i), Caunprim_rate_ext_on(:,i)] = calculate_Caprim_rate(Ca_bas(i), Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, 1);
end

% Calculate Steady State for the calculated internal calcium concentrations
for i = (1:length(Caunprim_rate_ext_on))
     [SS_normal_ext_on]= calculate_steady_state(par, Ca_bas(i), 1);
      SS_normal_ext_on(end+1) = SS_normal_ext_on(1)*(Caunprim_rate_ext_on(i)/Caprim_rate_ext_on(i));
      SS_total_ext_on{i} = SS_normal_ext_on/sum(SS_normal_ext_on);
end

% Calculate % occupied release site for the calculated internal calcium concentrations
for i = (1:length(SS_total_ext_on))
   percent_ext_on(:,i) = (sum(SS_total_ext_on{1,i})-SS_total_ext_on{1, i}(end))*100;
end


% Calculate priming rate and unpriming rate for the calculated internal calcium concentration 
for i = (1:length(Ca_bas))
    [Caprim_rate_ext_off(:,i), Caunprim_rate_ext_off(:,i)] = calculate_Caprim_rate(Ca_bas(i), Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, 2);
end

% Calculate Steady State for the calculated internal calcium concentrations
for i = (1:length(Caunprim_rate_ext_off))
     [SS_normal_ext_off]= calculate_steady_state(par, Ca_bas(i), 1);
      SS_normal_ext_off(end+1) = SS_normal_ext_off(1)*(Caunprim_rate_ext_off(i)/Caprim_rate_ext_off(i));
      SS_total_ext_off{i} = SS_normal_ext_off/sum(SS_normal_ext_off);
end

% Calculate % occupied release site for the calculated internal calcium concentrations
for i = (1:length(SS_total_ext_off))
   percent_ext_off(:,i) = (sum(SS_total_ext_off{1,i})-SS_total_ext_off{1, i}(end))*100;
end


%% plotting %occupied release sites (steady state) (Figure 2 C)

figure()
hold on
plot(CaExt_extr, percent_wt_ext)
plot(CaExt_extr, percent_ext_on)
ylim([0 100]);
xlim([0 6e-3])
ylabel('%occupied (steady state)')
xlabel('external calcium concentration')
set(gca, 'TickDir', 'out')
legend('WT','On','Location','best')
hold off

saveas(gcf,[save_folder ,date ,'_PercOccupied_CaExt'],'epsc')

%% plotting release probabilities

caext=[0.4 0.75 1.5 3 6]*1e-3;
for i=1:5
    a=CaExt_extr==caext(i);
    primatrestwt(i)=percent_wt_ext(find(a,1,"first"))*num_ves/100;
    primatreston(i)=percent_ext_on(find(a,1,"first"))*num_ves/100;
end

releaseprobwt1=fusedat10wt./primatrestwt;
releaseprobwt2=fusedat10wt./num_ves;
releaseprobon1=fusedat10on./primatreston;
releaseprobon2=fusedat10on./num_ves;

%(Figure S1 C)

figure()
hold on
plot(caext, releaseprobwt1,'-o');
plot(caext, releaseprobon1, '-o');
legend('WT','On','Location','best')
set(gca, 'TickDir', 'out')

xlim([0 6e-3])
ylim([0 1])
ylabel('P - /prim at rest' )
xlabel('[Ca^{2+}]_{ext} (mM)')
legend('wt','on','Location','best')
set(gca, 'TickDir', 'out')
box off

hold off

saveas(gcf,[save_folder ,date ,'_release_prob_prime_at_rest'],'epsc')

%(Figure 2 D)

figure()
hold on
plot(caext, releaseprobwt2, '-o');
plot(caext, releaseprobon2,'-o');
legend('WT','On','Location','best')
set(gca, 'TickDir', 'out')

xlim([0 6e-3])
ylim([0 1])
ylabel('P - nsites' )
xlabel('[Ca^{2+}]_{ext} (mM)')
legend('wt','on','Location','best')
set(gca, 'TickDir', 'out')
box off

hold off
saveas(gcf,[save_folder ,date ,'_release_prob_prime_numves'],'epsc')


%%
% Calculate priming rate and unpriming rate for the linear internal calcium concentration 
for i = (1:length(CaIntra))
    [Caprim_rate(:,i), Caunprim_rate(:,i)] = calculate_Caprim_rate(CaIntra(i), Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, 0);
end

% Calculate Steady State for the linear internal calcium concentrations
for i = (1:length(Caunprim_rate))
     [SS_normal]= calculate_steady_state(par, CaIntra(i), 1);
      SS_normal(end+1) = SS_normal(1)*(Caunprim_rate(i)/Caprim_rate(i));
      SS_total{i} = SS_normal/sum(SS_normal);
end
        
% Calculate % occupied release site for the linear internal calcium concentrations
for i = (1:length(SS_total))
   percent_wt(:,i) = (sum(SS_total{1,i})-SS_total{1, i}(end))*100;
end

% Calculate priming rate and unpriming rate for the calculated internal calcium concentration 
for i = (1:length(CaIntra))
    [Caprim_rate_on(:,i), Caunprim_rate_on(:,i)] = calculate_Caprim_rate(CaIntra(i), Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, 1);
end

% Calculate Steady State for the calculated internal calcium concentrations
for i = (1:length(Caunprim_rate_on))
     [SS_normal_on]= calculate_steady_state(par, CaIntra(i), 1);
      SS_normal_on(end+1) = SS_normal_on(1)*(Caunprim_rate_on(i)/Caprim_rate_on(i));
      SS_total_on{i} = SS_normal_on/sum(SS_normal_on);
end

% Calculate % occupied release site for the calculated internal calcium concentrations
for i = (1:length(SS_total_on))
   percent_on(:,i) = (sum(SS_total_on{1,i})-SS_total_on{1, i}(end))*100;
end


% Calculate priming rate and unpriming rate for the calculated internal calcium concentration 
for i = (1:length(CaIntra))
    [Caprim_rate_off(:,i), Caunprim_rate_off(:,i)] = calculate_Caprim_rate(CaIntra(i), Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, 2);
end

% Calculate Steady State for the calculated internal calcium concentrations
for i = (1:length(Caunprim_rate_off))
     [SS_normal_off]= calculate_steady_state(par, CaIntra(i), 1);
      SS_normal_off(end+1) = SS_normal_off(1)*(Caunprim_rate_off(i)/Caprim_rate_off(i));
      SS_total_off{i} = SS_normal_off/sum(SS_normal_off);
end

% Calculate % occupied release site for the calculated internal calcium concentrations
for i = (1:length(SS_total_off))
   percent_off(:,i) = (sum(SS_total_off{1,i})-SS_total_off{1, i}(end))*100;
end


%% Stable/unstable ratio as a function of the external calcium concentration
kD = (40.757e-9)^2;
K_M_current = 2.679*10^-3;
CaMax = 190e-9;
n = 2;

% figure()
% %subplot(2,2,3)
% hold on
% fplot (@(calc_ext) 100*((CaMax*calc_ext/(K_M_current+calc_ext))^n)/(kD+((CaMax*calc_ext/(K_M_current+calc_ext))^n)));
% ylim([0 100]);
% xlim([0 6e-3])
% 
% ylabel('steady state stable %')
% xlabel('external calcium concentration')
% set(gca, 'TickDir', 'out')
% box off

%% Stable/unstable ratio as a function of the internal calcium concentration, unpriming rate (Figure 2 B)
kD = (40.757e-9)^2;
K_M_current = 2.679*10^-3;
CaMax = 190e-9;
n = 2;


figure()
hold on
fplot (@(calc_int) 100*((calc_int)^n)/(kD+(calc_int)^n));
ylim([0 100]);
xlim([0 1.4e-07])
ylabel('steady state stable %')
xlabel('internal calcium concentration')
yyaxis right
plot (CaIntra, Caunprim_rate)
plot (CaIntra, Caunprim_rate_on)
ylabel('unpriming rate')
ylim([-50 300]);
legend('WT %unc13 bound to Ca-Cam', 'WT unpriming rate', 'Mutant unpriming rate')
set(gca, 'TickDir', 'out')
hold off
box off

saveas(gcf,[save_folder ,date ,'_StablePercent_UnprimRate_CaInt'],'epsc')


%% plotting internal Ca concentration related to external Ca concentration (Figure S1 A)
figure()
fplot (@(calc_ext) CaMax*calc_ext/(K_M_current+calc_ext));
ylabel('internal calcium concentration')
xlabel('external calcium concentration')
set(gca, 'TickDir', 'out')
xlim([0 6e-3])
ylim([0 1.5e-7])
box off
saveas(gcf,[save_folder ,date ,'_CaInt_CaExt'],'epsc')

%% plotting unpriming rates (Figure S1 B)
figure()
hold on
plot (CaExt_extr, Caunprim_rate_ext, 'LineWidt', 2)
plot (CaExt_extr, Caunprim_rate_ext_on, 'LineWidt', 2)
ylabel('unpriming rate')
xlabel('external calcium concentration')
legend('WT unpriming rate', 'Mutant unpriming rate')
set(gca, 'TickDir', 'out')
xlim([0 6e-3])
ylim([-50 300])
box off
saveas(gcf,[save_folder ,date ,'_UnprimRate_CaExt'],'epsc')
