
%% This script runs the bootstrapping analysis for line profiles obtained for both side and top view


%% Sideview alignment
folderpath = ''; % Set name to the path of the stacks here
%Within this folder, Image stacks should be generated in a subfolder called
%according to the genotype. 
genotype = 'Bdel';
Channel1 = 'BRP';
Channel2 = 'Unc13';
option = 'bilinear'; 
direction = 'Side';
larverange = [2, 38]; %Bdel
settings = [100]; %repetitions of bootstrapping

%% load lines - side view
%Rotate_and_align_images(folderpath, genotype, Channel1, Channel2, option,
%larverange); %Run this line, when images are not aligned yet.
[Bdel_v_BRP_scaled,Bdel_v_Unc13_scaled,Bdel_h_BRP_scaled,Bdel_h_Unc13_scaled, Bdel_BRP_images, Bdel_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

genotype = 'CaM';
larverange = [1, 30]; %Cam

%Rotate_and_align_images(folderpath, genotype, Channel1, Channel2, option,
%larverange); %Run this line, when images are not aligned yet.
[CaM_v_BRP_scaled,CaM_v_Unc13_scaled,CaM_h_BRP_scaled,CaM_h_Unc13_scaled, CaM_BRP_images, CaM_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

%% bootstrap lines - side view
lines1 = Bdel_v_BRP_scaled;
lines2 = CaM_v_BRP_scaled;
[lines1_mean, lines2_mean, Dkl_BRP] = bootstraplines(lines1,lines2, settings);

%Plot bootstrapped line profiles
figure()
plot([1:25]*20,lines1_mean', 'Color', [0.5, 0.5, 0.5, 0.4]);
hold on
plot([1:25]*20,lines2_mean','Color',[1 0 0, 0.2])
title('BRP')
xlim([0 25*20])
xlabel('Distance (nm)')
ylabel('Probability density')
box off
set(gca, 'Fontsize', 14)
mean(Dkl_BRP)

lines1 = Bdel_v_Unc13_scaled;
lines2 = CaM_v_Unc13_scaled;
[lines1_mean, lines2_mean, Dkl_Unc13] = bootstraplines(lines1,lines2, settings);
figure()
plot([1:25]*20,lines1_mean', 'Color', [0.5, 0.5, 0.5, 0.4])
hold on
plot([1:25]*20,lines2_mean','Color',[1 0 0, 0.2])
title('Unc13')

xlim([0 25*20])
xlabel('Distance (nm)')
ylabel('Probability density')
box off
set(gca, 'Fontsize', 14)

% KL divergence for BRP and Unc13 channel, in each channel Bdel is compared
% to CaM mutant
figure()
boxplot([Dkl_BRP,Dkl_Unc13],'Notch','on','Labels',{'BRP','Unc13'},'Whisker',1)
box off
set(gca, 'Fontsize', 14)
ylabel('Kullback-Leibler divergence')

[h, p, stats]=ttest2(Dkl_BRP, Dkl_Unc13)

%% Topview - load lines
folderpath = ;

genotype = 'Bdel';
larverange = [1, 47]; 
direction = 'Top'; 
[Bdel_v_BRP_scaled,Bdel_v_Unc13_scaled,Bdel_h_BRP_scaled,Bdel_h_Unc13_scaled, Bdel_BRP_images, Bdel_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

genotype = 'CaM';
larverange = [1, 42]; 

[CaM_v_BRP_scaled,CaM_v_Unc13_scaled,CaM_h_BRP_scaled,CaM_h_Unc13_scaled, CaM_BRP_images, CaM_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

%% Topview - Bootstrap line profiles
lines1 = Bdel_v_BRP_scaled;
lines2 = CaM_v_BRP_scaled;
[lines1_mean, lines2_mean, Dkl_BRP] = bootstraplines(lines1,lines2, settings);

%% Plot bootstrapped lines
figure()
plot([1:25]*20,lines1_mean', 'Color', [0.5, 0.5, 0.5, 0.4]);
hold on
plot([1:25]*20,lines2_mean','Color',[1 0 0, 0.2])
title('BRP')
xlim([0 25*20])
xlabel('Distance (nm)')
ylabel('Probability density')
%legend('Ctrl', 'Unc13CaM^w^r^w^r')
box off
set(gca, 'Fontsize', 14)
mean(Dkl_BRP)

lines1 = Bdel_v_Unc13_scaled;
lines2 = CaM_v_Unc13_scaled;
[lines1_mean, lines2_mean, Dkl_Unc13] = bootstraplines(lines1,lines2, settings);
figure()
plot([1:25]*20,lines1_mean', 'Color', [0.5, 0.5, 0.5, 0.4])
hold on
plot([1:25]*20,lines2_mean','Color',[1 0 0, 0.2])
title('Unc13')
%legend('Ctrl', 'Unc13CaM^w^r^w^r', 'box',  'off')
xlim([0 25*20])
xlabel('Distance (nm)')
ylabel('Probability density')
box off
set(gca, 'Fontsize', 14)
mean(Dkl_Unc13)

figure()
boxplot([Dkl_BRP,Dkl_Unc13],'Notch','on','Labels',{'BRP','Unc13'},'Whisker',1)
box off
set(gca, 'Fontsize', 14)
ylabel('Kullback-Leibler divergence')

[h, p, stats]=ttest2(Dkl_BRP, Dkl_Unc13)