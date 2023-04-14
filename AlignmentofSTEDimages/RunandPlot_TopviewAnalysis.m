%% This script runs the alignment of top view images, and plots average, substracted and scaled images, and line profiles.

folderpath = ''; % Set name to the path of the stacks here
%Within this folder, Image stacks should be generated in a subfolder called
%according to the genotype. 


genotype = 'Bdel';
Channel1 = 'BRP';
Channel2 = 'Unc13';
option = 'bilinear'; %option for rotation algorithm
larverange = [1, 47]; %Bdel

Align_topview(folderpath, genotype, Channel1, Channel2, larverange)

%extract data
L= length(imfinfo([folderpath, '/',genotype, '/Aligned_Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff']));
for i=1:L %larveno = larvenos
 %load images
Channel1_tiff = im2double(imread([folderpath, '/',genotype, '/Aligned_Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i));
Channel2_tiff = im2double(imread([folderpath, '/',genotype, '/Aligned_Top_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i) );
Bdel_BRP_images{i} = Channel1_tiff;
Bdel_Unc13_images{i} = Channel2_tiff;
end

genotype = 'CaM'; 
larverange = [1, 42];
Align_topview(folderpath, genotype, Channel1, Channel2, larverange)

%extract data
L= length(imfinfo([folderpath, '/',genotype, '/Aligned_Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff']));
for i=1:L %larveno = larvenos
 %load images
Channel1_tiff = im2double(imread([folderpath, '/',genotype, '/Aligned_Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i));
Channel2_tiff = im2double(imread([folderpath, '/',genotype, '/Aligned_Top_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i) );
CaM_BRP_images{i} = Channel1_tiff;
CaM_Unc13_images{i} = Channel2_tiff;
end

%% Plot average images

cat_img = (cat(3,CaM_Unc13_images{:}));
CaM_Unc13_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);

cat_img = (cat(3,Bdel_Unc13_images{:}));
Bdel_Unc13_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);

cat_img = (cat(3,CaM_BRP_images{:}));
CaM_BRP_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);

cat_img = (cat(3,Bdel_BRP_images{:}));
Bdel_BRP_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);

min_Cam_unc13 = min(min(CaM_Unc13_images_avg_scaled));
max_Cam_unc13 = max(max(CaM_Unc13_images_avg_scaled));
min_Bdel_unc13 = min(min(Bdel_Unc13_images_avg_scaled));
max_Bdel_unc13 = max(max(Bdel_Unc13_images_avg_scaled));
min_Cam_BRP = min(min(CaM_BRP_images_avg_scaled));
max_Cam_BRP = max(max(CaM_BRP_images_avg_scaled));
min_Bdel_BRP = min(min(Bdel_BRP_images_avg_scaled));
max_Bdel_BRP = max(max(Bdel_BRP_images_avg_scaled));

% plot images
figure()
heatmap(CaM_Unc13_images_avg_scaled)
caxis([min([min_Bdel_unc13, min_Cam_unc13]) max([max_Bdel_unc13, max_Cam_unc13])])
colormap(gray);
title('CaM unc13')
imwrite((CaM_Unc13_images_avg_scaled), [folderpath, '/Figure/Scaled images/Top_CaM_Unc13_images_avg.tiff']);
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Top_CaM_Unc13_images_avg.epsc']) %save heatmaps

figure()
heatmap(Bdel_Unc13_images_avg_scaled)
caxis([min([min_Bdel_unc13, min_Cam_unc13]) max([max_Bdel_unc13, max_Cam_unc13])])
colormap(gray);
title('Bdel unc13')
imwrite((Bdel_Unc13_images_avg_scaled), [folderpath, '/Figure/Scaled images/Top_Bdel_Unc13_images_avg.tiff']);
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Top_Bdel_Unc13_images_avg.epsc']) %save heatmaps

figure()
heatmap(Bdel_BRP_images_avg_scaled)
caxis([min([min_Bdel_BRP, min_Cam_BRP]) max([max_Bdel_BRP, max_Cam_BRP])])
colormap(gray);
title('Bdel BRP')
imwrite((Bdel_BRP_images_avg_scaled), [folderpath, '/Figure/Scaled images/Top_Bdel_BRP_images_avg.tiff']);
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Top_Bdel_BRP_images_avg.epsc']) %save heatmaps

figure()
heatmap(CaM_BRP_images_avg_scaled)
caxis([min([min_Bdel_BRP, min_Cam_BRP]) max([max_Bdel_BRP, max_Cam_BRP])])
colormap(gray);
title('CaM BRP')
imwrite((CaM_BRP_images_avg_scaled), [folderpath, '/Figure/Scaled images/Top_CaM_BRP_images_avg.tiff']);
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Top_CaM_BRP_images_avg.epsc']) %save heatmaps

%% Plot Substracted images
min_diff_Unc13 = min(min(CaM_Unc13_images_avg_scaled - Bdel_Unc13_images_avg_scaled));
max_diff_Unc13 = max(max(CaM_Unc13_images_avg_scaled - Bdel_Unc13_images_avg_scaled));
min_diff_BRP = min(min(CaM_BRP_images_avg_scaled - Bdel_BRP_images_avg_scaled));
max_diff_BRP = max(max(CaM_BRP_images_avg_scaled - Bdel_BRP_images_avg_scaled));

figure()
if x > max(max_diff_Unc13, max_diff_BRP) & y> abs(min(min_diff_Unc13, min_diff_BRP))
scalebar = round(x/(x+y)*256);
negColorMap = [zeros(1, 256-scalebar), linspace(0, 1, scalebar)];
posColorMap = [linspace(1, 0, 256-scalebar), zeros(1, scalebar)];
colorMap = [negColorMap; zeros(1, 256); posColorMap ]';
else
    disp('Change scaling in side view images - topview has a larger range')
    return
end

figure()
heatmap(CaM_Unc13_images_avg_scaled- Bdel_Unc13_images_avg_scaled)
colormap(colorMap);
caxis([-y x])
title('CaM - Bdel: unc13')
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Top_substracted_Unc13_images_avg.epsc']) %sav

figure()
heatmap(CaM_BRP_images_avg_scaled- Bdel_BRP_images_avg_scaled)
colormap(colorMap);
caxis([-y x])
title('CaM - Bdel: BRP')
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Top_substracted_BRP_images_avg.epsc'])
%% Top view lines

genotype = 'Bdel';
larverange = [1, 47]; 
direction = 'Top'; 
[Bdel_v_BRP_scaled,Bdel_v_Unc13_scaled,Bdel_h_BRP_scaled,Bdel_h_Unc13_scaled, Bdel_BRP_images, Bdel_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

genotype = 'CaM';
larverange = [1, 42]; 

[CaM_v_BRP_scaled,CaM_v_Unc13_scaled,CaM_h_BRP_scaled,CaM_h_Unc13_scaled, CaM_BRP_images, CaM_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

figure() %vertical line profiles BRP
plot([0:24]*20,sum(Bdel_v_BRP_scaled)/size(Bdel_v_BRP_scaled,1), 'k', 'LineWidth', 4); hold on
plot([0:24]*20,sum(CaM_v_BRP_scaled)/size(CaM_v_BRP_scaled,1),'r', 'LineWidth', 4);
xlim([0,24]*20)
set(gca,'TickDir','out');
box off
%plot((Bdel_v_BRP_scaled'), 'Color', [0.5, 0.5, 0.5); 
%plot((CaM_v_BRP_scaled'), 'Color', [1, 0.5,0]);
title('Vertical line - BRP')
legend('Bdel', 'CaM')
xlabel('nm')
ylabel('probability density')
%saveas(gcf, [folderpath, '/Figure/Verticalline_CaMvsBdel_BRP.eps']) %save heatmaps


figure() %horizontal line profiles BRP
plot([0:24]*20,sum(Bdel_h_BRP_scaled)/size(Bdel_h_BRP_scaled,1), 'k', 'LineWidth', 4); hold on %% ADJUSTrunning from end to beginning
plot([0:24]*20,sum(CaM_h_BRP_scaled)/size(CaM_h_BRP_scaled,1), 'r', 'LineWidth', 4)
xlim([0,24]*20)
set(gca,'TickDir','out');
box off
title('Horizontal line -BRP')
legend('Bdel', 'CaM')
xlabel('nm')
ylabel('probability density')
%saveas(gcf, [folderpath, '/Figure/Horizontalline_CaMvsBdel_BRP.eps']) %save heatmaps

figure() %vertical lines unc13
plot([0:24]*20,sum(Bdel_v_Unc13_scaled)/size(Bdel_v_Unc13_scaled, 1), 'k', 'LineWidth', 4); hold on
plot([0:24]*20,sum(CaM_v_Unc13_scaled)/size(CaM_v_Unc13_scaled,1), 'r', 'LineWidth', 4)
xlim([0,24]*20)
ylim([0, 0.08])
set(gca,'TickDir','out');
box off
title('Vertical line - Unc13')
legend('Bdel', 'CaM')
xlabel('nm')
ylabel('probability density')
%saveas(gcf, [folderpath, '/Figure/Verticalline_CaMvsBdel_Unc13.eps']) %save heatmaps

figure() %horizontal lines unc13
plot(sum(Bdel_h_Unc13_scaled)/size(Bdel_h_Unc13_scaled,1)); hold on
plot(sum(CaM_h_Unc13_scaled)/size(CaM_h_Unc13_scaled,1));
title('Horizontal line -Unc13')
legend('Bdel', 'CaM')

figure()
plot(sum(Bdel_v_BRP_scaled)/size(Bdel_v_BRP_scaled,1)); hold on
plot(sum(Bdel_v_Unc13_scaled)/size(Bdel_v_Unc13_scaled, 1));
legend('BRP', 'Unc13')
title('Bdel')

figure()
plot(sum(CaM_v_BRP_scaled)/size(CaM_v_BRP_scaled,1));
hold on
plot(sum(CaM_v_Unc13_scaled)/size(CaM_v_Unc13_scaled,1));
legend('BRP', 'Unc13')
title('CaM')


