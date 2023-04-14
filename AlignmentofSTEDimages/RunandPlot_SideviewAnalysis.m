
%% This script runs the alignment of side view images, and plots average, substracted and scaled images, and line profiles.

folderpath = ''; % Set name to the path of the stacks here
%Within this folder, Image stacks should be generated in a subfolder called
%according to the genotype. 

genotype = 'Bdel';
Channel1 = 'BRP'; 
Channel2 = 'Unc13';
option = 'bilinear'; %option for rotation algorithm
direction = 'Side'; %side or topview
larverange = [2, 38]; %Bdel, range for images

Rotate_and_align_images(folderpath, genotype, Channel1, Channel2, option, larverange);
[Bdel_v_BRP_scaled,Bdel_v_Unc13_scaled,Bdel_h_BRP_scaled,Bdel_h_Unc13_scaled, Bdel_BRP_images, Bdel_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

genotype = 'CaM';
larverange = [1, 30]; %Cam

Rotate_and_align_images(folderpath, genotype, Channel1, Channel2, option, larverange);
[CaM_v_BRP_scaled,CaM_v_Unc13_scaled,CaM_h_BRP_scaled,CaM_h_Unc13_scaled, CaM_BRP_images, CaM_Unc13_images]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, direction);

%% Sideview plotting scaled averages and substracted images

%compute min and max values in the heatmap, to ensure that all images have
%are displayed on the same scale
cat_img = (cat(3,CaM_Unc13_images{:}));
CaM_Unc13_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);
CaM_Unc13_images_avg_scaled = imrotate(CaM_Unc13_images_avg_scaled, 180, option, 'crop'); %turn image

cat_img = (cat(3,Bdel_Unc13_images{:}));
Bdel_Unc13_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);
Bdel_Unc13_images_avg_scaled = imrotate(Bdel_Unc13_images_avg_scaled, 180, option, 'crop'); %turn image

cat_img = (cat(3,CaM_BRP_images{:}));
CaM_BRP_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);
CaM_BRP_images_avg_scaled = imrotate(CaM_BRP_images_avg_scaled, 180, option, 'crop'); %turn image

cat_img = (cat(3,Bdel_BRP_images{:}));
Bdel_BRP_images_avg_scaled = mean(cat_img./(sum(sum(cat_img))),3);
Bdel_BRP_images_avg_scaled = imrotate(Bdel_BRP_images_avg_scaled, 180, option, 'crop');

min_Cam_unc13 = min(min(CaM_Unc13_images_avg_scaled));
max_Cam_unc13 = max(max(CaM_Unc13_images_avg_scaled));
min_Bdel_unc13 = min(min(Bdel_Unc13_images_avg_scaled));
max_Bdel_unc13 = max(max(Bdel_Unc13_images_avg_scaled));
min_Cam_BRP = min(min(CaM_BRP_images_avg_scaled));
max_Cam_BRP = max(max(CaM_BRP_images_avg_scaled));
min_Bdel_BRP = min(min(Bdel_BRP_images_avg_scaled));
max_Bdel_BRP = max(max(Bdel_BRP_images_avg_scaled));

%% Average rotated image in the Unc13 channel
figure() 
heatmap(CaM_Unc13_images_avg_scaled)
caxis([min([min_Bdel_unc13, min_Cam_unc13]) max([max_Bdel_unc13, max_Cam_unc13])])
colormap(gray);
title('CaM unc13')
imwrite((CaM_Unc13_images_avg_scaled), [folderpath, '/Figure/Scaled images/CaM_Unc13_images_avg.tiff']);
%saveas(gcf, [folderpath, '/Figure/Heatmap_Bdel_Unc13_images_avg.tiff']) %save heatmaps
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_CaM_Unc13_images_avg.eps']) %save heatmaps

figure()
heatmap(Bdel_Unc13_images_avg_scaled)
caxis([min([min_Bdel_unc13, min_Cam_unc13]) max([max_Bdel_unc13, max_Cam_unc13])])
colormap(gray);
title('Bdel unc13')
imwrite((Bdel_Unc13_images_avg_scaled), [folderpath, '/Figure/Scaled images/Bdel_Unc13_images_avg.tiff']);
%saveas(gcf, [folderpath, '/Figure/Heatmap_Bdel_Unc13_images_avg.tiff']) %save heatmaps
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Bdel_Unc13_images_avg.eps'])

%% Average rotated image in the BRP channel
figure()
heatmap(Bdel_BRP_images_avg_scaled)
caxis([min([min_Bdel_BRP, min_Cam_BRP]) max([max_Bdel_BRP, max_Cam_BRP])])
colormap(gray);
title('Bdel BRP')
imwrite((Bdel_BRP_images_avg_scaled), [folderpath, '/Figure/Scaled images/Bdel_BRP_images_avg.tiff']);
%saveas(gcf, [folderpath, '/Figure/Heatmap_Bdel_Unc13_images_avg.tiff']) %save heatmaps
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_Bdel_BRP_images_avg.eps'])

figure()
heatmap(CaM_BRP_images_avg_scaled)
caxis([min([min_Bdel_BRP, min_Cam_BRP]) max([max_Bdel_BRP, max_Cam_BRP])])
colormap(gray);
title('CaM BRP')
imwrite((CaM_BRP_images_avg_scaled), [folderpath, '/Figure/Scaled images/CaM_BRP_images_avg.tiff']);
%saveas(gcf, [folderpath, '/Figure/Heatmap_Bdel_Unc13_images_avg.tiff']) %save heatmaps
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_CaM_BRP_images_avg.eps'])

%% Substracted images 
min_diff_Unc13 = min(min(CaM_Unc13_images_avg_scaled - Bdel_Unc13_images_avg_scaled));
max_diff_Unc13 = max(max(CaM_Unc13_images_avg_scaled - Bdel_Unc13_images_avg_scaled));
min_diff_BRP = min(min(CaM_BRP_images_avg_scaled - Bdel_BRP_images_avg_scaled));
max_diff_BRP = max(max(CaM_BRP_images_avg_scaled - Bdel_BRP_images_avg_scaled));

x = max(max_diff_Unc13, max_diff_BRP);%max(CaM_Unc13_images_avg_scaled- Bdel_Unc13_images_avg_scaled);
y=abs(min(min_diff_Unc13, min_diff_BRP));%abs(min(CaM_Unc13_images_avg_scaled- Bdel_Unc13_images_avg_scaled));
scalebar = round(x/(x+y)*256);
negColorMap = [zeros(1, 256-scalebar), linspace(0, 1, scalebar)];
posColorMap = [linspace(1, 0, 256-scalebar), zeros(1, scalebar)];
colorMap = [negColorMap; zeros(1, 256); posColorMap]';

figure()
heatmap(CaM_Unc13_images_avg_scaled- Bdel_Unc13_images_avg_scaled)
colormap(colorMap);
caxis([min([min_diff_Unc13, min_diff_BRP]) max([max_diff_Unc13, max_diff_BRP])])
title('CaM - Bdel: unc13')
%imwrite((CaM_Unc13_images_avg- Bdel_Unc13_images_avg), [folderpath, '/Figure/CaM_Bdel_Unc13_images_avg.tiff']);
%saveas(gcf, [folderpath, '/Figure/Heatmap_CaM_Bdel_Unc13_images_avg.tiff']) %save heatmaps
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_CaM_Bdel_Unc13_images_avg.epsc']) %save heatmaps

figure()
heatmap(CaM_BRP_images_avg_scaled- Bdel_BRP_images_avg_scaled)
colormap(colorMap);
caxis([min([min_diff_Unc13, min_diff_BRP]) max([max_diff_Unc13, max_diff_BRP])])
title('CaM - Bdel: BRP')
%imwrite((CaM_BRP_images_avg- Bdel_BRP_images_avg), [folderpath, '/Figure/CaM_Bdel_BRP_images_avg.tiff']);
%saveas(gcf, [folderpath, '/Figure/Heatmap_CaM_Bdel_BRP_images_avg.tiff']) %save heatmaps
saveas(gcf, [folderpath, '/Figure/Scaled images/Heatmap_CaM_Bdel_BRP_images_avg.epsc']) %save heatmaps



%% Side view line profiles
figure()
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
saveas(gcf, [folderpath, '/Figure/Verticalline_CaMvsBdel_BRP.eps']) %save heatmaps

%
figure()
plot([0:24]*20,sum(Bdel_h_BRP_scaled)/size(Bdel_h_BRP_scaled,1), 'k', 'LineWidth', 4); hold on %% ADJUSTrunning from end to beginning
plot([0:24]*20,sum(CaM_h_BRP_scaled)/size(CaM_h_BRP_scaled,1), 'r', 'LineWidth', 4)
xlim([0,24]*20)
set(gca,'TickDir','out');
box off
title('Horizontal line -BRP')
legend('Bdel', 'CaM')
xlabel('nm')
ylabel('probability density')
saveas(gcf, [folderpath, '/Figure/Horizontalline_CaMvsBdel_BRP.eps']) %save heatmaps

%
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
saveas(gcf, [folderpath, '/Figure/Verticalline_CaMvsBdel_Unc13.eps']) %save heatmaps

