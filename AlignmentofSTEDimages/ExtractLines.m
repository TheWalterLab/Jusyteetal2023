function [line_v_1_scaled,line_v_2_scaled,line_h_1_scaled,line_h_2_scaled, image_1_cell, image_2_cell]=ExtractLines(folderpath, genotype, Channel1, Channel2, larverange, TopORSide)
% This function extracts the vertical and horizontal line going through the
% centre of the image.

%input
%folderpath: (string) the input folderpath contains the path to the folder with the
%images
%genotype: (string) name of analysed genotype, corresponding to the name of the
%folder with images. 
%Channel1, (string) name of imaged channel as in the image name
%Channel2, (string) name of imaged channel as in the image name
%larvenos = numbers of larve needed to be analysed.
%TopORSide = 'Side', when wanting to extract information from side view
%images and 'Top', when wanting to extract information from top view images

%output
%vertical and horizontal line profiles (indicated by v and h) in the first
%and second channel (indicated by 1 and 2)
%and collection of images from which these line profiles are extracted

%%

if strcmp(TopORSide, 'Side')
L= length(imfinfo([folderpath, '\',genotype, '\Rotated_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff']));
elseif strcmp(TopORSide, 'Top')
    L= length(imfinfo([folderpath, '\',genotype, '\Aligned_Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff']));
end

for i=1:L %larveno = larvenos


 %load images
if strcmp(TopORSide, 'Side')
Channel1_tiff = im2double(imread([folderpath, '/',genotype, '/Rotated_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i));
Channel2_tiff = im2double(imread([folderpath, '/',genotype, '/Rotated_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i) );
elseif strcmp(TopORSide, 'Top')
  Channel1_tiff = im2double(imread([folderpath, '/',genotype, '/Aligned_Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i));
Channel2_tiff = im2double(imread([folderpath, '/',genotype, '/Aligned_top_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tiff'],i) );
 end

size_image = size(Channel1_tiff);
line_v_1(i, :) = Channel1_tiff(:, round(size_image(2)/2));
line_v_2(i, :) = Channel2_tiff(:, round(size_image(2)/2));

if strcmp(TopORSide, 'Side')
[m, f] = max(line_v_1(i,:));
line_h_1(i, :) = Channel1_tiff(f,:);

[m, f] = max(line_v_2(i,:));
line_h_2(i, :) = Channel2_tiff(f,:);
elseif strcmp(TopORSide, 'Top')
line_h_1(i, :) = Channel1_tiff(round(size_image(2)/2),:);
line_h_2(i, :) = Channel2_tiff(round(size_image(2)/2),:);
end

line_v_1_scaled(i, :) = line_v_1(i,:)./sum(line_v_1(i,:));
line_v_2_scaled(i, :) = line_v_2(i,:)./sum(line_v_2(i,:));
line_h_1_scaled(i, :) = line_h_1(i,:)./sum(line_h_1(i,:));
line_h_2_scaled(i, :) = line_h_2(i,:)./sum(line_h_2(i,:));

image_1_cell{i} = Channel1_tiff;
image_2_cell{i} = Channel2_tiff;


end