function Align_topview(folderpath, genotype, Channel1, Channel2, larverange)
% This function aligns images based on a circle drawn in one of
% the channels
% The aligned images are saved in tiff files. 
%Alignment script for top view images

%Inputs:
%folderpath: (string) the input folderpath contains the path to the folder with the
%images
%genotype: (string) name of analysed genotype, corresponding to the name of the
%folder with images. 
%Channel1, (string) name of imaged channel as in the image name
%Channel2, (string) name of imaged channel as in the image name
%larvenos = numbers of larve needed to be analysed.

% genotype = 'Bdel';
% Channel1 = 'BRP';
% Channel2 = 'Unc13';

%larvenos = [2,3,5,8,10,12:16,18,20:22]; %Bdel
%larvenos = [1, 2,4,6:8,10:15]; %Cam

L= length(imfinfo([folderpath, '\',genotype, '\CombStack-Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif']));
%Check if file name already exists
if isfile([folderpath, '/',genotype,'/AlignedCircle_Top_', genotype,'_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'])
    disp('Rotation is already performed for this dataset')
    return
end

for i=1:L %larveno = larvenos
%%load images
Channel1_tiff = im2double(imread([folderpath, '/',genotype, '/CombStack-Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif'],i));
Line_tiff = im2double(imread([folderpath, '/',genotype, '/Circle_CombStack-Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif'],i) );
Channel2_tiff = im2double(imread([folderpath, '/',genotype, '/CombStack-Top_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif'],i) );


%% Adjust image
size_image = size(Line_tiff);

%Check if all images have the same size:
size_1 = size(Channel1_tiff);
size_2 = size(Channel2_tiff);
if size_1 ~= size_image | size_2 ~= size_image
    warning('Image sizes do not correspond')
    return
end

mid = round(size_image/2); %wanted midpoint of Circle

%find current midpoint of Circle
[x, y] = find(Line_tiff == 1);
rx = (max(x)-min(x))/2;
ry = (max(y)-min(y))/2;
if rx == ry
    midcurr = [min(x)+round(rx), min(y)+round(rx)];
elseif rx > ry
    midcurr = [min(x)+round(rx), round(mean(y(x == min(x))))];
elseif rx < ry
    midcurr = [round(mean(x(y==min(y)))),min(y)+round(ry)];
end


Line_new = Line_tiff;
Channel1_new = Channel1_tiff;
Channel2_new = Channel2_tiff;

if midcurr(1) <mid(1)
    
    %Adjust line plot
    Line_new = [ones(abs(midcurr(1)- mid(1)), size_image(1))*0; Line_new];
    Line_new(end-abs(midcurr(1)- mid(1))+1:end, : )=[];
    
    %Adjust Channel plots
    Channel1_new = [ones(abs(midcurr(1)- mid(1)), size_image(1))*0; Channel1_new];
    Channel1_new(end-abs(midcurr(1)- mid(1))+1:end, : )=[];
    
    Channel2_new = [ones(abs(midcurr(1)- mid(1)), size_image(1))*0; Channel2_new];
    Channel2_new(end-abs(midcurr(1)- mid(1))+1:end, : )=[];

    
elseif midcurr(1) > mid(1)
    
    %Adjust line plot    
    Line_new = [Line_new; ones( abs(midcurr(1)- mid(1)), size_image(1))*0];
    Line_new(1:abs(midcurr(1)- mid(1)), :) =[];
    
    %Adjust Channel plots    
    Channel1_new = [Channel1_new; ones( abs(midcurr(1)- mid(1)), size_image(1))*0];
    Channel1_new(1:abs(midcurr(1)- mid(1)), :) =[];
    
    Channel2_new = [Channel2_new; ones( abs(midcurr(1)- mid(1)), size_image(1))*0];
    Channel2_new(1:abs(midcurr(1)- mid(1)), :) =[];
    
end

if midcurr(2) <mid(2) 
    
    %Adjust line plot
    Line_new = [ones(size_image(2),abs(midcurr(2)- mid(2)))*0, Line_new];
    Line_new(:, end-abs(midcurr(2)- mid(2))+1:end) = [];
    
    %Adjust Channel plots
    Channel1_new = [ones(size_image(2),abs(midcurr(2)- mid(2)))*0, Channel1_new];
    Channel1_new(:, end-abs(midcurr(2)- mid(2))+1:end) = [];
    
    Channel2_new = [ones(size_image(2),abs(midcurr(2)- mid(2)))*0, Channel2_new];
    Channel2_new(:, end-abs(midcurr(2)- mid(2))+1:end) = [];
    
elseif midcurr(2) > mid(2)
    
    %Adjust line plot
    Line_new = [Line_new, ones(size_image(2), abs(midcurr(2)- mid(2)))*0];
    Line_new(:, 1:abs(midcurr(2)- mid(2))) =[];
    
    %Adjust Channel plots
    Channel1_new = [Channel1_new, ones(size_image(2), abs(midcurr(2)- mid(2)))*0];
    Channel1_new(:, 1:abs(midcurr(2)- mid(2))) =[];
    
    Channel2_new = [Channel2_new, ones(size_image(2), abs(midcurr(2)- mid(2)))*0];
    Channel2_new(:, 1:abs(midcurr(2)- mid(2))) =[];
    
end


%% saving
imwrite((Line_new), [folderpath, '/',genotype,'/AlignedCircle_Top_', genotype,'_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'], 'WriteMode', 'append');
imwrite((Channel1_new), [folderpath, '/',genotype,'/Aligned_Top_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'], 'WriteMode', 'append');
imwrite((Channel2_new), [folderpath, '/',genotype,'//Aligned_Top_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'], 'WriteMode', 'append');
