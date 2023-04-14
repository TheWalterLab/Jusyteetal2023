function Rotate_and_align_images(folderpath, genotype, Channel1, Channel2, option, larverange)
% This function rotates and aligns images based on a line drawn in one of
% the channels
% The rotated images are saved in tiff files. 

% !!! Script specifically implemented to align side view images !!! 

%Inputs:
%folderpath: (string) the input folderpath contains the path to the folder with the
%images
%genotype: (string) name of analysed genotype, corresponding to the name of the
%folder with images. 
%Channel1, (string) name of imaged channel as in the image name
%Channel2, (string) name of imaged channel as in the image name
%option - selected option for rotation of the images - 'nearest',
%'bilinear', 'bicubic'
%larvenos = numbers of larve needed to be analysed.

if isfile([folderpath, '/',genotype,'/RotatedLine_', genotype,'_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff']))
    disp('Rotation is already performed for this dataset')
    return
end

L= length(imfinfo([folderpath, '\',genotype, '\combStack-Side_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif']));

for i=1:L %larveno = larvenos
%%load images
Channel1_tiff = im2double(imread([folderpath, '/',genotype, '/combStack-Side_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif'],i));
Line_tiff = im2double(imread([folderpath, '/',genotype, '/Line_combStack-Side_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif'],i) );
%Line_tiff = im2double(imread([folderpath, '/',genotype, '/Lines/Line_Stack-Side_', genotype, '_' Channel1, '_Larve', num2str(larveno), '.tif']));
Channel2_tiff = im2double(imread([folderpath, '/',genotype, '/combStack-Side_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif'],i) );
%Channel2_tiff = im2double(imread([folderpath, '/',genotype, '/Stack-Side_', genotype, '_' Channel2, '_Larve', num2str(larveno), '.tif']));
Direction_Line_tiff =  im2double(imread([folderpath, '/',genotype, '/Direction_combStack-Side_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)) , '.tif'],i) );

%% Adjust image
size_image = size(Line_tiff);

%Check if all images have the same size:
size_1 = size(Channel1_tiff);
size_2 = size(Channel2_tiff);
if size_1 ~= size_image | size_2 ~= size_image
    warning('Image sizes do not correspond')
    return
end

mid = round(size_image/2); %wanted midpoint of line

%find current midpoint of line
[x, y] = find(Line_tiff == 1);
midcurr = [(round(x(1) + (x(end) - x(1))/2)), round(y(1)+(y(end) - y(1))/2)];

Line_new = Line_tiff;
Channel1_new = Channel1_tiff;
Channel2_new = Channel2_tiff;
DLine_new = Direction_Line_tiff; 

if midcurr(1) <mid(1)
    
    %Adjust line plot
    Line_new = [ones(abs(midcurr(1)- mid(1)), size_image(1))*0; Line_new];
    Line_new(end-abs(midcurr(1)- mid(1))+1:end, : )=[];
    
    %Adjust Channel plots
    Channel1_new = [ones(abs(midcurr(1)- mid(1)), size_image(1))*0; Channel1_new];
    Channel1_new(end-abs(midcurr(1)- mid(1))+1:end, : )=[];
    
    Channel2_new = [ones(abs(midcurr(1)- mid(1)), size_image(1))*0; Channel2_new];
    Channel2_new(end-abs(midcurr(1)- mid(1))+1:end, : )=[];
    
    %Adjust direction line
    DLine_new = [ones(abs(midcurr(1)- mid(1)), size_image(1))*0; DLine_new];
    DLine_new(end-abs(midcurr(1)- mid(1))+1:end, : )=[];
    
elseif midcurr(1) > mid(1)
    
    %Adjust line plot    
    Line_new = [Line_new; ones( abs(midcurr(1)- mid(1)), size_image(1))*0];
    Line_new(1:abs(midcurr(1)- mid(1)), :) =[];
    
    %Adjust Channel plots    
    Channel1_new = [Channel1_new; ones( abs(midcurr(1)- mid(1)), size_image(1))*0];
    Channel1_new(1:abs(midcurr(1)- mid(1)), :) =[];
    
    Channel2_new = [Channel2_new; ones( abs(midcurr(1)- mid(1)), size_image(1))*0];
    Channel2_new(1:abs(midcurr(1)- mid(1)), :) =[];
    
    %Adjust direction line
      DLine_new = [DLine_new; ones( abs(midcurr(1)- mid(1)), size_image(1))*0];
    DLine_new(1:abs(midcurr(1)- mid(1)), :) =[];
    
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
    
    DLine_new = [ones(size_image(2),abs(midcurr(2)- mid(2)))*0, DLine_new];
    DLine_new(:, end-abs(midcurr(2)- mid(2))+1:end) = [];
    
elseif midcurr(2) > mid(2)
    
    %Adjust line plot
    Line_new = [Line_new, ones(size_image(2), abs(midcurr(2)- mid(2)))*0];
    Line_new(:, 1:abs(midcurr(2)- mid(2))) =[];
    
    %Adjust Channel plots
    Channel1_new = [Channel1_new, ones(size_image(2), abs(midcurr(2)- mid(2)))*0];
    Channel1_new(:, 1:abs(midcurr(2)- mid(2))) =[];
    
    Channel2_new = [Channel2_new, ones(size_image(2), abs(midcurr(2)- mid(2)))*0];
    Channel2_new(:, 1:abs(midcurr(2)- mid(2))) =[];
    
    %adjust direction line
    DLine_new = [DLine_new, ones(size_image(2), abs(midcurr(2)- mid(2)))*0];
    DLine_new(:, 1:abs(midcurr(2)- mid(2))) =[];
    
end

[x, y] = find(Line_new == 1);
mid_new = [(round(x(1) + (x(end) - x(1))/2)), round(y(1)+(y(end) - y(1))/2)];
if mid_new ~= mid
    warning('Adjusted midpoint of line not equal to desired midpoint')
    return
end

%% Adjust angle
alpha = (atan((min(x)-max(x))/(min(y)-max(y))))*180/pi; %Compute angle

if x(1)>x(end)
    alpha = alpha *-1;
end

Line_new = imrotate(Line_new, alpha, option, 'crop');
Channel1_new = imrotate(Channel1_new, alpha, option, 'crop');
Channel2_new = imrotate(Channel2_new, alpha, option, 'crop');
DLine_new = imrotate(DLine_new, alpha, option, 'crop');


if sum(sum(DLine_new(1:mid(1)-1,:))) > sum(sum(DLine_new(mid(1)+1:size_1(1),:)))
    Line_new = imrotate(Line_new, 180, option, 'crop');
    Channel1_new = imrotate(Channel1_new, 180, option, 'crop');
    Channel2_new = imrotate(Channel2_new, 180, option, 'crop');
    DLine_new = imrotate(DLine_new, 180, option, 'crop');
end


[x, y] = find(Line_new >0);
mid_new = [(round(x(1) + (x(end) - x(1))/2)), round(y(1)+(y(end) - y(1))/2)];
if mid_new ~= mid
    warning('Adjusted midpoint of line after rotation not equal to desired midpoint')
    %return
end



%% saving
imwrite((Line_new), [folderpath, '/',genotype,'/RotatedLine_', genotype,'_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'], 'WriteMode', 'append');
imwrite((Channel1_new), [folderpath, '/',genotype,'/Rotated_', genotype, '_' Channel1, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'], 'WriteMode', 'append');
imwrite((Channel2_new), [folderpath, '/',genotype,'//Rotated_', genotype, '_' Channel2, '_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'], 'WriteMode', 'append');
imwrite((DLine_new), [folderpath, '/',genotype,'/Rotated_DirectionLine_', genotype,'_Larve', num2str(larverange(1)), '-',num2str(larverange(2)), '.tiff'], 'WriteMode', 'append');
end
