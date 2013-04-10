%%%%%%%%% Trimming Oval Template %%%%%%%%%%%%%%%

clc
clear all
close all
cd('/Users/dtaniguchi/Documents/Cuttlefish Project/Template Fitting')%setting the directory
directory = dir;

set(0,'defaultaxesfontsize',30);
set(0,'defaulttextfontsize',30);
set(0,'defaultlinelinewidth',2);

%% Reading in template, trimming, and finding edges 
    %Have to use the template because the images themselves don't always
    %have complete edges (i.e., there are often discontinuities)
template_img = imread('Oval_Template3.png');%<-------------------MAKE SURE HAVE CORRECT TEMPLATE

% Convert image to grayscale if colored
if size(template_img,3)==3
    template_img=rgb2gray(template_img);
end

% Find edges and make everything inside black
img_edge = edge(template_img,'prewitt'); %Edges of template
% [temp_row, temp_col] = find(template_img~=255); %finding non-white parts
% template_img_temp = template_img;
% template_img_temp(temp_row, temp_col) = 0; %Making all non-white parts black

% Trim template
rv = [];%matrix of the rows that don't contain edges
cv = [];%matrix of columns of the rotated image that don't contain edges
for ii = 1:size(template_img,1) % rows
    if isempty(find(img_edge(ii,:))) % indicates if there is an empty horizontal edge 
        rv = [rv ii];
    end
end
for ii = 1:size(template_img,2) % columns
    if isempty(find(img_edge(:,ii))) % indicates if there is an empty vertical edge 
        cv = [cv ii];
    end
end
template_img(rv(2:end-1),:) = [];%deletes the rows that didn't have any edges, but not all the way to edge or else can't find edges later
template_img(:,cv) = [];%deletes the columns that don't have any edges

imwrite(template_img,'Oval_Template3_Trimmed.png','png')