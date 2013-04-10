%%%%% Image Retexturizing%%%%%%%
%Use this code to fit a cuttlefish image to a cuttlefish-shaped template

%% Setting the stage

clc;
clear all;
close all;
cd('/Users/dtaniguchi/Documents/Cuttlefish Project/Template Fitting')

set(0,'defaultaxesfontsize',30);%Setting plot paramters
set(0,'defaulttextfontsize',30);
set(0,'defaultlinelinewidth',2);

%% Init parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic %Sets the timer
plot_flag = 1; %If this is 1, this script will plot images at the end of this code

%% Read and pre-process template image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in a template image
%template_img = imread('cuttlefish_template2.tif'); %Cuttlefish template
template_img = imread('Oval_Template2_Trimmed.png'); %Oval template

% Convert image to grayscale if colored
if size(template_img,3)==3
    template_img=rgb2gray(template_img);
end
%template_img = double(template_img);

% Pre-process template to ensure all colors inside the shape are
% black (i.e., value 0) and all colors outside the shape are white
% (i.e., 255)
template_edge = edge(template_img, 'prewitt');

% Go through each column and find the border of template_edge
for ii = 1:size(template_edge,2)
    template_edge_upper = find(template_edge(:,ii),1,'first');
    template_edge_lower = find(template_edge(:,ii),1,'last');
    for jj = 1:size(template_edge,1)%go through each row
        if jj > template_edge_upper && jj < template_edge_lower%indicates row inside the shape
            template_img(jj,ii) = 0;%make the pixel black
        else
            template_img(jj,ii) = 255;
        end
    end
end

%% Read and pre-process original image (orig_img) %%%%%%%%%%%%%%%%%%%%%%%%%
% Read in an original image
%orig_img=imread('Cuttlefish3.jpg'); %This image was cut out by hand in photoshop
orig_img = imread('1_Checkerboard_2_2012_11_16_13_50_3_Static_210_18.png'); %This image was cut out using the automated method

% Convert image to grayscale if colored
if size(orig_img,3)==3
    orig_img=rgb2gray(orig_img);
end

[r c] = size(orig_img);
orig_img = [128*ones(r,1), orig_img, 128*ones(r,1)];
orig_img = [128*ones(1,c+2); orig_img; 128*ones(1,c+2)];

% Make gray background color white, for automated segmented images
orig_img(find(orig_img ==128)) = 255;

% Pre-process the orig image so the cuttlefish head is facing left
rot_deg = image_rotation(orig_img);
orig_img_rotated = imrotate(orig_img,rot_deg);

% Pre-process orig_img_rotated to ensure all colors outside the shape are white
% (i.e., 255)
orig_rotated_edge = edge(orig_img_rotated, 'prewitt');
if rot_deg ~= 180 % there is no black background if the image is rotated multiples of 180 degrees 
    for ii = 1:size(orig_rotated_edge,2) % go through each column
        edge_upper = find(orig_rotated_edge(:,ii),1,'first');
        edge_lower = find(orig_rotated_edge(:,ii),1,'last');
        if ~isempty(edge_upper)
            for jj = 1:size(orig_rotated_edge,1)%go through each row
                if jj <= edge_upper+1 || jj >= edge_lower-1%indicates row outside the shape 
                    if orig_img_rotated(jj,ii) == 0 %if the pixel is black
                        orig_img_rotated(jj,ii) = 255;%make the pixel white
                    end
                end
            end
        else 
            orig_img_rotated(:,ii) = 255;
        end
    end
end

% Another way to pre-process orig_img_rotated to ensure all colors inside the shape are
% black (i.e., value 0) and all colors outside the shape are white
% (i.e., 255)
%orig_img_rotated(find(orig_img_rotated > 200)) = 0;%makes (almost) white into black.  Use 200 instead of 255 exactly because of roughness in the edges
%orig_img_rotated(find(orig_img_rotated == 0)) = 255;%makes all black into white

% Trim the orig_img_rotated of any empty vertical and horizontal edges
rv = [];%matrix of the rows that don't contain edges
cv = [];%matrix of columns of the rotated image that don't contain edges
for ii = 1:size(orig_rotated_edge,1) % rows
    if isempty(find(orig_rotated_edge(ii,:))) % indicates if there is an empty horizontal edge 
        rv = [rv ii];
    end
end
for ii = 1:size(orig_rotated_edge,2) % columns
    if isempty(find(orig_rotated_edge(:,ii))) % indicates if there is an empty vertical edge 
        cv = [cv ii];
    end
end
orig_img_rotated(rv,:) = [];%deletes the rows that didn't have any edges
orig_img_rotated(:,cv) = [];%deletes the columns that don't have any edges


%% Fitting processed image to the template

if ~isa(orig_img_rotated,'uint8') % if not already a type uint8
    orig_img_rotated = uint8(orig_img_rotated);%make images into uint8 format
end

% image_realign() tries to realign the template_img and orig_img_rotated such
% that the orig_img_rotated completely overlaps template_img
%retexturized_img -- if an overlap between the image and template is found, this is the cut out of the original image in the shape of the template
    %If no solution is found (flag ==1), this matrix is empty
%flag -- indicates if an overlap between the image and template is found (=0) or not found (=1)
%incdec -- indicates if the image needs to be increased (=1) or decreased (=0)
%offset_template -- tells the location where the template fits best in the image
[retexturized_img flag incdec offset_template] = image_realign(template_img, orig_img_rotated);

partial_coverage = 1; %initialization of the while loop below
    %(Within the loop, it's an indicator of whether or not there is a solution from the image_realign function)
  
orig_img_tmp = orig_img_rotated;%renaming the rotated image

while partial_coverage
    if incdec % if need to increase orig_img_tmp size       
        sprintf('Increasing orig image size');
        orig_img_tmp = imresize(orig_img_tmp, 1.05); %increase image by 105%
        [final_img partial_coverage incdec offset_template] = image_realign(template_img,orig_img_tmp);         
    else %if need to decrease orig_img_tmp size
        sprintf('Decreasing orig image size');
        orig_img_tmp = imresize(orig_img_tmp, 0.95);%decrease image by 95%
        [retexturized_img partial_coverage incdec2 offset_template2] = image_realign(template_img,orig_img_tmp);
        
        if ~isempty(offset_template2)
            offset_template = offset_template2;
        end
        
        if ~isempty(retexturized_img)
            final_img = retexturized_img;
        end
        
        if incdec2 == incdec % the orig2 may still be too big for template
            partial_coverage = 1; % keep the while loop alive
        else
            partial_coverage = 0; % optimal overlap found, terminate while loop
        end
    end
end

%% Post-processing of final_img by trimming it of any empty vertical and horizontal edges
final_img_edge = edge(offset_template, 'prewitt');%find edges of image fit to template 
rv = [];
cv = [];
for ii = 1:size(final_img_edge,1) % rows
    if isempty(find(final_img_edge(ii,:))) % empty horizontal edge exist
        rv = [rv ii];
    end
end
for ii = 1:size(final_img_edge,2) % columns
    if isempty(find(final_img_edge(:,ii))) % empty vertical edge exist
        cv = [cv ii];
    end
end
final_img(rv,:) = [];
final_img(:,cv) = [];

%% Plotting 
if plot_flag
    orig_img = uint8(orig_img);
    subplot(2,2,1)
    imshow(orig_img)
    title('Original')
    axis on
    subplot(2,2,2)
    imshow(template_img)
    title('Template')
    axis on
    subplot(2,2,3)
    imshow(orig_img_rotated)
    title('Original Rotated')
    axis on
    subplot(2,2,4)
    imshow(final_img)
    title('Template with Original (rotated) Texture')
    axis on
end
toc


