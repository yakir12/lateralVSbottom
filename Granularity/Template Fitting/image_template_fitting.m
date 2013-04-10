%%%%%%% Fitting Images to a Symmetric Template %%%%%%%%%%%%

%%%Use this code to fit the segmented (cut out of from their
%%%background) cuttlefish to a symmetric template (an oval)

%%% These images are selected from the folder of images that have been
%%% quality controlled

%%%This code is different than the fitting to the template for the
%%%Eigenfaces method because the template is symmetric rather than in the
%%%shape of the cuttlefish

%% Setting the stage
clc
clear all
close all
cd('/Users/dtaniguchi/Documents/Cuttlefish Project/Example Images, Quality Controlled')%setting the directory
%cd('/Users/dtaniguchi/Documents/Cuttlefish Project/Template Fitting')
directory = dir;

set(0,'defaultaxesfontsize',30);
set(0,'defaulttextfontsize',30);
set(0,'defaultlinelinewidth',2);

%% Reading in template and finding edges 
%     %Have to use the template because the images themselves don't always
%     %have complete edges (i.e., there are often discontinuities)
template_img = imread('Oval_Template2_Trimmed.png'); %Oval template%<-------------------MAKE SURE HAVE CORRECT TEMPLATE

% Convert image to grayscale if colored
if size(template_img,3)==3
    template_img=rgb2gray(template_img);
end
% 
% %img_edge = edge(oval_template,'prewitt');

%% Fitting the quality controlled images to the template
%Going through each subfoldfolder in quality controlled folder
for i=24%:size(directory,1);% Goes through each subfolder
    if ~isempty(strfind(directory(i).name,'Best'))%images I'm using are in folders that end in "Best"
        cd(directory(i).name)
        directory2 = dir;
        
        im_num = size(dir,1)-4;% number of images in the folder
        
        for j=1:im_num
            % Fitting each image to template
                %Images are the first argement 
                %Template, the second argument, is "Oval_Template2_Trimmed.png"
                %Last argument of image_retexturizing_funcition determines if want to plot fit image (1) or not (0)
            P = path;
            path(P, '/Users/dtaniguchi/Documents/Cuttlefish Project/Template Fitting');
            tic
            eval(['Cuttlefish_fit',num2str(j),' = image_retexturizing_function(directory2(',num2str(j+2),').name,template_img,0);'])
            %eval(['strcat({''fit_''},directory2(',num2str(j+3),').name) = image_retexturizing_function(directory2(',num2str(j+3),').name,template_img,0);'])
            toc
            % Save standardized images in mat file
            eval(['save Standardized_Oval_Images_',num2str(i),' Cuttlefish_fit* template_img'])

        end
        path(P);
        
        %Return to original directory
        cd('/Users/dtaniguchi/Documents/Cuttlefish Project/Example Images, Quality Controlled')%setting the directory
    end
end
%% Saving standardized images in training set folder
%save Standardized_Oval_Images Cuttlefish_fit* template_img
