%%%%%% Fast Fourier Transform and Filter %%%%%%%%%%%%%%%%%%

%%%%%%Use this code to create 6 octave-wide bandpass Gaussian filters and convolve
%%%%%%with standardized cuttlefish images to filter the images.  The energy from the filtered
%%%%%%images are then calculated and can be plotted.  

%%%Before can run this code, must also have run image_retexturizing.m to
%%%fit the image to a template

%% Setting the stage
clc
clear all
close all
cd('/Users/dtaniguchi/Documents/Cuttlefish Project/lateralVSbottom/Granularity')%setting the directory

set(0,'defaultaxesfontsize',30);
set(0,'defaulttextfontsize',30);
set(0,'defaultlinelinewidth',2);

plot_filt_flag = 1;%set to 1 if want to plot normalized fitlered images, otherwise set to 0
plot_energy_flag = 1; %set to 1 if want to plot the energy of the 6 filtered images, 0 o.w.

%% Read in and process image
%An assumption is made that the images have already been fitted to a template.  

orig_img = imread('Cuttlefish3.jpg');
if size(orig_img,3)==3%changing to gray scale
    orig_img=rgb2gray(orig_img);
end

% Finding the non-white parts of the background, averaging them, and replacing white parts with
%average
non_w=orig_img(find(orig_img<=200));%Finding (bascially) non-white values
mean_non_w=mean(non_w);
orig_img_bkgd = orig_img;
orig_img_bkgd(orig_img>200)=mean_non_w;

% Padding the image on all sizes by half the length of the longest
% dimension, to help with convolving below
[row_img col_img] = size(orig_img);
dim_high = max(row_img, col_img);
pad_img = [mean_non_w*(ones(floor((2*dim_high-row_img)/2),col_img));orig_img_bkgd;...
    mean_non_w*(ones(ceil((2*dim_high-row_img)/2),col_img))];%padding top and bottom

pad_img = [mean_non_w*ones(2*dim_high,floor((2*dim_high-col_img)/2)) pad_img ...
    mean_non_w*ones(2*dim_high,ceil((2*dim_high-col_img)/2))];%padding left and right

pad_img = double(pad_img);%changing to type double for use in convolution

% Finding edges of template for trimming after filtering
    %Using the template for edges rather than the image itself because
    %sometimes there are discontinuities in the edges of the image
template_img = imread('cuttlefish_template2.tif');
if size(template_img,3)==3
    template_img=rgb2gray(template_img);
end
img_edge = edge(template_img,'prewitt');


%% Creating the filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f1 f2] = freqspace(dim_high, 'meshgrid');%creates square matrix between -1 and 1 with length and width of dim_high
r = sqrt(f1.^2 + f2.^2);%finds distance of each point from center

%Defining matrices for filter heights
for m=1:6
    eval(['Hd',num2str(m),'=ones(dim_high);'])
end

% Defining frequency ranges for each filter
    % Frequency ranges are on an octave scale (meaning differ by 2^n)
Hd1((r<=0)|(r>2/64)) = 0;
Hd2((r<=2/64) | (r>4/64)) = 0;
Hd3((r<=4/64) | (r>8/64)) = 0;
Hd4((r<=8/64) | (r>16/64)) = 0;
Hd5((r<=16/64)| (r>32/64)) = 0;
Hd6((r<=32/64)| (r>1)) = 0;

% Creating the window for the filter
win = fspecial('gaussian',dim_high,0.1*dim_high);%create predefined 2-d filter, 
    %of size dim_high and with standard deviation of 10% of dim_high
win = win ./ max(win(:));  % Make the maximum window value be 1.

%Creating the filters
for n = 1:6
    eval(['filt',num2str(n),' = fwind2(Hd',num2str(n),',win);'])
end

%% Filtering the image, trimming, finding the energy, and normalizing (for plotting purposes)
energy_vec = NaN(1,6);% vector of energy values for filtered and trimmed images (only values in the template)
energy_trim_vec = NaN(1,6);% vector of energy of filtered and unpadded images (includes background)
energy_norm_vec = NaN(1,6);% vector of energy of filtered, trimmed, and normalized images

for mm = 1:6
    % Filtering
    eval(['filt_img = filter2(filt',num2str(mm),',pad_img,''same'');']) 
    
    % Trimming the padding
    trim_img_temp = filt_img(floor((2*dim_high-row_img)/2)+1:floor((2*dim_high-row_img)/2)+row_img,...
            floor((2*dim_high-col_img)/2)+1:floor((2*dim_high-col_img)/2)+col_img); % Trimming the padding
    energy_trim_vec(1,mm) = sum(sum(trim_img_temp.^2));
    trim_img = trim_img_temp;
        
    % Setting values outside template to zero
    for ii = 1:size(img_edge,2) % Going through columns
        edge_upper = find(img_edge(:,ii),1,'first');
        edge_lower = find(img_edge(:,ii),1,'last');
        
        for jj = 1:size(img_edge,1) % Going through rows
            if ~isempty(edge_upper) && ~isempty(edge_lower)
                if jj < edge_upper || jj > edge_lower
                    trim_img(jj,ii) = 0;
                end              
            end
        end
    end
    
    % Finding energy
    energy_vec(1,mm) = sum(sum(trim_img.^2));
    
    % Normalizing filtered image (so can see better when plotting)
    if plot_filt_flag
        if mm == 1; figure; end
        norm_filt_temp = (trim_img_temp - min(trim_img_temp(:)));% shifting the filtered values so they are all positive, to help with plotting in grayscale
        norm_filt_img = (norm_filt_temp./max(norm_filt_temp(:))).*255; % scaling values between 0 and 1, and then multiplying by 255 (grayscale is between 0 and 255)
        energy_norm_vec(1,mm) = sum(sum(norm_filt_img.^2));%energy of entire normalized image (including background)
        subplot(2,3,mm)
        imshow(uint8(norm_filt_img))
        eval(['title(''Filter ',num2str(mm),''')'])
    end
    
    clear filt_img trim_img norm_filt_temp norm_filt_img trim_img_temp
    
end

%% Plotting the energy 
% Plotting energy
if plot_energy_flag
    figure
    subplot(2,1,1)
    plot((1:6),energy_vec);
    xlabel('filter number')
    ylabel('energy of filtered image')
    title('Granularity Spectrum, Unnormalized, No Background')

    subplot(2,1,2)
    plot((1:6),energy_trim_vec./max(energy_trim_vec))
    xlabel('filter number')
    ylabel('energy/(max energy)')
    title('Granularity Spectrum, Normalized by Max, Includes Background')
end

%% Straight convolution method: Finding 2-D Guassian matrix (the ideal, isotropic filter) for convolving the image
% sigma = 2;%standard deviation of Gaussian kernel
% gaus_kern = NaN(dim_high,dim_high);%matrix for the Gaussian kernel for convolution
% for i = 1:dim_high
%     for j = 1:dim_high
%         gaus_kern(i,j) = (1/2/pi/sigma^2)*exp(-((floor(dim_high/2)-i).^2+(floor(dim_high/2)-j).^2)/2/sigma^2);
%     end
% end
% 
% %Filtering images and plotting
% 
% filt_img = conv2(pad_img, gaus_kern,'same');
%  
% figure
% imshow(uint8(filt_img))
% 


