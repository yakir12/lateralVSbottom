%%%%%%%%% Filtering Images Fit to Oval Template %%%%%%%%%%
%%%Use this code to filter the images (already standardized by fitting them to an
%%%oval template) with 6 bandpass filters

%% Setting the stage
clc
clear all
close all
cd('/Users/dtaniguchi/Documents/Cuttlefish Project/Example Images, Quality Controlled')%setting the directory

set(0,'defaultaxesfontsize',30);
set(0,'defaulttextfontsize',30);
set(0,'defaultlinelinewidth',2);

plot_filt_flag = 0;%set to 1 if want to plot normalized fitlered images, otherwise set to 0
plot_energy_flag = 1; %set to 1 if want to plot the energy of the 6 filtered images, 0 o.w.
num_im = 10; %Number of images to filter in each subfolder
oct = 128;% Reciprocal of min width of filter frequency

%% Loading the template and processing (for use in creating filter)
template_img = imread('Oval_Template2_Trimmed.png'); %Oval template%<-------------------MAKE SURE HAVE CORRECT TEMPLATE

% Convert image to grayscale if colored
if size(template_img,3)==3
    template_img=rgb2gray(template_img);
end

% Make double floating point
if isa(template_img,'uint8')
    template_img = double(template_img);
end

[template_row template_col] = size(template_img);
dim_high = max(template_row, template_col)+1;% Have to add 1 because template slightly smaller than standardized images<------

% Finding edges
img_edge_temp = edge(template_img,'prewitt');%Edges for use in finding energy later


% Padding template, to find edges for use in loop
if template_row < template_col 
    img_edge = [zeros(floor((template_col-template_row)/2),template_col); img_edge_temp; zeros(ceil((template_col-template_row)/2),template_col)];
else
    img_edge = [zeros(template_row,floor((template_row-template_col)/2)); img_edge_temp; zeros(template_row,ceil((template_row-template_col)/2))];
end

% Adding in extra row and col because template slightly smaller than fit images <-----------
[edge_row edge_col] = size(img_edge);
img_edge = [zeros(1,edge_col); img_edge];
img_edge = [zeros(edge_row + 1, 1) img_edge];

%% Creating filters
[f1 f2] = freqspace(dim_high, 'meshgrid');%creates square matrix between -1 and 1 with length and width of dim_high
%Have to add 1 above because template is slightly different size than image

r = sqrt(f1.^2 + f2.^2);%finds distance of each point from center

%Defining matrices for filter heights
for m=1:6
    eval(['Hd',num2str(m),'=ones(dim_high);'])
end

% Defining frequency ranges for each filter
% Frequency ranges are on an octave scale (meaning differ by 2^n)
Hd1((r<=0)|(r>1/oct)) = 0;
Hd2((r<=1/oct) | (r>3/oct)) = 0;
Hd3((r<=3/oct) | (r>7/oct)) = 0;
Hd4((r<=7/oct) | (r>15/oct)) = 0;
Hd5((r<=15/oct)| (r>31/oct)) = 0;
Hd6((r<=31/oct)| (r>63/oct)) = 0;

% Creating the window for the filter
win = fspecial('gaussian',dim_high,0.05*dim_high);%create predefined 2-d filter,
%of size dim_high and with standard deviation of 20% of dim_high
win = win ./ max(win(:));  % Make the maximum window value be 1.

%Creating the filters and fft'ing them
for n = 1:6
    eval(['filt',num2str(n),' = fwind2(Hd',num2str(n),',win);'])
    eval(['H',num2str(n),' = fft2(filt',num2str(n),');'])%fft the filter

end

%% Load images and create matrices
directory = dir;
for k=24%:size(directory,1);% Goes through each subfolder
    if ~isempty(strfind(directory(k).name,'Best'))%images I'm using are in folders that end in "Best"
        cd(directory(k).name)
        directory2 = dir;
        
        eval(['load Standardized_Oval_Images_',num2str(k)])

        energy_orig_mx = NaN(num_im,1);%energy in orginal images inside the template
        energy_orig_sq_mx = NaN(num_im,1);% energy with padding that makes image square
        energy_mx = NaN(num_im,6) ;%Filtered energy values just inside template (does NOT include backgd or padding).  Each row is an image, each column is from one filter
        energy_sq_mx = NaN(num_im,6);% Matrix of energy of filtered and unpadded/square images (includes background)
        energy_norm_mx = NaN(num_im,6);% Matrix of energy of filtered, trimmed, and normalized images
        energy_norm_mean_mx = NaN(num_im,6); %Matrix of mean pixel value of normalized images
        
        %% Loop going through each image
        
        %An assumption is made that the images have already been fitted to a template.
        for i=1:num_im
            
            %Read in and procces image
            eval(['orig_img = Cuttlefish_fit',num2str(i),';'])
            if size(orig_img,3)==3%changing to gray scale
                orig_img=rgb2gray(orig_img);
            end
            
            if isa(orig_img,'uint8')
                orig_img = double(orig_img);%change to double
            end
            
            % Finding the non-white parts of the background, averaging them, and replacing white parts with
            %average
            non_w=orig_img(find(orig_img<=200));%Finding (bascially) non-white values<-----------MAY NEED TO CHANGE TO 255
            mean_non_w=mean(non_w);
    
            energy_orig_mx(i,:) = sum(non_w.^2);% Finding energy inside template
                       
            orig_img_bkgd = orig_img;
            orig_img_bkgd(orig_img>200)=mean_non_w;%<-----------MAY NEED TO CHANGE TO 255

   %Making image square by padding with average value
   [row_img col_img] = size(orig_img);
    if row_img < col_img
        img_pad = [mean_non_w.*ones(floor((col_img-row_img)/2),col_img); orig_img_bkgd; mean_non_w.*ones(ceil((col_img-row_img)/2),col_img)];
    else
        img_pad = [mean_non_w.*ones(row_img,floor((row_img-col_img)/2)); orig_img_bkgd; mean_non_w.*ones(row_img,ceil((row_img-col_img)/2))];
    end                        
            
    % Energy of padded, square image
    energy_orig_sq_mx(i,1) = sum(img_pad(:).^2);
    
            %% Filtering the image, trimming, finding the energy, and normalizing (for plotting purposes)
            
            for mm = 1:6
         %Filtering image
        img_pad_freq = fft2(img_pad); %fft the image
        eval(['filt_out_freq = H',num2str(mm),'.*img_pad_freq;'])%filtering
        eval(['filt_out',num2str(mm),' = fftshift(ifft2(filt_out_freq));']) %inverse fft and putting mirror images back together
        eval(['filt_out_temp = filt_out',num2str(mm),';']) %To make easier to work with value
        
        energy_sq_mx(i,mm) = sum(filt_out_temp(:).^2);
        filt_out_temp2 = filt_out_temp;% for filtered image just inside template
        
                
                % Setting values outside template to zero
                for ii = 1:size(img_edge,2) % Going through columns
                    edge_upper = find(img_edge(:,ii),1,'first');
                    edge_lower = find(img_edge(:,ii),1,'last');
                    
                    for jj = 1:size(img_edge,1) % Going through rows
                        if ~isempty(edge_upper) && ~isempty(edge_lower)
                            if jj < edge_upper || jj > edge_lower %Outside the shape of the animal
                               filt_out_temp2(jj,ii) = 0;
                            end
                        end
                    end
                end
                
                % Finding energy
                energy_mx(i,mm) = sum(filt_out_temp2(:).^2);                
              
                % Normalizing filtered image (so can see better when plotting)
                norm_filt_temp = (filt_out_temp2 - min(filt_out_temp2(:)));% shifting the filtered values so they are all positive, to help with plotting in grayscale
                norm_filt_img = (norm_filt_temp./max(norm_filt_temp(:))).*255; % scaling values between 0 and 1, and then multiplying by 255 (grayscale is between 0 and 255)
                energy_norm_mx(i,mm) = sum(norm_filt_img(:).^2)/(size(norm_filt_img,1)*size(norm_filt_img,2)*255);%energy of entire normalized image (including background)
                energy_norm_mean_mx(i,mm) = mean(norm_filt_img(:))/(size(norm_filt_img,1)*size(norm_filt_img,2)*255); %mean pixel value, normalized by max value
                
                if plot_filt_flag
                    if mm == 1; 
                        figure(i*2-1); 
                        subplot(3,3,1)
                        imshow(uint8(orig_img))
                        
                        figure(i*2)
                        subplot(3,3,1)
                        imshow(uint8(orig_img))
                    end
                    
                    figure(i*2-1)
                    subplot(3,3,mm+1)                    
                    
                    imshow(uint8(filt_out_temp));
                    eval(['title(''Filter ',num2str(mm),''')'])
                    
                    figure(i*2)
                    subplot(3,3,mm+1)
                    imshow(uint8(norm_filt_img))
                    eval(['title(''Filter ',num2str(mm),''')'])
                end
                
                %clear filt_out_temp 
                
            end
            %clear filt_out_freq
            eval(['save Filtered_Info_',num2str(k),' energy_* filt_out*'])
            
        end
        %% Plotting the energy
        
        if plot_energy_flag
            figure
            plot((1:6),energy_mx)
            ylabel('energy/(specific max energy)')
            xlabel('filter number')
            title('Granularity Spectrum,Inside Template')

            figure
            plot((1:6),energy_mx./repmat(energy_orig_mx,1,6));%normalzing by max value (where image has all 255 values for each pixel)
            xlabel('filter number')
            ylabel('energy/(original energy)')
            title('Granularity Spectrum, Inside Template Divided by Original Energy')
                        
            figure
            plot((1:6),energy_sq_mx);
            xlabel('filter number')
            ylabel('energy of filtered image')
            title('Granularity Spectrum, Square Image')
            
            figure
            plot(1:6,energy_sq_mx./repmat(energy_orig_sq_mx,1,6))
            xlabel('filter number')
            ylabel('energy/(original energy)')
            title('Granularity Spectrum, Square Image Normalized by Original Energy')           
            
            figure
            plot(1:6,energy_norm_mx)
            xlabel('filter number')
            ylabel('normalized energy/max energy')
            title('Granularity Spectrum, Normalized to Make Good Picture, Includes Background')
            
            figure
            plot(1:6, energy_norm_mean_mx)
            xlabel('filter number')
            ylabel('normalized energy/max energy')
            title('Granularity Spectrum, Normalized to Make Good Picture, Includes Background')

            
        end
        
        %% Clearing matrices 
        %clear energy_trim_mx energy_norm_mx energy_mx
    end
end







