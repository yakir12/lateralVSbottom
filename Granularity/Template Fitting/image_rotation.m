function [rot_deg] = image_rotation(orig_img)
% This function takes as input an image, called orig_img, and uses a linear
% regression to determine the angle (rot_deg) at which it must be rotated to
% point the narrow end (i.e., cuttlefish head) to the left

%% Trim the orig_img of any empty vertical and horizontal edges
orig_img_edge = edge(orig_img, 'prewitt');

rv = [];%matrix for the empty rows
cv = [];%matrix for the empty columns
for ii = 1:size(orig_img_edge,1) % rows
    if isempty(find(orig_img_edge(ii,:))) % empty horizontal edge exist
        rv = [rv ii];
    end
end

for ii = 1:size(orig_img_edge,2) % columns
    if isempty(find(orig_img_edge(:,ii))) % empty vertical edge exist
        cv = [cv ii];
    end
end
orig_img(rv,:) = [];
orig_img(:,cv) = [];

%% Find the orig_img boundaries 
orig_img_edge = edge(orig_img, 'prewitt');

% Go through each column and find the border of orig_img_edge
orig_edge_upper = zeros(1,size(orig_img_edge,2));%matrix of indices of upper bound of edges
orig_edge_lower = zeros(1,size(orig_img_edge,2));%matrix of indices lower bound of edges
for ii = 1:size(orig_img_edge,2) 
        orig_edge_upper(ii) = find(orig_img_edge(:,ii),1,'first');
        orig_edge_lower(ii) = find(orig_img_edge(:,ii),1,'last');
end

%% Interlace the upper boundary points with lower boundary points
%Have to go through this step to make sure the edge boundaries are a
%function (meaning there is one unique y for each x) to do the regression
%below

orig_edge_upper = orig_edge_upper(1:2:end);
orig_edge_lower = orig_edge_lower(1:2:end);

if length(orig_edge_upper) > length(orig_edge_lower)
    orig_edge_lower = [orig_edge_lower 0];
elseif length(orig_edge_upper) < length(orig_edge_lower)
    orig_edge_upper = [orig_edge_upper 0];
end
    
%% Perform linear regression to find the rotational angle
orig_boundaries = [orig_edge_upper; orig_edge_lower];
orig_boundaries = reshape(orig_boundaries,2*size(orig_edge_upper,2),1);
b = regress(orig_boundaries,[ones(size(orig_boundaries,1),1) [1:size(orig_boundaries,1)]']);
rot_deg = (180/pi)*atan(b(2));

%% Rotate the image according to rot_deg, trim edges, and determine which end is narrow

% Rotating image
orig_rotated = imrotate(orig_img,rot_deg);%rotate image
orig_rotated(find(orig_rotated==0)) = 255;%change any black parts (which result from the rotation) to be white

% Trim the orig_rotated of any empty vertical and horizontal edges
orig_rotated_edge = edge(orig_rotated, 'prewitt');
rv = [];
cv = [];

for ii = 1:size(orig_rotated_edge,1) % rows
    if isempty(find(orig_rotated_edge(ii,:))) % empty horizontal edge exist
        rv = [rv ii];
    end
end
for ii = 1:size(orig_rotated_edge,2) % columns
    if isempty(find(orig_rotated_edge(:,ii))) % empty vertical edge exist
        cv = [cv ii];
    end
end
orig_rotated(rv,:) = [];
orig_rotated(:,cv) = [];

% Find the orig_img_rotated boundaries 
orig_rotated_edge = edge(orig_rotated, 'prewitt');

% Go through each column and find the border of orig_img_edge
orig_rotated_edge_upper = zeros(1,size(orig_rotated_edge,2));
orig_rotated_edge_lower = zeros(1,size(orig_rotated_edge,2));

for ii = 1:size(orig_rotated_edge,2) 
        orig_rotated_edge_upper(ii) = find(orig_rotated_edge(:,ii),1,'first');
        orig_rotated_edge_lower(ii) = find(orig_rotated_edge(:,ii),1,'last');
end

% Compare the "narrowness" of the boundaries on the left half of the image
% to those on the right half of the image to locate the head
    %If the narrow part (head) is on the right, rotate the image 180
    %degrees so head is on the left
ten_percent = round(.1*size(orig_rotated_edge,2));%To look at ten percent of the right and left parts of the image
left_narrowness = abs(orig_rotated_edge_upper(ten_percent) - orig_rotated_edge_lower(ten_percent));
right_narrowness = abs(orig_rotated_edge_upper(end-ten_percent) - orig_rotated_edge_lower(end-ten_percent));
if left_narrowness > right_narrowness % if the narrow part points to the right
    rot_deg = 180 + rot_deg;
end





