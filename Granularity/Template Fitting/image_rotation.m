function [rot_deg] = image_rotation(orig_img)
% This function takes as input an image, called orig_img, and rotates it in
% 5 degree increments until it is in the most horizontal position (i.e.,
% width/height ratio is the largest.  

wh_ratio = zeros(1,36); % width/height ratio

for i = 0:35 % To go through 180 degrees
    orig_img2 = imrotate(orig_img,5*i); % rotation by 5 degrees
    %orig_img2(find(orig_img2 > 200)) = 0; %Changes almost white values to black
    orig_img2(find(orig_img2 == 0)) = 255; %Change all the black to white because get black filled in when rotate
    
%% Trimming the orig_img2 of any empty vertical and horizontal edges
    orig_img_edge = edge(orig_img2, 'prewitt');
    rv = [];
    cv = [];
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
    orig_img2(rv,:) = [];
    orig_img2(:,cv) = [];
    
    wh_ratio(i+1) = size(orig_img2,2)/size(orig_img2,1);
end

%% Finding max ratio and use that for the degree of rotation
[val indx] = max(wh_ratio); 
rot_deg = 5*(indx - 1);

%% Rotating the image according to rot_deg and determine which end is narrow
orig_rotated = imrotate(orig_img,rot_deg);
orig_rotated(find(orig_rotated==0)) = 255;

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
orig_rotated_edge(rv,:) = [];
orig_rotated_edge(:,cv) = [];

% Go through each column and find the border of orig_img_edge
orig_rotated_edge_upper = zeros(1,size(orig_rotated_edge,2));
orig_rotated_edge_lower = zeros(1,size(orig_rotated_edge,2));
for ii = 1:size(orig_rotated_edge,2) 
        orig_rotated_edge_upper(ii) = find(orig_rotated_edge(:,ii),1,'first');
        orig_rotated_edge_lower(ii) = find(orig_rotated_edge(:,ii),1,'last');
end

% Compare the "narrowness" of the boundaries on the left half of the image
% to those on the right half of the image
ten_percent = round(.1*size(orig_rotated_edge,2));
left_narrowness = abs(orig_rotated_edge_upper(ten_percent) - orig_rotated_edge_lower(ten_percent));
right_narrowness = abs(orig_rotated_edge_upper(end-ten_percent) - orig_rotated_edge_lower(end-ten_percent));
if left_narrowness > right_narrowness % if the narrow part points to the right
    rot_deg = 180 + rot_deg;
end

