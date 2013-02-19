function [orig2 flag incdec offset_template_final] = image_realign(template,orig)
% Use this function to determine whether or not orig image can completely
% overlap the template image.  The template image is smaller than the orig
% image and starts a search at the center of orig image and moves
% circularly outward.  At each step, the amount of overlap between the two
% images is calculated.  If orig image completely overlaps the template
% image, then the search stops.  The output is orig2 which is the texture
% of orig in the shape of the template.  If the search ends with no
% possibility of complete overlap, then orig2 is returned as an empty
% value.
%
% Return cases:
% 1) Complete overlap found.  orig2 is the texture of orig in the shape of
% the template.  flag = 0, indicating an overlap is found.  incdec = 0,
% indicating that one should reduce the orig size and rerun this function
% to see if a tighter overlap can be found.
%
% 2) No overlap found.  orig2 = [] is returned.  flag = 1, indicating no
% overlap is found.  incdec = 1, indicating that one should increase the
% orig size and rerun this function.

if nargin<2
    error('-> Not enough input arguments');
end

if isequal(size(template), size(orig))%the return case 2, indicated above, 
        %because of the assumption that there is no overlap if the images are the same size
    orig2 = [];
    flag = 1;
    incdec = 1;
    offset_template = template;
else
    
    % Make sure all images in uint8 format
    if ~isa(template,'uint8')
        template = uint8(template);
    end
    
    if ~isa(orig,'uint8')
        orig = uint8(orig);
    end
    
    % These parameters make sure template starts in the middle
    % of the orig image and searches by moving circularly outward.  A
    % subset of these parameters also limit the search area to inside the
    % orig image
    [row col] = size(orig);
    [row2 col2] = size(template);
    center_x = floor(col/2);%finding mid point of orig image in x direction
    center_y = floor(row/2);%finding midpoint of orig imag in y direction
    offset_x = center_x - floor(col2/2);%difference between midpoints of the two images in x direction
    offset_y = center_y - floor(row2/2);%difference between midpoints of images in y direction
    
    % This is the core search algorithm where template is moved circularly
    % outward.  At each step, the amount of overlap between template and
    % orig is calculated.  If the orig image completely overlaps the
    % template, then search stops, and the output is the orig image texture
    % in the shape of the template.
    steps = 1;%increasing counter for how much the template moves in the circular motion in each direction
    sgn = 1;%indicates sign of offset to keep moving in a circular motion
    x_y = 0;%used in first if statement in jj loop below
    offset_x2 = offset_x;%amount of offset in the template
    offset_y2 = offset_y;
    x_steps = 0;%indicates the magnituede of how much to move the template
    y_steps = 0;
    flag = 1;%keeps the while loop below alive
    while offset_x2 >= 0 && offset_x2 < (col-col2) && offset_y2 >=0 && offset_y2 < (row-row2) && flag%to keep moving template until
            %there is no more room in the template image to move in either
            %the x or y direction, or there is overlap between the original
            %image and the template
        x_steps = 0;
        y_steps = 0;
        for ii = 1:2
            for jj = 1:steps
                if x_y == 1
                    x_steps = x_steps + 1;
                else
                    y_steps = y_steps + 1;
                end
                offset_template = 255*ones(row,col);%defining white area same size as orig
                offset_template( (1:row2)+offset_y2+sgn*y_steps, (1:col2)+offset_x2+sgn*x_steps) = template;%the new 
                    %template offset from its original position, with a white background
                offset_template = uint8(offset_template);
                
                
                % Determine if the orig fully covers template
                template_tmp = offset_template;
                template_tmp(find(orig > 240)) = 255;%find parts of orig that are basically white and make the corresponding
                    %locations in the temp template white
                    
                if isequal(template_tmp, offset_template)%basically a way to see if the orig image completely overlaps the template
                        %by comparing the offset template and the offset
                        %template turned white in certain parts
                    flag = 0; % orig fully covers template
                    offset_template_final = offset_template;
                end
            end
            x_y = 1 - x_y;%oscillates between 0 and 1 to determine if should go in x or y direction
        end
        offset_x2 = offset_x2 + sgn*x_steps;%determines how much to offset the template
        offset_y2 = offset_y2 + sgn*y_steps;
        sgn = -sgn;
        steps = steps + 1;
    end
    
    if flag == 0%if there is full overlap between the template and orig
        orig2 = orig;
        orig2(find(offset_template_final == 255)) = 255;%cut out shape of template inside orig image
        incdec = 0;
    else%there is no overlap possible, even after offsetting
        orig2 = [];
        offset_template_final = [];
        incdec = 1;
    end
end