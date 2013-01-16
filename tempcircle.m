close all
clear all
clc
%%
path = '/home/yg32/Documents/PostDoc/Holodeck/Code/data/';
fldr = dir(path);
kill = cellfun(@(x) strcmp(x(1),'.'),{fldr.name});
fldr(kill) = [];
n = length(fldr);
for i = 1%:n
    file = dir(fullfile(path,fldr(i).name,'*.jpg'));
    I = imread(fullfile(path,fldr(i).name,file(1).name));
    imshow(I)
    h = imellipse;
    wait(h);
    I = imcrop(I,h.getPosition);
    imshow(I)
    
    fid = fopen(fullfile(path,fldr(i).name,'details.txt'),'a');
    fprintf(fid,...
        [''],...
        h.getPosition);
    fclose(fid);
    
end