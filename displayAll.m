close all
clear all
clc
%%
sz = [7,4];
sz2 = [5,4];
path = '2processedata';
path2 = '32processedata';
fldr = dir(path);
killd = cellfun(@(x) strcmp(x(1),'.'),{fldr.name});
fldr(killd) = [];
n = length(fldr);
fid = fopen('table.csv','r');
% table = textscan(fid,'%s','Delimiter',','); %s %u %s %u %u %u %u %u %u %s %u %u
head = textscan(fid,repmat('%s',1,12),1,'Delimiter',','); 
table = textscan(fid,'%s %f %s %f %f %f %f %f %f %s %f %f','Delimiter',','); 
fclose(fid);
table{1} = cellfun(@(x) x(1:end-4),table{1},'uni',0);
%%
rows = sort(unique(table{10}));
cols = sort(unique(table{11}));
cols(isnan(cols)) = [];
dt = cat(2,table{4:9});
matlabpool(4)
parfor i = 1:n 
    files = dir(fullfile(path,fldr(i).name,'*.png')); 
    m = length(files);
    ind = cellfun(@(x) find(strcmp(x(1:end-4),table{1})),{files.name});
    cf = table{2}(ind(1));
    bknd = table{3}{ind(1)};
    dt0 = dt(ind(1),:);    
    dt0 = datestr(dt0,'mmm dd');
    imtitle = sprintf('Cuttlefish %u on background %s %s',cf,bknd,dt0);
    I = repmat({repmat({128*ones(240,420,'uint8')},sz)},sz2);
    for j = 1:m
        row = find(strcmp(table{10}{ind(j)},rows));
        col = find(table{11}(ind(j)) == cols);
        if isempty(col) 
            col = 1;
        end
        I{row,col}{table{12}(ind(j))} = imread(fullfile(path,fldr(i).name,files(j).name));
    end
    
    for j = 1:prod(sz2)
        I{j} = padarray(cell2mat(I{j}),[10,10],uint8(255));
    end
        
    I = cell2mat(I);
    imwrite(I,fullfile(path2,[imtitle,'.png']))
end
matlabpool close