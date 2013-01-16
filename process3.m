close all
clear all
clc
%%
% sz = [100,100];
% path = '/home/yg32/Documents/PostDoc/Holodeck/Code/processedata/';
sz = [100,98];
path = '/home/yg32/Documents/PostDoc/Holodeck/Code/processedata/';
fldr = dir(path);
killd = cellfun(@(x) strcmp(x(1),'.'),{fldr.name});
fldr(killd) = [];
n = 1;%length(fldr);
%%
for i = 1:n
    im = dir(fullfile(path,fldr(i).name,'*.jpg')); 
% cuttlefish bkgd year month day hour minute second screen replicate 
    t = cellfun(@(x) sscanf(strrep(x,'_',' '),'%*u %*s %*u %u %u %u %u %u %u %*s'),{im.name},'UniformOutput',false);
    t = cell2mat(t)';
    t = datenum(t);
    [~,ind] = sort(t);
    im = im(ind);
    m = length(im);
    I = zeros([sz,m],'uint8');
    for j = 1:m
        temp = imread(fullfile(path,fldr(i).name,im(j).name));
        I(:,:,j) = temp;%(:,:,2);
    end
%     I = reshape(I,prod(sz),m);
%     fid = fopen('kaka.txt');
%     fprintf(fid,repmat([repmat('%u ',1,m),'\n'],1,prod(sz)),I);
%     fclose(fid);

%     J = single(I);
%     kill = ~any(I,3);
%     bg = nan(sz);
%     bg(kill) = 0;
%     
%     for j = 1:m-1
%         d = (J(:,:,j)-J(:,:,j+1)).^2;
%         d(kill) = nan;
%         d = round(my_normab(d,0,99));
%         mask = d > 0;
%         tempI = I(:,:,j);
%         bg(mask) = tempI(mask);
%         if sum(isnan(bg)) == 0
%             break
%         end
%     end    
end
%%
for i = 40:sz(1)
    for j = 10:sz(2)
        plot(double(squeeze(I(i,j,:))),'ok'), ylim([0,255])        
        drawnow
        pause(.1)
    end
end
%%
sub2ind(sz,i,j)

% j
% 
% imagesc(bg)