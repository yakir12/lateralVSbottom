close all
clear all
clc
%%
sz = [100,100];
ab = [419,239];
maa = 45;
mia = 25;
maa = maa/sz(1);
mia = mia/sz(1);

stationary = dir('meanStationary/*.png');
stationary = cellfun(@(x) x(1:end-4),{stationary.name},'uni',0);

thdif= 1;
path1 = 'data';
path2 = 'processedata';
path3 = '2processedata';
fldr = dir(path2);
killd = cellfun(@(x) strcmp(x(1),'.'),{fldr.name});
fldr(killd) = [];
% killd = cellfun(@(x) any(strcmp(x,stationary)),{fldr.name});
% fldr(killd) = [];
n = length(fldr);
%%
for i = 1:n
    rmdir(fullfile(path3,fldr(i).name),'s');
    mkdir(fullfile(path3,fldr(i).name));
end

matlabpool(4)

parfor i = 1:n        
    files = dir(fullfile(path2,fldr(i).name,'*.png')); 
    I = imread(fullfile(path1,fldr(i).name,[files(1).name(1:end-3),'jpg']));
    sz1 = size(I);
    sz1(3) = [];
    c0 = fliplr(sz1)'/2;
    m = length(files);  
    I = zeros([sz,m],'uint8');
    for j = 1:m
        I(:,:,j) = imread(fullfile(path2,fldr(i).name,files(j).name));
    end
    II = reshape(I,[sz(1),sz(2)*m]);
    if any(strcmp(fldr(i).name,stationary))
        d = imread(fullfile('meanStationary',[fldr(i).name,'.png']));
        ac = regionprops(d,{'centroid','orientation'});
        both = ac.Centroid;
        row = both(2);
        col = both(1);        
        o = cat(1,ac.Orientation);
        
        c = [col;row]-c0;
        r = norm(c);
        theta1 = atan2(c(2),c(1));
        theta2 = theta1+deg2rad(o);
        col2 = r*cos(theta2)+c0(1);
        row2 = r*sin(theta2)+c0(2);
        rect = round([[col2,row2]-ab/2,ab]);
        
        for j = 1:m
            I = imadjust(rgb2gray(imread(fullfile(path1,fldr(i).name,[files(j).name(1:end-3),'jpg']))));            
            I2 = imrotate(I,-o,'crop');
            I3 = imcrop(I2,rect);
            if any(size(I3) < rect(3:4))
                I3 = padarray(I3,fliplr(rect(3:4)+1)-size(I3),0,'post');
            end
            imwrite(I3,fullfile(path3,fldr(i).name,files(j).name))
        end       
        %%
    else
    %%
        d = imabsdiff(II(:,1:sz(2)*(m-1)),II(:,sz(2)+1:end));
        d = reshape(d,[prod(sz),m-1]);    
        s = sum(d);
        [~,ind] = sort(s,'descend');
        keep = d(:,ind) < 30;%5;%thdif;
        keep1 = cumsum(keep,2);
        keep2 = cumsum(keep1,2);
        keep2(keep2 > 1) = 0;
        I1 = reshape(I(:,:,1:m-1),prod(sz),m-1);
        I1 = I1(:,ind);
        I1(~keep2) = 0;
        bg = sum(I1,2,'native');
        bg = reshape(bg,sz);    
        
%     imshow(bg,[])
%%
        BG = repmat(bg,[1,m]);
        d = imabsdiff(BG,II);
    %     d(repmat(isnan(bg),[1,m])) = 0;
        d = imadjust(d) > 128;
        d = imdilate(d,strel('disk',2));
    %     d = imfill(d,'holes');
        d = bwareaopen(d,100);
        d = bwconvhull(d,'objects');

        label = reshape(repmat(reshape(single(1:m),[1,1,m]),sz),sz(1),sz(2)*m);
        label(~d) = 0;
        ac = regionprops(label,{'centroid','orientation'});
        c = cat(1,ac.Centroid);
        o = cat(1,ac.Orientation);
        row = c(:,2);
        both = c(:,1);
        page = fix((both-1)/sz(2))+1;
        col = both-(page-1)*sz(2); 
        kill = isnan(row);
        o(kill) = NaN;

        Xnan = nan(m,3);
        Xnan(1:max(label(:)),:) = [row,col,o];    
        row = Xnan(:,1);
        col = Xnan(:,2);
        o = Xnan(:,3);   
        
        %%
        fact = mean(sz1./sz);
        c = fact*[col';row']-repmat(c0,1,m);
        r = sqrt(sum(c.^2));
        theta1 = atan2(c(2,:),c(1,:));
        theta2 = theta1+deg2rad(o');
        col2 = r.*cos(theta2)+c0(1);
        row2 = r.*sin(theta2)+c0(2);
        AB = repmat(ab,m,1);
        rect = round([[col2',row2']-AB/2,AB]);
    %%
        for j = 1:m
            if ~isnan(row(j))
                I = imadjust(rgb2gray(imread(fullfile(path1,fldr(i).name,[files(j).name(1:end-3),'jpg']))));
                I2 = imrotate(I,-o(j),'crop');
                I3 = imcrop(I2,rect(j,:));
                if any(size(I3) < fliplr(rect(j,3:4)+1))
                    I3 = padarray(I3,fliplr(rect(j,3:4)+1)-size(I3),0,'post');
                end
                imwrite(I3,fullfile(path3,fldr(i).name,files(j).name))
            end
        end    
    end
end
matlabpool close