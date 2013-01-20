close all
clear all
clc
%%
th = linspace(-pi,pi,100);
ntheta = 180;
thetad = linspace(-90,90,ntheta);
thetad(end) = [];
% shift90 = sum(cumsum(diff(thetad)) < 90);
thetar = deg2rad(thetad);
costheta = cos(thetar);
sintheta = sin(thetar);
path = '/home/yg32/Documents/PostDoc/Holodeck/Code/data/';
fldr = dir(path);
killd = cellfun(@(x) strcmp(x(1),'.'),{fldr.name});
fldr(killd) = [];
n = length(fldr);
smallpath = '/home/yg32/Documents/PostDoc/Holodeck/Code/processedData/';
%%
% for i = 2%:n
%     mkdir(fullfile(smallpath,fldr(i).name))
%     cmnd = ['mogrify -resize 100x100 -colorspace gray -auto-level -path ',smallpath,...
%         fldr(i).name,'/ ',...
%         fullfile(path,fldr(i).name,'*.jpg')];
%     [~,~] = unix(cmnd);
% end
%%
for i = 1%:n    
    file = dir(fullfile(smallpath,fldr(i).name,'*.jpg'));
    m = length(file);
    I = imread(fullfile(smallpath,fldr(i).name,file(1).name));
    sz = size(I);
    c = floor((fliplr(sz)'+1)/2);
    I = repmat({zeros(sz,'uint8')},1,m);
    for j = 1:m
        I{j} = imread(fullfile(smallpath,fldr(i).name,file(j).name));
    end
    bw = double(cat(3,I{:}));
    varbw = var(bw,[],3);
    nomove = sum(reshape(varbw,prod(sz),1));
    if nomove < 1e6
        d = varbw > 500;
        d = bwperim(d);
        [y,x] = find(d);
        ell = fit_ellipse(x,y);
        cos_phi = cos(ell.phi);
        sin_phi = sin(ell.phi);
        R = [ cos_phi sin_phi; -sin_phi cos_phi ];
        x = ell.X0+ell.a*cos(th);
        y = ell.Y0+ell.b*sin(th);
        xy = R*[x;y];
        d = poly2mask(xy(1,:),xy(2,:),sz(1),sz(2));
        d = repmat(d,[1,1,m]);
    else

    %     bg = std(bw,[],3);
        bg = trimmean(bw,90,3);
        d = (bw-repmat(bg,[1,1,m])).^2;
        killd = zeros(sz);
        killd(round(sz(1)/2),round(sz(2)/2)) = 1;
        killd = bwdist(killd);
        kill = false(sz);
        kill(killd > min(sz)/2*.85) = true;
        kill = repmat(kill,[1,1,m]);
        d(kill) = 0;
        d = reshape(d,sz(1),sz(2)*m);
        d = imclose(d,strel('disk',2));
        d =  d > 2000;%1500;
        d = imfill(d,'holes');
        d = bwareaopen(d,100);
        
        %%
        p = regionprops(d,{'Orientation','PixelIdxList','BoundingBox'});
        idx = cat(1,p.PixelIdxList);
        I2 = zeros(sz(1),sz(2)*m,'uint8');
        temp = cat(2,I{:});
        I2(idx) = temp(idx);        
        J = cellfun(@(x) imcrop(I2,x),{p.BoundingBox},'UniformOutput',false);
        for j = 1:m
            temp = imrotate(J{j},-p(j).Orientation);
            killcol = ~any(temp);
            killrow = ~any(temp,2); 
            temp(killrow,:) = [];
            temp(:,killcol) = [];
            J{j} = temp;
        end
        szs = cellfun(@size,J,'UniformOutput',false);
        szs = cat(1,szs{:});
        szs = max(szs);
        J = cellfun(@(x) padarray(x,szs-size(x),0,'post'),J,'UniformOutput',false);
        imdisp(J)
    end
%         %%
%         J = 
%         pix = {p.PixelList};
%         idx = cellfun(@length,pix);
%         pix = cat(1,pix{:});
%         row = pix(:,1);
%         all = pix(:,2);
%         page = fix((all-1)/sz(2))+1;
%         col = all-(page-1)*sz(2); 
%         ind = sub2ind(sz,row,col);
%         ind = mat2cell(ind,idx,1);
% %         row = mat2cell(row,idx,1);
% %         col = mat2cell(col,idx,1);
%         J = I;
%         for j = 1:m
%             J{j} = imrotate(I{j}(ind{j}),rad2deg(-p(j).Orientation));
%         end
%             
%         dfgsdfggsdfg
%         
%         %%
% 
% 
%     %     label = reshape(repmat(reshape(single(1:m),[1,1,m]),sz),sz(1),sz(2)*m);
%     %     label(~d) = 0;
%     %     ac = regionprops(label,{'Area','Centroid','PixelIdxList'});
%     %     area = cat(1,ac.Area);
%     %     c = cat(1,ac.Centroid);
%     %     row = c(:,2);
%     %     all = c(:,1);
%     % %     page = fix((all-1)/sz(2))+1;
%     % %     col = all-(page-1)*sz(2); 
%     % %     kill = isnan(row);
%     % %     area(kill) = NaN;
%     % %     Xnan = nan(m,3);
%     % %     Xnan(1:max(label(:)),:) = [row,col,area];    
%     % %     row = Xnan(:,1);
%     % %     col = Xnan(:,2);
%     % %     area = Xnan(:,3);
%     %     
%     %     d2 = bwselect(d,all,row,4);
% 
%         d = bwperim(d);
%         [y,x] = find(d);
%         z = fix((x-1)/sz(2))+1;
%         x = x-(z-1)*sz(2); 
%         ab = zeros(m,2);
%         ellid = zeros(m,3);
%         for j = 1:m
%             ell = fit_ellipse(x(z == j),y(z == j));
%             ab(j,:) = [ell.a,ell.b];
%             ellid(j,:) = [ell.phi,ell.X0,ell.Y0];
%         end
%         ab = median(ab);
%         d = reshape(d,[sz(1),sz(2),m]);
%         for j = 1:m
%             cos_phi = cos(ellid(j,1));
%             sin_phi = sin(ellid(j,1));
%             R = [ cos_phi sin_phi; -sin_phi cos_phi ];
%             x = ellid(j,2)+ab(1)*cos(th);
%             y = ellid(j,3)+ab(2)*sin(th);
%             xy = R*[x;y];
%             d(:,:,j) = poly2mask(xy(1,:),xy(2,:),sz(1),sz(2));
%         end
%             
%     %     d2 = reshape(d2,[sz(1),sz(2),m]);
% 
%     end
%     
%     d = reshape(d,[sz(1),sz(2)*m]);
%     d = regionprops(d,'Image');
%     %%
%     for j = 1:m
%         bw2 = imrotate(d(j).Image,rad2deg(-ellid(j,1)));
%         imshow(bw2)
%         sdfsddf
%     end
%     
%     
%     %%
%     vxcvxcvxxc
%     
%     for j = 1:m
% %         subplot(121)
%         imshow(d(:,:,j))
% %         subplot(122)
% %         imshow(d2(:,:,j))
%         drawnow
%         clf
%     end
% 
%         fgdfgdfg
%     %%
%     
%     p = regionprops(d,{'Area','Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation'});
%     
%     for j = 1:m
%         imshow(I{j},[])
%         hold on
%         plot(col(j),row(j),'or')
%         drawnow
%         clf
%     end
%     
%     gdfggfdfgdfg
%     %%
%     p = regionprops(d,{'Area','Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation'});
%     kill = [p.Area] < 200;
%     p(kill) = [];
%     
%     gsdfgdfg
%     
%     d = reshape(d,[sz(1),sz(2),m]);
    
% %     d = abs(bw-repmat(median(bw,3),[1,1,m]));
%     killd = zeros(sz);
%     killd(round(sz(1)/2),round(sz(2)/2)) = 1;
%     killd = bwdist(killd);
%     kill = false(sz);
%     kill(killd > min(sz)/2*.85) = true;
%     kill = repmat(kill,[1,1,m]);
%     d(kill) = 0;
%     tf = d > 20;
%     tf1 = reshape(tf,sz(1),sz(2)*m);
%     tf1 = bwmorph(tf1,'majority');
% %     tf1 = bwmorph(tf1,'open',10);
%     tf1 = bwmorph(tf1,'dilate');
%     tf1 = imfill(tf1,'holes');
%     tf1 = bwareaopen(tf1,80);
%     tf1 = bwmorph(tf1,'remove');
%     ind = find(tf1);
%     [row,col,page] = ind2sub([sz,m],ind);
%     p = zeros(m,5);
%     for j = 1:m
%         el = fit_ellipse(col(page == j),row(page == j));
%         p(j,:) = [el.a,el.b,el.phi,el.X0,el.Y0];
%     end
%     ab = trimmean(p(:,1:2),90,1);
%     a = ab(1);
%     b = ab(2);
%    
% % %     p(:,3) = rad2deg(p(:,3));
%     
%         theta_r         = linspace(0,2*pi);
%     for j = 1:m
%         cos_phi = cos(p(j,3));
%         sin_phi = sin(p(j,3));
%         X0 = p(j,4);
%         Y0 = p(j,5);
%         
%         R = [ cos_phi sin_phi; -sin_phi cos_phi ];
%     
%         % the axes
%         ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
%         horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
%         new_ver_line    = R*ver_line;
%         new_horz_line   = R*horz_line;
% 
%         % the ellipse
%         ellipse_x_r     = X0 + a*cos( theta_r );
%         ellipse_y_r     = Y0 + b*sin( theta_r );
%         rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
%         rotated_ellipse = round(rotated_ellipse);
%         rotated_ellipse = sub2ind(sz,rotated_ellipse(2,:),rotated_ellipse(1,:));
%         im = repmat(I{j},[1,1,3]);
%         im(rotated_ellipse) = 255;
%         
% %         im = imrotate(I{j},-p(j,3));
%         imshow(im)
%         drawnow
%     end
    
%     bw(~tf1) = 0;
% %     
%     figure
%     ax = axes;
%     for j = 1:m        
%         imshow(bw(:,:,j))
%         hold on
%         el = fit_ellipse(col(page == j),row(page == j),ax);
%         waitforbuttonpress
%     end
%     bw = std(double(rgb),[],3);
end
% %     I = mat2cell(bw,sz(1),sz(2)*ones(1,m));
% %     check = textscan(file(1).name(1:end-4),'%s','Delimiter','_');
% %     if strcmp(check{1}{2},'Checkerboard')
% %         for j = 1:m
% %             file(j).name(15) = [];
% %         end
% %     end
% %     t = cellfun(@(x) textscan(x(1:end-4),'%d%s%d%d%d%d%d%d%s%d%d%d','Delimiter','_'),{file.name},'UniformOutput',false); 
% %     t = cat(1,t{:});
% %     row = zeros(m,1);
% %     col = zeros(m,1);
% %     ti = zeros(m,1);
% %     for j = 1:m
% %         switch t{j,9}{:}
% %             case screens{1}
% %                 row(j) = find(checkersz == t{j,10});
% %                 col(j) = find(dists == t{j,11});
% %                 ti(j) = t{j,12};
% %             case screens{2}
% %                 row(j) = 4;
% %                 col(j) = find(dists == t{j,10});
% %                 ti(j) = t{j,11};
% %             case screens{3}
% %                 row(j) = 5;
% %                 col(j) = 1;
% %                 ti(j) = t{j,10};
% %         end
% %     end
% %     
% %     kill = ti ~= 20;
% %     I(kill) = [];
% %     t(kill,:) = [];
% %     row(kill) = [];
% %     col(kill) = [];
% %     ti(kill) = [];
% %     
% %     tiles = repmat({zeros(sz,'uint8')},5,4);
% %     
% %     for j = 1:length(I)
% %         tiles{row(j),col(j)} = I{j};
% %     end
% %       
% %     tiles = cell2mat(tiles);
% %     imwrite(tiles,['results/',t{1,2}{:},'_',num2str(t{1,1}),'.jpg'])
% %     
% %     
% % %     [trow,tcol] = ind2sub([5,5],ti);
% % %     row = (row-1)*5+trow;
% % %     col = (col-1)*5+tcol;
% % %     tiles = repmat({zeros(sz,'uint8')},5*5,5*4);
% % %     
% % %     for j = 1:m
% % %         tiles{row(j),col(j)} = I{j};
% % %     end
% % %       
% % %     tiles = cell2mat(tiles);
% % %     imwrite(tiles,[t{1,2}{:},'.jpg'])
% %     
% end
% %%
