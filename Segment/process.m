close all
clear all
clc
%%
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
% for i = 1%:n
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
    d = reshape(d,[sz(1),sz(2),m]);
    for j = 1:m
        im = d(:,:,j) > 1500;
% %         im = phantom([1,.7,.3,-.4,.4,33],200); sz = size(im); c = floor((fliplr(sz)'+1)/2);
        [r,x] = radon(im,thetad);
        h = max(r);
        
        [~,col] = max(h);
        cfang = thetar(col);
        
%         h = zscore(h);
%         f = fit(thetar(:),h(:),'gauss1');
%         [~,col] = min(abs(thetar-f.b1));
%         cfang = f.b1;
        
% % %         f = fit(thetar(:),h(:),'sin1');
% % %         theta1 = (pi/2-f.c1)/f.b1;
% % %         [~,col] = min(abs(thetar-theta1));

%         f = fit(x(:),r(:,col),'gauss1');
%         cfwi = 2*sqrt(2*log(2))*f.c1;
%         cfcen1 = f.b1;
        
        [~,row] = max(r(:,col));
        cfcen1 = x(row);
        
%         A2 = [costheta(col),sintheta(col);-sintheta(col),costheta(col)];
        
%         ax = [c(2),c(2);0,sz(1)];
%         ax1 = ax-repmat(c,1,2)+repmat([cfcen1;0],1,2);
%         ax2 = A2*ax1;
%         ax2 = ax2+repmat(c,1,2);        
%         imagesc(im) 
%         axis image
%         hold on
%         plot(ax(1,:),ax(2,:),'g')
%         plot(ax2(1,:),ax2(2,:),'y')
%         drawnow
%         waitforbuttonpress
    
        [~,col2] = min(abs(thetar-wrapToPi(thetar(col)+pi/2)));
        [~,row] = max(r(:,col2));
        cfcen2 = x(row);

%         f = fit(x(:),r(:,col2),'gauss1');
%         cfcen2 = f.b1;
%         cflen = 2*sqrt(2*log(2))*f.c1;
        
%         [~,row] = min(abs(x-f.b1));
%         A2 = [costheta(col2),sintheta(col2);-sintheta(col2),costheta(col2)];
%         
%         ax = [0,sz(2);c(1),c(1)];
%         ax1 = ax-repmat(c,1,2)+repmat([x(row);0],1,2);
%         ax2 = A2*ax1;
%         ax2 = ax2+repmat(c,1,2);      
%         
%         figure
%         imagesc(im) 
%         axis image
%         hold on
%         plot(ax(1,:),ax(2,:),'g')
%         plot(ax2(1,:),ax2(2,:),'y')
%         drawnow
        imagesc(im) 
        axis image
        set(gca,'ydir','normal')
        hold on
        p1 = c+[cos(-cfang);sin(-cfang)]*cfcen1;
        p2 = p1+[cos(-cfang-pi/2);sin(-cfang-pi/2)]*cfcen2;
        plot(c(1),c(2),'ow','markerfacecolor','r')
        plot(p1(1),p1(2),'ow','markerfacecolor','r')
        plot(p2(1),p2(2),'ow','markerfacecolor','g')
        drawnow
%         gfsdfgdfg
waitforbuttonpress
        clf
    end
    
    
    
    
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
