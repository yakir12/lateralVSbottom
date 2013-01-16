close all
clear all
clc
%%
ntheta = 180;
thetad = linspace(-90,90,ntheta);
thetad(end) = [];
thetar = deg2rad(thetad);
costheta = cos(thetar);
sintheta = sin(thetar);
fac = 10;
sz1 = 100;
sz2 = 1000;
path1 = 'data';
path2 = 'processedata';
fldr = dir(path2);
killd = cellfun(@(x) strcmp(x(1),'.'),{fldr.name});
fldr(killd) = [];
n = 1;%length(fldr);
%%
theta = linspace(0,2*pi,10);
figure('colormap',gray)
for i = 1:n
    im = dir(fullfile(path2,fldr(i).name,'*.png')); 
    killf = cellfun(@(x) strcmp(x(1),'m'),{im.name});
    im(killf) = [];
    m = length(im);  
    I = imread(fullfile(path2,fldr(i).name,im(1).name));
    sz = size(I);
    c = floor((fliplr(sz)'+1)/2);
    for j = 1:m
        I = imread(fullfile(path2,fldr(i).name,im(j).name));
        [r,x] = radon(I,thetad);
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
        imagesc(I) 
        axis image
        set(gca,'ydir','normal')
        hold on
        p1 = c+[cos(-cfang);sin(-cfang)]*cfcen1;
        p2 = p1+[cos(-cfang-pi/2);sin(-cfang-pi/2)]*cfcen2;
        plot(c(1),c(2),'ow','markerfacecolor','r')
        plot(p1(1),p1(2),'ow','markerfacecolor','r')
        plot(p2(1),p2(2),'ow','markerfacecolor','g')
        drawnow
        waitforbuttonpress
%         I1 = imread(fullfile(path2,fldr(i).name,im(j).name));
%         [y,x] = find(I1);
%         el = fit_ellipse(x,y);  
%         I = imread(fullfile(path1,fldr(i).name,[im(j).name(1:end-3),'jpg']));
%         def = [1,el.a*2/sz1-1,el.b*2/sz1-1,el.X0*2/sz1-1,el.Y0*2/sz1-1,rad2deg(el.phi)];
%         P = phantom(def,sz2);
%         subplot(2,2,1)
%         imshow(I,[])
%         subplot(2,2,2)
%         imshow(I1,[])
%         hold on
%         ellipse(el.a,el.b,el.phi,el.X0,el.Y0,'r')
%         subplot(2,2,3)
%         imshow(P,[])
%         
%         jsaghsjkadh
%         I = imrotate(I,rad2deg(el.phi));
%         rect = [el.X0-el.a/2,el.Y0-el.b/2,el.a,el.b];
%         I = imcrop(I,rect);
%         imshow(I)
%         drawnow
    end
end