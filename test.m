close all
clear all
clc
%%
paperwidth = 7; %inches
height = 9; %inches
cs = 5; %mm
%% convert to cm
paperwidth = paperwidth*2.54;
height = height*2.54;
cs = cs/10;
%%
width = paperwidth/3;
%%
h = figure('color','w','units','centimeters','PaperUnits','centimeters','PaperSize',[paperwidth,height],'PaperPosition',[0,0,paperwidth,height]);%,'Position',[0,0,paperwidth,height]);
axes('position',[0,0,1,1]);
patch([0,width,width,0],[0,0,height,height],.5*ones(1,3),'edgecolor','none')

nx = ceil(width/cs);
ny = ceil(height/cs);
x = linspace(width,2*width,nx);
x(x > 2*width) = [];
x(end+1) = 2*width;
y = linspace(0,height,ny);
y(y > height) = [];
y(end+1) = height;
for i = 1:nx-1
    for j = 1:ny-1
        patch([x(i),x(i+1),x(i+1),x(i)],[y(j),y(j),y(j+1),y(j+1)],(1+(-1)^i*(-1)^j)/2*ones(1,3),'edgecolor','none')        
    end
end

cs0 = .1; %in cm
nx = ceil(width/cs0);
ny = ceil(height/cs0);
x = linspace(2*width,3*width,nx);
x(x > 3*width) = [];
x(end+1) = 3*width;
y = linspace(0,height,ny);
y(y > height) = [];
y(end+1) = height;
col = rand(ny-1,nx-1);
for i = 1:nx-1
    for j = 1:ny-1
        patch([x(i),x(i+1),x(i+1),x(i)],[y(j),y(j),y(j+1),y(j+1)],col(j,i)*ones(1,3),'edgecolor','none')        
    end
end

axis image
axis off

print(h,'-dpdf','-r72','test.pdf')