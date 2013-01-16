close all
clear all
clc
delete('*.tiff')
delete('*.pdf')
%% parameters
viewd2 = 1219.2/2; %mm
checkersz = [2,4,7]; %mm
nchecker = length(checkersz);
pixelsz = .48; %mm (on work station precision 0.28194)
staticsz = 1; %mm
nviewd = 4;
viewd1 = linspace(10,viewd2,nviewd); %mm
papersz = [35,35]; %inch
paperdpi = 30; %dots per inch
screensz = [4*1920,1058]; %pixels
%% settings
paperres = sprintf('-r%d',paperdpi);
paperoptions = {'color','w','units','inches','PaperUnits','inches',...
    'PaperSize',papersz,...
    'PaperPosition',[0,0,papersz],'colormap',colormap(gray)};
close
%% background 1
h = figure(paperoptions{:});
axes('units','normalized','position',[0,0,1,1],'Visible','off')
patch([0,papersz(1),papersz(1),0],[0,0,papersz(2),papersz(2)],.5*ones(1,3),'EdgeColor','none')
axis image
print(h,'-dpdf',paperres,'Uniform.pdf')

I = 128*ones(fliplr(screensz),'uint8');
imwrite(I,'Uniform.tiff')
%% background2
for j = 1:nchecker
    h = figure(paperoptions{:});
    axes('units','normalized','position',[0,0,1,1],'Visible','off')
    background(papersz,'Checkerboard',checkersz(j))
    axis image
    print(h,'-dpdf',paperres,sprintf('Checkerboard_%d.pdf',checkersz(j)))
    for i = 1:nviewd
        cheackersz2 = viewd2/viewd1(i)*checkersz(j);
        n = round(cheackersz2/pixelsz);
        p = ceil(screensz(2)/n/2);
        q = ceil(screensz(1)/n/2);
        I = uint8(255*(checkerboard(n,p,q) > .5));
        I(screensz(2)+1:end,:) = [];
        I(:,screensz(1)+1:end) = [];
        imwrite(I,sprintf('Checkerboard_%d_%.0f.tiff',checkersz(j),viewd1(i)))
    end
end
%% background3
h = figure(paperoptions{:});    
axes('units','normalized','position',[0,0,1,1],'Visible','off')
background(papersz,'Static',staticsz)
axis image
print(h,'-dpdf',paperres,'Static.pdf')    
for i = 1:nviewd
    staticsz2 = viewd2/viewd1(i)*staticsz;
    n = staticsz2/pixelsz;
    p = round(screensz(2)/n);
    q = round(screensz(1)/n);
    I = uint8(255*(rand(p,q)));
    I = imresize(I,fliplr(screensz),'nearest','Method','box');
    imwrite(I,sprintf('Static_%.0f.tiff',viewd1(i)))
end