function exp1()
%% parameters
close all
path = 'data';
nrep = 1;
fadetime = 5; %in seconds
testtime = 5*60; %in seconds
cm = colormap(gray(256));
close
paperoptions = {'color','k','colormap',cm,'position',[1920*2+1,-32,1920*4,1058],'menubar','none'};
%% settings
[~,~] = unix('/usr/bin/killall PTPCamera');
imfiles = dir('Stimuli/*.tiff');
scrnames = cellfun(@(x) x(1:end-5),{imfiles.name},'UniformOutput',false);
n = length(imfiles);
info = imfinfo(['Stimuli/',imfiles(1).name]);
sz = [info.Width,info.Height];
fh = figure(paperoptions{:});
ah = axes('units','normalized','position',[0,0,1,1],'parent',fh,'color','k');
ih = image(zeros(fliplr(sz),'uint8'),'parent',ah);
axis image
axis off
bkgdfiles = dir('Stimuli/*.pdf');
bkgdnames = cellfun(@(x) x(1:end-4),{bkgdfiles.name},'UniformOutput',false);
data = [scrnames;repmat({zeros(fliplr(sz),'uint8')},1,n)];
I = struct(data{:});
for i = 1:n
    I.(scrnames{i}) = imread(strcat('Stimuli/',scrnames{i},'.tiff'));
end
%%
nfade = 40-2;%254;
afade = linspace(0+1/256,1-1/256,nfade);
x = erfinv(2*(afade-.5));
afade = [0,afade,1];
x = [x(1)-diff(x(1:2)),x,x(end)+diff(x(end-1:end))];
x = x-min(x);
x = x/max(x)*fadetime;
dfade = diff(x);
dfade = [dfade(1),dfade];
%%
S = {'1','2','3','4','5','6','7'};
[i,ok] = listdlg('ListString',S,'SelectionMode','single','Name',...
    'Cuttlefish','PromptString','Choose individual');
if ok
    cuttlefish = S{i};
else
    close(fh)
end
[i,ok] = listdlg('ListString',bkgdnames,'SelectionMode','single','Name',...
    'Bottom background','PromptString','Choose Background');
if ok
    bkgd = bkgdnames{i};
else
    close(fh)
end
ind = zeros(1,n*nrep);
ind(1:n) = randperm(n);
for i = 2:nrep
    done = false;
    while ~done
        ind2 = randperm(n);
        if ind2(1) ~= ind(end)
            ind((i-1)*n+1:i*n) = ind2;
            done = true;
        end
    end
end
order = scrnames(ind);
m = length(order);
time = datestr(now);
remarks = inputdlg('Remarks?','Remarks',1,{'All normal'});
if isempty(remarks)
    remarks = {'All normal'};
end

%folder name: cuttlefish bkgd year month day hour minute
fldr = strcat(cuttlefish,'_',bkgd,'_',sprintf('%d_%d_%d_%d_%d%*.0f',datevec(now)));
if exist(fullfile(path,fldr),'dir')
    movefile(fullfile(path,fldr),strcat(fullfile(path,fldr),'_OLD'))
end
mkdir(fullfile(path,fldr))

fid = fopen(fullfile(path,fldr,'details.txt'),'w');
txt = repmat('%s\t',1,m);
txt(end-1:end) = [];
fprintf(fid,...
    ['cuttlefish\t%s\nbackground\t%s\ntime\t%s\nremarks\t%s\nfolder\t%s\norder\t',...
    txt],...
    cuttlefish,bkgd,time,remarks{:},fldr,order{:});
fclose(fid);

%% crop
[~,~] = unix('/opt/local/bin/gphoto2 --capture-image-and-download --filename /Users/holodeck1/Documents/MATLAB/Code/temp.jpg --force-overwrite');
fh1 = figure;
drawnow
I1 = imread('temp.jpg');
imshow(I1,'InitialMagnification',48)
delete('temp.jpg')
[~,rect] = imcrop;
close(fh1)
%%

ok = questdlg('Start?','Start the displays','Yes','Abort','Yes'); 
if strcmp(ok,'Abort')    
    close(fh)
end

pause(600)

for i = 1:m
    set(ih,'CData',I.(order{i}))
    fadein()
    drawnow
    tic
    start = toc;
    finish = toc;
    inx = 0;
    while finish < start+testtime && inx <= 24 
        inx = inx+1;
        
        [~,~] = unix('/opt/local/bin/gphoto2 --capture-image-and-download --filename /Users/holodeck1/Documents/MATLAB/Code/temp.jpg');
        while ~exist('temp.jpg','file')
            pause(.5)
%             disp('kaka')
        end
        %image name: cuttlefish bkgd year month day hour minute second
        %screen replicate 
        filename = strcat(cuttlefish,'_',bkgd,'_',sprintf('%d_%d_%d_%d_%d_%.0f',datevec(now)),...
            '_',order{i},'_',sprintf('%d',inx),'.jpg');
        I1 = imread('temp.jpg');
        delete('temp.jpg')
        I1 = imcrop(I1,rect);
        imwrite(I1,fullfile(path,fldr,filename))
        finish = toc;
    end
    fadeout()
end

remarks = inputdlg('Finishing remarks?','Remarks',1,{'All normal'});
if isempty(remarks)
    remarks = {'All normal'};
end
fid = fopen(fullfile(path,fldr,'details.txt'),'a');
fprintf(fid,'\nfinishing remarks\t%s',remarks{:});
fclose(fid);

questdlg('Done','Remove cuttlefish','OK','Abort','OK');
close(fh)



%%%%%%%%%%%%%%%%%%%%%%%%%
function fadein()
    for ifade = 1:nfade
        set(fh,'Colormap',afade(ifade)*cm)
        drawnow
        pause(dfade(ifade))
    end
end
function fadeout()
    for ifade = nfade:-1:1
        set(fh,'Colormap',afade(ifade)*cm)
        drawnow
        pause(dfade(ifade))
    end
end

end