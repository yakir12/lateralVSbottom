close all
clear all
clc
%%
bbuff = 30;
sbuff = 10;
height = 420;
width = 240;
path = '2processedata';
fid = fopen('table.csv','r');
% table = textscan(fid,'%s','Delimiter',','); %s %u %s %u %u %u %u %u %u %s %u %u
head = textscan(fid,'%s',13,'Delimiter',','); 
head = head{1}';
table = textscan(fid,'%s %f %s %f %f %f %f %f %f %s %f %f %s','Delimiter',','); 
fclose(fid);
% table{1} = cellfun(@(x) [x(1:end-4),'.png'],table{1},'uni',0);
% table{13} = cellfun(@(x) regexprep(x,'data','2processedata'),table{13},'uni',0);
table{11}(isnan(table{11})) = 10;
%%
botu = sort(unique(table{3}));
scru = sort(unique(table{10}));
vdisu = sort(unique(table{11}));
cfu = sort(unique(table{2}));
botn = length(botu);
scrn = length(scru);
vdisn = length(vdisu);
cfn = length(cfu);
%%
I = repmat({repmat({repmat({128*ones(height,width,'uint8')},[1,cfn])},[vdisn,1])},[scrn,botn]);
for i = 1:botn
    for j = 1:scrn
        for k = 1:vdisn
            for l = 1:cfn
                tf = strcmp(botu{i},table{3}) & strcmp(scru{j},table{10}) & ...
                    vdisu(k) == table{11} & cfu(l) == table{2};
                ind = find(tf);
                if ~isempty(ind)
                    [~,this] = min(abs(ind-20));
                    ind = ind(this);
                    fldr = table{13}{ind};
                    I{j,i}{k}{l} = rot90(imread(fullfile(fldr,table{1}{ind})));
                end
            end
        end
    end
end

for i = 1:botn
    for j = 1:scrn
        for k = 1:vdisn
            I{j,i}{k} = padarray(cell2mat(I{j,i}{k}),[sbuff,sbuff],uint8(255));
        end
    end
end

for i = 1:botn
    for j = 1:scrn
        I{j,i} = padarray(cell2mat(I{j,i}),[bbuff,bbuff],uint8(0));
    end
end

for i = 1:botn
    I{end,i}(height+10+bbuff:end,:) = [];
end

I = cell2mat(I);
I(end:end+bbuff,:) = uint8(0);
% % I = flipud(I);
%%
figure('color','w')
imshow(I,'InitialMagnification',7)
hold on
set(gca,'ydir','normal')
for i = 1:botn
    text((2*bbuff+cfn*width+2*sbuff)*(i-.5),-400,botu{i},'HorizontalAlignment','center')
    for j = 1:cfn
        text((2*bbuff+cfn*width+2*sbuff)*(i-1)+bbuff+sbuff+(j-.5)*width,-150,num2str(cfu(j)),'HorizontalAlignment','center')
    end
end

text(size(I,2)/2,-700,'Bottom stimulus','HorizontalAlignment','center')

for i = 1:scrn-1
    text(-500,(2*bbuff+vdisn*(height+2*sbuff))*(i-.5),scru{i},'HorizontalAlignment','center','Rotation',90)
    for j = 1:vdisn
        text(-350,(2*bbuff+vdisn*(height+2*sbuff))*(i-1)+bbuff+sbuff+(j-.5)*height,-150,num2str(vdisu(j)),'VerticalAlignment','middle')
    end    
end
text(-500,(2*bbuff+vdisn*(height+2*sbuff))*(i)+bbuff+sbuff+height/2,scru{i+1},'HorizontalAlignment','center','Rotation',90)
text(-350,(2*bbuff+vdisn*(height+2*sbuff))*(i)+bbuff+sbuff+(1-.5)*height,-150,'-','VerticalAlignment','middle')

text(-800,size(I,1)/2,'Side stimulus','HorizontalAlignment','center','Rotation',90)

set(gcf,'position',[1108         870         784         749])

export_fig factor.pdf -pdf
