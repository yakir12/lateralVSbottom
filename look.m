close all
clear al
clc
%%
[~,~] = unix('/usr/bin/killall PTPCamera');
[~,~] = unix('/opt/local/bin/gphoto2 --capture-image-and-download --filename /Users/holodeck1/Documents/MATLAB/Code/temp.jpg --force-overwrite');
figure('position',[1921,1,1920,1007])
imshow('temp.jpg','InitialMagnification','fit')
delete('temp.jpg')