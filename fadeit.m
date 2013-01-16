function fadeit(a,nx,dx,h,cm)
for i = 1:nx
%     set(h,'AlphaData',a(i))
    set(h,'colormap',cm*a(i))
    drawnow
%     disp(i)
    pause(dx)
end