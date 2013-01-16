function background(papersz,kind,data)
axes('position',[0,0,1,1],'Visible','off');
switch kind
    case 'Checkerboard'
        cs = data*0.0393701; %inches
        nx = ceil(papersz(1)/cs);
        x = linspace(0,nx*cs,nx);
        x(x >= papersz(1)) = [];
        x(nx) = papersz(1);
        ny = ceil(papersz(2)/cs);
        y = linspace(0,ny*cs,ny);
        y(y >= papersz(2)) = [];
        y(ny) = papersz(2);
        [X,Y] = meshgrid(x,y);

        N = reshape(1:nx*ny,ny,nx);

        F = zeros(1,4);
        inx = 0;
        for i = 1:2:ny-1
            for j = 1:2:nx-1
                inx = inx+1;
                F(inx,:) = [N(i,j),N(i+1,j),N(i+1,j+1),N(i,j+1)];    
            end
        end
        for i = 2:2:ny-1
            for j = 2:2:nx-1
                inx = inx+1;
                F(inx,:) = [N(i,j),N(i+1,j),N(i+1,j+1),N(i,j+1)];    
            end
        end
        patch('Faces',F,'Vertices',[X(:),Y(:)],'FaceColor',1-get(gcf,'Color'),'EdgeColor','none')
    case 'Static'
        %%        
        cs = data*0.0393701; %inches
        nx = ceil(papersz(1)/cs);
        x = linspace(0,nx*cs,nx);
        x(x >= papersz(1)) = [];
        x(nx) = papersz(1);
        ny = ceil(papersz(2)/cs);
        y = linspace(0,ny*cs,ny);
        y(y >= papersz(2)) = [];
        y(ny) = papersz(2);
        [X,Y] = meshgrid(x,y);

        N = reshape(1:nx*ny,ny,nx);

        F = zeros(1,4);
        inx = 0;
        for i = 1:ny-1
            for j = 1:nx-1
                inx = inx+1;
                F(inx,:) = [N(i,j),N(i+1,j),N(i+1,j+1),N(i,j+1)];    
            end
        end
        patch('Faces',F,'Vertices',[X(:),Y(:)],'FaceVertexCData',rand(inx,1),'EdgeColor','none','FaceColor','flat')
        %%
end
axis off
axis image