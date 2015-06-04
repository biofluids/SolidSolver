function plotmesh(coords,nsd,nn,connect,nel,nen,color)
% Function to plot a mesh.  
f2D_3 = [1,2,3];
f2D_4 = [1,2,3,4];
hold on
if (nsd==2)
    for ele = 1:nel
        for i = 1:nen
            x(i,1:2) = coords(1:2,connect(i,ele));
        end
        scatter(x(:,1),x(:,2),'MarkerFaceColor','r');
        if (nen==3)
            patch('Vertices',x,'Faces',f2D_3,'FaceColor','none','EdgeColor',color);
        elseif(nen==4)
             patch('Vertices',x,'Faces',f2D_4,'FaceColor','none','EdgeColor',color);
        end
    end
end
axis equal
hold off
end