function plotmesh(coords,nsd,connect,nel,nen,color)
% Function to plot a mesh.  
f2D_3 = [1,2,3];
f2D_4 = [1,2,3,4];
f3D_4 = [[1,2,3];[1,4,2];[2,4,3];[3,4,1]];
f3D_8 = [[1,2,3,4];[5,8,7,6];[1,5,6,2];[2,3,7,6];[3,7,8,4];[4,8,5,1]];
hold on
if nsd==2
    for ele = 1:nel
        for i = 1:nen
            x(i,1:2) = coords(1:2,connect(i,ele));
        end
        scatter(x(:,1),x(:,2),'MarkerFaceColor',color);
        if (nen==3)
            patch('Vertices',x,'Faces',f2D_3,'FaceColor','none','EdgeColor',color);
        elseif(nen==4)
             patch('Vertices',x,'Faces',f2D_4,'FaceColor','none','EdgeColor',color);
        end
    end
elseif nsd==3
    for ele = 1:nel
        for i=1:nen
            x(i,1:3)=coords(1:3,connect(i,ele));
        end
        scatter3(x(:,1),x(:,2),x(:,3),'MarkerFaceColor',color);
        if nen==4 
            patch('Vertices',x,'Faces',f3D_4,'FaceColor','none','EdgeColor',color);
        elseif nen==8
            patch('Vertices',x,'Faces',f3D_8,'FaceColor','none','EdgeColor',color);
        end
    end
end         
axis equal
hold off
end