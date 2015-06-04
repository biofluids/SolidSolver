function [nn,nel,connect,coords]=meshgenerate(nsd,nen,a,b,h)
% a=length of the beam, b=thickness, h=meshsize
count = 0;
connect = zeros(nen,1);
if nsd == 2
    nelx = a/h;
    nely = b/h;
    nel = nelx*nely;
    nnx = nelx + 1;
    nny = nely + 1;
    nn = nnx*nny;
    count = 0;
    for i=1:nny
        for j=1:nnx
            count = count + 1;
            coords(1,count)=(j-1)*h;
            coords(2,count)=(i-1)*h;
        end
    end
    count = 0;
    if nen == 4   
        for i=1:nely
            temp1 = (i-1)*nnx + 1;
            temp2 = temp1 + nelx-1;
            for j=temp1:temp2
                count = count + 1;
                temp3 = [j;j+1;j+1+nnx;j+nnx];
                connect = [connect,temp3];
            end
        end
        connect = connect(:,2:end);    
    elseif nen == 3
        for i=1:nely
            temp1 = (i-1)*nnx + 1;
            temp2 = temp1+nelx-1;
            for j=temp1:temp2
                count=count+1;
                temp3=[j;j+1;j+1+nnx];
                connect = [connect,temp3];
                count=count+1;
                temp3=[j+1+nnx;j+nnx;j];
                connect = [connect,temp3];
            end
        end
        connect = connect(:,2:end);
        nel = 2*nel;                
    end
end
end
        
        