function [materialprops,nsd,gravity,ned,nn,coords,nel,nen,connect,ng,gnodes,nh,hnodes] = ...
    read_input(infile)
cellarray=textscan(infile,'%s');
cellno = 0;
materialprops = zeros(5,1);
for i=1:5
    cellno=cellno+2;
    materialprops(i)=str2num(cellarray{1}{cellno});
end
cellno = cellno + 2;
nsd = str2num(cellarray{1}{cellno});
cellno = cellno + 1;
gravity = zeros(nsd,1);
for i=1:nsd
    cellno=cellno+1;
    gravity(i)=str2num(cellarray{1}{cellno});
end
cellno = cellno + 2;
ned = str2num(cellarray{1}{cellno});
cellno = cellno + 2;
nn = str2num(cellarray{1}{cellno});
cellno = cellno + 1;
% coords
coords = zeros(nsd,nn);
for i = 1 : nn
     for j = 1 : nsd
       cellno = cellno + 1;
       coords(j,i) = str2num(cellarray{1}{cellno});
     end
end
cellno = cellno + 2;
nel = str2num(cellarray{1}{cellno});
cellno = cellno + 2;
nen = str2num(cellarray{1}{cellno});
cellno = cellno + 1;
% connect
connect = zeros(nel,nen);
for i = 1:nel
    for j = 1:nen
        cellno = cellno + 1;
        connect(i,j) = str2num(cellarray{1}{cellno});
    end
end
connect = transpose(connect);
cellno = cellno + 2;
ng = str2num(cellarray{1}{cellno});
cellno = cellno + 3;
% gnodes
gnodes = zeros(ng,3);
for i = 1:ng
    for j = 1:3
        cellno = cellno + 1;
        gnodes(i,j) = str2num(cellarray{1}{cellno});
    end
end
gnodes = transpose(gnodes);
cellno = cellno + 2;
nh = str2num(cellarray{1}{cellno});
% hnodes
cellno = cellno + 3;
hnodes = zeros(nh,4);
for i = 1:nh
    for j = 1:4
        cellno = cellno + 1;
        hnodes(i,j) = str2num(cellarray{1}{cellno});
    end
end
hnodes = transpose(hnodes);
end







