function write_geo(outfile,nsd,ned,nn,coords,nel,nen,connect,materialprops,w)
fprintf(outfile,'This is a geometry file\n'); % description line 1
fprintf(outfile,'This is a geometry file\n'); % description line2
fprintf(outfile,'node id given\n'); 
fprintf(outfile,'element id given\n');
fprintf(outfile,'part\n');
fprintf(outfile,'%10d\n',1); % part number
fprintf(outfile,'solid\n'); % description of part
fprintf(outfile,'coordinates\n');
fprintf(outfile,'%10d\n',nn); % number of nodes
% node ids
for i=1:nn
    fprintf(outfile,'%10d\n',i); 
end
% node coordinates: x;y;z
for i=1:nsd
    for j=1:nn
        row=(j-1)*nsd+i;
        coords(row)=coords(i,j)+w(row);
        fprintf(outfile,'%12.5e\n',coords(i,j));
    end
end
if nsd==2
    for i=1:nn
        fprintf(outfile,'%12.5e\n',0);
    end
end
% element type
if nsd==2
    if nen==3
        fprintf(outfile,'tria3\n');
    elseif nen==4
        fprintf(outfile,'quad4\n');
    end
elseif nsd==3
    if nen==4
        fprintf(outfile,'tetra4\n');
    elseif nen==8
        fprintf(outfile,'hexa8\n');
    end
end
% number of elements        
fprintf(outfile,'%10d\n',nel);
% element ids
for i=1:nel
    fprintf(outfile,'%10d\n',i);
end
% connectivity
for i=1:nel
    fprintf(outfile,'%10d',connect(:,i));
    fprintf(outfile,'\n');
end
end