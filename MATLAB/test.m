% this is used to validate Fortran code
infile = fopen('displacement.txt','r')
cellarray=textscan(infile,'%s');
cellno=1;
temp = zeros(nsd,nn);
for i=1:3
    for j=1:nn
        row=ned*(i-1)+j;
        temp(i,j) = str2num(cellarray{1}{cellno});
        cellno=cellno+1;
    end
end

for i=1:nn
    for j=1:ned
        row = ned*(i-1)+j;
        dofs(row) = temp(j,i);
    end
end
dofs = transpose(dofs);
Fint = internalforce(nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs);
K = Kint(nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs);
Fext = externalforce(nsd,ned,nn,nel,nen,no_bc2,materialprops,gravity,coords,connect,bc2);