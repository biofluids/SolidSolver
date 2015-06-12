function write_results(nsteps,nprint,dt,nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs,nt)
str1='/Users/jiecheng/Documents/SolidResults/solid.sig';
str2='/Users/jiecheng/Documents/SolidResults/solid.b';
str3='/Users/jiecheng/Documents/SolidResults/solid.d';
str4='/Users/jiecheng/Documents/SolidResults/solid.geo';
str5='/Users/jiecheng/Documents/SolidResults/solid.case';
str=num2str(nt,'%06d');
outfile1=fopen([str1,str],'w');
outfile2=fopen([str2,str],'w');
outfile3=fopen([str3,str],'w');
outfile4=fopen([str4,str],'w');
%% write case file
if nt==0
    outfile5=fopen(str5,'w');
    fprintf(outfile5,'FORMAT\n\n');
    fprintf(outfile5,'type:\tensight gold\n\n');
    fprintf(outfile5,'GEOMETRY\n\n');
    fprintf(outfile5,'model:\tsolid.geo******\tchange_coords_only\n\n');
    fprintf(outfile5,'VARIABLE\n\n');
    fprintf(outfile5,'vector per node:\tdisplacement\tsolid.d******\n');
    fprintf(outfile5,'tensor symm per element:\tstress\tsolid.sig******\n');
    fprintf(outfile5,'tensor symm per element:\tstrain\tsolid.b******\n\n');
    fprintf(outfile5,'TIME\n\n');
    fprintf(outfile5,'time set:\t1\n');
    fprintf(outfile5,'number of steps:\t%10d\n',nsteps/nprint+1);
    fprintf(outfile5,'filename start number:\t%10d\n',0);
    fprintf(outfile5,'filename increment:\t%10d\n',1);
    time=dt*[0:nprint:nsteps];
    fprintf(outfile5,'time values:\n');
    fprintf(outfile5,'%12.3f\n',time);
end
%% compute averaged Cauchy stress and left Cauchy-Green tensor for each element
elecoord = zeros(nsd,nen);
eledof = zeros(ned,nen);
stress=zeros(6,nen);
strain=zeros(6,nen);
% loop over elements
for ele = 1:nel
    % Extract coords of nodes, and dof for the element
    for a = 1:nen
        elecoord(:,a) = coords(:,connect(a,ele));
        for i = 1:ned
            eledof(i,a) = dofs(ned*(connect(a,ele)-1)+i);
        end
    end
    dxdxi = zeros(nsd,nsd);
    dNdx = zeros(nen,nsd);
    % set up quadrature points and weights
    npt = no_quad_pts(nsd,nen,0);
    xilist = quad_pts(nsd,nen,npt);
    % to evaluate average stress for each element 
    sum_sigma=zeros(6,1);
    sum_B=zeros(6,1);
    % looping over quadrature points
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nen,nsd,xi);
        N = sf(nen,nsd,xi);
        % compute the coords of the quadrature points
        x = elecoord*N;
        % set up the jacobian matrix
        dxdxi = elecoord*dNdxi;
        dxidx = inv(dxdxi);
        % computing dNdx
        dNdx = dNdxi*dxidx;
        % computing deformation gradient, F_ij = delta_ij + du_i/dx_j
        F = zeros(nsd,nsd);
        F = eye(nsd) + eledof*dNdx;
        % computing B = FF^T and J = det(F)
        B = F*transpose(F);
        J = det(F); 
        % computing Kirchhoff stress 
        sigma = Kirchhoffstress(ned,nsd,B,J,materialprops);
        % computing Cauchy stress
        sigma = sigma/J;
        % vectorize stress and strain
        if nsd==2
            sigma=[sigma(1,1);sigma(2,2);0;sigma(1,2);0;0];
            B=[B(1,1);B(2,2);0;B(1,2);0;0];
        else
            sigma=[sigma(1,1);sigma(2,2);sigma(3,3);sigma(1,2);sigma(1,3);sigma(2,3)];
            B=[B(1,1);B(2,2);B(3,3);B(1,2);B(1,3);B(2,3)];
        end
        % averaged stress in this element
        sum_sigma=sum_sigma+sigma;
        sum_B=sum_B+B;
    end
    % stress and strain for this element
    sigma=sum_sigma/npt;
    B=B/npt;
    stress(:,ele)=sigma;
    strain(:,ele)=B;                   
end
%% write stress, strain, and displacement to file
fprintf(outfile1,'This is a symm tensor per element file for stress\n');
fprintf(outfile1,'part\n');
fprintf(outfile1,'%10d\n',1);

fprintf(outfile2,'This is a symm tensor per element file for strain\n');
fprintf(outfile2,'part\n');
fprintf(outfile2,'%10d\n',1);

fprintf(outfile3,'This is a vector per node file for displacement\n');
fprintf(outfile3,'part\n');
fprintf(outfile3,'%10d\n',1);

% element type
if nsd==2
    if nen==3
        fprintf(outfile1,'tria3\n');
        fprintf(outfile2,'tria3\n');
    elseif nen==4
        fprintf(outfile1,'quad4\n');
        fprintf(outfile2,'quad4\n');      
    end
elseif nsd==3
    if nen==4
        fprintf(outfile1,'tetra4\n');
        fprintf(outfile2,'tetra4\n');      
    elseif nen==8
        fprintf(outfile1,'hexa8\n');
        fprintf(outfile2,'hexa8\n');
    end
end
% print stress, strain
for i=1:6
    for j=1:nel
        fprintf(outfile1,'%12.5e\n',stress(i,j));
        fprintf(outfile2,'%12.5e\n',strain(i,j));
    end
end
% print displacement
fprintf(outfile3,'coordinates\n');
for i=1:nsd
    for j=1:nn
        row=ned*(j-1)+i;
        fprintf(outfile3,'%12.5e\n',dofs(row));
    end
end
if nsd==2
    for i=1:nn
        fprintf(outfile3,'%12.5e\n',0);
    end
end
%% write geometry file
fprintf(outfile4,'This is a geometry file\n'); % description line 1
fprintf(outfile4,'This is a geometry file\n'); % description line2
fprintf(outfile4,'node id given\n'); 
fprintf(outfile4,'element id given\n');
fprintf(outfile4,'part\n');
fprintf(outfile4,'%10d\n',1); % part number
fprintf(outfile4,'solid\n'); % description of part
fprintf(outfile4,'coordinates\n');
fprintf(outfile4,'%10d\n',nn); % number of nodes
% node ids
for i=1:nn
    fprintf(outfile4,'%10d\n',i); 
end
% node coordinates: x;y;z
for i=1:nsd
    for j=1:nn
        row=(j-1)*nsd+i;
        coords(row)=coords(i,j)+dofs(row);
        fprintf(outfile4,'%12.5e\n',coords(i,j));
    end
end
if nsd==2
    for i=1:nn
        fprintf(outfile4,'%12.5e\n',0);
    end
end
% element type
if nsd==2
    if nen==3
        fprintf(outfile4,'tria3\n');
    elseif nen==4
        fprintf(outfile4,'quad4\n');
    end
elseif nsd==3
    if nen==4
        fprintf(outfile4,'tetra4\n');
    elseif nen==8
        fprintf(outfile4,'hexa8\n');
    end
end
% number of elements        
fprintf(outfile4,'%10d\n',nel);
% element ids
for i=1:nel
    fprintf(outfile4,'%10d\n',i);
end
% connectivity
for i=1:nel
    fprintf(outfile4,'%10d',connect(:,i));
    fprintf(outfile4,'\n');
end
%% end writing files
fclose(outfile1);
fclose(outfile2);
fclose(outfile3);
fclose(outfile4);
end


    
