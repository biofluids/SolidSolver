function write_variables(outfile,nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs)
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
    sum_sigma=zeros((nsd-1)*3,1);
    sum_B=zeros((nsd-1)*3,1);
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
%% write to file
fprintf(outfile,'This is a symm tensor per element file\n');
fprintf(outfile,'part\n');
fprintf(outfile,'%10d\n',1);
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
% print stress
for i=1:6
    for j=1:nel
        fprintf(outfile,'%12.5f\n',stress(i,j));
    end
end
end


    
