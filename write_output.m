function write_output(outfile,nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs)
fprintf(outfile,' Nodal Displacements: \n');
if nsd == 2
    fprintf(outfile,' Node      Coords         u1       u2 \n');
    for i = 1:nn
        fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f\n',...
                i,coords(1,i),coords(2,i),dofs(2*i-1),dofs(2*i));
    end
end
% print strain and stress at the quadrature points
fprintf(outfile,'\n\n Left Cauchy-Green deformation tensor and Cauchy stresses \n');
elecoord = zeros(nsd,nen);
eledof = zeros(ned,nen);
% loop over elements
for ele = 1:nel
    fprintf(outfile,' \n Element; %d ',ele);
    if nsd == 2
        fprintf(outfile,'  \n B_11    B_22    B_12     s_11         s_22      s_12 \n');
    end
    % Extract coords of nodes, and dof for the element
    for a = 1:nen
        for i = 1:nsd
            elecoord(i,a) = coords(i,connect(a,ele));
        end
        for i = 1:ned
            eledof(i,a) = dofs(ned*(connect(a,ele)-1)+i);
        end
    end
    dxdxi = zeros(nsd,nsd);
    dNdx = zeros(nen,nsd);
    % set up quadrature points and weights
    npt = no_quad_pts(nsd,nen,0);
    xilist = quad_pts(nsd,nen,npt);
    % looping over quadrature points
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nen,nsd,xi);
        N = sf(nen,nsd,xi);
        % compute the coords of the quadrature points
        for i = 1:nsd
            x(i) = 0.;
            for a = 1:nen
                x(i) = x(i) + elecoord(i,a)*N(a);
            end
        end
        % set up the jacobian matrix
        for i = 1:nsd
            for j = 1:nsd
                dxdxi(i,j) = 0;
                for a = 1:nen
                    dxdxi(i,j) = dxdxi(i,j) + dNdxi(a,j)*elecoord(i,a);
                end
            end
        end
        dxidx = inv(dxdxi);
        % computing dNdx
        for a = 1:nen
            for i = 1:nsd
                dNdx(a,i) = 0;
                for j = 1:nsd
                    dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
                end
            end
        end
        % computing deformation gradient, F_ij = delta_ij + du_i/dx_j
        F = zeros(nsd,nsd);
        for i = 1:nsd
            for j = 1:nsd
                if i == j
                    F(i,j) = 1.;
                end
                for a = 1:nen
                    F(i,j) = F(i,j) + dNdx(a,j)*eledof(i,a);
                end
            end
        end
        % computing B = FF^T and J = det(F)
        B = F*transpose(F);
        J = det(F); 
        % computing Kirchhoff stress 
        stress = Kirchhoffstress(ned,nsd,B,J,materialprops);
        % computing Cauchy stress
        stress = stress/J;
        if (nsd == 2) 
            fprintf(outfile,'%7.4f %7.4f %7.4f %9.4f %9.4f %9.4f\n', ...
            B(1,1),B(2,2),B(1,2),stress(1,1),stress(2,2),stress(1,2));
        end
    end
end
end


    
