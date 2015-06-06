function write_output(outfile,nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs)
fprintf(outfile,' Nodal Displacements: \n');
if ned == 2
    fprintf(outfile,' Node      Coords         u1       u2 \n');
    for i = 1:nn
        fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f\n',...
                i,coords(1,i),coords(2,i),dofs(2*i-1),dofs(2*i));
    end
elseif ned == 3
    fprintf(outfile,' Node            Coords            u1       u2       u3 \n');
    for i = 1:nn
    fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', ...
        i,coords(1,i),coords(2,i),coords(3,i),dofs(3*i-2),dofs(3*i-1),dofs(3*i));
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
        fprintf(outfile,'  \n int pt    Coords          B_11      B_22     B_12      s_11       s_22      s_12 \n');
    elseif nsd == 3
        fprintf(outfile,'\n int pt         Coords            B_11      B_22     B_33      B_12       B_13      B_23      s_11      s_22      s_33      s_12      s_13      s_23 \n');
    end
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
        stress = Kirchhoffstress(ned,nsd,B,J,materialprops);
        % computing Cauchy stress
        stress = stress/J;
        if nsd == 2
            fprintf(outfile,'%5d %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n', ...
                intpt,x(1),x(2),B(1,1),B(2,2),B(1,2),stress(1,1),stress(2,2),stress(1,2));
        elseif nsd == 3
            fprintf(outfile,'%5d %7.4f %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n',...
                intpt,x(1),x(2),x(3), ...
                B(1,1),B(2,2),B(3,3),B(1,2),B(1,3),B(2,3), ...
                stress(1,1),stress(2,2),stress(3,3),stress(1,2),stress(1,3),stress(2,3));
        end
    end
end
end


    
