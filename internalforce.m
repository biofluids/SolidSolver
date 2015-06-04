function Fint = internalforce(nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs)
% returns the global internal force
Fint = zeros(ned*nn,1);
elecoord = zeros(nsd,nen);
eledof = zeros(ned,nen);
%%
% Loop over elements
for ele = 1:nel
    % Extract coords of nodes, and dof for the element
    for a = 1:nen
        elecoord(:,a) = coords(:,connect(a,ele)); 
        for i = 1:ned
            eledof(i,a) = dofs(ned*(connect(a,ele)-1)+i);
        end
    end    
    %% compute the internal force of the element
    % 
    fint = zeros(ned*nen,1);
    % fully integration
    % set up quadrature points and weights, 0 means fully integration
    npt = no_quad_pts(nsd,nen,0);
    xilist = quad_pts(nsd,nen,npt);
    weights = quad_weights(nsd,nen,npt);
    % loop over quadrature points
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nen,nsd,xi);
        % set up the jacobian matrix
        dxdxi = zeros(nsd,nsd);
        dxdxi = elecoord*dNdxi;
        dxidx = inv(dxdxi);
        dt = det(dxdxi);
        % compute dNdx
        dNdx = zeros(nen,nsd);
        dNdx = dNdxi*dxidx;
        % compute deformation gradient, F_ij = delta_ij + du_i/dx_j
        F = zeros(nsd,nsd);
        F = eye(nsd) + eledof*dNdx;
        % compute left Cauchy-Green tensor, B = FF^T and J = det(F)
        B = F*transpose(F);
        J = det(F);
        % compute dNdy, in which y is the coord. after deformation
        dNdy = zeros(nen,nsd);
        Finv = inv(F);
        dNdy = dNdx*Finv;
        % compute the stress, Kirchhoffstress is used
        stress = Kirchhoffstress(ned,nsd,B,J,materialprops);
        % compute the element internal force
        for a = 1:nen
            for i = 1:nsd
                row = (a-1)*ned + i;
                for j = 1:nsd
                    fint(row) = fint(row) + stress(i,j)*dNdy(a,j)*weights(intpt)*dt;
                    fint(row) = fint(row) - stress(j,j)/nsd*dNdy(a,i)*weights(intpt)*dt;
                end
            end
        end       
    end
    % reduced integration
    % set up quadrature points and weights, 1 means reduced integration
    npt = no_quad_pts(nsd,nen,1);
    xilist = quad_pts(nsd,nen,npt);
    weights = quad_weights(nsd,nen,npt);
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nen,nsd,xi);
        % set up the jacobian matrix
        dxdxi = zeros(nsd,nsd);
        dxdxi = elecoord*dNdxi;
        dxidx = inv(dxdxi);
        dt = det(dxdxi);
        % compute dNdx
        dNdx=zeros(nen,nsd);
        dNdx = dNdxi*dxidx;
        % compute deformation gradient, F_ij = delta_ij + du_i/dx_j
        F = zeros(nsd,nsd);
        F = eye(nsd) + eledof*dNdx;
        % compute left Cauchy-Green tensor, B = FF^T and J = det(F)
        B = F*transpose(F);
        J = det(F);
        % compute dNdy, in which y is the coord. after deformation
        dNdy = zeros(nen,nsd);
        Finv = inv(F);
        dNdy = dNdx*Finv;   
        % compute the stress, Kirchhoffstress is used
        stress = Kirchhoffstress(ned,nsd,B,J,materialprops);
        % compute the element internal force
        for a = 1:nen
            for i = 1:ned
                row = (a-1)*ned + i;
                for j = 1:nsd
                    fint(row) = fint(row) + stress(j,j)/nsd*dNdy(a,i)*weights(intpt)*dt;
                end
            end
        end       
    end
    %% scatter the element internal force into the global internal force
    % 
    for a = 1:nen
        for i = 1:ned
            row = ned*(connect(a,ele)-1) + i;
            Fint(row) = Fint(row) + fint(ned*(a-1)+i);
        end
    end
end
end