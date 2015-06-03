function Fint = original_internalforce(nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs)
% returns the global internal force
Fint = zeros(ned*nn,1);
elecoord = zeros(nsd,nen);
eledof = zeros(ned,nen);
%%
% Loop over elements
for ele = 1:nel
    % Extract coords of nodes, and dof for the element
    for a = 1:nen
        for i = 1:nsd
            elecoord(i,a) = coords(i,connect(a,ele));
        end
        for i = 1:ned
            eledof(i,a) = dofs(ned*(connect(a,ele)-1)+i);
        end
    end    
    %% compute the internal force of the element
    % 
    %
    fint = zeros(ned*nen,1);
    % fully integration
    % set up quadrature points and weights
    npt = no_quad_pts(nsd,nen,0);
    xilist = quad_pts(nsd,nen,npt);
    weights = quad_weights(nsd,nen,npt);
    % loop over quadrature points
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nen,nsd,xi);
        N = sf(nen,nsd,xi);
        % set up the jacobian matrix
        dxdxi = zeros(nsd,nsd);
        for i = 1:nsd
            for j = 1:nsd
                for a = 1:nen
                    dxdxi(i,j) = dxdxi(i,j) + dNdxi(a,j)*elecoord(i,a);
                end
            end
        end
        dxidx = inv(dxdxi);
        dt = det(dxdxi);
        % compute dNdx
        dNdx=zeros(nen,nsd);
        for a = 1:nen
            for i = 1:nsd
                for j = 1:nsd
                    dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
                end
            end
        end
        % compute deformation gradient, F_ij = delta_ij + du_i/dx_j
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
        % compute left Cauchy-Green tensor, B = FF^T and J = det(F)
        B = F*transpose(F);
        J = det(F);
        % compute dNdy, in which y is the coord. after deformation
        dNdy = zeros(nen,nsd);
        Finv = inv(F);
        for a = 1:nen
            for i = 1:nsd
                for j = 1:nsd
                    dNdy(a,i) = dNdy(a,i) + dNdx(a,j)*Finv(j,i);
                end
            end
        end    
        % compute the stress, Kirchhoffstress is used
        stress = Kirchhoffstress(ned,nsd,B,J,materialprops);
        % the hydrostatic part, what about stress(3,3)?
        p = -trace(stress)/3.;
        %stress = stress + p*eye(nsd);
        % compute the element internal force
        for a = 1:nen
            for i = 1:ned
                row = (a-1)*ned + i;
                for j = 1:nsd
                    fint(row) = fint(row) + stress(i,j)*dNdy(a,j)*weights(intpt)*dt;
                end
            end
        end       
    end
    %% scatter the element internal force into the global internal force
    % 
    % 
    for a = 1:nen
        for i = 1:ned
            row = ned*(connect(a,ele)-1) + i;
            Fint(row) = Fint(row) + fint(ned*(a-1)+i);
        end
    end
end
end