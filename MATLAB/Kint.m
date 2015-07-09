function K = Kint(nsd,ned,nn,coords,nel,nen,connect,materialprops,dofs)
% returns the linearized derivative of the global internal force, K_int
% reduced integration used
K = zeros(ned*nn,ned*nn);
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
    %% compute the kint of the element
    %
    kint = zeros(ned*nen,ned*nen);
    % fully integration part
    % set up quadrature points and weights, 0 means fully integration
    npt = no_quad_pts(nsd,nen,0);
    xilist = quad_pts(nsd,nen,npt);
    weights = quad_weights(nsd,nen,npt);
    % looping over quadrature points
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nen,nsd,xi);
        % set up the jacobian matrix
        dxdxi = zeros(nsd,nsd);
        dxdxi = elecoord*dNdxi;
        dxidx = inv(dxdxi);
        dt = det(dxdxi);
        % computing dNdx
        dNdx = zeros(nen,nsd);
        dNdx = dNdxi*dxidx;
        % computing deformation gradient, F_ij = delta_ij + du_i/dx_j
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
        % compute the tensor C_ijkl
        C = materialstiffness(ned,nsd,B,J,materialprops);
        % compute the element kint
        for a = 1:nen
            for i = 1:ned
                for b = 1:nen
                    for k = 1:ned
                        row = (a-1)*ned + i;
                        col = (b-1)*ned + k;
                        for j = 1:nsd
                            for l = 1:nsd
                                kint(row,col) = kint(row,col) + ...
                                    C(i,j,k,l)*dNdy(a,j)*dNdy(b,l)*weights(intpt)*dt;
                                kint(row,col) = kint(row,col) - C(j,j,k,l)/nsd*dNdy(a,i)*dNdy(b,l)*weights(intpt)*dt;
                            end
                            kint(row,col) = kint(row,col) - ...
                                stress(i,j)*dNdy(a,k)*dNdy(b,j)*dt*weights(intpt);
                            kint(row,col) = kint(row,col) + stress(j,j)/nsd*dNdy(a,k)*dNdy(b,i)*dt*weights(intpt);
                        end                       
                    end
                end
            end
        end
    end
    % reduced integration part
    % set up quadrature points and weights, 1 means reduced integration
    npt = no_quad_pts(nsd,nen,1);
    xilist = quad_pts(nsd,nen,npt);
    weights = quad_weights(nsd,nen,npt);
    % looping over quadrature points
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nen,nsd,xi);
        % set up the jacobian matrix
        dxdxi = zeros(nsd,nsd);
        dxdxi = elecoord*dNdxi;
        dxidx = inv(dxdxi);
        dt = det(dxdxi);
        % computing dNdx
        dNdx = zeros(nen,nsd);
        dNdx = dNdxi*dxidx;
        % computing deformation gradient, F_ij = delta_ij + du_i/dx_j
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
        % compute the tensor C_ijkl
        C = materialstiffness(ned,nsd,B,J,materialprops);
        % compute the element kint
        for a = 1:nen
            for i = 1:ned
                for b = 1:nen
                    for k = 1:ned
                        row = (a-1)*ned + i;
                        col = (b-1)*ned + k;
                        for j = 1:nsd
                            for l = 1:nsd
                                kint(row,col) = kint(row,col) + ...
                                    C(j,j,k,l)/nsd*dNdy(a,i)*dNdy(b,l)*weights(intpt)*dt;
                            end
                            kint(row,col) = kint(row,col) - ...
                                stress(j,j)/nsd*dNdy(a,k)*dNdy(b,i)*dt*weights(intpt);
                        end                      
                    end
                end
            end
        end
    end
    %% scatter the element kint into the global Kint
    %
    for a = 1:nen
        for i = 1:ned
            for b = 1:nen
                for k = 1:ned
                    row = ned*(connect(a,ele)-1) + i;
                    col = ned*(connect(b,ele)-1) + k;
                    K(row,col) = K(row,col) + ...
                        kint(ned*(a-1)+i,ned*(b-1)+k);
                end
            end
        end
    end
end
end