function M = mass(nsd,ned,nn,coords,nel,nen,connect,materialprops)
% returns the global mass matrix
M = zeros(nn*ned,nn*ned);
elecoord = zeros(nsd,nen);
%%
% loop over elements
for ele = 1:nel
    % Extract coords of nodes, DOF for the current element
    for a = 1:nen
        elecoord(:,a) = coords(:,connect(a,ele));
    end
    %% compute the element mass matrix
    %
    %
    m = zeros(nen*ned,nen*ned);
    rho = materialprops(5);
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
        dxdxi = elecoord*dNdxi;
        dt = det(dxdxi);
        % Compute the element mass, no loop of j because M = delta_ij*M^~       
        for a = 1:nen
            for b = 1:nen
                for i = 1:ned
                    row = ned*(a-1)+i;
                    col = ned*(b-1)+i;
                    m(row,col) = m(row,col) + N(a)*N(b)*rho*dt*weights(intpt);  
                end
            end
        end
    end  
    %% scatter the element mass matrix into the global mass matrix
    %
    %
    for a = 1:nen
      for i = 1:ned
        for b = 1:nen
          for k = 1:ned
            row = ned*(connect(a,ele)-1)+i;
            col = ned*(connect(b,ele)-1)+k;
            M(row,col) = M(row,col) + m(ned*(a-1)+i,ned*(b-1)+k);
          end
        end
      end
    end
end
end