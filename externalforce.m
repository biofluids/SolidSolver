function Fext = externalforce(nsd,ned,nn,nel,nen,nh,materialprops,gravity,coords,connect,hnodes)
% returns the external force
Fext = zeros(ned*nn,1);
%% force due to traction
%%
% loop over the faces with prescribed traction
for n = 1:nh
    % extract the coords and traction of that face
    ele = hnodes(1,n);
    face = hnodes(2,n);
    nfacenodes=no_facenodes(nsd,nen);
    nodelist = facenodes(nsd,nen,face); 
    coord = zeros(nsd,nfacenodes);
    for a = 1:nfacenodes
        coord(:,a) = coords(:,connect(nodelist(a),ele));
    end
    traction = zeros(ned,1);
    for i = 1:nsd
        traction(i) = hnodes(i+2,n);
    end
    %% compute the force on the face
    %
    f = zeros(nfacenodes*ned,1);
    % set up quadrature points and weights
    npt = no_quad_pts(nsd-1,nfacenodes,0);
    xilist = quad_pts(nsd-1,nfacenodes,npt);
    weights = quad_weights(nsd-1,nfacenodes,npt);
    % loop over the quadrature points
    for intpt = 1:npt
        xi = xilist(:,intpt);
        dNdxi = sfder(nfacenodes,nsd-1,xi);
        N = sf(nfacenodes,nsd-1,xi);
        % set up the jacobian matrix
        dxdxi = zeros(nsd,nsd-1);
        dxdxi = coord*dNdxi;
        if nsd == 2
            dt = sqrt( (dxdxi(1,1))^2 + (dxdxi(2,1)^2) );
        elseif nsd == 3
            dt = sqrt( ((dxdxi(2,1)*dxdxi(3,2))-(dxdxi(2,2)*dxdxi(3,1)))^2 ...
                + ((dxdxi(1,1)*dxdxi(3,2))-(dxdxi(1,2)*dxdxi(3,1)))^2 ...
                + ((dxdxi(1,1)*dxdxi(2,2))-(dxdxi(1,2)*dxdxi(2,1)))^2 );
        end
        % compute the force
        for a = 1:nfacenodes
            for i = 1:nsd
                row = ned*(a-1)+i;
                f(row) = f(row) + traction(i)*N(a)*dt*weights(intpt);
            end
        end
    end       
    %% scatter the force into the global external force
    %
    for a = 1:nfacenodes
        for i = 1:ned
            row = (connect(nodelist(a),ele)-1)*ned + i;
            Fext(row) = Fext(row) + f((a-1)*ned+i);
        end
    end
end
%% force due to gravity
%%
elecoord = zeros(nsd,nen);
% loop over the elements
for ele = 1:nel
    % Extract coords of nodes, and dof for the element
    for a = 1:nen
        elecoord(:,a) = coords(:,connect(a,ele));
    end
    %% compute the force on this element
    %
    f = zeros(ned*nen,1);
    rho = materialprops(5);
    %  Set up integration points and weights
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
        % compute the force
        for a = 1:nen
            for i = 1:ned
                row = ned*(a-1)+i;
                f(row,1) = f(row,1) + gravity(i)*rho*N(a)*weights(intpt)*dt;
            end
        end
    end
    %% scatter the force into the global external force
    %
    for a = 1:nen
        for i = 1:ned
            row = ned*(connect(a,ele)-1)+i;
            Fext(row,1) = Fext(row,1) + f(ned*(a-1)+i,1);   
        end
    end
end
end 