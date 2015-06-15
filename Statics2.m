function Statics2(nsteps,dt,nprint,maxit,tol,relax,damp,nsd,ned,nen,materialprops,gravity,nn,coords,nel,connect,no_bc1,bc1,no_bc2,bc2)
%% MAIN FEM ANALYSIS PROCEDURE 
% the load is applied step by step
% Augmented Lagrangian Method is NOT used in this code
% reduced integration used in Fint and Kint

w = zeros(ned*nn,1);

% original values
write_results(nsteps,nprint,dt,nsd,ned,nn,coords,nel,nen,connect,materialprops,w,0);
Fext = externalforce(nsd,ned,nn,nel,nen,no_bc2,materialprops,gravity,coords,connect,bc2);
loadfactor = 0.;
increment=0.05;
step=1;
while (loadfactor < 1)
    if (loadfactor+increment>1)
        increment=1-loadfactor;
    end
    w1 = w;
    loadfactor = loadfactor + increment;
    err1 = 1.;
    err2 = 1.;
    nit = 0;
    fprintf('\n Step %f Load %f\n', step, loadfactor);
    
    while (((err1>tol)||(err2>tol)) && (nit<maxit))
        nit = nit + 1;
        Fint = internalforce(nsd,ned,nn,coords,nel,nen,connect,materialprops,w);
        A = Kint(nsd,ned,nn,coords,nel,nen,connect,materialprops,w);
        R = Fint - loadfactor*Fext;
        % fix the prescribed displacements
        for n = 1:no_bc1
            row = ned*(bc1(1,n)-1) + bc1(2,n);
            for col = 1:ned*nn
                A(row,col) = 0;
                A(col,row) = 0;
            end
            A(row,row) = 1.;
            R(row) = -bc1(3,n) + w(row); 
        end
        % solve for the correction
        dw = A\(-R);
        % check for convergence
        w = w + relax*dw;
        wnorm = dot(w,w);
        err1 = dot(dw,dw);
        err2 = dot(R,R);
        err1 = sqrt(err1/wnorm);
        err2 = sqrt(err2)/(ned*nn);
        fprintf('\n Iteration number %d Correction %8.3e Residual %8.3e tolerance %8.3e\n',nit,err1,err2,tol);
    end
    step=step+1;
    if nit == maxit
        w = w1;
        loadfactor = loadfactor - increment;
        increment=increment/2;
    elseif nit>8
        increment=increment*(1-0.1);
    elseif nit<8
        increment=increment*(1+0.1);
    end
end
% writing the output file
write_results(nsteps,nprint,dt,nsd,ned,nn,coords,nel,nen,connect,materialprops,w,1);
%% plot the original and deformed mesh
coords1 = zeros(ned,nn);
for i = 1:nn
    for j = 1:ned
        coords1(j,i) = coords(j,i) + w(ned*(i-1)+j);
    end
end
% plot the undeformed and deformed mesh
figure
plotmesh(coords,nsd,connect,nel,nen,'g');
hold on
plotmesh(coords1,nsd,connect,nel,nen,'r');
end