function Statics
% main function for static analysis
% reduced integration used in Fint and Kint
clear
clc
close all
%% read input file
outfile=fopen('output2.txt','w');
nsd=2;
ned=2;
nen=3;
a=3;
b=1;
h=0.1;
[nn,nel,connect,coords]=meshgenerate(nsd,nen,a,b,h);
%plotmesh(coords,nsd,nn,connect,nel,nen,'g');
gnodes = [31*(0:b/h)+1,31*(0:b/h)+1;ones(1,b/h+1),2*ones(1,b/h+1);zeros(1,2*(b/h+1))];
ng=size(gnodes,2);
hnodes = [2*(a/h)*(1:(b/h))-1;2*ones(1,b/h);zeros(1,b/h);0.005*ones(1,b/h);];
%hnodes = [(a/h)*(1:b/h);2*ones(1,b/h);zeros(1,b/h);0.005*ones(1,b/h)];
nh=size(hnodes,2);
gravity=[0;0];
materialprops=[2;1;15000;10;1];

%% MAIN FEM ANALYSIS PROCEDURE 
% the load is applied step by step
% Augmented Lagrangian Method is NOT used in this code

tol = 0.0001;
maxit = 30;
relax = 1.;
nsteps = 10;

w = zeros(ned*nn,1);
for step = 1:nsteps
    loadfactor = step/nsteps;
    err1 = 1.;
    err2 = 1.;
    nit = 0;
    fprintf('\n Step %f Load %f\n', step, loadfactor);
    Fext = loadfactor*externalforce(nsd,ned,nn,nel,nen,nh,materialprops,gravity,coords,connect,hnodes);
    while (((err1>tol)||(err2>tol)) && (nit<maxit))
        nit = nit + 1;
        Fint = internalforce(nsd,ned,nn,coords,nel,nen,connect,materialprops,w);
        A = Kint(nsd,ned,nn,coords,nel,nen,connect,materialprops,w);
        R = Fint - Fext;
        % fix the prescribed displacements
        for n = 1:ng
            row = ned*(gnodes(1,n)-1) + gnodes(2,n);
            for col = 1:ned*nn
                A(row,col) = 0;
                A(col,row) = 0;
            end
            A(row,row) = 1.;
            R(row) = -gnodes(3,n) + w(row); 
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
    % writing the output file
    fprintf(outfile,'==================================================\n');
    fprintf(outfile,'Step %f Load %f\n',step,loadfactor);
    write_output(outfile,nsd,ned,nn,coords,nel,nen,connect,materialprops,w);
end
%% plot the original and deformed mesh
coords1 = zeros(ned,nn);
for i = 1:nn
    for j = 1:ned
        coords1(j,i) = coords(j,i) + w(ned*(i-1)+j);
    end
end
% plot the undeformed and deformed mesh
figure
plotmesh(coords,nsd,nn,connect,nel,nen,'g');
hold on
plotmesh(coords1,nsd,nn,connect,nel,nen,'r');
fclose(outfile);
end