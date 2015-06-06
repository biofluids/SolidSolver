function Statics
% main function for static analysis
% reduced integration used in Fint and Kint
clear
clc
close all
%% read input file
outfile=fopen('output2d4.txt','w');
%{
nsd=2;
ned=2;
nen=4;
a=6;
b=1;
h=0.2;
[nn,nel,connect,coords]=meshgenerator(nsd,nen,a,b,0,h);
%plotmesh(coords,nsd,connect,nel,nen,'g');
gnodes = [31*(0:b/h)+1,31*(0:b/h)+1;ones(1,b/h+1),2*ones(1,b/h+1);zeros(1,2*(b/h+1))];
ng=size(gnodes,2);
%hnodes = [2*(a/h)*(1:(b/h))-1;2*ones(1,b/h);zeros(1,b/h);0.005*ones(1,b/h);];
hnodes = [(a/h)*(1:b/h);2*ones(1,b/h);zeros(1,b/h);0.005*ones(1,b/h)];
nh=size(hnodes,2);
gravity=[0.;0.];
materialprops=[2;1;15000;1000;1];
%}
%3d hexatedral
%{
x=[0:3,0:3,0:3,0:3,0:3,0:3];
y=[-1*ones(1,4),zeros(1,4),ones(1,4),-1*ones(1,4),zeros(1,4),ones(1,4)];
z=[zeros(1,12),ones(1,12)];
coords=[transpose(x),transpose(y),transpose(z)];
coords=transpose(coords);
connect=[1,2,3,5,6,7;2,3,4,6,7,8;6,7,8,10,11,12;...
    5,6,7,9,10,11;13,14,15,17,18,19;14,15,16,18,19,20;...
    18,19,20,22,23,24;17,18,19,21,22,23];
nsd=3;
ned=3;
nen=8;
nn=24;
nel=6;
ng=9;
gnodes=[1,1,1,13,13,5,17,21,9;1,2,3,1,2,1,1,1,1;zeros(1,9)];
nh=2;
hnodes=[3,6;4,4;3.,3.;0.,0.;0.,0.];
gravity=[0;0;0];
materialprops=[2;1;0;10;1];
%}
% 3d tetradetral
%
coords=[0,1,1,0,0,1,1,0;0,0,1,1,0,0,1,1;0,0,0,0,1,1,1,1];
connect=[5,5,1,4,8,8;4,8,4,3,7,7;1,4,2,2,4,3;6,6,6,6,6,6];
nsd=3;
ned=3;
nen=4;
nel=6;
nn=8;
ng=7;
nh=2;
gnodes=[1,1,1,5,5,8,4;1,2,3,1,2,1,1;zeros(1,7)];
hnodes=[4,6;3,3;3,3;0,0;0,0];
gravity=[0;0;0];
materialprops=[2;1;0;10;1];
%
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
plotmesh(coords,nsd,connect,nel,nen,'g');
hold on
plotmesh(coords1,nsd,connect,nel,nen,'r');
fclose(outfile);
end