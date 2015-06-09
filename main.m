clear
clc
close all
infile1=fopen('input.txt','r');
infile2=fopen('hexa_test.inp','r');
outfile=fopen('output_hexa.txt','w');
[nsd,ned,nen,materialprops,gravity,nn,coords,nel,connect,boundary]=read_input(infile1,infile2);
%plotmesh(coords,nsd,connect,nel,nen,'r')
%% Boundary Conditions
bc1=zeros(3,1);
bc2=zeros(2+nsd,1);
subset1=[1];
subset2=[2];
% essential boundary conditions
for i=1:length(subset1)
    k=subset1(i);
    fprintf('Input the 3-digit essential boundary conditions for surface %d\n',k);
    fprintf('In 3D, 101 means fixing x and z displacements\n')
    prompt='In 2D, the unit digit 1 means plain strain\n';
    flag=input(prompt);
    if flag>=100
        temp=[boundary(k).nodes;ones(1,length(boundary(k).nodes));zeros(1,length(boundary(k).nodes))];
        bc1=[bc1,temp];
        flag=flag-100;
    end
    if flag>=10
        temp=[boundary(k).nodes;2*ones(1,length(boundary(k).nodes));zeros(1,length(boundary(k).nodes))];
        bc1=[bc1,temp];
        flag=flag-10;
    end
    if nsd==3
        if flag==1
            temp=[boundary(k).nodes;3*ones(1,length(boundary(k).nodes));zeros(1,length(boundary(k).nodes))];
            bc1=[bc1,temp];
        end
    else
        planestrain=flag;
    end
end
bc1=bc1(:,2:end);
% natural boundary conditions
for i=1:length(subset2)
    k=subset2(i);
    fprintf('Input the traction vector for surface %d\n',k);
    prompt='For example: [0.5,0,0] or [0.5,0]\n';
    flag=input(prompt);
    if nsd==3
        temp=[boundary(k).face2ele;flag(1)*ones(1,size(boundary(k).face2ele,2));...
            flag(2)*ones(1,size(boundary(k).face2ele,2));flag(3)*ones(1,size(boundary(k).face2ele,2))];
        bc2=[bc2,temp];
    else
        temp=[boundary(k).face2ele;flag(1)*ones(1,size(boundary(k).face2ele,2));...
            flag(2)*ones(1,size(boundary(k).face2ele,2))];
        bc2=[bc2,temp];
    end
end
bc2=bc2(:,2:end);
no_bc1=size(bc1,2);
no_bc2=size(bc2,2);
%% Call the solver
Statics(nsd,ned,nen,materialprops,gravity,nn,coords,nel,connect,no_bc1,bc1,no_bc2,bc2,outfile)
    
    
    