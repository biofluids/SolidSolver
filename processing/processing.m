function processing
clear
clc
close all
infile=fopen('hexa_test.inp','r');
nsd=3;
nen=8;
subset1=[1]; % input the id of sets that you want to apply essential bc on
subset2=[2]; % input the id of sets that you want to apply natural bc on
essential=111; % 110 represents fixing x and y; in 2d, the unit digit should be 1
natural=[0,0.1,0];
%
%
%%
content=fgets(infile);
% number of nodes and coords
while(strncmp('*Node',content,5)==0)
    content=fgets(infile);
end
content=fgets(infile);
while(strncmp('*',content,1)==0)
    temp=str2num(content);
    nn=temp(1);
    coords(:,nn)=transpose(temp(2:(nsd+1)));
    content=fgets(infile);
end
% number of elements and connect
while(strncmp('*Element',content,8)==0)
    content=fgets(infile);
end
content=fgets(infile);
while(strncmp('*',content,1)==0)
    temp=str2num(content);
    nel=temp(1);
    connect(:,nel)=transpose(temp(2:(1+nen)));
    content=fgets(infile);
end
% number of boundaries
nb=0;
while ~feof(infile)
    if strncmp('*Nset, nset=Set-',content,16)==1 
        % new set detected
        nb=nb+1;
        % name of the new set
        boundary(nb).name=['Set-' num2str(nb)];
        content=strtrim(content); % remove the leading and trailing spaces
        fix=content((length(content)-7):end); % last 8 characters
        % nodes in the new set
        if strncmp('generate',fix,8) % structural mesh
            content=fgets(infile);
            temp=str2num(content);
            % nodes in the new set
            boundary(nb).nodes=temp(1):temp(3):temp(2); 
            content=fgets(infile);
        else % unstructural mesh
            content=fgets(infile);
            boundary(nb).nodes=0;
            while(strncmp('*',content,1)==0)
                boundary(nb).nodes=[boundary(nb).nodes,str2num(content)];
                content=fgets(infile);
            end
            boundary(nb).nodes=boundary(nb).nodes(2:end);
        end
        % elements in the new set
        content=strtrim(content); % remove the leading and trailing spaces
        fix=content((length(content)-7):end); % last 8 characters
        if strncmp('generate',fix,8) % structural mesh
            content=fgets(infile);
            temp=str2num(content);
            % nodes in the new set
            boundary(nb).elements=temp(1):temp(3):temp(2); 
            content=fgets(infile);
        else % unstructural mesh
            content=fgets(infile);
            boundary(nb).elements=0;
            while(strncmp('*',content,1)==0)
                boundary(nb).elements=[boundary(nb).elements,str2num(content)];
                content=fgets(infile);
            end
            boundary(nb).elements=boundary(nb).elements(2:end);
        end   
    else
        content=fgets(infile);
    end
end
% face2ele
for n=1:nb
    boundary(n).face2ele=zeros(2,1);
    for i=1:length(boundary(n).elements)
        ele=boundary(n).elements(i);
        node=connect(:,ele);
        for face=1:(no_faces(nsd,nen))
            list=facenodes(nsd,nen,face);
            temp=transpose(node(list));
            if isequal(ismember(temp,boundary(n).nodes),ones(1,length(temp)))
                boundary(n).face2ele=[boundary(n).face2ele,[ele;face]];
            end
        end
    end
    boundary(n).face2ele=boundary(n).face2ele(:,2:end);
end
fclose(infile);
%
%
%%
bc1=zeros(3,1);
bc2=zeros(2+nsd,1);
% essential boundary conditions
for i=1:length(subset1)
    k=subset1(i);
    %fprintf('Input the 3-digit essential boundary conditions for surface %d\n',k);
    %fprintf('In 3D, 101 means fixing x and z displacements\n')
    %prompt='In 2D, the unit digit 1 means plain strain\n';
    %fix=input(prompt);
    if essential>=100
        temp=[boundary(k).nodes;ones(1,length(boundary(k).nodes));zeros(1,length(boundary(k).nodes))];
        bc1=[bc1,temp];
        essential=essential-100;
    end
    if essential>=10
        temp=[boundary(k).nodes;2*ones(1,length(boundary(k).nodes));zeros(1,length(boundary(k).nodes))];
        bc1=[bc1,temp];
        essential=essential-10;
    end
    if nsd==3
        if essential==1
            temp=[boundary(k).nodes;3*ones(1,length(boundary(k).nodes));zeros(1,length(boundary(k).nodes))];
            bc1=[bc1,temp];
        end
    else
        planestrain=essential;
    end
end
bc1=bc1(:,2:end);
% natural boundary conditions
for i=1:length(subset2)
    k=subset2(i);
    %fprintf('Input the traction vector for surface %d\n',k);
    %prompt='For example: [0.5,0,0] or [0.5,0]\n';
    %traction=input(prompt);
    if nsd==3
        temp=[boundary(k).face2ele;natural(1)*ones(1,size(boundary(k).face2ele,2));...
            natural(2)*ones(1,size(boundary(k).face2ele,2));natural(3)*ones(1,size(boundary(k).face2ele,2))];
        bc2=[bc2,temp];
    else
        temp=[boundary(k).face2ele;natural(1)*ones(1,size(boundary(k).face2ele,2));...
            natural(2)*ones(1,size(boundary(k).face2ele,2))];
        bc2=[bc2,temp];
    end
end
bc2=bc2(:,2:end);
no_bc1=size(bc1,2);
no_bc2=size(bc2,2);
%%
% coords
outfile1=fopen('coords.txt','w');
fprintf(outfile1,'%10d\n',nn);
for i=1:nn
    fprintf(outfile1,'%12.8f\t%12.8f\t%12.8f\n',coords(:,i));
end
fclose(outfile1);
% connect
outfile2=fopen('connect.txt','w');
fprintf(outfile2,'%10d\t%10d\n',nel,nen);
for i=1:nel
    for j=1:nen
        fprintf(outfile2,'%10d\t',connect(j,i));
    end
    fprintf(outfile2,'\n');
end
fclose(outfile2);
% bc1
outfile3=fopen('bc1.txt','w');
fprintf(outfile3,'%10d\n',no_bc1);
for i=1:no_bc1
    for j=1:3
        fprintf(outfile3,'%10d\t',bc1(j,i));
    end
    fprintf(outfile3,'\n');
end
fclose(outfile3);
% bc2
outfile4=fopen('bc2.txt','w');
fprintf(outfile4,'%10d\n',no_bc2);
for i=1:no_bc2
    for j=1:5
        fprintf(outfile4,'%12.8f\t',bc2(j,i));
    end
    fprintf(outfile4,'\n');
end
fclose(outfile4);

end
%% ========================= functions ====================================
%
function list=facenodes(nsd,nen,face)
list = zeros(no_facenodes(nsd,nen),1);
i3 = [2,3,1];
i4 = [2,3,4,1];
if nsd==2
    if nen==3
        list = [face;i3(face)];
    elseif nen==4
        list = [face;i4(face)];
    end
elseif nsd==3
    if nen==4
        if face==1
            list = [1,2,3];
        elseif face==2
            list = [1,4,2];
        elseif face==3
            list = [2,4,3];
        elseif face==4
            list = [3,4,1];
        end
    elseif nen==8
        if face==1
            list = [1,2,3,4];
        elseif face==2
            list = [5,8,7,6];
        elseif face==3
            list = [1,5,6,2];
        elseif face==4
            list = [2,3,7,6];
        elseif face==5
            list = [3,7,8,4];
        elseif face==6
            list = [4,8,5,1];
        end
    end
end
end

function n=no_facenodes(nsd,nen)
if nsd==2
    n=2;
elseif nsd==3
    if nen==4
        n=3;
    elseif nen==8
        n=4;
    end
end
end

function n=no_faces(nsd,nen)
if nsd==2
    if nen==3
        n=3;
    elseif nen==4
        n=4;
    end
elseif nsd==3
    if nen==4
        n=4;
    elseif nen==8
        n=6;
    end
end
end