function [simu_type,tol,maxit,relax,nsteps,dt,nprint,damp,nsd,ned,nen,materialprops,gravity,nn,coords,nel,connect,boundary]=read_input(infile1,infile2)
%% read input file

cellarray=textscan(infile1,'%s');

% simulation parameters

cellno=6;
simu_type=str2num(cellarray{1}{cellno});
cellno=cellno+2;
tol=str2num(cellarray{1}{cellno});
cellno=cellno+6;
maxit=str2num(cellarray{1}{cellno});
cellno=cellno+2;
relax=str2num(cellarray{1}{cellno});
cellno=cellno+7;
nsteps=str2num(cellarray{1}{cellno});
cellno=cellno+3;
dt=str2num(cellarray{1}{cellno});
cellno=cellno+3;
nprint=str2num(cellarray{1}{cellno});
cellno=cellno+2;
damp=str2num(cellarray{1}{cellno});
if simu_type==0
    nsteps=1;
    dt=1.;
    nprint=1; 
end

% mesh parameters

cellno=cellno+5;
nsd=str2num(cellarray{1}{cellno});
ned=nsd;
cellno=cellno+7;
nen=str2num(cellarray{1}{cellno});

% material parameters

materialprops=zeros(5,1);
cellno=cellno+6;
materialprops(1)=str2num(cellarray{1}{cellno});
cellno=cellno+2;
for i=1:4
    cellno=cellno+2;
    materialprops(i+1)=str2num(cellarray{1}{cellno});
end
gravity=zeros(3,1);
cellno=cellno+1;
for i=1:3
    cellno=cellno+1;
    gravity(i)=str2num(cellarray{1}{cellno});
end
fclose(infile1);

%% read mesh file

content=fgets(infile2);
% number of nodes and coords
while(strncmp('*Node',content,5)==0)
    content=fgets(infile2);
end
content=fgets(infile2);
while(strncmp('*',content,1)==0)
    temp=str2num(content);
    nn=temp(1);
    coords(:,nn)=transpose(temp(2:(nsd+1)));
    content=fgets(infile2);
end
% number of elements and connect
while(strncmp('*Element',content,8)==0)
    content=fgets(infile2);
end
content=fgets(infile2);
while(strncmp('*',content,1)==0)
    temp=str2num(content);
    nel=temp(1);
    connect(:,nel)=transpose(temp(2:(1+nen)));
    content=fgets(infile2);
end
% number of boundaries
nb=0;
while ~feof(infile2)
    if strncmp('*Nset, nset=Set-',content,16)==1 
        % new set detected
        nb=nb+1;
        % name of the new set
        boundary(nb).name=['Set-' num2str(nb)];
        content=strtrim(content); % remove the leading and trailing spaces
        flag=content((length(content)-7):end); % last 8 characters
        % nodes in the new set
        if strncmp('generate',flag,8) % structural mesh
            content=fgets(infile2);
            temp=str2num(content);
            % nodes in the new set
            boundary(nb).nodes=temp(1):temp(3):temp(2); 
            content=fgets(infile2);
        else % unstructural mesh
            content=fgets(infile2);
            boundary(nb).nodes=0;
            while(strncmp('*',content,1)==0)
                boundary(nb).nodes=[boundary(nb).nodes,str2num(content)];
                content=fgets(infile2);
            end
            boundary(nb).nodes=boundary(nb).nodes(2:end);
        end
        % elements in the new set
        content=strtrim(content); % remove the leading and trailing spaces
        flag=content((length(content)-7):end); % last 8 characters
        if strncmp('generate',flag,8) % structural mesh
            content=fgets(infile2);
            temp=str2num(content);
            % nodes in the new set
            boundary(nb).elements=temp(1):temp(3):temp(2); 
            content=fgets(infile2);
        else % unstructural mesh
            content=fgets(infile2);
            boundary(nb).elements=0;
            while(strncmp('*',content,1)==0)
                boundary(nb).elements=[boundary(nb).elements,str2num(content)];
                content=fgets(infile2);
            end
            boundary(nb).elements=boundary(nb).elements(2:end);
        end   
    else
        content=fgets(infile2);
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
fclose(infile2);
end
    