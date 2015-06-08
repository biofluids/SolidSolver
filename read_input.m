function [nn,coords,nel,connect,boundary]=read_input(nsd,nen,infile)
content=fgets(infile);
%% number of nodes and coords
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
%% number of elements and connect
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
%% number of boundaries
nb=0;
while ~feof(infile)
    if strncmp('*Nset, nset=Set-',content,16)==1 
        % new set detected
        nb=nb+1;
        % name of the new set
        boundary(nb).name=['Set-' num2str(nb)];
        content=strtrim(content); % remove the leading and trailing spaces
        flag=content((length(content)-7):end); % last 8 characters
        % nodes in the new set
        if strncmp('generate',flag,8) % structural mesh
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
        flag=content((length(content)-7):end); % last 8 characters
        if strncmp('generate',flag,8) % structural mesh
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
%% face2ele
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
end
    