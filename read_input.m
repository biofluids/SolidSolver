function [nn,coords,nel,connect,set1,set2nodes,set2eles]=read_input(nsd,nen,infile)
content=fgets(infile);
while(strncmp('*Node',content,5)==0)
    content=fgets(infile);
end
content=fgets(infile);
% coords
while(strncmp('*',content,1)==0)
    temp=str2num(content);
    nn=temp(1);
    coords(:,nn)=transpose(temp(2:(nsd+1)));
    content=fgets(infile);
end
% connect
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
%
% set1: nodes with essential bc
set1=0;
%
while(strncmp('*Nset, nset=Set-1',content,17)==0)
    content=fgets(infile);
end
content=strtrim(content); % remove the trailing spaces
content=content((length(content)-7):end); % last 8 characters
if strncmp('generate',content,8) % structural mesh
    content=fgets(infile);
    temp=str2num(content);
    set1=temp(1):temp(2);
else % unstructural mesh
    while(strncmp('*',content,1)==0)
        content=fgets(infile);
        set1=[set1,str2num(content)];
    end
    set1=set1(2:end);
end
% set2: nodes with natural bc
% set2nodes
set2nodes=0;
while(strncmp('*Nset, nset=Set-2',content,17)==0)
    content=fgets(infile);
end
content=strtrim(content); % remove the trailing spaces
content=content((length(content)-7):end); % last 8 characters
if strncmp('generate',content,8) % structural mesh
    content=fgets(infile);
    temp=str2num(content);
    set2nodes=temp(1):temp(2);
else % unstructural mesh
    while(strncmp('*',content,1)==0)
        content=fgets(infile);
        set2nodes=[set2nodes,str2num(content)];
    end
    set2nodes=set2nodes(2:end);
end
% set2elements
set2eles=0;
while(strncmp('*Elset, elset=Set-2',content,18)==0)
    content=fgets(infile);
end
content=strtrim(content); % remove the trailing spaces
content=content((length(content)-7):end); % last 8 characters
if strncmp('generate',content,8) % structural mesh
    content=fgets(infile);
    temp=str2num(content);
    set2eles=temp(1):temp(2);
else % unstructural mesh
    while(strncmp('*',content,1)==0)
        content=fgets(infile);
        set2eles=[set2eles,str2num(content)];
    end
    set2eles=set2eles(2:end);
end
end


    

