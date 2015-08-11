%function processing
nsd=3;
nen=8;
filename='vessel.inp';
infile=fopen(filename,'r');
%
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
%% set: name,node,element,dof
while(strncmp('*End Instance',content,13)==0)
    content=fgets(infile);
end
% count and read sets
nset=0;
while ~feof(infile)
    if strncmp('*Nset, nset=Set-',content,16)==1 
        nset=nset+1;
        set(nset).name=['Set-' num2str(nset)];
        content=strtrim(content); % remove the leading and trailing spaces
        flag=content((length(content)-7):end);
        if strncmp('generate',flag,8) % structured mesh
            content=fgets(infile);
            temp=str2num(content);
            set(nset).node=temp(1):temp(3):temp(2); 
            content=fgets(infile);
        else % unstructured mesh
            content=fgets(infile);
            set(nset).node=0;
            while(strncmp('*',content,1)==0)
                set(nset).node=[set(nset).node,str2num(content)];
                content=fgets(infile);
            end
            set(nset).node=set(nset).node(2:end);
        end
        content=strtrim(content); % remove the leading and trailing spaces
        flag=content((length(content)-7):end); 
        if strncmp('generate',flag,8) % structured mesh
            content=fgets(infile);
            temp=str2num(content);
            set(nset).element=temp(1):temp(3):temp(2); 
            content=fgets(infile);
        else % unstructured mesh
            content=fgets(infile);
            set(nset).element=0;
            while(strncmp('*',content,1)==0)
                set(nset).element=[set(nset).element,str2num(content)];
                content=fgets(infile);
            end
            set(nset).element=set(nset).element(2:end);
        end   
    elseif(strncmp('*Elset, elset=_Surf-',content,20)==1)
        break
    else    
        content=fgets(infile);
    end
end
fclose(infile);
infile=fopen(filename,'r');
%% surface: name, element, face, pressure, traction
while ~feof(infile)
    content=fgets(infile);
    if (strncmp('*End Instance',content,13)==1)
        break
    end
end
nsf=0;
while ~feof(infile)
    % Two possibilities: new surface, or new element face of a existing surface
    if strncmp('*Elset, elset=_Surf-',content,20)
        % new surface
        nsf=nsf+1;
        surface(nsf).id=[str2num(content(21)),str2num(content(24))];
        content=strtrim(content); % remove the leading and trailing spaces
        flag=content((length(content)-7):end);
        if strncmp('generate',flag,8) % structured mesh
            content=fgets(infile);
            temp=str2num(content);
            surface(nsf).number=(temp(2)-temp(1))/temp(3)+1;
            surface(nsf).element=temp(1):temp(3):temp(2); 
            content=fgets(infile);
        else % unstructured mesh
            content=fgets(infile);
            surface(nsf).element=0;
            while(strncmp('*',content,1)==0)
                surface(nsf).element=[surface(nsf).element,str2num(content)];
                content=fgets(infile);
            end
            surface(nsf).element=surface(nsf).element(2:end);
            surface(nsf).number=length(surface(nsf).element);
        end
        surface(nsf).face=surface(nsf).id(2)*ones(1,surface(nsf).number);
    elseif (strncmp('*End Assembly',content,13)==1)
        break
    else
        content=fgets(infile);
    end
end 
% form the supersurface
nsuper=0;
for i=1:nsf
    if surface(i).id(1)>nsuper
        nsuper=nsuper+1;
        super(nsuper).name=['Surf-' num2str(nsuper)];
        super(nsuper).element=0;
        super(nsuper).face=0;
        super(nsuper).element=[super(nsuper).element,surface(i).element];
        super(nsuper).face=[super(nsuper).face,surface(i).face];
    else
        id=surface(i).id(1);
        super(id).element=[super(id).element,surface(i).element];
        super(id).face=[super(id).face,surface(i).face];
    end
end
for i=1:nsuper
    super(i).element=super(i).element(2:end);
    super(i).face=super(i).face(2:end);
    super(i).pressure=1;
end
%% boundary conditions
isbc=0;
for i=1:nset
    set(i).dof=0;
end
while ~feof(infile)
    content=fgets(infile);
    if (strncmp('** BOUNDARY CONDITIONS',content,22)==1)
        isbc=1;
        content=fgets(infile);
        while (strncmp('** -',content,4)==0) 
            content=fgets(infile);
            if strncmp('*Boundary',content,9)
                while (strncmp('**',content,2)==0)
                    content=fgets(infile);
                    for i=1:nset
                        if strncmp(set(i).name,content,length(set(i).name))
                            if (str2num(content(length(set(i).name)+3)) == 1)
                                set(i).dof=set(i).dof+100;
                            elseif (str2num(content(length(set(i).name)+3)) == 2)
                                set(i).dof=set(i).dof+10;
                            elseif (str2num(content(length(set(i).name)+3)) == 3)
                                set(i).dof=set(i).dof+1;
                            end
                            set(i).value=str2num(content(length(set(i).name)+6));
                        end
                    end          
                end
            end
        end
        break
    end
end
fclose(infile);
%% pressure load
isload=0;
infile=fopen(filename,'r');
while ~feof(infile)
    content=fgets(infile);
    if strncmp('*Dsload',content,7)
        isload=1;
        content=fgets(infile);
        for i=1:nsuper
            if strncmp(super(i).name,content,length(super(i).name))
                super(i).pressure=str2num(content((length(super(i).name)+6):end));
            end
        end
    elseif strncmp('** OUTPUT',content,9)
        break
    end
end
fclose(infile);
%% write file: coords
outfile1=fopen('coords.txt','w');
fprintf(outfile1,'%10d\t%10d\n',nsd,nn);
if (nsd==3)
    for i=1:nn
        fprintf(outfile1,'%12.8f\t%12.8f\t%12.8f\n',coords(:,i));
    end
else
    for i=1:nn
        fprintf(outfile1,'%12.8f\t%12.8f\n',coords(:,i));
    end
end
fclose(outfile1);
%% write file: connect
outfile2=fopen('connect.txt','w');
fprintf(outfile2,'%10d\t%10d\n',nel,nen);
for i=1:nel
    for j=1:nen
        fprintf(outfile2,'%10d\t',connect(j,i));
    end
    fprintf(outfile2,'\n');
end
fclose(outfile2);
%% write file: bc
outfile3=fopen('bc.txt','w');
if (isbc==1)
    bc=[0;0;0.0];
    for i=1:nset
        if (set(i).dof>=100) % x is fixed
            bc=[bc,[set(i).node;ones(1,length(set(i).node));...
                set(i).value*ones(1,length(set(i).node))]];
            set(i).dof=set(i).dof-100;
        end
        if (set(i).dof>=10) % y is fixed
            bc=[bc,[set(i).node;2*ones(1,length(set(i).node));...
                set(i).value*ones(1,length(set(i).node))]];
            set(i).dof=set(i).dof-10;
        end
        if (set(i).dof==1) % z is fixed
            bc=[bc,[set(i).node;3*ones(1,length(set(i).node));...
                set(i).value*ones(1,length(set(i).node))]];
            set(i).dof=set(i).dof-1;
        end 
    end
    bc=bc(:,2:end);
    fprintf(outfile3,'%10d\n',size(bc,2));
    for i=1:size(bc,2)
        for j=1:2
            fprintf(outfile3,'%10d\t',bc(j,i));
        end
        fprintf(outfile3,'%12.8f',bc(3,i));
        fprintf(outfile3,'\n');
    end
else
    fprintf(outfile3,'%10d\n',0);
end
fclose(outfile3);
%% write file: load
outfile4=fopen('load.txt','w');
if (isload == 1)
    load=zeros(2+nsd,1);
    no_load=0;
    for i=1:nsuper
        no_load=no_load+length(super(i).element);
        super(i).norm=zeros(nsd,length(super(i).element));
        for j=1:length(super(i).element)
            face=super(i).face(j);
            temp=facenodes(nsd,nen,face);
            list=connect(temp,super(i).element(j));
            elecoord=coords(:,list);
            if nsd==3
                norm=cross(elecoord(:,2)-elecoord(:,1),elecoord(:,3)-elecoord(:,2));
            else
                a=elecoord(:,2)-elecoord(:,1);
                norm=[-a(2),a(1)];  
            end
            norm=norm/sqrt(dot(norm,norm));
            super(i).norm(:,j)=transpose(norm);
        end  
    end
    fprintf(outfile4,'%10d\n',no_load);
    for i=1:nsuper
        for j=1:length(super(i).element)
            fprintf(outfile4,'%10d\t',super(i).element(j),super(i).face(j));
            fprintf(outfile4,'%12.8f\t',super(i).pressure*super(i).norm(:,j));
            fprintf(outfile4,'\n');            
        end
    end    
else
    fprintf(outfile4,'%10d\n',0);
end
fclose(outfile4);
%end
%{
%% auxilary functions
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
            list = [2,6,7,3];
        elseif face==5
            list = [3,7,8,4];
        elseif face==6
            list = [4,8,5,1];
        end
    end
end
end
%}