function [nn,coords,nel,connect]=read_input(nsd,nen,infile)
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
end