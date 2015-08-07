infile=fopen('load.txt','r');
content=fgets(infile);
i=0;
while ~feof(infile)
    content=fgets(infile);
    i=i+1;
    temp=str2num(content);
    a(i,:)=temp(1,3:4);
end
plot(a(:,1),a(:,2),'*')