clear
clc
close all
infile=fopen('hexa.inp','r');
%infile=fopen('tetra.inp','r');
nsd=3;
nen=8;
[nn,coords,nel,connect,ng,gnodes,nh,hnodes]=read_input(3,8,infile);
%[nn,coords,nel,connect,ng,gnodes,nh,hnodes]=read_input(3,4,infile);
plotmesh(coords,nsd,connect,nel,nen,'r')
fclose(infile)