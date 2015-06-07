clear
clc
close all
infile=fopen('hexa.inp','r');
%infile=fopen('tetra.inp','r');
[nn,coords,nel,connect,set1,set2nodes,set2eles]=read_input(3,8,infile);
%[nn,coords,nel,connect,set1,set2nodes,set2eles]=read_input(3,4,infile);