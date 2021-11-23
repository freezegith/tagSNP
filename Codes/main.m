clc;
clear;
close all;

%adjacent_matrix = load('.\dataset\karate.dat');
%[Z, H] = FastNewman(adjacent_matrix);



adjacent_matrix=load('LD.txt');
[m,n]=size(adjacent_matrix);              %行数、列数
%{
for i=1:m
   for j=1:n
       if adjacent_matrix(i,j)>0.8
            adjacent_matrix(i,j) = 1;
       else adjacent_matrix(i,j) = 0;
       end
   end
end
%}
    [m,n]=FastNewman(adjacent_matrix);



