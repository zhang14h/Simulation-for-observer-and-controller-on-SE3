clc;
clear;
close all;
[x,y,z] = meshgrid(-10:.5:10,-10:.5:10,-10:.5:10);
a = 5;
b = 3;
c = 2;

d = x.^2/a^2 + y.^2/b^2 + z.^2/c^2;

surf(x,y,z,d)
colorbar