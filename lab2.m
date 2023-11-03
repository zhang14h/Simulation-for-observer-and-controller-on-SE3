clc;
clear;
close all;
num = [1600];
den = [1 0 0.83];
g  = tf(num,den);
numc = [0.04674 1];
denc = [0.0041 1];
gc = tf(numc,denc);
[Gm,PM,Wcg,Wcp] = margin(g);
margin(g);
figure()
margin(gc*g);
