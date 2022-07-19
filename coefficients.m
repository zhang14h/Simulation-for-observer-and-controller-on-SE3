clear;
clc;
close all; 

mw=0.089; %Mass of weight (kg)
mr=0.069; %Mass of rod (kg)
yw= 0.32; %distance from rod pivot center to weight center
yr= 0.42; %distance from rod pivot center to end of rod
Lr=0.43; %Total length of rod (m)
c1=0.02; %Friction coefficient for base
Rw=0.050; %Radius of round weight (m), only used for detailed expressions of Jy and Jz_bar
g=9.81; %Acceleration of gravity (m/s^2)
J1=0.0166; %Inertia of lower disk including Pendulum base hardware
Rh=0.25; %Distance from disk rotation axis to hub of pendulum (m)
c2 = 0.01;


lr = yr - 1/2*Lr;
lw = yw;
m = mr + mw;
lcg = (mr*lr + mw*lw)/m;

Jx = 1/12*mr*Lr^2 + 1.4*mw*Rw^2;
Jy = 1/4*mw*Rw^2;
Jz = 1/12*mr*Lr^2 + 1/2*mw*Rw^2;

J1_bar = J1+m*Rh^2;
Jx_bar = Jx+m*lcg^2;
Jz_bar = Jz+m*lcg^2;

p = Jz_bar*(J1_bar+Jy) - (m*Rh*lcg)^2;

A = [0       1             0                      0;
     0 -c1*Jz_bar/p   -m^2*lcg^2*Rh*g/p           0;
     0       0             0                      1;
     0 c1*m*Rh*lcg/p  +m*lcg*g*(J1_bar+Jy)/p       0]
 
B = [0;
     Jz_bar/p;
     0;
     -m*Rh*lcg/p]
 
C = [1 0 0 0;
     0 1 0 0;
     0 0 1 0
     0 0 0 1];
 
D = [0;0;0;0];


[b,a] = ss2tf(A,B,C,D,1);
tf1= tf(b(1,:),a)
tf2= tf(b(2,:),a)

pole = [-1 -2 -3 -4];
K = place(A,B,pole);
eig(A-B*K)