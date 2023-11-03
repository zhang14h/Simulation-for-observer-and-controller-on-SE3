mw=0.089   %Mass of weight (kg)
mr=0.069    %Mass of rod (kg)
Lr=0.43; %Total length of rod (m)
c1=0.02; %Friction coefficient for base
Rw=0.050; %Radius of round weight (m), only used for detailed expressions of Jy and Jz_bar
g=9.81; %Acceleration of gravity (m/s^2)
yr=.42;
ym=.32;
J1=0.0166; %Inertia of lower disk including Pendulum base hardware
Rh=0.25; %Distance from disk rotation axis to hub of pendulum (m)

lw = ym;
m=mr+mw;
lr = yr-1/2*Lr;
lcg=(mr*lr+mw*lw)/m;
Jx=1/12*mr*Lr^2+1/2*mw*Rw^2;
Jy=1/4*mw*Rw^2;
Jz=1/12*mr*Lr^2+1/2*mw*Rw^2;

J1_bar = J1+m*Rh^2;
Jx_bar = Jx+m*lcg^2;
Jz_bar = Jz+m*lcg^2;


p=Jz_bar*(J1_bar+Jy)-(m*Rh*lcg)^2;




A=[0 1 0 0; 
   0 -c1*Jz_bar/p -m^2*lcg^2*Rh*g/p 0; 
   0 0 0 1; 
   0 c1*m*Rh*lcg/p m*lcg*g*(J1_bar+Jy)/p 0]

B=1/p*[0; 
    Jz_bar; 
    0; 
    -m*Rh*lcg]*16

C=[1 0 0 0; 
   0 0 1 0];

D=[0;
   0];


[b,a] = ss2tf(A,B,C,D,1);
theta1= tf(b(1,:),a)
theta2= tf(b(2,:),a)


Q = C'*C;
R = 1;
K = lqr(A,B,Q,R)

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'r'};
outputs = {'x'; 'phi'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
r =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with LQR Control')



