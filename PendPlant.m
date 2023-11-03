
%PendPlant.m
% This script builds the numerical transfer function and state space models for the ECP
% pendulum accessory.  It is compatible with plants created by use of the trajectoru in
% conjunction with ECP models 205, 210, 220, and 750.  Inverted and noninverted dynamic 
% models may be generated.

%%ENTER SYSTEM TYPE 
%(1= MODEL 205 [TORSION], 2= MODEL 210 [LINEAR], 3= MODEL 220 [INDUSTRIAL], 4= MODEL 750[GYRO])
systype=1


%% INITIALIZE
% The following user-specified values vary with plant configuration
% Inverted v. noninverted operation
INV=1 % Set =+1 for inverted operation, =-1 for non-inverted

%ADJUSTABLE DYNAMIC PARAMETERS
yw=INV*.32;	%distance from rod pivot center to weight center
yr=INV*.42;	%distance from rod pivot center to end of rod

%Ksys=1

%FIXED PARAMETERS
mw=0.089;	%Mass of weight (kg)
mr=.069;%;	%Mass of rod (kg)
Lr=INV*.43; %Total length of rod (m)
c1=.02; %Friction coefficient for base
Rw=.050; %Radius of round weight (m), only used for detailed expressions of Jy and Jz_bar
grav=9.81; %Acceleration of gravity (m/s^2)

if systype==1
	J1=.0166; %Inertia of lower disk including Pendulum base hardware
	Rh=0.25; %Distance from disk rotation axis to hub of pendulum (m)
	Ksys_prime=.00675; %"Ksys" or "Khw" without encoder gains
	Kenc1=2546; %Counts / radian
	Kenc4=2608; %Counts / radian
end
if systype==2
	m1=1.22 %m1=2.44 %Mass of base including carriage, equivalent motor/drive inertia, and Pendulum base hardware
	Ksys_prime=.0565; %"Ksys" or "Khw" without encoder gains
	%Ksys_prime=5.027; %"Ksys" or "Khw" without encoder gains
	Kenc1=2546*89 %Counts / m (includes encoder pulley gain)
	Kenc4=2608; %Counts / radian
	test=1
end
if systype==3
	J1=TBD %Inertia of drive disk including Pendulum base hardware
	Rh=0.25 %Distance from drive disk rotation axis to hub of pendulum (m)
end
if systype==4
	J1=TBD %Inertia of all components rotating about the vertical axis including Pendulum base hardware
	Rh=0.25 %Distance from yoke rotation axis to hub of pendulum (m)
end

% Derived Properties
m=mr+mw;
lcgr=yr-Lr/2;
lcg=(mr*lcgr+mw*yw)/m;

%For all but Model 210, Build plant model
if systype~=2
	%SOLVE FOR INERTIA PARAMETERS
	Rh=0.25; %Distance from disk rotation axis to hub of pendulum (m)
	J1_bar=J1+m*Rh^2;
	Jz_bar=mr*(Lr^2/12+lcgr^2)+mw*yw^2;
	Jy=0;
	%The following may be used if Rw is to be accounted for (generally this may be neglected)
	%Jz_bar=mr*(Lr^2/12+lcgr^2)+mw*(1/2*Rw^2+yw^2);
	%Jy=mw*Rw^2/4;
	
	%Transfer Function Expressions
	N1=Kenc1*Ksys_prime*[Jz_bar 0 -m*grav*lcg];
	N2=-[Kenc4*Ksys_prime*m*Rh*lcg 0];
	d_den=[(J1_bar+Jy)*Jz_bar-(m*Rh*lcg)^2 c1*Jz_bar -(J1_bar+Jy)*m*grav*lcg -c1*m*grav*lcg];
	N1=N1/d_den(1),N2=N2/d_den(1),d_den=d_den/d_den(1); %Makes d_den monic
	D1=[d_den 0]
	D2=d_den
	%State Space realization
	a_param=1/((J1_bar+Jy)*Jz_bar-(m*Rh*lcg)^2);
	A1=[0 1 0 0];
	A3=[0 0 0 1];
	A2=a_param*[0 -c1*Jz_bar -m^2*lcg^2*Rh*grav*Kenc1/Kenc4 0];
	A4=a_param*[0 c1*m*Rh*lcg*Kenc4/Kenc1 (J1_bar+Jy)*m*grav*lcg 0];
	A=[A1;A2;A3;A4]
	B=Ksys_prime*a_param*[0;Kenc1*Jz_bar;0;-Kenc4*m*Rh*lcg]
	
end
%If Model 210, Build plant model
if systype==2
	Jz_bar=mr*(Lr^2/12+lcgr^2)+mw*yw^2;
	Jy=0;
	%The following may be used if Rw is to be accounted for (generally this may be neglected)
	%Jz_bar=mr*(Lr^2/12+lcgr^2)+mw*(1/2*Rw^2+yw^2);
	%Jy=mw*Rw^2/4;
	
	%Transfer Function Expressions
	N1=Kenc1*Ksys_prime*[Jz_bar 0 -m*grav*lcg];
	%N2=-[Kenc4*Ksys_prime*m*lcg 0];
	N2=[Kenc4*Ksys_prime*m*lcg 0];
	d_den=[(m1+m)*Jz_bar-(m*lcg)^2 c1*Jz_bar -(m1+m)*m*grav*lcg -c1*m*grav*lcg];
	N1=N1/d_den(1),N2=N2/d_den(1),d_den=d_den/d_den(1); %Makes d_den monic
	D1=[d_den 0]
	D2=d_den
	%State Space realization
	a_param=1/((m1+m)*Jz_bar-(m*lcg)^2);
	A1=[0 1 0 0];
	A3=[0 0 0 1];
	%A2=a_param*[0 -c1*Jz_bar -m^2*lcg^2*grav*Kenc1/Kenc4 0];
	A2=a_param*[0 -c1*Jz_bar m^2*lcg^2*grav*Kenc1/Kenc4 0];
	%A4=a_param*[0 c1*m*lcg*Kenc4/Kenc1 (m1+m)*m*grav*lcg 0]
	A4=a_param*[0 -c1*m*lcg*Kenc4/Kenc1 (m1+m)*m*grav*lcg 0];
	A=[A1;A2;A3;A4]
	%B=Ksys_prime*a_param*[0;Kenc1*Jz_bar;0;-Kenc4*m*lcg]
	B=Ksys_prime*a_param*[0;Kenc1*Jz_bar;0;Kenc4*m*lcg] 
	
end

