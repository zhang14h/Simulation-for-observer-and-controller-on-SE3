%InvPendPolePlace.m
%  The following generates a family of pole placement controllers and plots the
%  resulting closed loop poles and step and singular value responses.
%  You must first have generated the A,& B,matrices.


Cstr=[1 0 0 0]


Pcl1=[-1-1*i;-1+1*i;-2;-1.5]*2*pi
Pcl2=1.25*Pcl1
Pcl3=[-1.2;-1.4;-1.6;-1.8]*2*pi
Pcl4=1.25*Pcl3

Kpp1=place(A,B,Pcl1)
Kpp2=place(A,B,Pcl2)
Kpp3=place(A,B,Pcl3)
Kpp4=place(A,B,Pcl4)

Pcl1=eig(A-B*Kpp1)/2/pi
Pcl2=eig(A-B*Kpp2)/2/pi
Pcl3=eig(A-B*Kpp3)/2/pi
Pcl4=eig(A-B*Kpp4)/2/pi

% Plot Roots (should be same as prescribed roots)
plot(Pcl1,'x')
hold;
plot(Pcl2,'o');
plot(Pcl3,'*');
plot(Pcl4,'+');
title('Pole Place Controller, Closed Loop Poles For Various Control Effort Weights')
xlabel('Real  (Hz)')
ylabel('Imaginary  (Hz)')

grid;
hold;
pause
axis normal
% User may define other labels & axis limits, frequencies in Hz.
% Stike any key to continue

%Generate Step Responses at theta1 & theat2 to commanded step

CC=[[1 0 0 0];[0 0 1 0]];
steppp1=step([A-B*Kpp1],B*Kpp1,CC,[0 0],1,t);
steppp2=step([A-B*Kpp2],B*Kpp2,CC,[0 0],1,t);
steppp3=step([A-B*Kpp3],B*Kpp3,CC,[0 0],1,t);
steppp4=step([A-B*Kpp4],B*Kpp4,CC,[0 0],1,t);

%Example Step Response Plot
plot(t,[steppp1 steppp2 steppp3 steppp4]),grid,axis([0,5,-.4,1.2]),
xlabel('Time (s)'),ylabel('Amplitude')
pause

%Singular value plots
syslqrcl1=ss([A-B*Kpp1],B*Klqr1(:,1),CC,0);
syslqrcl2=ss([A-B*Kpp2],B*Klqr2(:,1),CC,0);
syslqrcl3=ss([A-B*Kpp3],B*Klqr3(:,1),CC,0);
%Example SV response
bode(syslqrcl3)
