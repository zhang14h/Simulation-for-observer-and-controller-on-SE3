% gain for each angle
clear; 
clc;

%signal decay value for each angle in dB (vertical and horizontal)
Gtx1H = [-19; -12.3; -9; -10; -10.2; -9.5; -10; -9; -6.3; -5.5; -4; -3.2; 
        -2.5; -1.1; -0.5; -1; -0.3; 0; -0.4; -1; -2; -2.5; -2.7; -3.2;
        -4.1; -4.3; -5.1;-6.9; -8.5; -10.7; -11.9; -12.5; -16.9].';
Gtx1V = [-19; -17; -9.8; -7; -7.5; -6.8; -5.3; -5.8; -4.2; -3.8; -3; -2.2;
         -1.6; -0.7; -1; -0.6; -0.2; -0.4; -0.5; -1.1 ; -1.3; -1.2; -1.8; 
         -3; -4.1; -4.3; -8.3; -8.1; -8.2; -8.1; -8.7; -10.8; -13.8].';    

%sampled angle     
Atx1H = [-80:5:80]*pi/180;
Atx1V = [-80:5:80]*pi/180;

%signal decay gain
Rtx1H = 10.^(Gtx1H/20);
Rtx1V = 10.^(Gtx1V/20);
Rtx1 = Rtx1V.'*Rtx1H*9;

% for i = 1:33
%     boundH(1,i) = Rtx1H(i)*cos(Atx1H(i));
%     boundH(2,i) = Rtx1H(i)*sin(Atx1H(i));
%     boundV(1,i) = Rtx1V(i)*cos(Atx1V(i));
%     boundV(2,i) = Rtx1V(i)*sin(Atx1V(i));
% end;

%counter
temp = 1;

% find the furthest distance for the Radar at each angle (combine H and V)
for i = 1:33
    for j = 1:33
       bound(1,temp) = Rtx1(j,i)*cos(Atx1H(i))*cos(Atx1V(j));
       bound(2,temp) = Rtx1(j,i)*cos(Atx1H(i))*sin(Atx1V(j));
       bound(3,temp) = Rtx1(j,i)*sin(Atx1H(i));
       bound(4,temp) = Atx1H(i)/pi*180;
       bound(5,temp) = Atx1V(j)/pi*180;
       temp = temp + 1;
    end;
end;

%plot the shape

% plot(boundH(1,:), boundH(2,:));
% figure()
% plot(boundV(1,:), boundV(2,:));


shp = alphaShape(bound(1,:).',bound(2,:).',bound(3,:).',2);
plot(shp,'FaceColor','red')
% fill3(bound(1,:),bound(2,:),bound(3,:),'g');
% K = convhull(bound(1,:),bound(2,:),bound(3,:));
% plot3(bound(1,K),bound(2,K),bound(3,K));