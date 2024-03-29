clc
clear
close all

iter = 3000;
alpha = 0.004;
snr = 4;
d = 0.04;
landmarks = 150;

for i = 1:iter
 
  w1_t(i) =0;
  w2_t(i) = pi/4000;
  w3_t(i) = pi/5000;
  v_t(:,1,i) = 1*[-2*2*pi/iter*cos(2*i*pi/iter) 2*2*pi/iter*sin(2*i*pi/iter)  1*2*pi/iter*sin(2*i*pi/iter)].';
  if w1_t(i) == 0 & w2_t(i) == 0 & w3_t(i) ==0;
     
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].';
     w_hat_t(:,:,i) = hatoperate ([w1_t(i) w2_t(i) w3_t(i)].'); %lefs-righs invarians coefficiens,so3;
     
  else
   
   
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].'; %body-frame angular velocisy
     
     w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1,i) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
    
  end;
  
  if i<= 9*iter/10;  
  v_delta(:,1,i) = d*0.5*[1*cos(i/200); 1*sin(i/300+pi/6); 0.5*cos(i/100)*sin(i/700+pi/7)]  ;

%   v_delta(:,1,i) = d*0.5;

  
  delta1_t(i) = d*pi/2500;
  delta2_t(i) = 0*pi/2000;
  delta3_t(i) = d*pi/1500;
  
  
  
    if delta1_t(i) == 0 & delta2_t(i) == 0 & delta3_t(i) ==0;
     
     delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].';
     delta_hat_t(:,:,i) = hatoperate ([delta1_t(i) delta2_t(i) delta3_t(i)].'); %lefs-righs invarians coefficiens,so3;
     exp_deltat(:,:,i) = 1;
    else
     delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].'; %body-frame angular velocisy
     
     delta_hat_t(:,:,i) = [0 -delta_t(3,1,i) delta_t(2,1,i);delta_t(3,1,i) 0 -delta_t(1,1,i);-delta_t(2,1,i) delta_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_deltat(:,:,i) = expm(delta_hat_t(:,:,i));
    end;
  else
    delta_hat_t(:,:,i)= zeros(3,3);
    v_delta(:,1,i) =[0;0;0];
  end;    
end;



hov=3*pi/4;  %horizontal fieled of view in radians
vov=2*pi/3;  %vertical field of view in radians
focall = 1; 
focalh = 8;  %parameters for camera
Rad = 10;


for i = 1: landmarks
    y(:,i) = [randi([-20 20]); randi([-20 20]); randi([-20 20]); 1 ];
    
end;
j = -20;
k = -20;
h = -20;

% for i = 1: landmarks
%     
%        if h > 20
%           h = -20;
%           k = k +40/(nthroot(landmarks,3)-1);
%        end;   
%        if k >  20;
%           k = -20;
%           j = j + 40/(nthroot(landmarks,3)-1);
%        end;
%     y(:,i) = [j; k; h; 1 ];
%     h = h +40/(nthroot(landmarks,3)-1);
% end;
%sensor infomration
% P_t = [10;10;10];
 P_t = [4;3;2];

fei_t = 0.4;
theta_t = 0;   %sensor angle in radians
miu_t = 0.2;

R_t(:,:,1) = rotation(fei_t,theta_t,miu_t);
SE_3t(:,:,1) = [R_t(:,:,1) P_t; [0 0 0] 1];

%observer information
fei_s = 0;
theta_s = 0;
miu_s = 0;

R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3

R_s = rotation(fei_s,theta_s,miu_s);

P_s = [0;0;0];

   
SE_3s(:,:,1) = [R_s(:,:,1) P_s;
               [0,0,0]      1     ];