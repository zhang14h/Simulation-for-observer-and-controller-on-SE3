clc
clear
close all
a0 =1;

for i = 1:200
  w1_t(i) = pi/90*(2*rand(1)-1);
  w2_t(i) = pi/90*(2*rand(1)-1);
  w3_t(i) = pi/90*(2*rand(1)-1);


  if w1_t(i) == 0 & w2_t(i) == 0 & w3_t(i) ==0;
     exp_t(:,:,i) = 1;
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].';
     w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
  else
   
   
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].'; %body-frame angular velocisy
     
     w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1,i) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_t(:,:,i) = expm(w_hat_t(:,:,i));
  end;
end;

fei_t= 0.6;
theta_t= -0.8 ;
miu_t= -0.5;
rpy = [fei_t;
       theta_t;
       miu_t];



R_t(:,:,1) = rotation(fei_t,theta_t,miu_t); %SO3 


px_t = 3;
py_t = 3;
pz_t = 3;
P_t = [px_t;
       py_t;
       pz_t];
         





fei_s = 0;
theta_s =0;
miu_s = 0;

R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3

R_s = rotation(fei_s,theta_s,miu_s);
px_s = 0;
py_s = 0;
pz_s = 0;
P_s = [px_s;
       py_s;
       pz_s];

for i = 1:200;
   
%     dF(i) = sqrt(6-trace(2*transpose(R_t(:,:,i))*R_s(:,:,i)));
%     oF(i) = 1/(1-dF(i)/sqrt(8));
%     eF(i) = (sqrt(transpose(P_t-P_s)*(P_t-P_s)));
%     dV(i) = oF(i) * eF(i);
    dV(i) = real(visiondistance(R_t(:,:,i),R_s(:,:,i),P_t,P_s));
    R_t(:,:,i+1) = exp_t(:,:,i)*R_t(:,:,i);
    error(:,:,i) = R_s(:,:,i).'*(R_t(:,:,i)+0.02*randi(1,3));
    proj= 0.5*1/2*((error(:,:,i) - error(:,:,i).'));
    A = R_t(:,:,i)*w_t(:,1,i)+0.04*rand(1);
    A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0]+proj;
    if A == 0
      exp_s = 1;
    else   
%        a_s = sqrt(A(2,1)^2+A(3,1)^2+A(3,2)^2);
%        A = A/a_s;
%        exp_s = [1 0 0;0 1 0;0 0 1] + A*sin(a_s)+ A^2*(1-cos(a_s));
     exp_s = expm(A);
    end;
    R_s(:,:,i+1) = exp_s*R_s(:,:,i);
    [roll_s1(i),roll_s2(i),roll_s3(i)] = derotation(R_s(:,:,i));
%     
    [roll_t1(i),roll_t2(i),roll_t3(i)] = derotation(R_t(:,:,i));
%     
end 

plot(roll_t1/pi*180,'--')
hold on;
plot(roll_s1/pi*180)
figure;
plot(roll_t2/pi*180,'--')
hold on;
plot(roll_s2/pi*180)
figure;
plot(roll_t3/pi*180,'--');
hold on;
plot(roll_s3/pi*360);
figure 
plot (dV);

