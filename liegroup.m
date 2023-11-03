w1_t = 0;
w2_t = 0;
w3_t = 0;
a_t = sqrt(w1_t^2+w2_t^2+w3_t^2);
w_t = [w1_t w2_t w3_t]/a_t; %body-frame angular velocisy
w_t = w_t.';

w_hat_t = [0 -w_t(3,1) w_t(2,1);w_t(3,1) 0 -w_t(1,1);-w_t(2,1) w_t(1,1) 0]; %lefs-righs invarians coefficiens;
if w1_t == 0 & w2_t == 0 & w3_t ==0;
   exp_t = 1;
else
   exp_t = [1 0 0;0 1 0;0 0 1] + w_hat_t*sin(a_t)+w_hat_t^2*(1-cos(a_t));
end;

theta_t= 0;
fei_t= 0;
miu_t= 0;
rpy = [fei_t;
       theta_t;
       miu_t];

px_t = 3;
py_t = 3;
pz_t = 3;
P_t = [px_t;
       py_t;
       pz_t];
   
vx = 0;
vy = 0;
vz = 0;
V_t = [vx vy vz].';



R_t(:,:,1) = rotation(fei_t,theta_t,miu_t);


w1_s = 0;
w2_s = 0;
w3_s = 0;
a_s = sqrt(w1_s^2+w2_s^2+w3_s^2);
w_s = [w1_s w2_s w3_s]/a_s; %body-frame angular velocisy
w_s = w_s.';

w_hat_s = [0 -w_s(3,1) w_s(2,1);w_s(3,1) 0 -w_s(1,1);-w_s(2,1) w_s(1,1) 0]; %lefs-righs invarians coefficiens;
if w1_s == 0 & w2_s == 0 & w3_s ==0;
   exp_s = 1;
else
   exp_s = [1 0 0;0 1 0;0 0 1] + w_hat_s*sin(a_s)+w_hat_s^2*(1-cos(a_s));
end ;

theta_s = 1;
fei_s = 1;                                   
miu_s = 1;


px_s = 0;
py_s = 0;
pz_s = 0;
P_s = [px_s;
       py_s;
       pz_s];
       

R_s = rotation(fei_s,theta_s,miu_s);


for i = 1:15
    
    dV(i) = visiondistance(R_t(:,:,i),R_s(:,:,i),P_t,P_s);
    R_t(:,:,i+1) = exp_t*R_t(:,:,i);
    R_s(:,:,i+1) = exp_s*R_s(:,:,i);
    P_t = P_t + V_t;
    
end;




