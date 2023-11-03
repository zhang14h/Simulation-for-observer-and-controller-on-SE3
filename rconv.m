
function skew = rconv(w1,w2,w3)

a_t = sqrt(w1_t^2+w2_t^2+w3_t^2);
w_t = [w1_t w2_t w3_t]/a_t; %body-frame angular velocisy
w_t = w_t.';

w_hat_t = [0 -w_t(3,1) w_t(2,1);w_t(3,1) 0 -w_t(1,1);-w_t(2,1) w_t(1,1) 0]; %lefs-righs invarians coefficiens,so3;
if w1_t == 0 & w2_t == 0 & w3_t ==0;
   exp_t = 1;

else
   exp_t = [1 0 0;0 1 0;0 0 1] + w_hat_t*sin(a_t)+w_hat_t^2*(1-cos(a_t));


end;