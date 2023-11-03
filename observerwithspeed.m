clc
clear
close all
a0 =1;
iter =60;
for i = 1:iter
  w1_t(i) = 3*pi/90*(2*rand(1)-1);
  w2_t(i) = 3*pi/90*(2*rand(1)-1);
  w3_t(i) = 3*pi/90*(2*rand(1)-1);
  v_t(:,1,i) = 1*[2*rand(1)-1,2*rand(1)-1,2*rand(1)-1].';
  v_delta(:,1,i) = 0*[1*cos(1*i),1*cos(i/2),1*cos(2*i)] ;
  v_o(:,1,i) = 0*[1*cos(i/3),1*sin(i/3),1*cos(i/3)] ;
  delta1_t(i) = 0*pi/90;
  delta2_t(i) = 0*pi/90;
  delta3_t(i) = 0*pi/90;
  mo1_t(i) = 0*pi/90*sin(i/2);
  mo2_t(i) = 0*pi/90*sin(i/4);
  mo3_t(i) = 0*pi/90*sin(3*i);
  
  if w1_t(i) == 0 & w2_t(i) == 0 & w3_t(i) ==0;
     
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].';
     w_hat_t(:,:,i) = hatoperate ([w1_t(i) w2_t(i) w3_t(i)].'); %lefs-righs invarians coefficiens,so3;
     exp_t(:,:,i) = 1;
  else
   
   
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].'; %body-frame angular velocisy
     
     w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1,i) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_t(:,:,i) = exp(w_hat_t(:,:,i));
  end;
  if delta1_t(i) == 0 & delta2_t(i) == 0 & delta3_t(i) ==0;
     
     delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].';
     delta_hat_t(:,:,i) = hatoperate ([delta1_t(i) delta2_t(i) delta3_t(i)].'); %lefs-righs invarians coefficiens,so3;
     exp_deltat(:,:,i) = 1;
  else
   
   
     delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].'; %body-frame angular velocisy
     
     delta_hat_t(:,:,i) = [0 -delta_t(3,1,i) delta_t(2,1,i);delta_t(3,1,i) 0 -delta_t(1,1,i);-delta_t(2,1,i) delta_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_deltat(:,:,i) = expm(delta_hat_t(:,:,i));
  end;
  if mo1_t(i) == 0 & mo2_t(i) == 0 & mo3_t(i) ==0;
     
     mo_t(:,1,i) = [mo1_t(i) mo2_t(i) mo3_t(i)].';
     mo_hat_t(:,:,i) = hatoperate ([mo1_t(i) mo2_t(i) mo3_t(i)].'); %lefs-righs invarians coefficiens,so3;
     exp_mot(:,:,i) = [1 0 0; 0 1 0;0 0 1];
  else
   
   
     mo_t(:,1,i) = [mo1_t(i) mo2_t(i) mo3_t(i)].'; %body-frame angular velocisy
     
     mo_hat_t(:,:,i) = [0 -mo_t(3,1,i) mo_t(2,1,i);mo_t(3,1,i) 0 -mo_t(1,1,i);-mo_t(2,1,i) mo_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_mot(:,:,i) = expm(mo_hat_t(:,:,i));
  end;
     SE3_mo(:,:,i) =expm([ mo_hat_t(:,:,i) v_o(:,1,i) ; [ 0 0 0] 0]);
end;

fei_t= 0;
theta_t= 0  ;
miu_t= 0;
rpy = [fei_t;
       theta_t;
       miu_t];



R_t(:,:,1) = rotation(fei_t,theta_t,miu_t); %SO3 


px_t = 0;
py_t = 0;
pz_t = 0;
P_t= [px_t;
      py_t;
      pz_t];
   

         
SE3_t(:,:,1) = [R_t(:,:,1) P_t;
               [0,0,0]            1     ];


fei_s = 1.5;
theta_s =1.5;
miu_s = 1.5;

R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3

R_s = rotation(fei_s,theta_s,miu_s);
px_s = 5;
py_s = 5;
pz_s = 5;
P_s= [px_s;
       py_s;
       pz_s];

v_s(:,1) = [0,0,0].';
   
SE3_s(:,:,1) = [R_s(:,:,1) P_s;
               [0,0,0]      1     ];


 
xlim([-50 50]);
ylim([-50,50]);
zlim([-50,50]);
hold on
for i = 1:iter;
   
   
%     dV(i) = real(visiondistance(R_t(:,:,i),R_s(:,:,i),P_t,P_s));
%     R_t(:,:,i+1) = exp_t(:,:,i)*R_t(:,:,i);
%     
%     error(:,:,i) = R_s(:,:,i).'*(R_t(:,:,i)+0.02*randi(1,3));
%     
%     proj= 0.5*1/2*((error(:,:,i) - error(:,:,i).'));
%     A = w_t(:,1,i)+0.04*rand(1);
%     A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0]+proj;
%     if A == 0
%       exp_s = 1;
%     else   
% %        a_s = sqrt(A(2,1)^2+A(3,1)^2+A(3,2)^2);
% %        A = A/a_s;
% %        exp_s = [1 0 0;0 1 0;0 0 1] + A*sin(a_s)+ A^2*(1-cos(a_s));
%      exp_s = expm(A);
%     end;
%     
%     R_s(:,:,i+1) = exp_s*R_s(:,:,i);
%     
    

% end of so3 estimation;    
    
    se3_t(:,:,i) =  [w_hat_t(:,:,i)+delta_hat_t(:,:,i) v_t(:,1,i)+v_delta(:,1,i);
                  [0,0,0]        0   ];
    
    SE3_t(:,:,i+1) = expm(se3_t(:,:,i))*SE3_t(:,:,i);
    TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)]; 
    SE3_tm(:,:,i) = SE3_mo(:,:,i)* SE3_t(:,:,i);
    
    TRtm(:,:,i) = [SE3_tm(1,1,i) SE3_tm(1,2,i) SE3_tm(1,3,i);
                 SE3_tm(2,1,i) SE3_tm(2,2,i) SE3_tm(2,3,i);
                 SE3_tm(3,1,i) SE3_tm(3,2,i) SE3_tm(3,3,i)];
             
    TPtm(:,i) = [SE3_tm(1,4,i) SE3_tm(2,4,i) SE3_tm(3,4,i)];
    
    TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
                 SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
                 SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
             
    TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];  
    dv1(i) = real(visiondistance(TRtm(:,:,i),TRt(:,:,i),TPtm(:,i),TPt(:,i)));
    dv2(i) = real(visiondistance(TRtm(:,:,i),TRs(:,:,i),TPtm(:,i),TPs(:,i)));
    dv3(i) = real(visiondistance(TRt(:,:,i),TRs(:,:,i),TPt(:,i),TPs(:,i)));
    error3(:,:,i) =  TRs(:,:,i).'*(TRtm(:,:,i));
    
    A = w_t(:,1,i);
    A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0];
    omega(:,:,i) = A+ 0.5*1/2*((error3(:,:,i) - error3(:,:,i).'));
    
%     A = w_t(:,1,i)+0.04*rand(1);
%     A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0];
    
    
     v_s(:,i) = v_t(:,i)-(TPs(:,i)-TRs(:,:,i)*TRtm(:,:,i).'*TPtm(:,i));
    se3_s(:,:,i) =  [omega(:,:,i) v_s(:,i);
                     [0,0,0]        0   ] ; 
    SE3_s(:,:,i+1) = expm(se3_s(:,:,i))*SE3_s(:,:,i);
    [roll_s1(i),roll_s2(i),roll_s3(i)] = derotation(TRs(:,:,i));
    
    
    roll_s(:,i) =  [roll_s1(i),roll_s2(i),roll_s3(i)];
%     
    [roll_t1(i),roll_t2(i),roll_t3(i)] = derotation(TRt(:,:,i));
    roll_t(:,i) = [roll_t1(i),roll_t2(i),roll_t3(i)];
    
    
end
k1 = 0;
k2 = 0;
k3 = 0;
for i=1:iter
k1 = k1+(dv1(i));  
k2 = k2+(dv2(i));
k3 = k3+(dv3(i));
end;;


% subplot(3,1,1)
% plot (dv1)
% txt = {'dV = ' ,k1};
% text(iter/2,10,txt);
% subplot(3,1,2)
% plot (dv2)
% txt = {'dV = ' ,k2};
% text(iter/2,10,txt);
  subplot(2,1,1)
  plot (dv3);
  title('Observer proposed in equation (8)');
  grid
  xlabel('Time')
  ylabel('Estimation error')
%   txt = {'dV = ' ,k3};
%   text(iter/2,10,txt);
%   subplot(2,1,2)
  
  
%   xlabel('Time')
%   ylabel('Estimation error')
%   grid
% %   txt = {'dV = ' ,k2};
% %   text(iter/2,10,txt);
%   title('Our proposed observer k_0=1.5');
% %   subplot(2,1,1)
% %   plot (dv3);
% %   xlabel('Time')
% %   ylabel('Estimation error')
% %   grid
% %   title('Our proposed observer k_0=1');
%   subplot(3,1,3)
%   plot (dv3);
%   txt = {'dV = ' ,k3};
% text(iter/2,10,txt);
  