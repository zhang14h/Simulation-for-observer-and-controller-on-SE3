clc
clear
close all
a0 =1;
iter = 300;


for i = 1:iter
    if rem(i,100) <=40;
       q(i) = 0.05*pi/90;
    else
       q(i) = 0;
    end;
%   w1_t(i) = 3*pi/90*(2*rand(1)-1);
%   w2_t(i) = 3*pi/90*(2*rand(1)-1);
%   w3_t(i) = 3*pi/90*(2*rand(1)-1);
%  
%   v_delta(:,1,i) = 1*[1*cos(1*i),1*cos(i/2),1*cos(2*i)] ;
%   
%   delta1_t(i) = 0.3*pi/90;
%   delta2_t(i) = 0.3*pi/90;
%   delta3_t(i) = 0.3*pi/90;
%  
%     
  v_delta(:,1,i) = 0.08*[1*cos(i/3),1*cos(i/4),1*cos(i/6)] ;
  
  delta1_t(i) = 0.06*cos(i/8)+10*q(i);
  delta2_t(i) = 0*pi/90;
  delta3_t(i) = 0*pi/90;
% % %   
% % %   
  
%   v_delta(:,1,i) = 1*[1*cos(1*i),1*cos(1.5*i),1*cos(1*i)] ;
%   
%   delta1_t(i) = 70*q(i);
%   delta2_t(i) = 0.01*pi/90;
%   delta3_t(i) = 0.03*pi/90;
% %   
% % % % %   
%   v_delta(:,1,i) = 0*[1*cos(1*i),1*cos(i/2),1*cos(2*i)] ;
% 
%   delta1_t(i) = 0*pi/90;
%   delta2_t(i) = 0*pi/90;
%   delta3_t(i) = 0*pi/90;
%   

  
%   if w1_t(i) == 0 & w2_t(i) == 0 & w3_t(i) ==0;
%      
%      w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].';
%      w_hat_t(:,:,i) = hatoperate ([w1_t(i) w2_t(i) w3_t(i)].'); %lefs-righs invarians coefficiens,so3;
%      exp_t(:,:,i) = 1;
%   else
%    
%    
%      w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].'; %body-frame angular velocisy
%      
%      w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1,i) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
%      exp_t(:,:,i) = exp(w_hat_t(:,:,i));
%   end;
%   if delta1_t(i) == 0 & delta2_t(i) == 0 & delta3_t(i) ==0;
%      
%      delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].';
%      delta_hat_t(:,:,i) = hatoperate ([delta1_t(i) delta2_t(i) delta3_t(i)].'); %lefs-righs invarians coefficiens,so3;
%      exp_deltat(:,:,i) = 1;
%   else
%    
%    
%      delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].'; %body-frame angular velocisy
%      
%      delta_hat_t(:,:,i) = [0 -delta_t(3,1,i) delta_t(2,1,i);delta_t(3,1,i) 0 -delta_t(1,1,i);-delta_t(2,1,i) delta_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
%      exp_deltat(:,:,i) = expm(delta_hat_t(:,:,i));
%   end;
    delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].';  
     delta_hat_t(:,:,i) = [0 -delta_t(3,1,i) delta_t(2,1,i);delta_t(3,1,i) 0 -delta_t(1,1,i);-delta_t(2,1,i) delta_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_deltat(:,:,i) = expm(delta_hat_t(:,:,i));
  
end;
% fei_t= 0.2;
% theta_t= 0;
% miu_t= 0.4;
fei_t= 1*0.99*pi
theta_t= 0*0.98*pi/3;
miu_t= 1*0.97*pi;
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


% SE3_t(:,:,1) = [0.1313  -0.0230    0.9911 2; 0.3377   0.9410   -0.0230 2; -0.9320    0.3377    0.1313 2; 0 0 0 1];  

fei_s = 1*0.99*pi;
theta_s = 0;
miu_s = 1*0.97*pi;

R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3

R_s = rotation(fei_s,theta_s,miu_s);
px_s = 0;
py_s = 0;
pz_s = 0;
P_s= [px_s;
       py_s;
       pz_s];

v_s(:,1) = [0,0,0].';
   
SE3_s(:,:,1) = [R_s(:,:,1) P_s;
               [0,0,0]      1     ];
         
%       

R = rotation( 0, 0 , 0 );
P = [5 4 1 ];
T = [R P.'; [ 0 0 0] 1];


y(:,1) = [2;2;2;1]; y(:,2) = [4; 1; 1; 1]; y(:,3) = [1; 2;3;1]; y(:,4) = [3;2;1;1]; y(:,5) = [-1;3;1;1]; y(:,6) = [2;2;3;1]; 
y(:,7) = [2; 1; 2; 1]; y(:,8) = [3;4;5;1] ; y(:,9) = [ 5;3;1;1] ; y(:,10) = [7;8;9;1]; y(:,11) = [2;5;7;1]; y(:,12) = [3;6;8;1];
y(:,13) = [1;5;9;1]; y(:,14) = [3;7;8;1];

% y(:,1) = [2;2;2;1]; y(:,2) = [4;4;4;1]; y(:,3) = [1; 2;3;1]; y(:,4) = [3;2;1;1]; y(:,5) = [-1;3;1;1]; y(:,6) = [2;2;3;1]; 
% y(:,7) = [2; 1; 2; 1]; y(:,8) = [3;4;5;1] ; y(:,9) = [ 5;3;1;1] ; y(:,10) = [7;8;9;1]; y(:,11) = [2;5;7;1]; y(:,12) = [3;6;8;1];
% y(:,13) = [1;5;9;1]; y(:,14) = [3;7;8;1];

% y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1]; y(:,5) = [1.7;3.4;5.1;1]; y(:,6) = [1.9;3.8;5.7;1]; 
% y(:,7) = [2.2;4.4;6.6;1]; y(:,8) = [3;5;4;1] ;  y(:,9) = [ 5;3;1;1] ; y(:,10) = [7;8;9;1]; y(:,11) = [2;5;7;1]; y(:,12) = [3;6;8;1];
% y(:,13) = [1;5;9;1]; y(:,14) = [3;7;8;1];

%  y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1]; y(:,5) = [1.7;3.4;5.1;1]; y(:,6) = [1.9;3.8;5.7;1]; y(:,7) = [2.2;4.4;6.6;1]; 
% y(:,8) = [3;5;4;1] ;  y(:,9) = [3.6;6;4.8;1] ; y(:,10) = [3.9;6.5;5.2;1]; y(:,11) = [7.2;13;10.4;1]; y(:,12) = [5.1;8.5;6.8;1];
% y(:,13) = [5.7;9.5;7.2;1]; y(:,14) = [3;7;8;1];
% 
% y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1]; y(:,5) = [1.7;3.4;5.1;1]; y(:,6) = [3.6;6;4.8;1];  y(:,7) = [3;5;4;1]; 

% y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.5;3.5;5;1]; y(:,4) = [3;7;8;1]; y(:,5) = [3;5;4;1];  y(:,6) = [3.6;6;4.8;1] ;y(:,7) = [7.2;13;10.4;1];

% y(:,1) = [1;2;3;1]; y(:,2) = [2;7;1;1]; y(:,3) = [1.5;3.5;5;1]; y(:,4) = [3;7;8;1]; y(:,5) = [3;5;4;1];  y(:,6) = [5;4;1;1] ;y(:,7) = [7.2;13;10.4;1];

% y(:,1) = [1;2;3;1]; y(:,2) = [2;7;1;1]; y(:,3) = [3;7;8;1]; y(:,4) = [3;5;4;1];
% y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1];
% y(:,1) = [1;2;3;1];  y(:,2) = [2;7;1;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1];

Ut(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
Ui(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
figure()
for i = 1:iter;
%     se3_t(:,:,i) =  [w_hat_t(:,:,i)+delta_hat_t(:,:,i) v_t(:,1,i)+v_delta(:,1,i);
%                   [0,0,0]        0   ];
    TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)]; 
   
    
    
    
    error6(:,:,i) = SE3_t(:,:,i)*inv(T);
    
    
    omega(:,:,i) =   1*projection(1/4*((error6(:,:,i) - error6(:,:,i).')));
    
%     SE3_t(:,:,70) =[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 0];
%     v_s(:,i) = v_t(:,1,i)+0.1*rand(1)+ 0.3.*(TPt(:,i)+0.03*rand(1)-TPs(:,i));
%     Ui(:,:,i) =  -inv(SE3_s(:,:,i))*(omega(:,:,i))*SE3_s(:,:,i);
%     Ui(:,:,i) =  -inv(T)*(omega(:,:,i))*(T);
    Ut(:,:,i) =  -inv(T)*(omega(:,:,i))*(T)
%      Ut(:,:,i) =  Ui(:,:,i)-[delta_hat_t(:,:,i) v_delta(:,1,i);[0,0,0]    0 ];
     
%     U(:,:,i) =  -(omega(:,:,i))-[delta_hat_t(:,:,i) v_delta(:,1,i);[0,0,0]    0 ]-omega4(:,:,i);
    SE3_t(:,:,i+1) = SE3_t(:,:,i)*expm(Ut(:,:,i));
    
    
    dv2(i) = real(visiondistance(R,TRt(:,:,i),P.',TPt(:,i)));
   
%     dv1(i) = norm(error4(:,:,i)-eye(4),'fro'); 
%       dv2(i) =norm(error6(:,:,i)-eye(4),'fro'); 
%       dv3(i) =norm(error7(:,:,i)-eye(4),'fro');
    plot3([TPt(1,i),TPt(1,i)+5*TRt(1,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(1,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(1,3,i)],'color',[10 10 10]/255);
    hold on 
    xlim([-20 20]);
    ylim([-20,20]);
    zlim([-20,20]);
    
    plot3([TPt(1,i),TPt(1,i)+5*TRt(2,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(2,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(2,3,i)],'color',[100 100 100]/255);
    plot3([TPt(1,i),TPt(1,i)+5*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(3,3,i)],'color',[200 200 200]/255);
    plot3(y(1,:),y(2,:),y(3,:),'o','Color','k','MarkerSize',5)
    hold off  
    pause(0.01)

end;
k1 = 0;
k2 = 0;
k3 = 0;
for i=1:iter

k2 = k2+(dv2(i));

end;
% for i= 1:iter
% 
% trn3(i) = trace(projection(inv(error5(:,:,i))-inv(error5(:,:,i)).')*projection(inv(error5(:,:,i))-inv(error5(:,:,i)).').');
% trn4(i) = trace(inv(error5(:,:,i)));
% trn5(i)= trace(projection((inv(error5(:,:,i))*SE3_s(:,:,i)*inv(T))-(inv(error5(:,:,i))*SE3_s(:,:,i)*inv(T)).')*projection((SE3_s(:,:,i)*inv(T)-(SE3_s(:,:,i)*inv(T)).')).');
% 
% end
% display (k1);
% figure()
%  
%   plot (dv1);
%   title('observer and proposed controlled output');
%   grid
figure(1)  
  p1 = plot (dv2,'-r','linewidth',1); m1 = "Ec with Kc = 0";
%   title('attitude and position control performance with Kc = 0');
%  
%   txt = {'dV = ' ,k2};
%   text(iter/2,10,txt);
  grid
  hold on
  p2 = plot (dv3,'-.g','linewidth',1); m2 = "Eo with Ko = 0";
  title('Closed-loop Performance ');
  xlabel('Time');
  ylabel('Error');
  legend([p1,p2],[m1,m2]);
% 
% figure(2)  
%   p2 = plot (dv3,'-.g','linewidth',1);
%   title('observer in [24] with Ko=0');
%   xlabel('Time');
%   ylabel('Estimation error');
% %   txt = {' = ' ,k3};
% %   text(iter/2,3,txt);
%   grid
%   hold on