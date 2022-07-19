% clc
% % clear
% % close all
% a0 =1;
% iter =500;
% % for i = 1:iter
% %   w1_t(i) = pi/90*(2*rand(1)-1);
% %   w2_t(i) = pi/90*(2*rand(1)-1);
% %   w3_t(i) = pi/90*(2*rand(1)-1);
% %   v_t(:,1,i) = [2*rand(1)-1,2*rand(1)-1,2*rand(1)-1].';
% % 
% %   delta1_t(i) = 0.5*pi/90*(2*rand(1)-1)*sin(i);
% %   delta2_t(i) = 0.5*pi/90*(2*rand(1)-1)*sin(i);
% %   delta3_t(i) = 0.5*pi/90*(2*rand(1)-1)*sin(i);
% %   if w1_t(i) == 0 & w2_t(i) == 0 & w3_t(i) ==0;
% %      
% %      w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].';
% %      w_hat_t(:,:,i) = hatoperate ([w1_t(i) w2_t(i) w3_t(i)].'); %lefs-righs invarians coefficiens,so3;
% %      exp_t(:,:,i) = 1;
% %   else
% %    
% %    
% %      w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].'; %body-frame angular velocisy
% %      
% %      w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1,i) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
% %      exp_t(:,:,i) = exp(w_hat_t(:,:,i));
% %   end;
% %   if delta1_t(i) == 0 & delta2_t(i) == 0 & delta3_t(i) ==0;
% %      
% %      delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].';
% %      delta_hat_t(:,:,i) = hatoperate ([delta1_t(i) delta2_t(i) delta3_t(i)].'); %lefs-righs invarians coefficiens,so3;
% %      exp_deltat(:,:,i) = 1;
% %   else
% %    
% %    
% %      delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].'; %body-frame angular velocisy
% %      
% %      delta_hat_t(:,:,i) = [0 -delta_t(3,1,i) delta_t(2,1,i);delta_t(3,1,i) 0 -delta_t(1,1,i);-delta_t(2,1,i) delta_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
% %      exp_deltat(:,:,i) = exp(delta_hat_t(:,:,i));
% %   end;
% % end;
% 
% fei_t= 0.5;
% theta_t= -0.2 ;
% miu_t= -0.3;
% rpy = [fei_t;
%        theta_t;
%        miu_t];
% 
% 
% 
% R_t(:,:,1) = rotation(fei_t,theta_t,miu_t); %SO3 
% 
% 
% px_t =0;
% py_t = 0;
% pz_t = 0;
% P_t= [px_t;
%       py_t;
%       pz_t];
%    
% 
%          
% SE3_t(:,:,1) = [R_t(:,:,1) P_t;
%                [0,0,0]            1     ];
% 
% 
% fei_s = 0;
% theta_s =0;
% miu_s = 0;

% R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3
% 
% R_s = rotation(fei_s,theta_s,miu_s);
% px_s = 2;
% py_s = 2;
% pz_s = 2;
% P_s= [px_s;
%        py_s;
%        pz_s];
% 
% v_s(:,1) = [0,0,0].';
%    
% SE3_s(:,:,1) = [R_s(:,:,1) P_s;
%                [0,0,0]      1     ];
% 
% y(:,1) = [2;2;2; 1]; y(:,2) = [1; 1; 1; 0]; y(:,3) = [1; 2;3;1]; y(:,4) = [3;2;1;0]; y(:,5) = [-1;3;1;1]; y(:,6) = [2;2;3;0];
% y(:,7) = [2; 1; 2; 1]; y(:,8) = [1; 0; 1; 0];
% Q = 0.01*[3 1 2 1; 1 1 1 1; 2 1 1 1; 1 1 1 1];
% Q=  Q*Q.'+eye(4);
Q = 0.001 * [1 2 3 1; 4 3 1 1; 0 1 0 1; 1 2 7 3];
Q=  Q+eye(4);
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
    
    SE3_t(:,:,i+1) = SE3_t(:,:,i)*expm(se3_t(:,:,i));
    TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)];
    
    TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
                 SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
                 SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
             
    TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];  
%   
    
   
    error5(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
    dv2(i) = norm(error5(:,:,i)-eye(4),'fro');
    omega1(:,:,i)=zeros(4,4);
    omega2(:,:,i)=zeros(4,4);
    for p = 1:14
       omega1(:,:,i) = omega1(:,:,i)+0.01*(-error5(:,:,i)*y(:,p)-SE3_s(:,:,i)*deltay(:,p,i)+y(:,p))*transpose(y(:,p));
      omega2(:,:,i) = omega2(:,:,i) +0.006*(-error5(:,:,i)*y(:,p)-SE3_s(:,:,i)*deltay(:,p,i)+y(:,p))*transpose(y(:,p))/norm(y(:,p),2);
    end;
    
%     A = w_t(:,1,i)+0.04*rand(1);
%     A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0];
    omega(:,:,i) =   projection(1/4*((omega1(:,:,i) - omega1(:,:,i).')));
    
    
%     v_s(:,i) = v_t(:,1,i)+0.1*rand(1)+ 0.3.*(TPt(:,i)+0.03*rand(1)-TPs(:,i));
%     se3_s(:,:,i) =  [w_hat_t(:,:,i) v_t(:,1,i);
%                   [0,0,0]        0   ]-inv(SE3_s(:,:,i))*omega(:,:,i)*SE3_s(:,:,i)-inv(SE3_s(:,:,i))*projection (0.5*[1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1]*1/4*((omega2(:,:,i) - omega2(:,:,i).')))*SE3_s(:,:,i); 
    se3_s(:,:,i) =  [w_hat_t(:,:,i) v_t(:,1,i);
                  [0,0,0]        0   ]-inv(SE3_s(:,:,i))*omega(:,:,i)*SE3_s(:,:,i)-inv(SE3_s(:,:,i))*projection (1/4*(Q*omega2(:,:,i) - omega2(:,:,i).'*Q.'))*SE3_s(:,:,i); 
    SE3_s(:,:,i+1) = SE3_s(:,:,i)*expm(se3_s(:,:,i));
    [roll_s1(i),roll_s2(i),roll_s3(i)] = derotation(TRs(:,:,i));
    roll_s(:,i) =  [roll_s1(i),roll_s2(i),roll_s3(i)];
%     
    [roll_t1(i),roll_t2(i),roll_t3(i)] = derotation(TRt(:,:,i));
    roll_t(:,i) = [roll_t1(i),roll_t2(i),roll_t3(i)];
    
%     plot3([TPs(1,i),TPs(1,i)+10*TRs(1,1,i)],[TPs(2,i),TPs(2,i)+10*TRs(1,2,i)],[TPs(3,i),TPs(3,i)+10*TRs(1,3,i)],'r');
%     hold on 
%     xlim([-20 20]);
%     ylim([-20,20]);
%     zlim([-20,20]);
%     plot3([TPs(1,i),TPs(1,i)+10*TRs(2,1,i)],[TPs(2,i),TPs(2,i)+10*TRs(2,2,i)],[TPs(3,i),TPs(3,i)+10*TRs(2,3,i)],'g');
%     plot3([TPs(1,i),TPs(1,i)+10*TRs(3,1,i)],[TPs(2,i),TPs(2,i)+10*TRs(3,2,i)],[TPs(3,i),TPs(3,i)+10*TRs(3,3,i)],'b');
%     plot3([TPt(1,i),TPt(1,i)+10*TRt(1,1,i)],[TPt(2,i),TPt(2,i)+10*TRt(1,2,i)],[TPt(3,i),TPt(3,i)+10*TRt(1,3,i)],'color',[10 10 10]/255);
%     plot3([TPt(1,i),TPt(1,i)+10*TRt(2,1,i)],[TPt(2,i),TPt(2,i)+10*TRt(2,2,i)],[TPt(3,i),TPt(3,i)+10*TRt(2,3,i)],'color',[100 100 100]/255);
%     plot3([TPt(1,i),TPt(1,i)+10*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+10*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+10*TRt(3,3,i)],'color',[200 200 200]/255);
%     plot3(y(1,:),y(2,:),y(3,:),'o','Color','k','MarkerSize',5)
%     hold off  
%     pause(0.01)
%     
end 
k2 = 0;
for i=1:iter
k2 = k2+(dv2(i));
end;

figure()
 
  plot (dv2);
  title('observer and proposed controlled output');
  grid
  

% figure()
% subplot(3,1,1);
%   plot(roll_t1/pi*180,'--')
%   hold on;
%   plot(roll_s1/pi*180)
%   ylim('auto')
%   subplot(3,1,2);
%   plot(roll_t2/pi*180,'--')
%   hold on;
%   plot(roll_s2/pi*180)
%   ylim('auto')
%   subplot(3,1,3)
%   plot(roll_t3/pi*180,'--');
%   hold on;
%   plot(roll_s3/pi*180);
%   ylim('auto')
% figure()
%   subplot(3,1,1);
% plot(TPs(1,:),'--');
% hold on 
% plot(TPt(1,:));
% ylim('auto')
% subplot(3,1,2);
% plot(TPs(2,:),'--');
% hold on 
% plot(TPt(2,:));
% ylim('auto')
% subplot(3,1,3);
% plot(TPs(3,:),'--');
% hold on 
% plot(TPt(3,:));
% ylim('auto')
