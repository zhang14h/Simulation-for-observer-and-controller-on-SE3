% clc
% clear
% close all
% a0 =1;
% iter = 300;
% 
% 
% for i = 1:iter
% %   w1_t(i) = 3*pi/90*(2*rand(1)-1);
% %   w2_t(i) = 3*pi/90*(2*rand(1)-1);
% %   w3_t(i) = 3*pi/90*(2*rand(1)-1);
%  
%   v_delta(:,1,i) = 1*[1*cos(1*i),1*cos(i/2),1*cos(2*i)] ;
%   v_o(:,1,i) = 1*[1*cos(i/3),1*sin(i/3),1*cos(i/3)] ;
%   delta1_t(i) = 0*pi/90;
%   delta2_t(i) = 0*pi/90;
%   delta3_t(i) = 0*pi/90;
%   mo1_t(i) = 0.5*pi/90*sin(i/2);
%   mo2_t(i) = 0.5*pi/90*sin(i/4);
%   mo3_t(i) = 0.5*pi/90*sin(3*i);
%   
% %   if w1_t(i) == 0 & w2_t(i) == 0 & w3_t(i) ==0;
% %      
% %      w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].';
% %      w_hat_t(:,:,i) = hatoperate ([w1_t(i) w2_t(i) w3_t(i)].'); %lefs-righs invarians coefficiens,so3;
% %      exp_t(:,:,i) = 1;
% %   else
% %    ugbyhnjm
% %    
% %      w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].'; %body-frame angular velocisy
% %      
% %      w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1,i) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
% %      exp_t(:,:,i) = exp(w_hat_t(:,:,i));
% %   end;
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
%   if mo1_t(i) == 0 & mo2_t(i) == 0 & mo3_t(i) ==0;
%      
%      mo_t(:,1,i) = [mo1_t(i) mo2_t(i) mo3_t(i)].';
%      mo_hat_t(:,:,i) = hatoperate ([mo1_t(i) mo2_t(i) mo3_t(i)].'); %lefs-righs invarians coefficiens,so3;
%      exp_mot(:,:,i) = [1 0 0; 0 1 0;0 0 1];
%   else
%    
%    
%      mo_t(:,1,i) = [mo1_t(i) mo2_t(i) mo3_t(i)].'; %body-frame angular velocisy
%      
%      mo_hat_t(:,:,i) = [0 -mo_t(3,1,i) mo_t(2,1,i);mo_t(3,1,i) 0 -mo_t(1,1,i);-mo_t(2,1,i) mo_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
%      exp_mot(:,:,i) = expm(mo_hat_t(:,:,i));
%   end;
%      SE3_mo(:,:,i) =expm([ mo_hat_t(:,:,i) v_o(:,1,i) ; [ 0 0 0] 0]);
% end;
% fei_t= 0;
% theta_t= 0 ;
% miu_t= 0;
% rpy = [fei_t;
%        theta_t;
%        miu_t];
% 
% R_t(:,:,1) = rotation(fei_t,theta_t,miu_t); %SO3 
% 
% 
% px_t = 0;
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
% 
% R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3
% 
% R_s = rotation(fei_s,theta_s,miu_s);
% px_s = 1;
% py_s = 2;
% pz_s = 1;
% P_s= [px_s;
%        py_s;
%        pz_s];
% 
% v_s(:,1) = [0,0,0].';
%    
% SE3_s(:,:,1) = [R_s(:,:,1) P_s;
%                [0,0,0]      1     ];
%            
%            
% R = rotation( pi/8, pi/7 ,-pi/5 );
% P = [10 10 7];
% T = [R P.'; [ 0 0 0] 1];
% 
% y(:,1) = [2;2;2;1]; y(:,2) = [1; 1; 1; 0]; y(:,3) = [1; 2;3;1]; y(:,4) = [3;2;1;0]; y(:,5) = [-1;3;1;1]; y(:,6) = [2;2;3;0]; 
% y(:,7) = [2; 1; 2; 1]; y(:,8) = [1; 0; 1; 0];
% 
% 
Ut(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
Ui(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
SE3_s(:,:,1) = [R_s(:,:,1) P_s;
               [0,0,0]      1     ];
SE3_t(:,:,1) = [R_t(:,:,1) P_t;
                [0,0,0]            1     ];

Ut(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
Ui(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
% Q = 0.00*[3 1 2 1; 1 1 1 1; 2 1 1 1; 1 1 1 1];
% Q=  Q*Q.'+0.5*eye(4);
Q = 0.01 * [1 2 3 1; 4 3 1 1; 2 1 2 1; 1 2 1 3];
Q=  Q+1*eye(4);

% figure
for i = 1:iter;
%     se3_t(:,:,i) =  [w_hat_t(:,:,i)+delta_hat_t(:,:,i) v_t(:,1,i)+v_delta(:,1,i);
%                   [0,0,0]        0   ];
   TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)]; 
   
    TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
                 SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
                 SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
    TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];
%     dv7(i) = real(visiondistance(R,TRs(:,:,i),P.',TPs(:,i)));
%     dv8(i) = real(visiondistance(R,TRt(:,:,i),P.',TPt(:,i)));
%     dv9(i) = real(visiondistance(TRs(:,:,i),TRt(:,:,i),TPs(:,i),TPt(:,i)));
    
    error4(:,:,i) =  SE3_s(:,:,i)*inv(T);
    error6(:,:,i) = SE3_t(:,:,i)*inv(T);
    error7(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
    omega1(:,:,i)=zeros(4,4);
    omega2(:,:,i)=zeros(4,4);
    omega3(:,:,i)=zeros(4,4);
%     omegak(:,:,i)=zeros(4,4);
%     omega_ki(:,:,i)=zeros(4,4);
    
    omega(:,:,i) =  1*projection(1/4*((error4(:,:,i) - error4(:,:,i).')));
%     omega4(:,:,i) =  0.9*inv(SE3_s(:,:,i))*projection(1/4*((error4(:,:,i) - eye(4))-(error4(:,:,i) - eye(4)).'))*SE3_s(:,:,i);
   
%     SE3_t(:,:,70) =[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 0];
%     v_s(:,i) = v_t(:,1,i)+0.1*rand(1)+ 0.3.*(TPt(:,i)+0.03*rand(1)-TPs(:,i));
%      Ui(:,:,i) =  -inv(SE3_s(:,:,i))*(omega(:,:,i))*SE3_s(:,:,i);
    Ui(:,:,i) =  -inv(T)*(omega(:,:,i))*(T);
    Ut(:,:,i) =  Ui(:,:,i)+[delta_hat_t(:,:,i) v_delta(:,1,i);[0,0,0]    0 ];;
%     U(:,:,i) =  -(omega(:,:,i))-[delta_hat_t(:,:,i) v_delta(:,1,i);[0,0,0]    0 ]-omega4(:,:,i);
    SE3_t(:,:,i+1) = SE3_t(:,:,i)*expm(Ut(:,:,i));
   error5(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
    for p = 1:lengthy
      
        omega2(:,:,i) = omega2(:,:,i)+0.01*(-error5(:,:,i)*y(:,p)-SE3_s(:,:,i)*deltay(:,p,i)+y(:,p))*transpose(y(:,p));
       omega3(:,:,i) = omega3(:,:,i)+0.006*(-error5(:,:,i)*y(:,p)-SE3_s(:,:,i)*deltay(:,p,i)+y(:,p))*transpose(y(:,p))/norm(y(:,p),2);
    end;
    omega_t(:,:,i) =   projection(1/4*((omega2(:,:,i) - omega2(:,:,i).')));
   se3_s = Ui(:,:,i)-inv(SE3_s(:,:,i))*omega_t(:,:,i)*SE3_s(:,:,i)-inv(SE3_s(:,:,i))*projection (1/4*(Q*omega3(:,:,i) - omega3(:,:,i).'*Q.'))*SE3_s(:,:,i); 
%     se3_s = Ui(:,:,i)-inv(SE3_s(:,:,i))*omega_t(:,:,i)*SE3_s(:,:,i);
  
    SE3_s(:,:,i+1) = SE3_s(:,:,i)*expm(se3_s);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     error9(:,:,i) = SE3_s(:,:,i+1)*inv(SE3_tm(:,:,i));
%     for p = 1:14
%        omegak(:,:,i) = omegak(:,:,i)+0.01 *(error9(:,:,i)*y(:,p)-y(:,p))*transpose(y(:,p));
%     end;
%     omega_ki(:,:,i)=  projection(1/4*((omegak(:,:,i) - omegak(:,:,i).')));
%     se3_s = Ui(:,:,i)-inv(SE3_s(:,:,i))*omega_ki(:,:,i)*SE3_s(:,:,i);
%     SE3_s(:,:,i+1) = SE3_s(:,:,i)*expm(se3_s);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [roll_s1(i),roll_s2(i),roll_s3(i)] = derotation(TRs(:,:,i));
    roll_s(:,i) =  [roll_s1(i),roll_s2(i),roll_s3(i)];
%     dv7(i) = real(visiondistance(R,TRs(:,:,i),P.',TPs(:,i)));
%     dv8(i) = real(visiondistance(R,TRt(:,:,i),P.',TPt(:,i)));
%     dv9(i) = real(visiondistance(TRs(:,:,i),TRt(:,:,i),TPs(:,i),TPt(:,i)));
    dv7(i) = norm(error4(:,:,i)-eye(4),'fro');
    dv8(i) = norm(error6(:,:,i)-eye(4),'fro');
    dv9(i) = norm(error7(:,:,i)-eye(4),'fro');
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(1,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(1,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(1,3,i)],'--r');
%     hold on 
%     xlim([-20 20]);
%     ylim([-20,20]);
%     zlim([-20,20]);
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(2,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(2,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(2,3,i)],'--g');
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(3,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(3,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(3,3,i)],'--b');
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(1,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(1,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(1,3,i)],'color',[10 10 10]/255);
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(2,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(2,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(2,3,i)],'color',[100 100 100]/255);
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(3,3,i)],'color',[200 200 200]/255);
%     hold off  
%     pause(0.01)

end;
k7 = 0;
k8 = 0;
k9 = 0;
for i=1:iter
k7 = k7+(dv7(i));
k8 = k8+(dv8(i));
k9 = k9+(dv9(i));
end;
for i= 1:iter

trdc3(i) = trace(projection(inv(error5(:,:,i))-inv(error5(:,:,i)).')*projection(inv(error5(:,:,i))-inv(error5(:,:,i)).').');
trdc4(i) = trace(inv(error5(:,:,i)));
trdc5(i)= trace(projection((inv(error5(:,:,i))*SE3_s(:,:,i)*inv(T))-(inv(error5(:,:,i))*SE3_s(:,:,i)*inv(T)).')*projection((SE3_s(:,:,i)*inv(T)-(SE3_s(:,:,i)*inv(T)).')).');

end
% display (k1);
% figure()
%  
%   plot (dv1);
%   title('observer and proposed controlled output');
%   grid
% figure()  
%   p1 = plot (dv2,'-r','linewidth',1); m1 = " controller with no damping";
%   hold on;
%   p2 = plot (dv5,'--g','linewidth',1); m2 = "controller with k = 1.5";
%   hold on;
%   p5 = plot (dv8,':b','linewidth',1); m5 = "controller with k = 1.5";
%  legend([p1,p2,p5],[m1,m2,m5])
% %   legend([p1,p2],[m1,m2])
%   title('attitude and position control performance');
%   xlabel('Time');
%   ylabel('Tracking error');
% %   txt = {'dV = ' ,k2};
% %   text(iter/2,10,txt);
%   grid
% 
% figure()  
%   p3 = plot (dv3,'-r','linewidth',1); m3 = "observer in (8)";
%   hold on;
%   p4 = plot (dv6,'--g','linewidth',1); m4 = "observer with k = 2";
%   hold on;
%   p6 = plot (dv9,':b','linewidth',1); m6 = " observer with k = 4";
%   legend([p3,p4,p6],[m3,m4,m6]);
%   title('Performance of the observer in (8) and the proposed observers');
%   xlabel('Time');
%   ylabel('Estimation error');
% %   txt = {' = ' ,k3};
% %   text(iter/2,3,txt);
%   grid
%   