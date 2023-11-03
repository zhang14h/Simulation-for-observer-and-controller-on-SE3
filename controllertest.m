clc
clear
close all
a0 =1;
iter = 500;

% fei_t= 0.4;
% theta_t= 0.3;
% miu_t= 0.5;
fei_t= 0;
theta_t= 1*pi;
miu_t= 0.3*pi;
rpy = [fei_t;
       theta_t;
       miu_t];

R_t(:,:,1) = rotation(fei_t,theta_t,miu_t); %SO3 


px_t = 3;
py_t = 4;
pz_t = 1;
P_t= [px_t;
      py_t;
      pz_t];
   

         
SE3_t(:,:,1) = [R_t(:,:,1) P_t;
                [0,0,0]            1     ];

R = rotation( 0 , 0 , 0 );
P = [0 0 0];
T = [R P.'; [ 0 0 0] 1];

% Q = 0.01 * [1 2 0 0; 0 3 1 0; 0 0 5 4; 0 0 0 10];
Q = [1 1 1 1; 2 3 1 1; 3 1 2 1; 1 1 1 3];
Q=  0.7*eye(4)+ 0*[0 1 0 0 ; 0 0 1 0; 0 0 0 0; 0 0 0 0];



for i = 1:iter
%%   
%%%%%%%%%%%%%   Calculating Avoidance Algorithm
%     U_avoid(:,:,i) = zeros(4,4);
%     for j = 1:lengthavoid
%         p_obstacle(:,j) = inv(SE3_t(:,:,i))*y_obstacle(:,j); 
%         omega_avoid(:,:,j) =  p_obstacle(:,j)*(p_obstacle(:,j)-[0;0;0;1]).';
%         a_avoid = 1*exp(1 - 1*norm(p_obstacle(:,j)-[0;0;0;1],'fro'));
% %         a_avoid = 0;
%         U_avoid(:,:,i) = U_avoid(:,:,i) + a_avoid*projection(1/4*((omega_avoid(:,:,j) - omega_avoid(:,:,j).'))); 
%     end;
%%    
%%%%%%%%%%%%%   Error Calculation
    TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)]; 
   
%     TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
%                  SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
%                  SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
%     TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];
    
    
%     error4(:,:,i) =  SE3_s(:,:,i)*inv(T);
    error6(:,:,i) = SE3_t(:,:,i)*inv(T);
%     error7(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
    
%     dv1(i) = norm(T-SE3_s(:,:,i),'fro');  %T - hatX
%     dv2(i) = norm(T-SE3_t(:,:,i),'fro');  %T - X
%     dv3(i) = norm(SE3_s(:,:,i)-SE3_t(:,:,i),'fro');
%     dv1(i) = norm(SE3_s(:,:,i)*inv(T)-eye(4),'fro');  %T - hatX
    dv3(i) = norm(SE3_t(:,:,i)*inv(T)-eye(4),'fro');  %T - X
%     dv3(i) = norm(SE3_s(:,:,i)*inv(SE3_t(:,:,i))-eye(4),'fro');
    omega1(:,:,i)=zeros(4,4);
    omega2(:,:,i)=zeros(4,4);
    omega3(:,:,i)=zeros(4,4);
    
    omega(:,:,i) =   round(1*projection(1/4*((error6(:,:,i) - error6(:,:,i).'))),2);
    omega4(:,:,i) =  round(inv(T)*projection(1/4*(Q*(error6(:,:,i) - eye(4))-(error6(:,:,i) - eye(4)).'*Q.'))*T,2);

%%   
%%%%%%%%%%%%%   Measurement Noise
%     deltay(:,:,i) = zeros(4,landmarks);  
%% 
%%%%%%%%%%%%%   Ui = U_trajecotry + U_avoid
    Ui(:,:,i) = -0.1*inv(T)*(omega(:,:,i))*(T)-0.1*omega4(:,:,i);
%%%%%%%%%%%%%   Ut = Ui + modling uncertianty      
    Ut(:,:,i) =  Ui(:,:,i);
%%   
%%%%%%%%%%%%%   System Trajectory
    SE3_t(:,:,i+1) = SE3_t(:,:,i)*expm(Ut(:,:,i));
%%   
%%%%%%%%%%%%%   Observer Convergence Trajectory
%     error5(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
%     for p = 1:lengthy
%         cframept(:,p) =  inv(SE3_t(:,:,i)) *  [y(1,p) y(2,p) y(3,p) 1].';
%          deltay(:,p,i) = 0*[0.3*cos(i*p/2) 0.2*sin(i*p/5) 0.1*sin(i*p/4) 0].';
%          omega2(:,:,i) = omega2(:,:,i)+0.001*(SE3_s(:,:,i)*(cframept(:,p)+deltay(:,p,i))-[y(1,p) y(2,p) y(3,p) 1].')*transpose([y(1,p) y(2,p) y(3,p) 1].');
% %        omega2(:,:,i) = omega2(:,:,i)+0.01*inv(SE3_s(:,:,i)).'*(inv(SE3_s(:,:,i))*y(:,p)-inv(SE3_t(:,:,i))*y(:,p))*transpose(y(:,p));
%     end;
%     omega_t(:,:,i) =   projection(1/4*((omega2(:,:,i) - omega2(:,:,i).')));
%     se3_s = Ui(:,:,i)-inv(SE3_s(:,:,i))*omega_t(:,:,i)*SE3_s(:,:,i); 
%   
%     SE3_s(:,:,i+1) = SE3_s(:,:,i)*expm(se3_s);
%%    
%%%%%%%%%%%%%   Two Measures    
%     dv1(i) = real(visiondistance(R,TRs(:,:,i),P.',TPs(:,i)));
%     dv2(i) = real(visiondistance(R,TRt(:,:,i),P.',TPt(:,i)));
%     dv3(i) = real(visiondistance(TRs(:,:,i),TRt(:,:,i),TPs(:,i),TPt(:,i)));
        %X - hatX 
%%   
%%%%%%%%%%%%%   Real Time Plot
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(1,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(1,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(1,3,i)],'--r');
%     hold on 
%     
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(2,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(2,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(2,3,i)],'--g');
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(3,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(3,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(3,3,i)],'--b');
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(1,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(1,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(1,3,i)],'--r');
%     hold on
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(2,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(2,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(2,3,i)],'--g');
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(3,3,i)],'--b');
% %     plot3(y(1,:),y(2,:),y(3,:),'o','Color','k','MarkerSize',5);
%     
%     plot3([T(1,4),T(1,4)+4*T(1,1)],[T(2,4),T(2,4)+4*T(1,2)],[T(3,4),T(3,4)+4*T(1,3)],'k','linewidth',7);
%      
%     
% plot3([T(1,4),T(1,4)+4*T(2,1)],[T(2,4),T(2,4)+4*T(2,2)],[T(3,4),T(3,4)+4*T(2,3)],'k','linewidth',7);
% plot3([T(1,4),T(1,4)+4*T(3,1)],[T(2,4),T(2,4)+4*T(3,2)],[T(3,4),T(3,4)+4*T(3,3)],'k','linewidth',7);
%     hold off  
%     pause(0.01)
%%
end;


L2_2 = sqrt(sum(dv3*dv3.')/iter);
Linf_2 = max(dv3);

% L2_3 = sqrt(sum(dv3*dv3.')/iter);
% Linf_3 = max(dv3);

% for i= 1:iter
% 
% trn3(i) = trace(projection(inv(error5(:,:,i))-inv(error5(:,:,i)).')*projection(inv(error5(:,:,i))-inv(error5(:,:,i)).').');
% trn4(i) = trace(inv(error5(:,:,i)));
% trn5(i)= trace(projection((inv(error5(:,:,i))*SE3_s(:,:,i)*inv(T))-(inv(error5(:,:,i))*SE3_s(:,:,i)*inv(T)).')*projection((SE3_s(:,:,i)*inv(T)-(SE3_s(:,:,i)*inv(T)).')).');
% 
% end

%%
%%%%%%%%%%%%%   Trajectory and Obstacle Plot
figure();
for i = 1:iter
    px1(i) = SE3_t(1,4,i);
    py1(i) = SE3_t(2,4,i);
    pz1(i) = SE3_t(3,4,i);
%     hpx1(i) = SE3_s(1,4,i);
%     hpy1(i) = SE3_s(2,4,i);
%     hpz1(i) = SE3_s(3,4,i);
end;    
 plot3(px1,py1,pz1,'b');  
 

